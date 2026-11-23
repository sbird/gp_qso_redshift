function baseline = load_civ_injection_baseline(chunkDir, testMask)
%LOAD_CIV_INJECTION_BASELINE Read original results in filtered-test order.
files=dir(fullfile(chunkDir,'processed_qsos_tst_posterior_corrected_S_*_E_*.mat'));
assert(~isempty(files),'No original processed chunks in %s.',chunkDir);
ranges=zeros(numel(files),2);
for k=1:numel(files)
    tok=regexp(files(k).name,'_S_(\d+)_E_(\d+)\.mat$','tokens','once');
    assert(~isempty(tok),'Invalid filename %s.',files(k).name);
    ranges(k,:)=[str2double(tok{1}) str2double(tok{2})];
end
[ranges,order]=sortrows(ranges,1); files=files(order);
assert(all(ranges(:,1)>=1 & ranges(:,2)>=ranges(:,1)) && ...
    all(ranges(2:end,1)>ranges(1:end-1,2)),'Invalid or overlapping original chunks.');
assert(islogical(testMask)&&isvector(testMask),'testMask must be logical.');
parent=find(testMask(:)); n=numel(parent);
baseline=struct('parent_index',parent,'processed',false(n,1),'valid',false(n,1), ...
    'zmin',nan(n,1),'zmax',nan(n,1),'chunk_ranges',ranges,'chunk_files',{{files.name}});
fields={'p_c4','p_c4L1','p_no_c4','map_z_c4L2','map_z_c4L1'};
for k=1:numel(fields), baseline.(fields{k})=nan(n,7); end
for k=1:numel(files)
    path=fullfile(files(k).folder,files(k).name);
    C=load(path,'test_ind','min_z_c4s','max_z_c4s',fields{:});
    assert(all(isfield(C,[{'test_ind','min_z_c4s','max_z_c4s'},fields])), ...
        'Required variable missing in %s.',path);
    assert(islogical(C.test_ind)&&isequal(C.test_ind(:),testMask(:)), ...
        'Chunk test_ind differs from current science mask in %s.',path);
    rows=(ranges(k,1):ranges(k,2)).'; count=numel(rows);
    assert(rows(end)<=n && numel(C.min_z_c4s)==count && numel(C.max_z_c4s)==count, ...
        'Invalid original chunk dimensions in %s.',path);
    for j=1:numel(fields)
        x=C.(fields{j}); assert(ismatrix(x)&&size(x,1)==count&&size(x,2)==7, ...
            'Invalid %s dimensions in %s.',fields{j},path);
        baseline.(fields{j})(rows,:)=x;
    end
    pn=C.p_no_c4; p1=C.p_c4L1; p2=C.p_c4;
    allFinite=isfinite(pn)&isfinite(p1)&isfinite(p2);
    allAbsent=isnan(pn)&isnan(p1)&isnan(p2);
    validP=allFinite & pn>=0 & pn<=1 & p1>=0 & p1<=1 & p2>=0 & p2<=1 & abs(pn+p1+p2-1)<1e-6;
    validZ=(isfinite(C.map_z_c4L1)&C.map_z_c4L1>0 & ...
        isfinite(C.map_z_c4L2)&C.map_z_c4L2>0) | allAbsent;
    baseline.valid(rows)=validP(:,1)&all((validP&validZ)|allAbsent,2);
    baseline.zmin(rows)=C.min_z_c4s(:); baseline.zmax(rows)=C.max_z_c4s(:);
    baseline.processed(rows)=true;
end
if ranges(1,1)~=1 || any(ranges(2:end,1)>ranges(1:end-1,2)+1) || ranges(end,2)<n
    warning('Only supplied original chunks can be used for injections.');
end
end
