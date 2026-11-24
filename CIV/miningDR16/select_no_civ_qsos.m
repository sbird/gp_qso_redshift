%% Select candidate no-C IV sightlines from processed DR16 chunks.
% This is a probabilistic selection, not proof that a spectrum is C IV-free.
% Default: first search has p_no_c4 >= 0.95 and no saved search has
% p_c4 >= 0.95. The first requirement also excludes high-probability
% singlets because the three model probabilities sum to one.
% Run from the gpCIVIC directory; outputs are a MAT file and a CSV table.

chunkDir = fullfile(pwd,'processedFiles');
chunkPattern = 'processed_qsos_tst_posterior_corrected_S_*_E_*.mat';
catalogFile = '/home/reza/civic/REW-Work/catalog.mat';
outputDir = fullfile(pwd,'no_civ_selection');
minNullProbability = 0.95;
doubletThreshold = 0.95;

assert(isfinite(minNullProbability) && minNullProbability>=0 && minNullProbability<=1 && ...
    isfinite(doubletThreshold) && doubletThreshold>0 && doubletThreshold<=1, ...
    'Probability thresholds must lie in [0,1] and doubletThreshold must be positive.');
files=dir(fullfile(chunkDir,chunkPattern));
assert(~isempty(files),'No processed chunks found in %s.',chunkDir);
ranges=zeros(numel(files),2);
for k=1:numel(files)
    tok=regexp(files(k).name,'_S_(\d+)_E_(\d+)\.mat$','tokens','once');
    assert(~isempty(tok),'Invalid chunk name: %s',files(k).name);
    ranges(k,:)=[str2double(tok{1}) str2double(tok{2})];
end
[ranges,order]=sortrows(ranges,1); files=files(order);
assert(all(ranges(:,1)>=1 & ranges(:,2)>=ranges(:,1)) && ...
    all(ranges(2:end,1)>ranges(1:end-1,2)),'Invalid or overlapping chunk ranges.');
if ranges(1,1)~=1 || any(ranges(2:end,1)>ranges(1:end-1,2)+1)
    warning('Chunk coverage has gaps. Unprocessed QSOs are not selected.');
end
assert(isfile(catalogFile),'Cannot find catalog: %s',catalogFile);
cat=load(catalogFile,'all_QSO_ID_dr16');
assert(isfield(cat,'all_QSO_ID_dr16'),'Catalog is missing all_QSO_ID_dr16.');
qsoIds=cat.all_QSO_ID_dr16;
assert(iscellstr(qsoIds) || isstring(qsoIds), ...
    'Expected all_QSO_ID_dr16 to contain QSO names as a cell array or string array.');
qsoIds=string(qsoIds(:));
referenceTest=[]; results=cell(numel(files),1); nProcessed=0; nInvalid=0;
for k=1:numel(files)
    path=fullfile(files(k).folder,files(k).name);
    C=load(path,'test_ind','p_no_c4','p_c4L1','p_c4');
    assert(all(isfield(C,{'test_ind','p_no_c4','p_c4L1','p_c4'})), ...
        'Missing required variables in %s.',path);
    assert(islogical(C.test_ind) && isvector(C.test_ind) && numel(C.test_ind)==numel(qsoIds), ...
        'Invalid test_ind or catalog length in %s.',path);
    if isempty(referenceTest), referenceTest=C.test_ind(:); else
        assert(isequal(referenceTest,C.test_ind(:)),'test_ind differs in %s.',path);
    end
    parentRows=find(referenceTest); first=ranges(k,1); last=ranges(k,2);
    assert(last<=numel(parentRows),'Chunk range exceeds test sample: %s',path);
    expected=last-first+1;
    assert(ismatrix(C.p_no_c4) && size(C.p_no_c4,1)==expected && size(C.p_no_c4,2)>=1 && ...
        isequal(size(C.p_no_c4),size(C.p_c4L1),size(C.p_c4)), ...
        'Probability dimensions do not match chunk range in %s.',path);
    pn=C.p_no_c4; p1=C.p_c4L1; p2=C.p_c4;
    finiteAll=isfinite(pn)&isfinite(p1)&isfinite(p2);
    absent=isnan(pn)&isnan(p1)&isnan(p2);
    goodValues=finiteAll & pn>=0 & pn<=1 & p1>=0 & p1<=1 & p2>=0 & p2<=1 & ...
        abs(pn+p1+p2-1)<1e-6;
    validRow=goodValues(:,1) & all(goodValues | absent,2);
    % Stopped searches are NaN; they must never be treated as detections.
    detected=any(finiteAll & p2>=doubletThreshold,2);
    keep=validRow & pn(:,1)>=minNullProbability & ~detected;
    local=find(keep); testIndex=first+local-1; parentIndex=parentRows(testIndex);
    results{k}=table(testIndex(:),parentIndex(:),qsoIds(parentIndex),pn(local,1),p1(local,1),p2(local,1), ...
        'VariableNames',{'test_index','parent_index','qso_id','p_no_c4_first','p_c4L1_first','p_c4_first'});
    nProcessed=nProcessed+expected; nInvalid=nInvalid+nnz(~validRow);
    fprintf('%s: %d selected / %d rows (%d invalid probability rows excluded).\n', ...
        files(k).name,numel(local),expected,nnz(~validRow));
end
noCIV=vertcat(results{:});
selected_parent_mask=false(numel(qsoIds),1); selected_parent_mask(noCIV.parent_index)=true;
selected_test_indices=noCIV.test_index; selected_parent_indices=noCIV.parent_index;
selected_qso_ids=noCIV.qso_id; chunkRanges=ranges; chunkFiles={files.name};
if ~isfolder(outputDir), mkdir(outputDir); end
stem=sprintf('no_civ_pnull_%0.3f_pdoublet_%0.3f',minNullProbability,doubletThreshold);
matFile=fullfile(outputDir,[stem '.mat']); csvFile=fullfile(outputDir,[stem '.csv']);
assert(~isfile(matFile)&&~isfile(csvFile),'Output already exists; move it or choose different thresholds.');
save(matFile,'noCIV','selected_parent_mask','selected_test_indices','selected_parent_indices', ...
    'selected_qso_ids','minNullProbability','doubletThreshold','chunkRanges','chunkFiles', ...
    'nProcessed','nInvalid','catalogFile','-v7.3');
writetable(noCIV,csvFile);
fprintf('Selected %d / %d processed QSOs; excluded %d invalid rows.\n',height(noCIV),nProcessed,nInvalid);
fprintf('Saved %s\nSaved %s\n',matFile,csvFile);
