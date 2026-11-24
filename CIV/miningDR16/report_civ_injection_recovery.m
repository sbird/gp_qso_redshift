function report_civ_injection_recovery(truth,baselineSelected,processedFile,outputDir,matchKms)
%REPORT_CIV_INJECTION_RECOVERY Report detection at the injected redshift.
% Recovery uses L2 MAP redshift within matchKms and p_c4 > threshold. This
% is an outcome criterion; posterior REW bias and original detections are
% recorded but do not define recovery. All valid processed trials enter
% the denominator. Empty bins are NaN, not zero completeness.
assert(isscalar(matchKms)&&isfinite(matchKms)&&matchKms>0,'Invalid match velocity.');
R=load(processedFile,'p_c4','p_c4L1','p_no_c4','map_z_c4L2','test_ind', ...
    'posterior_REW_1548_draws');
need={'p_c4','p_no_c4','map_z_c4L2','test_ind','posterior_REW_1548_draws'};
assert(all(isfield(R,need)),'Missing processed recovery fields.');
n=height(truth); assert(isequal(find(R.test_ind(:)),truth.parent_index(:)), ...
    'Processed spectra do not align with injection truth.');
assert(size(R.p_c4,1)==n && isequal(size(R.p_c4),size(R.p_no_c4),size(R.map_z_c4L2)), ...
    'Invalid recovery array dimensions.');
probability_fields={'p_c4','p_no_c4'};
if isfield(R,'p_c4L1')
    assert(isequal(size(R.p_c4),size(R.p_c4L1)), ...
        'Invalid singlet posterior dimensions in legacy recovery output.');
    probability_fields={'p_c4','p_c4L1','p_no_c4'};
end
assert(isequal(baselineSelected.parent_index(:),truth.parent_index(:)), 'Baseline row order mismatch.');
assert(size(R.posterior_REW_1548_draws,1)==n && ...
    size(R.posterior_REW_1548_draws,3)==size(R.p_c4,2), 'Invalid REW draw dimensions.');
probabilityTolerance=1e-12; % allow only floating-point boundary roundoff
probability_sum=zeros(size(R.p_c4));
valid_probabilities=true(size(R.p_c4));
finiteAll=true(size(R.p_c4));
absent=true(size(R.p_c4));
for probabilityField=probability_fields
    probability=R.(probabilityField{1});
    probability_sum=probability_sum+probability;
    valid_probabilities=valid_probabilities & ...
        probability>=-probabilityTolerance & probability<=1+probabilityTolerance;
    finiteAll=finiteAll & isfinite(probability);
    absent=absent & isnan(probability);
end
valid=finiteAll(:,1) & valid_probabilities(:,1) & ...
    abs(probability_sum(:,1)-1)<1e-6;
assert(all(valid),'%d injections have invalid first-search posteriors; inspect processing before reporting.',nnz(~valid));
good=finiteAll & valid_probabilities & abs(probability_sum-1)<1e-6;
assert(all(good(:)|absent(:)),'Invalid later-search probabilities; inspect processing before reporting.');
% Treat boundary roundoff as exact 0 or 1 in downstream recovery decisions.
for probabilityField=probability_fields
    x=R.(probabilityField{1});
    x(x<0 & x>=-probabilityTolerance)=0;
    x(x>1 & x<=1+probabilityTolerance)=1;
    R.(probabilityField{1})=x;
end
assert(all(isfinite(R.map_z_c4L2(good)) & R.map_z_c4L2(good)>0), ...
    'A completed search has an invalid L2 MAP redshift.');
thresholds=[0.80 0.85 0.90 0.95]; nt=numel(thresholds);
matchedP=nan(n,1); matchedZ=nan(n,1); matchedSearch=nan(n,1); matchDv=nan(n,1);
matchedREW=nan(n,1);
baselineP=nan(n,1); baselineDv=nan(n,1);
for i=1:n
    dv=299792.458*abs(R.map_z_c4L2(i,:)-truth.z_abs(i))/(1+truth.z_abs(i));
    ok=isfinite(R.p_c4(i,:))&R.p_c4(i,:)>=0&R.p_c4(i,:)<=1&isfinite(dv)&dv<=matchKms;
    if any(ok)
        indices=find(ok); [matchedP(i),k]=max(R.p_c4(i,indices)); j=indices(k);
        matchedZ(i)=R.map_z_c4L2(i,j); matchedSearch(i)=j; matchDv(i)=dv(j);
        draws=R.posterior_REW_1548_draws(i,:,j);
        if all(isfinite(draws)), matchedREW(i)=median(draws); end
    end
    d0=299792.458*abs(baselineSelected.map_z_c4L2(i,:)-truth.z_abs(i))/(1+truth.z_abs(i));
    ok=isfinite(d0)&d0<=matchKms&isfinite(baselineSelected.p_c4(i,:));
    if any(ok), ids=find(ok); [baselineP(i),k]=max(baselineSelected.p_c4(i,ids)); baselineDv(i)=d0(ids(k)); end
end
recovered=false(n,nt); results=truth;
results.matched_p_c4=matchedP; results.matched_z=matchedZ; results.matched_search=matchedSearch;
results.match_dv_kms=matchDv; results.original_match_p_c4=baselineP; results.original_match_dv_kms=baselineDv;
results.matched_REW_1548_median=matchedREW;
results.REW_1548_median_minus_injected=matchedREW-truth.rew1548_A;
for t=1:nt
    recovered(:,t)=matchedP>thresholds(t);
    results.(sprintf('recovered_p%02d',round(100*thresholds(t))))=recovered(:,t);
end
% Bin edges include all parameter values of the configured prior.
zEdges=[0 1.4 1.8 2.2 2.6 3 3.5 4 4.5 5 6 Inf];
nEdges=[12.5 13 13.5 14 14.5 15 15.5 16.1];
wEdges=[0 .1 .2 .3 .4 .6 .9 1.2 1.5 2 3 5 Inf];
snrEdges=[2 3 5 10 20 50 Inf];
iz=discretize(truth.z_abs,zEdges); in=discretize(truth.logN,nEdges); iw=discretize(truth.rew1548_A,wEdges);
isnr=discretize(truth.snr,snrEdges);
assert(all(~isnan(iz)&~isnan(in)&~isnan(iw)&~isnan(isnr)),'Injection truth outside completeness bins.');
dims=[numel(zEdges)-1 numel(nEdges)-1 numel(wEdges)-1]; sub=[iz in iw];
inserted=accumarray(sub,1,dims,@sum,0); recoveredCounts=zeros([dims nt]); completeness=nan([dims nt]);
snrDims=[numel(zEdges)-1 numel(wEdges)-1 numel(snrEdges)-1]; snrSub=[iz iw isnr];
insertedZWSNR=accumarray(snrSub,1,snrDims,@sum,0);
recoveredZWSNR=zeros([snrDims nt]); completenessZWSNR=nan([snrDims nt]);
insertedZN=sum(inserted,3); insertedZW=reshape(sum(inserted,2),dims(1),dims(3));
completenessZN=nan(dims(1),dims(2),nt); completenessZW=nan(dims(1),dims(3),nt);
for t=1:nt
    rc=zeros(dims); if any(recovered(:,t)), rc=accumarray(sub(recovered(:,t),:),1,dims,@sum,0); end
    rs=zeros(snrDims); if any(recovered(:,t)), rs=accumarray(snrSub(recovered(:,t),:),1,snrDims,@sum,0); end
    recoveredZWSNR(:,:,:,t)=rs; xs=nan(snrDims); ok=insertedZWSNR>0;
    xs(ok)=rs(ok)./insertedZWSNR(ok); completenessZWSNR(:,:,:,t)=xs;
    recoveredCounts(:,:,:,t)=rc; x=nan(dims); good=inserted>0; x(good)=rc(good)./inserted(good); completeness(:,:,:,t)=x;
    rz=sum(rc,3); x=nan(size(insertedZN)); good=insertedZN>0; x(good)=rz(good)./insertedZN(good); completenessZN(:,:,t)=x;
    rw=reshape(sum(rc,2),dims(1),dims(3)); x=nan(size(insertedZW)); good=insertedZW>0; x(good)=rw(good)./insertedZW(good); completenessZW(:,:,t)=x;
end
save(fullfile(outputDir,'recovery_results.mat'),'results','thresholds','matchKms','zEdges','nEdges','wEdges','snrEdges', ...
    'inserted','recoveredCounts','completeness','insertedZN','insertedZW','completenessZN','completenessZW', ...
    'insertedZWSNR','recoveredZWSNR','completenessZWSNR','-v7.3');
% Tables cannot write cell arrays of variable-sized interval matrices to CSV.
flat=removevars(results,{'allowed_intervals','blocked_intervals','baseline_doublet_z','baseline_singlet_z'});
writetable(flat,fullfile(outputDir,'recovery_per_injection.csv'));
summary=table(thresholds(:),repmat(n,nt,1),sum(recovered,1).',mean(recovered,1).', ...
    'VariableNames',{'p_threshold','injected','recovered','recovery_fraction'});
writetable(summary,fullfile(outputDir,'recovery_summary.csv'));
fig=figure('Visible','off','Color','w','Position',[100 100 1400 1500]); tiledlayout(fig,nt,2,'TileSpacing','compact');
zLabels=arrayfun(@(a,b)sprintf('%.1f-%.1f',a,b),zEdges(1:end-1),zEdges(2:end),'UniformOutput',false);
nLabels=arrayfun(@(a,b)sprintf('%.1f-%.1f',a,b),nEdges(1:end-1),nEdges(2:end),'UniformOutput',false);
wLabels=arrayfun(@(a,b)sprintf('%.1f-%.1f',a,b),wEdges(1:end-1),wEdges(2:end),'UniformOutput',false);
for t=1:nt
    ax=nexttile; h=imagesc(ax,completenessZN(:,:,t).'); h.AlphaData=isfinite(completenessZN(:,:,t).'); axis(ax,'xy'); clim(ax,[0 1]); colorbar(ax);
    xticks(ax,1:dims(1)); xticklabels(ax,zLabels); xtickangle(ax,45); yticks(ax,1:dims(2)); yticklabels(ax,nLabels);
    xlabel(ax,'Injected absorber redshift'); ylabel(ax,'log N (cm^{-2})'); title(ax,sprintf('p > %.2f',thresholds(t)));
    ax=nexttile; h=imagesc(ax,completenessZW(:,:,t).'); h.AlphaData=isfinite(completenessZW(:,:,t).'); axis(ax,'xy'); clim(ax,[0 1]); colorbar(ax);
    xticks(ax,1:dims(1)); xticklabels(ax,zLabels); xtickangle(ax,45); yticks(ax,1:dims(3)); yticklabels(ax,wLabels);
    xlabel(ax,'Injected absorber redshift'); ylabel(ax,'W_{1548} (Angstrom)'); title(ax,sprintf('p > %.2f',thresholds(t)));
end
exportgraphics(fig,fullfile(outputDir,'completeness_surfaces.png'),'Resolution',180); close(fig);
fprintf('Recovery summary (match %.0f km/s):\n',matchKms); disp(summary);
fprintf('These projections reflect the sampled sightlines and injected parameter distribution.\n');
end
