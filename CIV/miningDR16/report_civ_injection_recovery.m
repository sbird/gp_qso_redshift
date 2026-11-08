function report_civ_injection_recovery(truth, processedFile, outputDir, matchKms)
%REPORT_CIV_INJECTION_RECOVERY Match L2 MAP redshifts and summarize recovery.
% One truth row per processed row, with a 300 km/s default matching radius.
% All injections enter the denominator, including poor-coverage trials.
% Counts are stratified by true z, log10(N/cm^-2) and W_r,1548 (Angstrom).
% The 2-D projections marginalize over this pilot's injection distribution.
if nargin<4, matchKms=300; end
assert(isfinite(matchKms)&&isscalar(matchKms)&&matchKms>0,'Invalid matching tolerance.');
assert(isfolder(outputDir)&&isfile(processedFile),'Missing injection output or processed file.');
required={'z_abs','logN','rew1548_A','parent_index'};
assert(istable(truth)&&all(ismember(required,truth.Properties.VariableNames)), ...
    'Truth table lacks required fields.');
R=load(processedFile,'p_c4','map_z_c4L2','test_ind');
assert(all(isfield(R,{'p_c4','map_z_c4L2','test_ind'})), 'Processed file lacks recovery variables.');
assert(isequal(size(R.p_c4),size(R.map_z_c4L2)) && size(R.p_c4,1)==height(truth), ...
    'Processed rows do not align with truth table.');
assert(isequal(find(R.test_ind(:)),truth.parent_index(:)), ...
    'Processed parent order does not match truth table.');
assert(all(isfinite(truth.z_abs)&truth.z_abs>=0&isfinite(truth.logN)& ...
    isfinite(truth.rew1548_A)&truth.rew1548_A>=0), 'Invalid injection truth values.');
thresholds=[0.80 0.85 0.90 0.95]; n=height(truth); nt=numel(thresholds);
matchedP=nan(n,1); matchedZ=nan(n,1); matchedSearch=nan(n,1); matchedDv=nan(n,1);
for i=1:n
    p=R.p_c4(i,:); z=R.map_z_c4L2(i,:);
    dv=299792.458*abs(z-truth.z_abs(i))/(1+truth.z_abs(i));
    valid=isfinite(p)&p>=0&p<=1&isfinite(z)&dv<=matchKms;
    if any(valid)
        slots=find(valid); [matchedP(i),j]=max(p(slots)); slot=slots(j);
        matchedZ(i)=z(slot); matchedSearch(i)=slot; matchedDv(i)=dv(slot);
    end
end
recovered=false(n,nt);
for t=1:nt, recovered(:,t)=matchedP>thresholds(t); end
results=truth;
results.matched_search=matchedSearch; results.matched_z=matchedZ;
results.match_dv_kms=matchedDv; results.matched_p_c4=matchedP;
for t=1:nt
    results.(sprintf('recovered_p%02d',round(100*thresholds(t))))=recovered(:,t);
end
% Broad pilot bins; empty cells remain NaN and counts are saved for auditing.
zEdges=[0 1.4 1.8 2.2 2.6 3 3.5 4 4.5 5 6 Inf];
logNEdges=12.5:0.4:16.1;
rewEdges=[0 .1 .2 .3 .4 .6 .9 1.2 1.5 2 3 5 Inf];
iz=discretize(truth.z_abs,zEdges); in=discretize(truth.logN,logNEdges);
iw=discretize(truth.rew1548_A,rewEdges);
assert(all(~isnan(iz)&~isnan(in)&~isnan(iw)), 'A truth value falls outside completeness bin edges.');
dims=[numel(zEdges)-1,numel(logNEdges)-1,numel(rewEdges)-1];
sub=[iz in iw]; inserted=accumarray(sub,1,dims,@sum,0);
recoveredCounts=zeros([dims nt]); completeness=nan([dims nt]);
insertedZN=sum(inserted,3); insertedZW=reshape(sum(inserted,2),dims(1),dims(3));
recoveredZN=zeros(dims(1),dims(2),nt); recoveredZW=zeros(dims(1),dims(3),nt);
completenessZN=nan(size(recoveredZN)); completenessZW=nan(size(recoveredZW));
for t=1:nt
    if any(recovered(:,t))
        rc=accumarray(sub(recovered(:,t),:),1,dims,@sum,0);
    else
        rc=zeros(dims);
    end
    recoveredCounts(:,:,:,t)=rc;
    temp=nan(dims); use=inserted>0; temp(use)=rc(use)./inserted(use); completeness(:,:,:,t)=temp;
    rzn=sum(rc,3); rzw=reshape(sum(rc,2),dims(1),dims(3));
    recoveredZN(:,:,t)=rzn; recoveredZW(:,:,t)=rzw;
    temp=nan(size(insertedZN)); use=insertedZN>0; temp(use)=rzn(use)./insertedZN(use); completenessZN(:,:,t)=temp;
    temp=nan(size(insertedZW)); use=insertedZW>0; temp(use)=rzw(use)./insertedZW(use); completenessZW(:,:,t)=temp;
end
save(fullfile(outputDir,'recovery_results.mat'),'results','thresholds','matchKms', ...
    'zEdges','logNEdges','rewEdges','inserted','recoveredCounts','completeness', ...
    'insertedZN','insertedZW','recoveredZN','recoveredZW','completenessZN','completenessZW', ...
    'processedFile','-v7.3');
writetable(results,fullfile(outputDir,'recovery_per_injection.csv'));
summary=table(thresholds(:),repmat(n,nt,1),sum(recovered,1).',mean(recovered,1).', ...
    'VariableNames',{'p_threshold','injected','recovered','recovery_fraction'});
writetable(summary,fullfile(outputDir,'recovery_summary.csv'));
% Plot cells by bin index so the open-ended upper bins remain visible.
fig=figure('Visible','off','Color','w','Position',[100 100 1400 1500]);
tiledlayout(fig,nt,2,'TileSpacing','compact');
zLabels=compose('%.1f-%.1f',zEdges(1:end-1),zEdges(2:end));
nLabels=compose('%.1f-%.1f',logNEdges(1:end-1),logNEdges(2:end));
wLabels=compose('%.1f-%.1f',rewEdges(1:end-1),rewEdges(2:end));
for t=1:nt
    ax=nexttile; imagesc(ax,completenessZN(:,:,t).'); axis(ax,'xy'); clim(ax,[0 1]); colorbar(ax);
    xticks(ax,1:dims(1)); xticklabels(ax,zLabels); xtickangle(ax,45);
    yticks(ax,1:dims(2)); yticklabels(ax,nLabels);
    xlabel(ax,'Injected z'); ylabel(ax,'Injected log_{10} N (cm^{-2})');
    title(ax,sprintf('p > %.2f: z vs log N',thresholds(t)));
    ax=nexttile; imagesc(ax,completenessZW(:,:,t).'); axis(ax,'xy'); clim(ax,[0 1]); colorbar(ax);
    xticks(ax,1:dims(1)); xticklabels(ax,zLabels); xtickangle(ax,45);
    yticks(ax,1:dims(3)); yticklabels(ax,wLabels);
    xlabel(ax,'Injected z'); ylabel(ax,'Injected W_{1548} (Angstrom)');
    title(ax,sprintf('p > %.2f: z vs REW',thresholds(t)));
end
exportgraphics(fig,fullfile(outputDir,'completeness_surfaces.png'),'Resolution',180);
close(fig);
fprintf('Recovery summary (tolerance %.0f km/s):\n',matchKms); disp(summary);
fprintf('Saved per-injection results, counts, completeness and plots in %s\n',outputDir);
fprintf('Empty bins are NaN; projections depend on the injection parameter distribution.\n');
end
