function report_civ_zrefine_pilot(truth,processedFile,outputDir,matchKms)
%REPORT_CIV_ZREFINE_PILOT Compare broad and conditional local injection fits.
% This reports only the search matching the known injected redshift. It does
% not change p_CIV, and is a diagnostic rather than a selection correction.
R=load(processedFile,'p_c4','map_z_c4L2','test_ind','refine_status', ...
    'refine_z_min','refine_z_max','refine_ess','refine_max_logL', ...
    'refine_z_draws','refine_REW_1548_draws','posterior_REW_1548_draws');
assert(isequal(find(R.test_ind(:)),truth.parent_index(:)),'Truth and processing rows differ.');
n=height(truth); slot=nan(n,1); p=nan(n,1); broadREW=nan(n,1); localREW=nan(n,1);
localESS=nan(n,1); localLo=nan(n,1); localHi=nan(n,1); truthInside=false(n,1);
localMedianZ=nan(n,1); localEdgeFraction=nan(n,1); status=zeros(n,1);
for i=1:n
    v=299792.458*abs(R.map_z_c4L2(i,:)-truth.z_abs(i))/(1+truth.z_abs(i));
    valid=isfinite(v)&v<=matchKms&isfinite(R.p_c4(i,:));
    if ~any(valid), continue; end
    ids=find(valid); [p(i),j]=max(R.p_c4(i,ids)); slot(i)=ids(j); k=slot(i);
    draws=R.posterior_REW_1548_draws(i,:,k);
    if all(isfinite(draws)), broadREW(i)=median(draws); end
    status(i)=R.refine_status(i,k);
    if status(i)==1
        localLo(i)=R.refine_z_min(i,k); localHi(i)=R.refine_z_max(i,k);
        localESS(i)=R.refine_ess(i,k); truthInside(i)=truth.z_abs(i)>=localLo(i)&&truth.z_abs(i)<=localHi(i);
        draws=R.refine_REW_1548_draws(i,:,k); zs=R.refine_z_draws(i,:,k);
        assert(all(isfinite(draws))&&all(isfinite(zs)),'Invalid refined draws for parent %d.',truth.parent_index(i));
        localREW(i)=median(draws); localMedianZ(i)=median(zs);
        width=localHi(i)-localLo(i);
        localEdgeFraction(i)=mean(zs<=localLo(i)+.05*width | zs>=localHi(i)-.05*width);
    end
end
comparison=table(truth.parent_index,truth.z_abs,truth.rew1548_A,slot,p,broadREW,localREW, ...
    localREW-truth.rew1548_A,localESS,localLo,localHi,truthInside,localMedianZ,localEdgeFraction,status, ...
    'VariableNames',{'parent_index','injected_z','injected_REW_1548','matched_search','broad_p_c4', ...
    'broad_REW_median','local_REW_median','local_REW_minus_injected','local_ESS','local_z_min', ...
    'local_z_max','truth_inside_local_window','local_z_median','fraction_draws_near_window_edges','local_status'});
writetable(comparison,fullfile(outputDir,'zrefine_pilot_comparison.csv'));
save(fullfile(outputDir,'zrefine_pilot_comparison.mat'),'comparison','matchKms','-v7.3');
matched=isfinite(slot); success=status==1;
fprintf('Matched broad redshift: %d/%d; local refinement succeeded: %d/%d.\n',nnz(matched),n,nnz(success),n);
fprintf('Injected redshift outside local window: %d/%d successful local fits.\n',nnz(success&~truthInside),nnz(success));
if any(success)
    fprintf('Median local ESS: %.3g; below 2: %d/%d.\n',median(localESS(success)),nnz(success&localESS<2),nnz(success));
    fprintf('Median local REW minus injection: %.4f Angstrom.\n',median(localREW(success)-truth.rew1548_A(success)));
end
fprintf('Saved %s\n',fullfile(outputDir,'zrefine_pilot_comparison.csv'));
end
