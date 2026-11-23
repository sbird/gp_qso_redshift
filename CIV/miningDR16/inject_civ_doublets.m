function [injected,truth,info] = inject_civ_doublets(preloaded,zQSO,baseline,nInject,varargin)
%INJECT_CIV_DOUBLETS Inject into real spectra away from original detections.
% The input parent arrays and baseline must correspond to the same catalog.
% One doublet is injected per spectrum. Known L1 and L2 MAP features with
% p > DetectionThreshold are excluded by ExclusionKms from both injected
% lines; originals are not masked. Unprocessed/invalid rows are ineligible.
% Default selection is approximately balanced across QSO-z/SNR bins. Output
% truth includes sampling probability and the full allowed/blocked intervals.
% Noise is retained; a running median estimates the underlying signal for an
% additive absorption decrement. This is an approximation, not a clean
% reconstruction of quasar continua or blended absorbers.
p=inputParser;
p.addParameter('Seed',42,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x==fix(x));
p.addParameter('DetectionThreshold',0.75,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('ExclusionKms',1000,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<299792.458);
p.addParameter('SNRThreshold',3,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
p.addParameter('MinZAbs',-Inf,@(x)isnumeric(x)&&isscalar(x)&&~isnan(x));
p.addParameter('MinREW1548',-Inf,@(x)isnumeric(x)&&isscalar(x)&&~isnan(x));
p.addParameter('ZEdges',[0 1.7 2 2.5 3 4 5 Inf],@(x)isnumeric(x)&&isvector(x)&&all(diff(x)>0));
p.addParameter('SNREdges',[2 3 5 10 20 50 Inf],@(x)isnumeric(x)&&isvector(x)&&all(diff(x)>0));
p.addParameter('ContinuumPixels',101,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=3&&x==fix(x));
p.parse(varargin{:}); opts=p.Results;
assert(isnumeric(nInject)&&isscalar(nInject)&&isfinite(nInject)&&nInject>=1&&nInject==fix(nInject),'Invalid injection count.');
assert(isnumeric(zQSO)&&isvector(zQSO),'zQSO must be a vector.'); zQSO=zQSO(:); nParent=numel(zQSO);
fields={'all_wavelengths','all_flux','all_noise_variance','all_pixel_mask','all_sigma_pixel'};
assert(isstruct(preloaded)&&isscalar(preloaded)&&all(isfield(preloaded,fields)),'Missing preloaded fields.');
for k=1:numel(fields)
    assert(iscell(preloaded.(fields{k}))&&numel(preloaded.(fields{k}))==nParent,'Invalid %s.',fields{k});
end
assert(exist('voigt_iP_transition','file')==3,'Compile voigt_iP_transition MEX first.');
% Declare script-supplied names in this function workspace. In particular,
% MATLAB otherwise resolves `width` as the built-in table width function.
width=[]; nAVG=[]; min_lambda=[]; max_lambda=[]; pixel_spacing=[];
fit_min_log_nciv=[]; fit_max_log_nciv=[]; min_sigma=[]; max_sigma=[];
civ_1548_wavelength=[]; civ_1550_wavelength=[];
set_parameters_dr16
assert(isscalar(width)&&width==3,'MEX padding assumes width 3.');
oldRng=rng; cleanup=onCleanup(@()rng(oldRng)); %#ok<NASGU>
rng(opts.Seed,'twister');
parent=baseline.parent_index(:); n=numel(parent);
assert(all(parent>=1&parent<=nParent)&&numel(unique(parent))==n,'Invalid baseline parent mapping.');
snr=nan(n,1); allowed=cell(n,1); blocked=cell(n,1); lengthAllowed=zeros(n,1);
knownDoublet=cell(n,1); knownSinglet=cell(n,1);
for i=find(baseline.processed & baseline.valid).'
    q=parent(i); w=preloaded.all_wavelengths{q}(:); f=preloaded.all_flux{q}(:);
    v=preloaded.all_noise_variance{q}(:); m=preloaded.all_pixel_mask{q}(:); d=preloaded.all_sigma_pixel{q}(:);
    if ~isfinite(zQSO(q))||numel(w)<3||any([numel(f) numel(v) numel(m) numel(d)]~=numel(w))|| ...
        any(~isfinite(w))||any(diff(w)<=0)||~isreal(v)||any(v<0), continue; end
    snr(i)=median(f./sqrt(v));
    if ~isfinite(snr(i))||snr(i)<=opts.SNRThreshold, continue; end
    good=w/(1+zQSO(q))>=min_lambda & w/(1+zQSO(q))<=max_lambda & d>0;
    if nnz(good)<3||nnz(good & ~logical(m))<3||any(~isfinite(d(good))), continue; end
    zlo=baseline.zmin(i); zhi=baseline.zmax(i);
    if ~isfinite(zlo)||~isfinite(zhi)||zlo<0||zhi<=zlo, continue; end
    sel2=isfinite(baseline.p_c4(i,:))&baseline.p_c4(i,:)>opts.DetectionThreshold;
    sel1=isfinite(baseline.p_c4L1(i,:))&baseline.p_c4L1(i,:)>opts.DetectionThreshold;
    zd=baseline.map_z_c4L2(i,sel2); zs=baseline.map_z_c4L1(i,sel1);
    [allowed{i},blocked{i}]=civ_allowed_injection_intervals(zlo,zhi,zd,zs,opts.ExclusionKms);
    if ~isempty(allowed{i})
        allowed{i}(:,1)=max(allowed{i}(:,1),opts.MinZAbs);
        allowed{i}=allowed{i}(allowed{i}(:,2)>allowed{i}(:,1),:);
    end
    knownDoublet{i}=zd; knownSinglet{i}=zs;
    lengthAllowed(i)=sum(allowed{i}(:,2)-allowed{i}(:,1));
end
zb=discretize(zQSO(parent),opts.ZEdges); sb=discretize(snr,opts.SNREdges);
valid=baseline.processed&baseline.valid&isfinite(snr)&lengthAllowed>0&~isnan(zb)&~isnan(sb);
candidates=find(valid); assert(numel(candidates)>=nInject,'Only %d eligible spectra; %d requested.',numel(candidates),nInject);
strata=sub2ind([numel(opts.ZEdges)-1,numel(opts.SNREdges)-1],zb(candidates),sb(candidates));
labels=unique(strata); groups=cell(numel(labels),1); cursor=ones(numel(labels),1);
for k=1:numel(labels), x=candidates(strata==labels(k)); groups{k}=x(randperm(numel(x))); end
selected=zeros(nInject,1); chosenStratum=zeros(nInject,1); count=0;
while count<nInject
    for k=randperm(numel(labels))
        if cursor(k)<=numel(groups{k})
            count=count+1; selected(count)=groups{k}(cursor(k)); chosenStratum(count)=k; cursor(k)=cursor(k)+1;
            if count==nInject, break; end
        end
    end
end
stratumAvailable=cellfun(@numel,groups); stratumSelected=cursor-1;
samplingProbability=stratumSelected(chosenStratum)./stratumAvailable(chosenStratum);
for k=1:numel(fields), injected.(fields{k})=preloaded.(fields{k})(parent(selected)); end
injected.parent_index=parent(selected); injected.z_qso=zQSO(parent(selected));
zAbs=zeros(nInject,1);
for q=1:nInject
    intervals=allowed{selected(q)}; widths=intervals(:,2)-intervals(:,1);
    u=rand*sum(widths); j=find(u<cumsum(widths),1); before=sum(widths(1:j-1));
    zAbs(q)=intervals(j,1)+u-before;
end
logN=nan(nInject,1); sigma=nan(nInject,1);
rew1548=nan(nInject,1); rew1550=nan(nInject,1); cov1548=zeros(nInject,1); cov1550=zeros(nInject,1);
for q=1:nInject
    w=double(injected.all_wavelengths{q}(:)); f=injected.all_flux{q}; originalSize=size(f); f=double(f(:));
    v=double(injected.all_noise_variance{q}(:)); mask=logical(injected.all_pixel_mask{q}(:));
    dispix=double(injected.all_sigma_pixel{q}(:)); zq=injected.z_qso(q);
    model=w/(1+zq)>=min_lambda & w/(1+zq)<=max_lambda & dispix>0;
    wm=w(model); dm=dispix(model); factor=nAVG+1; pad=width*factor;
    wf=finer(wm,nAVG).'; df=finer(dm,nAVG).';
    wp=[logspace(log10(min(wm))-width*pixel_spacing,log10(min(wm))-pixel_spacing/factor,pad).'; ...
        wf; logspace(log10(max(wm))+pixel_spacing/factor,log10(max(wm))+width*pixel_spacing,pad).'];
    dp=[repmat(dm(1),pad,1);df;repmat(dm(end),pad,1)];
    % Measure the intrinsic REW on a line-centered grid. The QSO modeling
    % window can truncate a line near its edge, but that truncation should
    % affect recovery, not the injected absorber's true equivalent width.
    ew_rest_grid=linspace(civ_1548_wavelength-40,civ_1548_wavelength+40,4001).';
    ew_obs_core=ew_rest_grid*(1+zAbs(q)); ew_step=ew_obs_core(2)-ew_obs_core(1);
    ew_pad=width*factor;
    ew_wavelengths=[ew_obs_core(1)-ew_step*(ew_pad:-1:1).';ew_obs_core; ...
        ew_obs_core(end)+ew_step*(1:ew_pad).'];
    ew_wdisp=nan(size(ew_wavelengths));
    % For a targeted high-REW run, draw from the existing logN/sigma priors
    % until the intrinsic 1548-A REW meets the requested floor.
    accepted=false; maxCandidateREW=-Inf;
    for attempt=1:10000
        logN(q)=fit_min_log_nciv+(fit_max_log_nciv-fit_min_log_nciv)*rand;
        sigma(q)=min_sigma+(max_sigma-min_sigma)*rand;
        a1=voigt_iP_transition(wp,zAbs(q),10^logN(q),1,sigma(q),dp,factor);
        ew_profile=voigt_iP_transition(ew_wavelengths,zAbs(q),10^logN(q),1, ...
            sigma(q),ew_wdisp,factor);
        candidateREW=trapz(ew_rest_grid,1-ew_profile);
        if isfinite(candidateREW), maxCandidateREW=max(maxCandidateREW,candidateREW); end
        if isfinite(candidateREW)&&all(isfinite(ew_profile))&&candidateREW>opts.MinREW1548
            accepted=true; break
        end
    end
    assert(accepted,['Could not draw intrinsic REW1548 > %.4g A after 10000 trials ', ...
        'for parent row %d at z_abs=%.5f (maximum trial REW %.4g A).'], ...
        opts.MinREW1548,injected.parent_index(q),zAbs(q),maxCandidateREW);
    profile=voigt_iP_transition(wp,zAbs(q),10^logN(q),3,sigma(q),dp,factor);
    a2=voigt_iP_transition(wp,zAbs(q),10^logN(q),2,sigma(q),dp,factor);
    assert(numel(profile)==numel(wf)&&all(isfinite(profile))&&all(isfinite(a1))&&all(isfinite(a2)), ...
        'Invalid Voigt profile for parent row %d.',injected.parent_index(q));
    rew1548(q)=candidateREW; rew1550(q)=trapz(wf,1-a2)/(1+zAbs(q));
    good=~mask&isfinite(f)&isfinite(v)&v>0; continuumInput=f; continuumInput(~good)=NaN;
    continuum=movmedian(continuumInput,opts.ContinuumPixels,'omitnan');
    nativeProfile=profile(1:factor:numel(wf)); modelRows=find(model);
    affectedLocal=~mask(model)&abs(1-nativeProfile)>1e-10; affected=modelRows(affectedLocal);
    assert(all(isfinite(continuum(affected)) & continuum(affected)>0), ...
        'Invalid continuum in parent row %d; inspect before proceeding.',injected.parent_index(q));
    f(affected)=f(affected)-continuum(affected).*(1-nativeProfile(affectedLocal));
    injected.all_flux{q}=reshape(f,originalSize);
    dz=(300/299792.458)*(1+zAbs(q));
    cov1548(q)=nnz(good&model&abs(w/civ_1548_wavelength-1-zAbs(q))<=dz);
    cov1550(q)=nnz(good&model&abs(w/civ_1550_wavelength-1-zAbs(q))<=dz);
end
truth=table(selected,parent(selected),zQSO(parent(selected)),snr(selected),zb(selected),sb(selected), ...
    baseline.zmin(selected),baseline.zmax(selected),lengthAllowed(selected),samplingProbability, ...
    zAbs,logN,sigma,rew1548,rew1550,cov1548,cov1550,allowed(selected),blocked(selected), ...
    knownDoublet(selected),knownSinglet(selected), ...
    'VariableNames',{'test_index','parent_index','z_qso','snr','z_stratum','snr_stratum', ...
    'search_z_min','search_z_max','allowed_z_length','spectrum_selection_probability', ...
    'z_abs','logN','sigma_cm_s','rew1548_A','rew1550_A','coverage1548_pixels','coverage1550_pixels', ...
    'allowed_intervals','blocked_intervals','baseline_doublet_z','baseline_singlet_z'});
info=struct('options',opts,'n_processed',nnz(baseline.processed),'n_valid_baseline',nnz(baseline.valid), ...
    'n_eligible',numel(candidates),'stratum_labels',labels,'stratum_available',stratumAvailable, ...
    'stratum_selected',stratumSelected,'z_distribution','uniform in allowed redshift length per selected sightline', ...
    'parameter_distribution','uniform logN and sigma','noise_model','additive continuum estimate; noise unchanged', ...
    'profile_model','voigt_iP_transition, native WDISP','nAVG',nAVG, ...
    'chunk_ranges',baseline.chunk_ranges,'chunk_files',{baseline.chunk_files});
end
