function [injected, truth, info] = inject_civ_doublets(preloaded, zQSO, nInject, varargin)
%INJECT_CIV_DOUBLETS Inject one simulated C IV doublet per selected spectrum.
%
% [injected, truth, info] = inject_civ_doublets(preloaded,zQSO,nInject)
% preloaded is a scalar struct containing full parent-catalog cell arrays:
% all_wavelengths, all_flux, all_noise_variance, all_pixel_mask and
% all_sigma_pixel (WDISP, in native pixels). zQSO is the corresponding
% parent-catalog quasar redshift vector. The original inputs are unchanged.
%
% Name-value options:
%   'Eligible'       logical parent-catalog mask (default: all rows)
%   'Seed'           nonnegative integer random seed (default: 1)
%   'ZEdges'         QSO redshift stratification edges (default: [1.7 2 2.5 3 4 Inf])
%   'SNREdges'       SNR stratification edges (default: [3 5 10 20 Inf])
%   'ContinuumPixels' running-median window in native pixels (default: 101)
%
% Sampling is without replacement, approximately balanced across occupied
% zQSO/SNR cells. SNR is median(flux./sqrt(noise_variance)), as in drawer.
% log10(N/cm^-2), sigma (cm/s) and absorber redshift are drawn uniformly
% over the ranges defined in set_parameters_dr16. The search-window limits
% follow process_qsos_dr16_par_fixed. No doublet-coverage cut is imposed;
% masked or poorly covered injections remain in the truth table.
%
% The injected signal is flux - continuum_estimate*(1-transmission), with
% a running median as continuum estimate. Existing noise and masks remain
% unchanged. This approximates an additive signal on noisy data; it is not
% an exact simulation of continuum structure, blending or photon noise.
% The MEX profile uses the same WDISP and fine-grid convention as the GP
% processor. Returned REWs integrate the convolved profiles over the modeled
% wavelength range, so edge truncation is possible and must be inspected.
%
% Output fields in injected are selected cell arrays, in truth row order.
% parent_index maps each selected row to the original parent catalog.
% This function does not run recovery or overwrite any data on disk.

p=inputParser;
p.addParameter('Eligible',true(numel(zQSO),1),@(x)islogical(x)&&isvector(x));
p.addParameter('Seed',1,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x==fix(x));
p.addParameter('ZEdges',[1.7 2 2.5 3 4 Inf],@(x)isnumeric(x)&&isvector(x)&&all(diff(x)>0));
p.addParameter('SNREdges',[3 5 10 20 Inf],@(x)isnumeric(x)&&isvector(x)&&all(diff(x)>0));
p.addParameter('ContinuumPixels',101,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=3&&x==fix(x));
p.parse(varargin{:}); opts=p.Results;
assert(isnumeric(zQSO)&&isvector(zQSO),'zQSO must be a numeric vector.');
assert(isnumeric(nInject)&&isscalar(nInject)&&isfinite(nInject)&&nInject>=1&&nInject==fix(nInject),'nInject must be a positive integer.');
zQSO=zQSO(:); nParent=numel(zQSO);
assert(numel(opts.Eligible)==nParent,'Eligible must match zQSO.');
fields={'all_wavelengths','all_flux','all_noise_variance','all_pixel_mask','all_sigma_pixel'};
assert(isstruct(preloaded)&&isscalar(preloaded)&&all(isfield(preloaded,fields)),'preloaded is missing required fields.');
for k=1:numel(fields)
    assert(iscell(preloaded.(fields{k}))&&numel(preloaded.(fields{k}))==nParent, ...
        '%s must be a cell array matching zQSO.',fields{k});
end
assert(exist('voigt_iP_transition','file')==3,'Compile voigt_iP_transition MEX before injecting.');
% The parameter file is a script; keep its variables in this function workspace.
set_parameters_dr16
assert(width==3,'MEX convolution assumes width=3.');
previousRng=rng; cleanup=onCleanup(@() rng(previousRng)); %#ok<NASGU>
rng(opts.Seed,'twister');

% Calculate the exact drawer SNR and processor search windows on eligible rows.
snr=nan(nParent,1); zlo=nan(nParent,1); zhi=nan(nParent,1);
for idx=find(opts.Eligible(:)&isfinite(zQSO)).'
    w=preloaded.all_wavelengths{idx}(:); f=preloaded.all_flux{idx}(:);
    v=preloaded.all_noise_variance{idx}(:); mask=preloaded.all_pixel_mask{idx}(:);
    dispix=preloaded.all_sigma_pixel{idx}(:);
    if numel(w)<3 || any([numel(f),numel(v),numel(mask),numel(dispix)]~=numel(w)) || ...
            any(~isfinite(w)) || any(diff(w)<=0), continue; end
    if ~isreal(v) || any(v<0), continue; end
    snr(idx)=median(f./sqrt(v));
    model=w/(1+zQSO(idx))>=min_lambda & w/(1+zQSO(idx))<=max_lambda & dispix>0;
    usable=model & ~logical(mask);
    if nnz(model)<3 || nnz(usable)<3 || any(~isfinite(dispix(model))), continue; end
    zlo(idx)=min_z_c4(w(usable),zQSO(idx));
    zhi(idx)=max_z_c4(zQSO(idx),max_z_cut);
end
zb=discretize(zQSO,opts.ZEdges); sb=discretize(snr,opts.SNREdges);
valid=opts.Eligible(:)&isfinite(zlo)&isfinite(zhi)&zhi>zlo&~isnan(zb)&~isnan(sb);
candidates=find(valid);
assert(numel(candidates)>=nInject,'Only %d eligible spectra; %d requested.',numel(candidates),nInject);
% Randomize within each stratum, then draw round-robin across nonempty strata.
stratum=sub2ind([numel(opts.ZEdges)-1,numel(opts.SNREdges)-1],zb(candidates),sb(candidates));
labels=unique(stratum); groups=cell(numel(labels),1); cursor=ones(numel(labels),1);
for k=1:numel(labels)
    x=candidates(stratum==labels(k)); groups{k}=x(randperm(numel(x)));
end
selected=zeros(nInject,1); count=0;
while count<nInject
    order=randperm(numel(labels));
    for k=order
        if cursor(k)<=numel(groups{k})
            count=count+1; selected(count)=groups{k}(cursor(k)); cursor(k)=cursor(k)+1;
            if count==nInject, break; end
        end
    end
end
for k=1:numel(fields), injected.(fields{k})=preloaded.(fields{k})(selected); end
injected.parent_index=selected; injected.z_qso=zQSO(selected);
zAbs=zlo(selected)+(zhi(selected)-zlo(selected)).*rand(nInject,1);
logN=fit_min_log_nciv+(fit_max_log_nciv-fit_min_log_nciv).*rand(nInject,1);
sigma=min_sigma+(max_sigma-min_sigma).*rand(nInject,1);
rew1548=nan(nInject,1); rew1550=nan(nInject,1); coverage1548=zeros(nInject,1); coverage1550=zeros(nInject,1);
for q=1:nInject
    w=double(injected.all_wavelengths{q}(:)); f=injected.all_flux{q}; fshape=size(f); f=double(f(:));
    v=double(injected.all_noise_variance{q}(:)); mask=logical(injected.all_pixel_mask{q}(:));
    dp=double(injected.all_sigma_pixel{q}(:));
    model=w/(1+zQSO(selected(q)))>=min_lambda & w/(1+zQSO(selected(q)))<=max_lambda & dp>0;
    wm=w(model); dm=dp(model); factor=nAVG+1; pad=width*factor;
    wf=finer(wm,nAVG).'; df=finer(dm,nAVG).';
    wp=[logspace(log10(min(wm))-width*pixel_spacing,log10(min(wm))-pixel_spacing/factor,pad).'; ...
        wf; logspace(log10(max(wm))+pixel_spacing/factor,log10(max(wm))+width*pixel_spacing,pad).'];
    dsp=[repmat(dm(1),pad,1); df; repmat(dm(end),pad,1)];
    native=1:factor:numel(wf);
    profile=voigt_iP_transition(wp,zAbs(q),10^logN(q),3,sigma(q),dsp,factor);
    a1=voigt_iP_transition(wp,zAbs(q),10^logN(q),1,sigma(q),dsp,factor);
    a2=voigt_iP_transition(wp,zAbs(q),10^logN(q),2,sigma(q),dsp,factor);
    assert(numel(profile)==numel(wf)&&all(isfinite(profile))&&all(isfinite(a1))&&all(isfinite(a2)), ...
        'Nonfinite or misaligned Voigt profile for parent row %d.',selected(q));
    rew1548(q)=trapz(wf,1-a1)/(1+zAbs(q)); rew1550(q)=trapz(wf,1-a2)/(1+zAbs(q));
    good=~mask & isfinite(f) & isfinite(v) & v>0;
    continuumInput=f; continuumInput(~good)=NaN;
    continuum=movmedian(continuumInput,opts.ContinuumPixels,'omitnan');
    modelIndices=find(model); nativeProfile=profile(native);
    affectedLocal=~mask(model) & abs(1-nativeProfile)>1e-10;
    affected=modelIndices(affectedLocal);
    assert(all(isfinite(continuum(affected)) & continuum(affected)>0), ...
        'Invalid local continuum estimate in parent row %d.',selected(q));
    f(affected)=f(affected)-continuum(affected).*(1-nativeProfile(affectedLocal));
    injected.all_flux{q}=reshape(f,fshape);
    dv=300/299792.458;
    coverage1548(q)=nnz(good & model & abs(w/civ_1548_wavelength-1-zAbs(q))<=dv*(1+zAbs(q)));
    coverage1550(q)=nnz(good & model & abs(w/civ_1550_wavelength-1-zAbs(q))<=dv*(1+zAbs(q)));
end
truth=table(selected,zQSO(selected),snr(selected),zb(selected),sb(selected),zlo(selected),zhi(selected), ...
    zAbs,logN,sigma,rew1548,rew1550,coverage1548,coverage1550, ...
    'VariableNames',{'parent_index','z_qso','snr','z_stratum','snr_stratum','search_z_min','search_z_max', ...
    'z_abs','logN','sigma_cm_s','rew1548_A','rew1550_A','coverage1548_pixels','coverage1550_pixels'});
info=struct('options',opts,'sampling','uniform logN, sigma and search-window redshift', ...
    'eligible_spectra',nnz(valid),'occupied_strata',numel(labels),'nAVG',nAVG, ...
    'noise_model','observed noise retained; additive local-continuum decrement', ...
    'profile_model','voigt_iP_transition with native-pixel WDISP');
end
