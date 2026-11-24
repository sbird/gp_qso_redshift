function inspect_civ_injection_recovery(runDir,parentIds)
%INSPECT_CIV_INJECTION_RECOVERY Plot saved injection trials without reprocessing.
% parentIds is a vector of parent-catalog row numbers; empty selects the
% three pilot REW outliers (if present), two weak misses and two recoveries.
% Original and injected flux are shown alongside a reconstructed GP posterior
% L2 model at the matched search (or first search if no match). The latter
% requires the original trained GP model and processor settings. This is a
% reconstruction, not a model vector saved by the processor.
if nargin<2, parentIds=[]; end
required={'recovery_results.mat','injection_truth.mat','injected_spectra.mat','processed_injections.mat'};
for k=1:numel(required), assert(isfile(fullfile(runDir,required{k})),'Missing %s in %s.',required{k},runDir); end
originalFile=fullfile(pwd,'data','dr16','processed','preloaded_qsos_dr16.mat');
modelFile=fullfile(pwd,'data','dr7','processed','learned_model-C13-full.mat');
assert(isfile(originalFile),'Original preloaded spectra not found: %s.',originalFile);
assert(isfile(modelFile),'Trained GP model not found: %s.',modelFile);
assert(exist('voigt_iP_transition','file')==3,'Compiled Voigt MEX is required.');
A=load(fullfile(runDir,'recovery_results.mat'),'results'); T=A.results;
B=load(fullfile(runDir,'injection_truth.mat'),'injection_info');
S=load(fullfile(runDir,'injected_spectra.mat'),'injected_spectra'); S=S.injected_spectra;
R=load(fullfile(runDir,'processed_injections.mat'),'p_c4','p_c4L1','p_no_c4', ...
    'map_z_c4L2','map_N_c4L2','map_sigma_c4L2','map_z_c4L1','test_ind');
G=load(modelFile,'rest_wavelengths','mu','M');
assert(istable(T)&&all(ismember({'parent_index','z_qso','z_abs','rew1548_A','snr', ...
    'matched_search','matched_p_c4','matched_REW_1548_median','recovered_p95'},T.Properties.VariableNames)), ...
    'Recovery table missing required fields.');
assert(isequal(find(R.test_ind(:)),T.parent_index(:)), 'Processor and recovery row ordering differ.');
assert(numel(unique(S.parent_index))==numel(S.parent_index)&&height(T)==numel(S.parent_index), ...
    'Injected parent rows are not unique or complete.');
if isempty(parentIds)
    parentIds=intersect([633763 721937 650760],T.parent_index,'stable');
    missed=find(~logical(T.recovered_p95));
    [~,o]=sort(T.rew1548_A(missed)); missed=missed(o(1:min(2,numel(o))));
    recovered=find(logical(T.recovered_p95));
    [~,o]=sort(abs(T.matched_REW_1548_median(recovered)-T.rew1548_A(recovered)),'ascend');
    recovered=recovered(o(1:min(2,numel(o))));
    parentIds=unique([parentIds(:);T.parent_index(missed);T.parent_index(recovered)],'stable');
end
assert(isnumeric(parentIds)&&isvector(parentIds)&&~isempty(parentIds),'Provide at least one parent ID.');
[found,rows]=ismember(parentIds(:),T.parent_index); assert(all(found),'A requested parent row is not in this run.');
% MATLAB script-generated names must be declared in this function workspace.
width=[]; nAVG=[]; min_lambda=[]; max_lambda=[]; pixel_spacing=[];
civ_1548_wavelength=[]; civ_1550_wavelength=[]; set_parameters_dr16
assert(width==3,'Voigt MEX assumes padding width 3.');
assert(isfield(B,'injection_info')&&isfield(B.injection_info,'nAVG')&& ...
    B.injection_info.nAVG==nAVG,'Current nAVG differs from the injection run.');
% The saved results do not include dv_mask; 250 km/s was used by the launcher.
dvMask=250; fprintf('Reconstructing previous-search masks with dv_mask=%g km/s.\n',dvMask);
muI=griddedInterpolant(G.rest_wavelengths,G.mu,'linear');
MI=griddedInterpolant({G.rest_wavelengths,1:size(G.M,2)},G.M,'linear');
P=matfile(originalFile); sz=size(P,'all_flux');
assert(numel(sz)==2 && (sz(1)==1 || sz(2)==1),'Expected vector all_flux in original preloaded file.');
outputDir=fullfile(runDir,['inspection_' datestr(now,'yyyymmdd_HHMMSS')]);
assert(~isfolder(outputDir),'Inspection output directory already exists.'); mkdir(outputDir);
for rr=1:numel(rows)
    i=rows(rr); parent=T.parent_index(i); s=find(S.parent_index==parent,1);
    assert(~isempty(s),'Injected spectrum not found for parent %d.',parent);
    if sz(1)==1, originalCell=P.all_flux(1,parent); else, originalCell=P.all_flux(parent,1); end
    original=double(originalCell{1}(:)); w=double(S.all_wavelengths{s}(:));
    flux=double(S.all_flux{s}(:)); variance=double(S.all_noise_variance{s}(:));
    mask=logical(S.all_pixel_mask{s}(:)); dispix=double(S.all_sigma_pixel{s}(:));
    assert(numel(original)==numel(w)&&numel(flux)==numel(w),'Original and injection spectrum lengths differ.');
    matched=isfinite(T.matched_search(i));
    if matched, slot=T.matched_search(i); else, slot=1; end
    assert(slot>=1&&slot<=size(R.p_c4,2)&&slot==fix(slot),'Invalid matched search.');
    model=nan(size(w)); modelNote='reconstructed posterior L2 model';
    try
        % Reproduce the processor's fine grid, prior-search masks and GP fit.
        modelRegion=w/(1+T.z_qso(i))>=min_lambda & w/(1+T.z_qso(i))<=max_lambda & dispix>0;
        wm=w(modelRegion); dm=dispix(modelRegion); fineFactor=nAVG+1; pad=width*fineFactor;
        wf=finer(wm,nAVG).'; df=finer(dm,nAVG).';
        wp=[logspace(log10(min(wm))-width*pixel_spacing,log10(min(wm))-pixel_spacing/fineFactor,pad).'; ...
            wf;logspace(log10(max(wm))+pixel_spacing/fineFactor,log10(max(wm))+width*pixel_spacing,pad).'];
        dp=[repmat(dm(1),pad,1);df;repmat(dm(end),pad,1)];
        profile=voigt_iP_transition(wp,R.map_z_c4L2(i,slot),10^R.map_N_c4L2(i,slot),3, ...
            R.map_sigma_c4L2(i,slot),dp,fineFactor);
        native=profile(1:fineFactor:numel(wf)); usable=~mask(modelRegion);
        fullRows=find(modelRegion); fitRows=fullRows(usable); fitW=w(fitRows);
        A2=native(usable); rest=fitW/(1+T.z_qso(i));
        mu=muI(rest).*A2; M=MI({rest,1:size(G.M,2)}).*A2;
        surviving=true(numel(fitRows),1); z1548=fitW/civ_1548_wavelength-1; z1550=fitW/civ_1550_wavelength-1;
        for prev=1:slot-1
            if R.p_c4(i,prev)>0.75
                z=R.map_z_c4L2(i,prev); dz=dvMask/299792.458*(1+z);
                surviving=surviving & abs(z1548-z)>dz & abs(z1550-z)>dz;
            elseif R.p_c4L1(i,prev)>0.75
                z=R.map_z_c4L1(i,prev); dz=dvMask/299792.458*(1+z);
                surviving=surviving & abs(z1548-z)>dz;
            else
                error('Cannot reproduce previous-search masking at search %d.',prev);
            end
        end
        fit=mean_posterior_continuum(flux(fitRows(surviving)),mu(surviving),M(surviving,:),variance(fitRows(surviving)));
        model(fitRows(surviving))=fit;
    catch ME
        modelNote=['model reconstruction failed: ' ME.message];
        warning('Parent %d: %s',parent,modelNote);
    end
    fig=figure('Visible','off','Color','w','Position',[100 100 1250 720]);
    tiledlayout(fig,2,1,'TileSpacing','compact');
    centers=[civ_1548_wavelength civ_1550_wavelength]*(1+T.z_abs(i));
    lineNames=[1548 1550];
    for line=1:2
        ax=nexttile; vel=299792.458*(w/centers(line)-1);
        near=abs(vel)<=2500;
        stairs(ax,vel(near),original(near),'Color',[.65 .65 .65],'DisplayName','Original flux'); hold(ax,'on');
        stairs(ax,vel(near),flux(near),'Color',[0 .3 .7],'DisplayName','Injected flux');
        plot(ax,vel(near),model(near),'Color',[.85 .25 .05],'LineWidth',1.4,'DisplayName',modelNote);
        xline(ax,0,'--k','Injected center','HandleVisibility','off');
        xlabel(ax,sprintf('Velocity relative to injected %.0f (km/s)',lineNames(line)));
        ylabel(ax,'Normalized flux'); grid(ax,'on'); xlim(ax,[-2500 2500]);
        if line==1, legend(ax,'Interpreter','none','Location','best'); end
    end
    if matched
        status=sprintf('matched search %d, p=%.4f, recovered median W=%.3f A',slot,T.matched_p_c4(i),T.matched_REW_1548_median(i));
    else
        status=sprintf('no match; showing first-search L2 model (p=%.4f)',R.p_c4(i,1));
    end
    sgtitle(fig,sprintf('Parent %d | z_{inj}=%.5f | SNR=%.1f | injected W=%.3f A | %s', ...
        parent,T.z_abs(i),T.snr(i),T.rew1548_A(i),status),'Interpreter','tex');
    exportgraphics(fig,fullfile(outputDir,sprintf('parent_%d.png',parent)),'Resolution',180); close(fig);
    fprintf('Saved parent %d: %s\n',parent,status);
end
writetable(T(rows,{'parent_index','z_qso','snr','z_abs','rew1548_A','matched_p_c4', ...
    'matched_search','matched_REW_1548_median','recovered_p95'}),fullfile(outputDir,'inspected_cases.csv'));
fprintf('Inspection plots saved in %s\n',outputDir);
end
