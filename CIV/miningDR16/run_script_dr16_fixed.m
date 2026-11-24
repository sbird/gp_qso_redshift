
% Batch arguments are supplied by submit_injection_recovery.sh as environment
% variables. Without them, the interactive settings below remain in control.
injection_batch_start_text = getenv('CIV_INJECTION_START');
injection_batch_end_text = getenv('CIV_INJECTION_END');
injection_batch_mode = ~isempty(injection_batch_start_text) || ~isempty(injection_batch_end_text);
if injection_batch_mode
    assert(~isempty(injection_batch_start_text) && ~isempty(injection_batch_end_text), ...
        'Both CIV_INJECTION_START and CIV_INJECTION_END are required.');
    injection_range_start = str2double(injection_batch_start_text);
    injection_range_end = str2double(injection_batch_end_text);
    assert(isfinite(injection_range_start) && isfinite(injection_range_end) && ...
        injection_range_start >= 1 && injection_range_end >= injection_range_start && ...
        injection_range_start == fix(injection_range_start) && injection_range_end == fix(injection_range_end), ...
        'Injection range must consist of positive integer test-set indices.');
end
% Set to 1 to plot an existing injection run without rerunning the pipeline.
recovery_inspector = 0;
if injection_batch_mode, recovery_inspector = 0; end
recovery_inspection_dir = fullfile(pwd,'injection_recovery','pilot_seed_42');
% Empty selects three pilot outliers plus representative misses and recoveries.
recovery_inspection_parent_ids = [713100 256723 597816 410674 577519 ...
                                  221575 330121 79637 342685 584946];
if recovery_inspector == 1
    inspect_civ_injection_recovery(recovery_inspection_dir, ...
        recovery_inspection_parent_ids);
    return
elseif recovery_inspector ~= 0
    error('recovery_inspector must be 0 or 1.');
end

fprintf('Setting parameters ...\n')
% Injection pilot: use original detections to avoid known features.
% Each run gets a separate directory; change the seed or directory to rerun.
%injector = 0;
injection_count = 5000;
injection_seed = 42;
injection_chunk_dir = fullfile(pwd,'data','dr16','processed','latestProcessedMATfiles');
injection_detection_threshold = 0.75; % original L1 or L2 posterior
injection_exclusion_kms = 1000;         % line-to-line separation
injection_output_dir = fullfile(pwd,'injection_recovery',sprintf('pilot_seed_%d',injection_seed));
injection_match_kms = 300;
injection_min_zabs = -Inf;
injection_min_rew1548 = -Inf;
snr_threshold = 2;
injection_snr_threshold = snr_threshold;
injection_resume = false;
% Processing is chunked so the script can be launched interactively or by a
% batch scheduler.  Adjust these three values for a batch job; `num_quasars`
% is clamped to the number of eligible preloaded spectra below.
% ind_S = 1;                   % plotting run: filtered test-set range
% num_quasars = 5000;           % covers all nine plotted_qso_ind entries
% cores = 6;
hpcc = 1;
% cores = 32;
cataloging  = 0;
preloading  = 0;
sampling    = 0;
plotting    = 0;
processing  = 1;
merging     = 0;
EWer        = 0;
pltP        = 0;
CredInt     = 0;
dv_mask     = 250;
voigtPrep = 0;
Averaging =1; 
maskType = 1;
priorType = 1;
saving=1;
if injection_batch_mode
    % Batch execution selects original filtered-test rows, not parent rows.
    injector = 1; plotting = 0; hpcc = 1;
    injection_series_text = getenv('CIV_INJECTION_SERIES');
    if isempty(injection_series_text), injection_series = 1;
    else, injection_series = str2double(injection_series_text); end
    assert(isfinite(injection_series) && injection_series >= 1 && injection_series == fix(injection_series), ...
        'CIV_INJECTION_SERIES must be a positive integer.');
    injection_seed = injection_range_start + injection_series - 1;
    assert(injection_seed <= 2^32-1,'Injection seed exceeds MATLAB twister range.');
    injection_count_text = getenv('CIV_INJECTION_COUNT');
    if ~isempty(injection_count_text), injection_count = str2double(injection_count_text); end
    assert(isfinite(injection_count) && injection_count >= 1 && injection_count == fix(injection_count), ...
        'CIV_INJECTION_COUNT must be a positive integer.');
    assert(injection_count <= injection_range_end-injection_range_start+1, ...
        'Requested injections exceed the number of test rows in this range.');
    batch_cores_text = getenv('SLURM_CPUS_PER_TASK');
    if ~isempty(batch_cores_text), cores = str2double(batch_cores_text); end
    local_cluster = parcluster('local');
    assert(isfinite(cores) && cores >= 1 && cores == fix(cores), ...
        'Requested worker count must be a positive integer.');
    if cores > local_cluster.NumWorkers
        warning('run_script_dr16:WorkerCap', ...
            'Requested %d workers; local profile allows %d. Capping pool to %d.', ...
            cores,local_cluster.NumWorkers,local_cluster.NumWorkers);
        cores = local_cluster.NumWorkers;
    end
    injection_output_dir = fullfile(pwd,'injection_recovery', ...
        sprintf('batch_old_rule_S_%d_E_%d_seed_%d', ...
        injection_range_start,injection_range_end,injection_seed));
    injection_min_zabs_text = getenv('CIV_INJECTION_MIN_ZABS');
    if isempty(injection_min_zabs_text), injection_min_zabs = -Inf;
    else, injection_min_zabs = str2double(injection_min_zabs_text); end
    injection_min_rew_text = getenv('CIV_INJECTION_MIN_REW1548');
    if isempty(injection_min_rew_text), injection_min_rew1548 = -Inf;
    else, injection_min_rew1548 = str2double(injection_min_rew_text); end
    injection_snr_text = getenv('CIV_INJECTION_SNR_THRESHOLD');
    if isempty(injection_snr_text), injection_snr_threshold = snr_threshold;
    else, injection_snr_threshold = str2double(injection_snr_text); end
    assert(~isnan(injection_min_zabs) && ~isnan(injection_min_rew1548) && ...
        isfinite(injection_snr_threshold), 'Injection targeting values must be numeric.');
    injection_resume_text = getenv('CIV_INJECTION_RESUME');
    assert(isempty(injection_resume_text) || strcmp(injection_resume_text,'1'), ...
        'CIV_INJECTION_RESUME must be empty or 1.');
    injection_resume = strcmp(injection_resume_text,'1');
    if isfinite(injection_min_zabs) || isfinite(injection_min_rew1548) || injection_snr_threshold~=snr_threshold
        injection_output_dir = [injection_output_dir sprintf('_zabs_gt_%g_rew_gt_%g_snr_gt_%g', ...
            injection_min_zabs,injection_min_rew1548,injection_snr_threshold)];
    end
    if injection_series > 1
        injection_output_dir = [injection_output_dir sprintf('_series_%d',injection_series)];
    end
    % Use a MEX compiled for this cluster; concurrent jobs must not rebuild it.
    voigtPrep = 0;
    addpath(genpath(fullfile(pwd,'minFunc_2012')));
    assert(exist('voigt_iP_transition','file')==3, ...
        'Compile voigt_iP_transition for the HPCC before submitting injection batches.');
    fprintf('Injection batch: filtered test rows %d-%d; %d injections; series %d; seed %d.\n', ...
        injection_range_start,injection_range_end,injection_count,injection_series,injection_seed);
end
set_parameters_dr16
if injector == 1 && ~singlet_model_enabled
    injection_output_dir = [injection_output_dir '_doublet_only'];
end
if doublet_coverage_enabled
    coverage_gate_label = 'ON';
else
    coverage_gate_label = 'OFF';
end
fprintf('Doublet coverage gate: %s; processor sample name: %s.\n', ...
    coverage_gate_label, testing_set_name);
if singlet_model_enabled
    singlet_model_label = 'ON';
else
    singlet_model_label = 'OFF';
end
fprintf('Singlet model competition: %s.\n', singlet_model_label);
if injector == 1
    injection_output_tag = '';
    if ~doublet_coverage_enabled, injection_output_tag = [injection_output_tag '_gateoff']; end
    if ~singlet_model_enabled, injection_output_tag = [injection_output_tag '_doublet_only']; end
    injection_processor_file = fullfile(injection_output_dir, ...
        ['processed_injections' injection_output_tag '.mat']);
end
training_set_name
fprintf('Building catalogs ...\n')
if cataloging == 1
    build_catalog_dr16
end
variables_to_load= {'all_plate_dr16', 'all_mjd_dr16', 'all_fiber_dr16', ...
'all_QSO_ID_dr16', 'all_RA_dr16', 'all_DEC_dr16', 'all_zqso_dr16','all_BAL_PROB'};
load('data/dr16/processed/catalogDR16.mat', variables_to_load{:});


fprintf('preparing voigt.c ...\n')
if voigtPrep == 1 
    cd minFunc_2012
    addpath(genpath(pwd))
    mexAll
    cd ..
    
    if hpcc 
        mex voigt_iP_transition.c -lcerf -I/rhome/rmona003/bin/include/ -L/rhome/rmona003/bin/lib64
    else
        mex voigt_iP_transition.c -lcerf
    end

      
end



catDR7 = load(sprintf('%s/catalog', processed_directory(releasePrior)));
filter_flagsDR7 = load('data/dr7/processed/filter_flags');
prior_ind = (filter_flagsDR7.filter_flags==0); 
all_z_civ_C13 = catDR7.all_z_civ1;
% all_REW_1548_DR7 = catDR7.all_EW1;
% all_REW_1550_DR7 = catDR7.all_EW2;


fprintf('Learning model ...\n')

variables_to_load = { 'max_noise_variance', ...
                   'minFunc_options', 'rest_wavelengths', 'mu', ...
                    'initial_M', 'M',  'log_likelihood', ...
                    };

load('data/dr7/processed/learned_model-C13-full.mat', variables_to_load{:});

fprintf('Generating samples for integrating out parameters in the model...\n')
% variables_to_load = {'offset_z_samples', 'offset_sigma_samples',...
%                      'log_nciv_samples', 'nciv_samples'};
if sampling==1
    generate_c4_samples
else
    load(sprintf('data/dr16/processed/civ_samples_N-%d-%d-Sigma-%d-%d-Num-%d.mat',...
        floor(fit_min_log_nciv*100), floor(fit_max_log_nciv*100), floor(min_sigma/1e5), floor(max_sigma/1e5), ...
        num_C4_samples));
end
    
fprintf(sprintf('%d Samples are generated\n', num_C4_samples));
% load preprocessed QSOs
fprintf('Preloading QSOs ...\n')
if preloading == 1
    preload_qsos_dr16_fixed
end
load('data/dr16/processed/preloaded_qsos_dr16.mat');


fprintf('preparing testing and prior indices ...\n')
test_ind = (filter_flags==0);

if injector == 1
    assert(processing==1 && Averaging==1 && saving==1 && plotting==0 && ...
        merging==0 && EWer==0 && pltP==0 && CredInt==0, ...
        'Injection mode requires fixed processing, saving on, and optional downstream steps off.');
    if injection_resume
        assert(isfolder(injection_output_dir) && ...
            isfile(fullfile(injection_output_dir,'injected_spectra.mat')) && ...
            isfile(fullfile(injection_output_dir,'injection_truth.mat')) && ...
            ~isfile(injection_processor_file), ...
            'Resume requires saved injection inputs and no completed processor output in %s.', ...
            injection_output_dir);
        saved_truth=load(fullfile(injection_output_dir,'injection_truth.mat'));
        saved_spectra=load(fullfile(injection_output_dir,'injected_spectra.mat'),'injected_spectra');
        assert(all(isfield(saved_truth,{'injection_truth','injection_info','baselineSelected', ...
            'injection_seed','injection_chunk_dir'})) && isfield(saved_spectra,'injected_spectra'), ...
            'Saved injection inputs are incomplete in %s.',injection_output_dir);
        assert(saved_truth.injection_seed==injection_seed && ...
            strcmp(saved_truth.injection_chunk_dir,injection_chunk_dir) && ...
            saved_truth.injection_info.options.MinZAbs==injection_min_zabs && ...
            saved_truth.injection_info.options.MinREW1548==injection_min_rew1548 && ...
            saved_truth.injection_info.options.SNRThreshold==injection_snr_threshold, ...
            'Saved injection settings do not match this resume request.');
        if injection_batch_mode
            assert(saved_truth.injection_range_start==injection_range_start && ...
                saved_truth.injection_range_end==injection_range_end && ...
                saved_truth.injection_series==injection_series, ...
                'Saved injection batch range/series does not match this resume request.');
        end
        injection_truth=saved_truth.injection_truth;
        injection_info=saved_truth.injection_info;
        baselineSelected=saved_truth.baselineSelected;
        injected_spectra=saved_spectra.injected_spectra;
        assert(injection_count <= height(injection_truth), ...
            'Requested %d injections, but the saved batch contains only %d.', ...
            injection_count,height(injection_truth));
        if injection_count < height(injection_truth)
            selected_truth=injection_truth(1:injection_count,:);
            [found,spectrum_rows]=ismember(selected_truth.parent_index, ...
                injected_spectra.parent_index);
            assert(all(found) && numel(unique(spectrum_rows))==injection_count, ...
                'Saved spectra do not align uniquely with the requested truth rows.');
            for field_name={'all_wavelengths','all_flux','all_noise_variance', ...
                    'all_pixel_mask','all_sigma_pixel'}
                injected_spectra.(field_name{1})=injected_spectra.(field_name{1})(spectrum_rows);
            end
            injected_spectra.parent_index=injected_spectra.parent_index(spectrum_rows);
            injected_spectra.z_qso=injected_spectra.z_qso(spectrum_rows);
            injection_truth=selected_truth;
            baseline_fields=fieldnames(baselineSelected);
            for field_index=1:numel(baseline_fields)
                value=baselineSelected.(baseline_fields{field_index});
                if size(value,1)==height(saved_truth.injection_truth)
                    baselineSelected.(baseline_fields{field_index})=value(1:injection_count,:);
                end
            end
        end
        fprintf('Resuming %d of %d saved injections from %s; processor output %s.\n', ...
            height(injection_truth),height(saved_truth.injection_truth), ...
            injection_output_dir,injection_processor_file);
    else
        assert(~isfolder(injection_output_dir),'Injection output directory exists: %s',injection_output_dir);
        baseline=load_civ_injection_baseline(injection_chunk_dir,test_ind);
    end
    if ~injection_resume && injection_batch_mode
        assert(injection_range_end <= numel(baseline.parent_index), ...
            'Requested test row %d exceeds the %d test-set rows.', ...
            injection_range_end,numel(baseline.parent_index));
        batch_rows = false(size(baseline.processed));
        batch_rows(injection_range_start:injection_range_end) = true;
        baseline.processed = baseline.processed & batch_rows;
        fprintf('Original processed rows in requested range: %d.\n',nnz(baseline.processed));
        previous_parents = zeros(0,1);
        for previous_series = 1:injection_series-1
            previous_seed = injection_range_start + previous_series - 1;
            previous_dir = fullfile(pwd,'injection_recovery', ...
                sprintf('batch_old_rule_S_%d_E_%d_seed_%d', ...
                injection_range_start,injection_range_end,previous_seed));
            if isfinite(injection_min_zabs) || isfinite(injection_min_rew1548) || injection_snr_threshold~=snr_threshold
                previous_dir = [previous_dir sprintf('_zabs_gt_%g_rew_gt_%g_snr_gt_%g', ...
                    injection_min_zabs,injection_min_rew1548,injection_snr_threshold)];
            end
            if previous_series > 1
                previous_dir = [previous_dir sprintf('_series_%d',previous_series)];
            end
            if ~singlet_model_enabled
                previous_dir = [previous_dir '_doublet_only'];
            end
            previous_file = fullfile(previous_dir,'injection_truth.mat');
            assert(isfile(previous_file), ...
                'Previous series truth file missing: %s. Run earlier series first.',previous_file);
            prior = load(previous_file,'injection_truth','injection_range_start','injection_range_end');
            assert(isfield(prior,'injection_truth') && istable(prior.injection_truth) && ...
                isfield(prior,'injection_range_start') && isfield(prior,'injection_range_end') && ...
                prior.injection_range_start == injection_range_start && prior.injection_range_end == injection_range_end && ...
                all(ismember({'parent_index','test_index'},prior.injection_truth.Properties.VariableNames)), ...
                'Invalid or mismatched prior batch truth: %s.',previous_file);
            rows = prior.injection_truth.test_index(:);
            parents = prior.injection_truth.parent_index(:);
            assert(all(isfinite(rows) & rows == fix(rows) & rows >= injection_range_start & rows <= injection_range_end) && ...
                isequal(baseline.parent_index(rows),parents), ...
                'Prior batch does not match the current test-set mapping: %s.',previous_file);
            previous_parents = [previous_parents; parents]; %#ok<AGROW>
        end
        assert(numel(unique(previous_parents)) == numel(previous_parents), ...
            'Earlier series already contain repeated sightlines in this range.');
        baseline.processed(ismember(baseline.parent_index,previous_parents)) = false;
        fprintf('Excluded %d sightlines used in previous series.\n',numel(previous_parents));
    end
    if ~injection_resume
    preload_for_injection=struct('all_wavelengths',{all_wavelengths},'all_flux',{all_flux}, ...
        'all_noise_variance',{all_noise_variance},'all_pixel_mask',{all_pixel_mask}, ...
        'all_sigma_pixel',{all_sigma_pixel});
    [injected_spectra,injection_truth,injection_info]=inject_civ_doublets( ...
        preload_for_injection, all_zqso_dr16, baseline,injection_count,'Seed',injection_seed, ...
        'DetectionThreshold',injection_detection_threshold,'ExclusionKms',injection_exclusion_kms, ...
        'SNRThreshold',injection_snr_threshold,'MinZAbs',injection_min_zabs, ...
        'MinREW1548',injection_min_rew1548);
    clear preload_for_injection
    % The processor filters parent arrays in ascending parent-row order.
    [parent_order,spectrum_order]=sort(injected_spectra.parent_index);
    injection_truth=injection_truth(spectrum_order,:);
    assert(isequal(parent_order(:),injection_truth.parent_index(:)),'Injection row alignment failed.');
    baseline_rows=injection_truth.test_index;
    baselineSelected=struct('parent_index',baseline.parent_index(baseline_rows), ...
        'p_c4',baseline.p_c4(baseline_rows,:),'p_c4L1',baseline.p_c4L1(baseline_rows,:), ...
        'map_z_c4L2',baseline.map_z_c4L2(baseline_rows,:), ...
        'map_z_c4L1',baseline.map_z_c4L1(baseline_rows,:));
    clear baseline baseline_rows
    ind_S=1; num_quasars=injection_count;
    mkdir(injection_output_dir);
    save(fullfile(injection_output_dir,'injected_spectra.mat'),'injected_spectra','-v7.3');
    save(fullfile(injection_output_dir,'injection_truth.mat'),'injection_truth','injection_info', ...
        'baselineSelected','injection_chunk_dir','injection_seed','injection_match_kms', ...
        'injection_detection_threshold','injection_exclusion_kms','-v7.3');
    if injection_batch_mode
        save(fullfile(injection_output_dir,'injection_truth.mat'), ...
            'injection_range_start','injection_range_end','injection_series','-append');
    end
    rng(injection_seed,'twister');
    end
    % The processor filters parent arrays in ascending parent-row order.
    [parent_order,spectrum_order]=sort(injected_spectra.parent_index);
    assert(isequal(parent_order(:),injection_truth.parent_index(:)),'Injection row alignment failed.');
    for injection_row=1:numel(parent_order)
        all_flux{parent_order(injection_row)}=injected_spectra.all_flux{spectrum_order(injection_row)};
    end
    test_ind=false(size(test_ind)); test_ind(parent_order)=true;
    ind_S=1; num_quasars=height(injection_truth);
    clear injected_spectra
elseif injector ~= 0
    error('injector must be 0 or 1.');
end

% `process_qsos_dr16_par` indexes the filtered test set.  Validate and clamp
% the requested chunk here so its saved filename always records the actual
% [start, end] range and the final chunk cannot run past the catalog.
all_num_quasars = nnz(test_ind);
if ind_S < 1 || ind_S > all_num_quasars
    error('run_script_dr16:InvalidStartIndex', ...
        'ind_S (%d) must be between 1 and %d eligible QSOs.', ...
        ind_S, all_num_quasars);
end
if num_quasars < 1 || floor(num_quasars) ~= num_quasars
    error('run_script_dr16:InvalidChunkSize', ...
        'num_quasars must be a positive integer.');
end
if cores < 1 || floor(cores) ~= cores
    error('run_script_dr16:InvalidCoreCount', ...
        'cores must be a positive integer.');
end
num_quasars = min(num_quasars, all_num_quasars - ind_S + 1);
ind_E = ind_S + num_quasars - 1;


if processing==1
    fprintf('processing QSO: %d to %d\n\n', ind_S, ind_E);
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', cores);
    elseif pool.NumWorkers ~= cores
        warning('run_script_dr16:ExistingPool', ...
            ['Reusing the existing pool with %d workers (requested %d).'], ...
            pool.NumWorkers, cores);
    end
    if Averaging==1 && singlet_model_enabled
        process_qsos_dr16_par_fixed
    elseif Averaging==1
        process_qsos_dr16_par_doublet_only
    else
        assert(singlet_model_enabled, ...
            'The doublet-only processor currently requires Averaging=1.');
        process_qsos_dr16_par_nofine
    end
    if injector==1
        report_civ_injection_recovery(injection_truth, baselineSelected, injection_processor_file, ...
            injection_output_dir, injection_match_kms);
    end
    
end

if merging==1
    mergeProcessedDR16

end

if EWer==1
    EW_dr16_voigt
end

if pltP  ==1
   pltPost;
end

if CredInt==1
    CI_sort
end
