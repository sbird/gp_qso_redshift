%% L2 effective sample size in processed C IV chunks
% From the gpCIVIC root, run: check_chunk_ess
% Change pCut to select a different probability threshold.
chunkDir = 'processedFiles';
files = dir(fullfile(chunkDir, ...
    'processed_qsos_tst_posterior_corrected_S_*_E_*.mat'));
assert(~isempty(files), 'No processed chunks found in %s.', chunkDir);
ranges = zeros(numel(files), 2);
for k = 1:numel(files)
    tok = regexp(files(k).name, '_S_(\d+)_E_(\d+)\.mat$', 'tokens', 'once');
    assert(~isempty(tok), 'Invalid chunk name: %s', files(k).name);
    ranges(k,:) = [str2double(tok{1}), str2double(tok{2})];
end
[ranges, order] = sortrows(ranges, 1); files = files(order);
assert(all(ranges(:,1) >= 1 & ranges(:,2) >= ranges(:,1)) && ...
    all(ranges(2:end,1) > ranges(1:end-1,2)), ...
    'Invalid or overlapping chunk ranges.');
if ranges(1,1) ~= 1 || any(ranges(2:end,1) > ranges(1:end-1,2) + 1)
    warning('Chunk ranges have gaps; results cover only supplied chunks.');
end

pCut = 0.90;
essLimits = [2 5 10 20 50];
blockSize = 100;  % QSOs loaded at a time

nSearch = 0;
nInvalid = 0;
searchCounts = zeros(size(essLimits));
qsoCounts = zeros(size(essLimits));

% for k = 1:numel(files)
for k = 1:1
    file = fullfile(files(k).folder, files(k).name);
    m = matfile(file);
    p = m.p_c4;
    [nQSO, nSlots] = size(p);
    qsoLowESS = false(nQSO, numel(essLimits));

    for slot = 1:nSlots
        for first = 1:blockSize:nQSO
            rows = first:min(first + blockSize - 1, nQSO);
            selected = isfinite(p(rows, slot)) & p(rows, slot) > pCut;
            if ~any(selected), continue; end

            % Each row contains one search's proposal log likelihoods.
            L = m.sample_log_likelihoods_c4L2(rows, :, slot);

            % -Inf is allowed: it represents a zero-weight proposal.
            valid = ~any(isnan(L) | L == Inf, 2) & any(isfinite(L), 2);
            nInvalid = nInvalid + nnz(selected & ~valid);
            use = selected & valid;
            if ~any(use), continue; end

            L = L(use, :);
            weights = exp(L - max(L, [], 2));
            ess = sum(weights, 2).^2 ./ sum(weights.^2, 2);
            selectedRows = rows(use);
            nSearch = nSearch + numel(ess);

            for j = 1:numel(essLimits)
                low = ess < essLimits(j);
                searchCounts(j) = searchCounts(j) + nnz(low);
                qsoLowESS(selectedRows(low), j) = true;
            end
        end
    end

    qsoCounts = qsoCounts + sum(qsoLowESS, 1);
    fprintf('Checked %s\n', files(k).name);
end

fprintf('\nL2 searches with p_CIV > %.2f and valid weights: %d\n', pCut, nSearch);
fprintf('Selected searches excluded for invalid weights: %d\n', nInvalid);
if nSearch > 0
    for j = 1:numel(essLimits)
        fprintf('ESS < %-3g: %6d searches (%5.1f%%); %6d QSOs\n', ...
            essLimits(j), searchCounts(j), ...
            100 * searchCounts(j) / nSearch, qsoCounts(j));
    end
end
