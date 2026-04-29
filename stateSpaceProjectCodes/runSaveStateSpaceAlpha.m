% runSaveStateSpaceData
% Master runner for state-space analysis. Implements a TWO-PASS architecture:
%
% PASS 1: Extract electrode-time matrices for all subjects, run parallel
%         analysis to find individual embedding dimensions (k).
%
% GLOBAL K: Compute K_std = median(k) of Controls during EO1 baseline.
%           This ensures all subjects' metrics are in the same dimensional
%           space for fair comparison.
%
% PASS 2: Project all subjects to K_std dimensions, compute dynamical
%         metrics (velocity, volume, radius) and TNN intrinsic dimension.

%% ========================= Configuration ================================
folderSourceString = 'C:\AIML\patagonia';
gridType = 'EEG';

% Analysis parameters
badEyeCondition = 'ep';
badTrialVersion = 'v8';
stRange = [0.25 1.25];         % 1-second stimulus window
gammaRange = [8 13];            % Bandpass filter for alpha
targetFs = 50;                  % Downsample target (Nyquist-safe for 13 Hz)
nPermutations = 100;            % Parallel analysis permutations (use 1000 for final run)

% Protocols: EO1 (baseline) vs M1 (meditation) — both eyes-open
protocolsToAnalyze = {'EO1','EC1', 'M1'};
numProtocols = length(protocolsToAnalyze);

saveFolderName = fullfile(fileparts(mfilename('fullpath')), 'savedDataAlpha_AllStates');
if ~exist(saveFolderName, 'dir'); mkdir(saveFolderName); end

%% ========================= Load Subject Info ============================
% Get paired subject list (30 meditator-control pairs)
pairedSubjectNameList = getPairedSubjectsBK1; % 30 x 2 cell array
numPairs = size(pairedSubjectNameList, 1);

% Get experiment dates for all subjects
[allSubjectNames, expDateList] = getDemographicDetails('BK1');

fprintf('=== State Space Analysis: %d pairs, %d protocols ===\n', numPairs, numProtocols);

%% ========================= PASS 1: Parallel Analysis ====================
% Determine individual k values for all subjects and find global K_std

fprintf('\n--- PASS 1: Computing embedding dimensions ---\n');

% Storage for Pass 1 results
% allResults(pairIdx, groupIdx, protocolIdx) — groupIdx: 1=meditator, 2=control
allResults = struct();

for pairIdx = 1:numPairs
    for grpIdx = 1:2
        subjectName = pairedSubjectNameList{pairIdx, grpIdx};
        expDate = expDateList{strcmp(subjectName, allSubjectNames)};

        for pIdx = 1:numProtocols
            protocolName = protocolsToAnalyze{pIdx};
            fprintf('  Pass 1: %s (%s) [pair %d/%d, grp %d]\n', ...
                subjectName, protocolName, pairIdx, numPairs, grpIdx);

            % Build gamma-filtered, downsampled, trial-concatenated matrix
            [dataMatrix, goodElecList, ~, numGoodTrials] = ...
                getElectrodeTimeMatrix(folderSourceString, subjectName, ...
                gridType, expDate, protocolName, badEyeCondition, ...
                badTrialVersion, stRange, gammaRange, targetFs);

            if isempty(dataMatrix)
                warning('Skipping %s %s: insufficient data', subjectName, protocolName);
                allResults(pairIdx, grpIdx, pIdx).k = NaN;
                allResults(pairIdx, grpIdx, pIdx).valid = false;
                continue;
            end

            % Run parallel analysis
            [k, eigenvalues, threshold95, pcCoeffs, pcScores] = ...
                computeParallelAnalysis(dataMatrix, nPermutations);

            % Store results for Pass 2
            allResults(pairIdx, grpIdx, pIdx).subjectName = subjectName;
            allResults(pairIdx, grpIdx, pIdx).expDate = expDate;
            allResults(pairIdx, grpIdx, pIdx).protocolName = protocolName;
            allResults(pairIdx, grpIdx, pIdx).k = k;
            allResults(pairIdx, grpIdx, pIdx).eigenvalues = eigenvalues;
            allResults(pairIdx, grpIdx, pIdx).threshold95 = threshold95;
            allResults(pairIdx, grpIdx, pIdx).pcCoeffs = pcCoeffs;
            allResults(pairIdx, grpIdx, pIdx).dataMatrix = dataMatrix;
            allResults(pairIdx, grpIdx, pIdx).goodElecList = goodElecList;
            allResults(pairIdx, grpIdx, pIdx).numGoodTrials = numGoodTrials;
            allResults(pairIdx, grpIdx, pIdx).valid = true;

            fprintf('    k = %d, %d elecs, %d trials\n', k, length(goodElecList), numGoodTrials);
        end
    end
end

%% ========================= Compute Global K_std =========================
% K_std = median k of Controls during EO1 (baseline)
% Controls are grpIdx = 2, EO1 is pIdx = 1

controlEO1_k = zeros(numPairs, 1);
for pairIdx = 1:numPairs
    if allResults(pairIdx, 2, 1).valid
        controlEO1_k(pairIdx) = allResults(pairIdx, 2, 1).k;
    else
        controlEO1_k(pairIdx) = NaN;
    end
end

K_std = round(nanmedian(controlEO1_k));
fprintf('\n=== Global K_std = %d (median of Control-EO1 k values) ===\n', K_std);
fprintf('  Control EO1 k range: [%d, %d], mean: %.1f\n', ...
    min(controlEO1_k), max(controlEO1_k), nanmean(controlEO1_k));

%% ========================= PASS 2: Metrics & TNN =======================
% Project all subjects to K_std dimensions, compute metrics and TNN

fprintf('\n--- PASS 2: Computing metrics in K_std=%d dimensional space ---\n', K_std);

for pairIdx = 1:numPairs
    for grpIdx = 1:2
        for pIdx = 1:numProtocols
            if ~allResults(pairIdx, grpIdx, pIdx).valid
                continue;
            end

            subjectName = allResults(pairIdx, grpIdx, pIdx).subjectName;
            protocolName = allResults(pairIdx, grpIdx, pIdx).protocolName;
            fprintf('  Pass 2: %s (%s)\n', subjectName, protocolName);

            dataMatrix = allResults(pairIdx, grpIdx, pIdx).dataMatrix;
            pcCoeffs = allResults(pairIdx, grpIdx, pIdx).pcCoeffs;

            % Project to exactly K_std dimensions
            centeredData = dataMatrix - mean(dataMatrix, 2);
            projected = centeredData' * pcCoeffs(:, 1:K_std); % T_total x K_std

            % Compute dynamical metrics in standardized space
            % Pass full eigenvalues (for PR) and numGoodTrials (for
            % within-trial velocity — avoids concatenation-boundary spikes)
            eigenvalues = allResults(pairIdx, grpIdx, pIdx).eigenvalues;
            numGoodTrials_r = allResults(pairIdx, grpIdx, pIdx).numGoodTrials;
            metrics = computeStateSpaceMetrics(projected, K_std, eigenvalues, numGoodTrials_r);

            % Compute TNN on FULL high-dimensional data (not k-truncated)
            [tnnD, muValues] = computeTNNDimension(dataMatrix');

            % Store in results
            allResults(pairIdx, grpIdx, pIdx).projected = projected;
            allResults(pairIdx, grpIdx, pIdx).metrics = metrics;
            allResults(pairIdx, grpIdx, pIdx).tnnD = tnnD;
            allResults(pairIdx, grpIdx, pIdx).muValues = muValues;

            fprintf('    Velocity=%.4f, logVol=%.2f, Radius=%.4f, PR=%.2f, TNN_d=%.2f\n', ...
                metrics.velocity, metrics.logVolume, metrics.radius, metrics.PR, tnnD);
        end
    end
end

%% ========================= Save Results =================================

fprintf('\n--- Saving results ---\n');

% Save per-subject files
for pairIdx = 1:numPairs
    for grpIdx = 1:2
        for pIdx = 1:numProtocols
            if ~allResults(pairIdx, grpIdx, pIdx).valid
                continue;
            end

            subjectName = allResults(pairIdx, grpIdx, pIdx).subjectName;
            protocolName = allResults(pairIdx, grpIdx, pIdx).protocolName;

            % Extract results for saving (exclude large dataMatrix to save space)
            resultData = rmfield(allResults(pairIdx, grpIdx, pIdx), 'dataMatrix');
            resultData.K_std = K_std;

            fileName = fullfile(saveFolderName, ...
                [subjectName '_' protocolName '_stateSpace.mat']);
            save(fileName, '-struct', 'resultData');
        end
    end
end

% Save summary file with all results aggregated
summaryFile = fullfile(saveFolderName, 'stateSpaceSummary.mat');
save(summaryFile, 'K_std', 'controlEO1_k', 'pairedSubjectNameList', ...
    'protocolsToAnalyze', 'gammaRange', 'targetFs', 'stRange', ...
    'nPermutations', 'badEyeCondition', 'badTrialVersion');

% Build summary table
fprintf('\n=== Summary Table ===\n');
fprintf('%-10s %-6s %5s %8s %10s %8s %6s %8s\n', ...
    'Subject', 'Proto', 'k', 'Velocity', 'logVolume', 'Radius', 'PR', 'TNN_d');
fprintf('%s\n', repmat('-', 1, 70));

for pairIdx = 1:numPairs
    for grpIdx = 1:2
        for pIdx = 1:numProtocols
            if ~allResults(pairIdx, grpIdx, pIdx).valid; continue; end
            r = allResults(pairIdx, grpIdx, pIdx);
            grpStr = {'Med', 'Con'};
            fprintf('%-10s %-6s %5d %8.4f %10.2f %8.4f %6.2f %8.2f  [%s]\n', ...
                r.subjectName, r.protocolName, r.k, ...
                r.metrics.velocity, r.metrics.logVolume, r.metrics.radius, ...
                r.metrics.PR, r.tnnD, grpStr{grpIdx});
        end
    end
end

fprintf('\nDone. Results saved to: %s\n', saveFolderName);
