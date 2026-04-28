function displayStateSpaceAllSubjects(pairedSubjectNameList, protocolsToAnalyze, ...
    savedDataFolder, displaySettings)
% displayStateSpaceAllSubjects - Load saved state-space results and generate
% all figures and statistical tables.
%
% INPUTS:
%   pairedSubjectNameList - 30x2 cell array (col 1=meditator, col 2=control)
%   protocolsToAnalyze    - cell array (e.g., {'EO1', 'M1'})
%   savedDataFolder       - path to savedData folder
%   displaySettings       - struct with display options

if ~exist('displaySettings', 'var'); displaySettings = struct(); end
if ~isfield(displaySettings, 'parametricTest');  displaySettings.parametricTest = 0; end
if ~isfield(displaySettings, 'pairedDataFlag');  displaySettings.pairedDataFlag = 1; end

% Default colors for Meditator (magenta) and Control (blue) — matching existing code
colorMed = [1 0 1];   % magenta
colorCon = [0 0 1];   % blue
violinColors = [{colorMed} {colorCon}];

numPairs = size(pairedSubjectNameList, 1);
numProtocols = length(protocolsToAnalyze);

%% ========================= Load All Results =============================
% Arrays: (numPairs x numProtocols) for each group
kMed   = NaN(numPairs, numProtocols);
kCon   = NaN(numPairs, numProtocols);
velMed = NaN(numPairs, numProtocols);
velCon = NaN(numPairs, numProtocols);
volMed = NaN(numPairs, numProtocols);
volCon = NaN(numPairs, numProtocols);
radMed = NaN(numPairs, numProtocols);
radCon = NaN(numPairs, numProtocols);
prMed  = NaN(numPairs, numProtocols);   % Participation Ratio
prCon  = NaN(numPairs, numProtocols);
varMed = NaN(numPairs, numProtocols);   % Total Variance (power proxy for ANCOVA)
varCon = NaN(numPairs, numProtocols);
tnnMed = NaN(numPairs, numProtocols);
tnnCon = NaN(numPairs, numProtocols);

% Store projected data for trajectory plots
projectedStore = cell(numPairs, 2, numProtocols); % (pair, grp, protocol)

K_std = NaN;

for pairIdx = 1:numPairs
    for grpIdx = 1:2
        subjectName = pairedSubjectNameList{pairIdx, grpIdx};
        for pIdx = 1:numProtocols
            protocolName = protocolsToAnalyze{pIdx};
            fileName = fullfile(savedDataFolder, ...
                [subjectName '_' protocolName '_stateSpace.mat']);

            if ~exist(fileName, 'file')
                warning('Missing: %s', fileName);
                continue;
            end

            r = load(fileName);
            if isnan(K_std) && isfield(r, 'K_std')
                K_std = r.K_std;
            end

            if grpIdx == 1 % Meditator
                kMed(pairIdx, pIdx)   = r.k;
                velMed(pairIdx, pIdx) = r.metrics.velocity;
                volMed(pairIdx, pIdx) = r.metrics.logVolume;
                radMed(pairIdx, pIdx) = r.metrics.radius;
                tnnMed(pairIdx, pIdx) = r.tnnD;
                if isfield(r.metrics, 'PR')
                    prMed(pairIdx, pIdx) = r.metrics.PR;
                end
                if isfield(r.metrics, 'totalVariance')
                    varMed(pairIdx, pIdx) = r.metrics.totalVariance;
                end
            else % Control
                kCon(pairIdx, pIdx)   = r.k;
                velCon(pairIdx, pIdx) = r.metrics.velocity;
                volCon(pairIdx, pIdx) = r.metrics.logVolume;
                radCon(pairIdx, pIdx) = r.metrics.radius;
                tnnCon(pairIdx, pIdx) = r.tnnD;
                if isfield(r.metrics, 'PR')
                    prCon(pairIdx, pIdx) = r.metrics.PR;
                end
                if isfield(r.metrics, 'totalVariance')
                    varCon(pairIdx, pIdx) = r.metrics.totalVariance;
                end
            end

            if isfield(r, 'projected')
                projectedStore{pairIdx, grpIdx, pIdx} = r.projected;
            end
        end
    end
end

fprintf('K_std = %d\n', K_std);

%% ========================= Figure 1: Scree Plots =======================
% Show 4 representative subjects (2 Med, 2 Con) with parallel analysis

figure('Name', 'Figure 1: Scree + Parallel Analysis', 'Position', [100 100 1200 800]);
repPairs = [1 min(15, numPairs)]; % pick 2 pairs for representative subjects
repPairs = unique(repPairs);
plotIdx = 0;
for rp = 1:length(repPairs)
    pairIdx = repPairs(rp);
    for grpIdx = 1:2
        plotIdx = plotIdx + 1;
        subjectName = pairedSubjectNameList{pairIdx, grpIdx};

        % Load EO1 result for scree
        fileName = fullfile(savedDataFolder, ...
            [subjectName '_EO1_stateSpace.mat']);
        if exist(fileName, 'file')
            r = load(fileName);
            ax = subplot(2, 2, plotIdx);
            grpLabels = {'Meditator', 'Control'};
            plotScreeWithPA(r.eigenvalues, r.threshold95, r.k, K_std, ...
                subjectName, ['EO1 - ' grpLabels{grpIdx}], ax);
        end
    end
end
sgtitle(sprintf('Scree Plots with Parallel Analysis (K_{std} = %d)', K_std));

%% ========================= Figure 2: Embedding Dimension k ==============
figure('Name', 'Figure 2: Embedding Dimension k', 'Position', [100 100 800 400]);

minPointsForViolin = 3; % need at least 3 points for ksdensity
for pIdx = 1:numProtocols
    subplot(1, numProtocols, pIdx);
    validIdx = ~isnan(kMed(:,pIdx)) & ~isnan(kCon(:,pIdx));
    data = {kMed(validIdx, pIdx), kCon(validIdx, pIdx)};
    ds = displaySettings;
    ds.xTickLabels = {'Med', 'Con'};
    ds.commonYLim = 0;
    allVals = [data{1}; data{2}];
    margin = max(1, 0.1 * (max(allVals) - min(allVals) + eps));
    ds.setYLim = [min(allVals)-margin, max(allVals)+margin];
    if sum(validIdx) >= minPointsForViolin
        displayViolinPlot(data, violinColors, 1, 1, 1, displaySettings.pairedDataFlag, ds);
    else
        bar([nanmean(data{1}) nanmean(data{2})]);
        set(gca, 'XTickLabel', {'Med', 'Con'});
    end
    yline(K_std, 'm--', ['K_{std}=' num2str(K_std)], 'LineWidth', 1.5);
    title(protocolsToAnalyze{pIdx});
    ylabel('Embedding Dimension k');

    % Paired Wilcoxon test
    validIdx = ~isnan(kMed(:,pIdx)) & ~isnan(kCon(:,pIdx));
    if sum(validIdx) >= 5
        pVal = signrank(kMed(validIdx,pIdx), kCon(validIdx,pIdx));
    else
        pVal = NaN;
    end
    fprintf('k [%s]: Med=%.1f+/-%.1f, Con=%.1f+/-%.1f, p=%.4f\n', ...
        protocolsToAnalyze{pIdx}, nanmean(kMed(:,pIdx)), nanstd(kMed(:,pIdx))/sqrt(sum(validIdx)), ...
        nanmean(kCon(:,pIdx)), nanstd(kCon(:,pIdx))/sqrt(sum(validIdx)), pVal);
end
sgtitle('Embedding Dimension (k) by Group');

%% ========================= Figure 3: Dynamical Metrics (2x3) ===========
figure('Name', 'Figure 3: Dynamical Metrics', 'Position', [50 50 1400 700]);

metricNames = {'Velocity', 'Log Volume', 'Radius'};
medData = {velMed, volMed, radMed};
conData = {velCon, volCon, radCon};
yLabels = {'Trajectory Velocity', 'Log State Space Volume', 'Excursion Radius'};

for pIdx = 1:numProtocols
    for mIdx = 1:3
        subplot(numProtocols, 3, (pIdx-1)*3 + mIdx);
        validIdx = ~isnan(medData{mIdx}(:,pIdx)) & ~isnan(conData{mIdx}(:,pIdx));
        data = {medData{mIdx}(validIdx, pIdx), conData{mIdx}(validIdx, pIdx)};
        ds = displaySettings;
        ds.xTickLabels = {'Med', 'Con'};
        ds.commonYLim = 0;
        allVals = [data{1}; data{2}];
        margin = 0.1 * (max(allVals) - min(allVals) + eps);
        ds.setYLim = [min(allVals)-margin, max(allVals)+margin];
        if sum(validIdx) >= minPointsForViolin && ~isempty(data{1}) && ~isempty(data{2})
            displayViolinPlot(data, violinColors, 1, 1, 1, displaySettings.pairedDataFlag, ds);
        elseif ~isempty(data{1}) && ~isempty(data{2})
            bar([nanmean(data{1}) nanmean(data{2})]);
            set(gca, 'XTickLabel', {'Med', 'Con'});
        end

        if pIdx == 1; title(metricNames{mIdx}); end
        if mIdx == 1; ylabel(sprintf('%s\n%s', protocolsToAnalyze{pIdx}, yLabels{mIdx}));
        else; ylabel(yLabels{mIdx}); end

        % Statistical test
        if sum(validIdx) > 5
            pVal = signrank(medData{mIdx}(validIdx,pIdx), conData{mIdx}(validIdx,pIdx));
            sigStr = getSigString(pVal);
            text(1.5, max(ylim)*0.95, sprintf('p=%.3f %s', pVal, sigStr), ...
                'HorizontalAlignment', 'center', 'FontSize', 9);
        end
    end
end
sgtitle(sprintf('Dynamical Metrics in K_{std}=%d Dimensional Space', K_std));

%% ========================= Figure 4: TNN Intrinsic Dimension ============
figure('Name', 'Figure 4: TNN Intrinsic Dimension', 'Position', [100 100 1000 400]);

for pIdx = 1:numProtocols
    subplot(1, numProtocols, pIdx);
    validIdx = ~isnan(tnnMed(:,pIdx)) & ~isnan(tnnCon(:,pIdx));
    data = {tnnMed(validIdx, pIdx), tnnCon(validIdx, pIdx)};
    ds = displaySettings;
    ds.xTickLabels = {'Med', 'Con'};
    ds.commonYLim = 0;
    allVals = [data{1}; data{2}];
    margin = 0.1 * (max(allVals) - min(allVals) + eps);
    ds.setYLim = [min(allVals)-margin, max(allVals)+margin];
    if sum(validIdx) >= minPointsForViolin && ~isempty(data{1}) && ~isempty(data{2})
        displayViolinPlot(data, violinColors, 1, 1, 1, displaySettings.pairedDataFlag, ds);
    elseif ~isempty(data{1}) && ~isempty(data{2})
        bar([nanmean(data{1}) nanmean(data{2})]);
        set(gca, 'XTickLabel', {'Med', 'Con'});
    end
    title(protocolsToAnalyze{pIdx});
    ylabel('TNN Intrinsic Dimension (d)');

    if sum(validIdx) > 5
        pVal = signrank(tnnMed(validIdx,pIdx), tnnCon(validIdx,pIdx));
        fprintf('TNN [%s]: Med=%.2f+/-%.2f, Con=%.2f+/-%.2f, p=%.4f\n', ...
            protocolsToAnalyze{pIdx}, nanmean(tnnMed(:,pIdx)), nanstd(tnnMed(:,pIdx))/sqrt(sum(validIdx)), ...
            nanmean(tnnCon(:,pIdx)), nanstd(tnnCon(:,pIdx))/sqrt(sum(validIdx)), pVal);

        % Intra-group std (reliability metric)
        fprintf('  Intra-group std(d): Med=%.2f, Con=%.2f\n', ...
            nanstd(tnnMed(:,pIdx)), nanstd(tnnCon(:,pIdx)));
    end
end
sgtitle('TNN Intrinsic Dimension by Group');

%% ========================= Figure 5: Participation Ratio ===============
% PR is amplitude-independent: it captures how many PCA dimensions share
% variance equally, without being inflated by total signal power.
figure('Name', 'Figure 5: Participation Ratio', 'Position', [100 100 800 400]);

for pIdx = 1:numProtocols
    subplot(1, numProtocols, pIdx);
    validIdx = ~isnan(prMed(:,pIdx)) & ~isnan(prCon(:,pIdx));
    data = {prMed(validIdx, pIdx), prCon(validIdx, pIdx)};
    ds = displaySettings;
    ds.xTickLabels = {'Med', 'Con'};
    ds.commonYLim = 0;
    allVals = [data{1}; data{2}];
    margin = 0.1 * (max(allVals) - min(allVals) + eps);
    ds.setYLim = [max(0, min(allVals)-margin), max(allVals)+margin];
    if sum(validIdx) >= minPointsForViolin && ~isempty(data{1}) && ~isempty(data{2})
        displayViolinPlot(data, violinColors, 1, 1, 1, displaySettings.pairedDataFlag, ds);
    elseif ~isempty(data{1}) && ~isempty(data{2})
        bar([nanmean(data{1}) nanmean(data{2})]);
        set(gca, 'XTickLabel', {'Med', 'Con'});
    end
    title(protocolsToAnalyze{pIdx});
    ylabel('Participation Ratio (PR)');

    if sum(validIdx) > 5
        pVal = signrank(prMed(validIdx,pIdx), prCon(validIdx,pIdx));
        sigStr = getSigString(pVal);
        text(1.5, max(ylim)*0.95, sprintf('p=%.3f %s', pVal, sigStr), ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
        fprintf('PR [%s]: Med=%.2f+/-%.2f, Con=%.2f+/-%.2f, p=%.4f\n', ...
            protocolsToAnalyze{pIdx}, nanmean(prMed(:,pIdx)), nanstd(prMed(:,pIdx))/sqrt(sum(validIdx)), ...
            nanmean(prCon(:,pIdx)), nanstd(prCon(:,pIdx))/sqrt(sum(validIdx)), pVal);
    end
end
sgtitle(sprintf('Participation Ratio (PR) in K_{std}=%d Dimensional Space', K_std));

%% ========================= ANCOVA: Volume ~ Group + TotalVariance =======
% Volume is partially confounded by total signal power (meditators have
% higher gamma power → mathematically larger ellipsoid). ANCOVA disentangles
% the Group effect from the TotalVariance (power) covariate.
fprintf('\n=== ANCOVA: LogVolume ~ Group + TotalVariance ===\n');
fprintf('%-6s %-20s %10s %10s %12s\n', 'Proto', 'Predictor', 'Coeff', 'SE', 'p-value');
fprintf('%s\n', repmat('-', 1, 62));

for pIdx = 1:numProtocols
    validIdx = ~isnan(volMed(:,pIdx)) & ~isnan(volCon(:,pIdx)) & ...
               ~isnan(varMed(:,pIdx)) & ~isnan(varCon(:,pIdx));
    nValid = sum(validIdx);
    if nValid < 6
        fprintf('%s: insufficient data for ANCOVA (n=%d)\n', protocolsToAnalyze{pIdx}, nValid);
        continue;
    end

    % Stack meditator (group=1) and control (group=0) observations
    volVec  = [volMed(validIdx, pIdx); volCon(validIdx, pIdx)];
    grpVec  = [ones(nValid, 1);        zeros(nValid, 1)];
    varVec  = [varMed(validIdx, pIdx); varCon(validIdx, pIdx)];

    tbl = table(volVec, grpVec, varVec, ...
        'VariableNames', {'LogVolume', 'Group', 'TotalVariance'});
    lm  = fitlm(tbl, 'LogVolume ~ Group + TotalVariance');

    % Extract coefficients table (rows: Intercept, Group, TotalVariance)
    ct = lm.Coefficients;
    for rowName = {'Group', 'TotalVariance'}
        rn = rowName{1};
        if ismember(rn, ct.Properties.RowNames)
            fprintf('%-6s %-20s %10.4f %10.4f %12.4f\n', ...
                protocolsToAnalyze{pIdx}, rn, ct{rn,'Estimate'}, ...
                ct{rn,'SE'}, ct{rn,'pValue'});
        end
    end
    fprintf('  R² = %.3f,  Adjusted R² = %.3f\n', lm.Rsquared.Ordinary, lm.Rsquared.Adjusted);
    fprintf('  Interpretation: Group p-value tests whether Meditator > Control\n');
    fprintf('  volume remains significant AFTER controlling for total power.\n\n');
end

%% ========================= Figure 6: 3D Trajectories ===================
figure('Name', 'Figure 6: 3D Trajectory Gallery', 'Position', [50 50 1200 900]);

% Find first pair where all 4 panels have data
repPair = 1;
for tryPair = 1:numPairs
    allPresent = true;
    for gi = 1:2
        for pi = 1:numProtocols
            if isempty(projectedStore{tryPair, gi, pi})
                allPresent = false;
            end
        end
    end
    if allPresent
        repPair = tryPair;
        break;
    end
end

grpLabels = {'Meditator', 'Control'};
plotIdx = 0;
for grpIdx = 1:2
    for pIdx = 1:numProtocols
        plotIdx = plotIdx + 1;
        ax = subplot(2, numProtocols, plotIdx);
        subjectName = pairedSubjectNameList{repPair, grpIdx};
        projected = projectedStore{repPair, grpIdx, pIdx};
        if ~isempty(projected)
            plotTrajectory3D(projected, subjectName, ...
                protocolsToAnalyze{pIdx}, grpLabels{grpIdx}, ax);
        else
            title(sprintf('%s - %s (no data)', subjectName, protocolsToAnalyze{pIdx}));
        end
    end
end
sgtitle('Neural State Trajectories (PC1-PC3)');

%% ========================= Within-Subject EO1 vs M1 ====================
if numProtocols >= 2
    fprintf('\n=== Within-Subject Condition Comparison (EO1 vs M1) ===\n');
    fprintf('%-15s %-12s %-12s %-10s %-10s\n', 'Metric', 'Group', 'EO1 mean', 'M1 mean', 'p-value');
    fprintf('%s\n', repmat('-', 1, 60));

    metricPairs = {velMed, velCon, 'Velocity'; volMed, volCon, 'LogVolume'; ...
                   radMed, radCon, 'Radius'; prMed, prCon, 'PR'; tnnMed, tnnCon, 'TNN_d'};

    for mRow = 1:size(metricPairs, 1)
        for grpIdx = 1:2
            grpLabelsShort = {'Med', 'Con'};
            if grpIdx == 1
                d = metricPairs{mRow, 1};
            else
                d = metricPairs{mRow, 2};
            end
            validIdx = ~isnan(d(:,1)) & ~isnan(d(:,2));
            if sum(validIdx) > 5
                [pVal, ~] = signrank(d(validIdx,1), d(validIdx,2));
                fprintf('%-15s %-12s %-12.4f %-10.4f %-10.4f\n', ...
                    metricPairs{mRow, 3}, grpLabelsShort{grpIdx}, ...
                    nanmean(d(:,1)), nanmean(d(:,2)), pVal);
            end
        end
    end
end

end

%% ========================= Helper Functions =============================
function sigStr = getSigString(pVal)
    if pVal < 0.005
        sigStr = '***';
    elseif pVal < 0.01
        sigStr = '**';
    elseif pVal < 0.05
        sigStr = '*';
    else
        sigStr = 'N.S.';
    end
end
