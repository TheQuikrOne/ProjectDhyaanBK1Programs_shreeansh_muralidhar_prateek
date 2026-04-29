function generateExpandedReport(dataFolder, protocols, medState, bandLabel)
% generateExpandedReport - Print findings.md-style stats for an expanded dataset.
% INPUTS:
%   dataFolder - subfolder under stateSpaceProjectCodes (e.g. 'savedDataAlpha_AllStates')
%   protocols  - cell array, e.g. {'EO1','EC1','M1'} (the meditation state must be among them)
%   medState   - the meditation-state protocol name for between-group tests
%   bandLabel  - 'Alpha' or 'Gamma' (for header)

addpath(genpath('C:\AIML\patagonia\nsp_project'));

savedDataFolder = fullfile(fileparts(mfilename('fullpath')), dataFolder);
pairedSubjectNameList = getPairedSubjectsBK1;
numPairs = size(pairedSubjectNameList,1);
nP = length(protocols);

K_std = NaN;
data = struct();
metricNames = {'k','velocity','logVolume','radius','PR','TNN_d','totalVariance'};
for m = 1:length(metricNames)
    data.(metricNames{m}) = struct('Med', NaN(numPairs, nP), 'Con', NaN(numPairs, nP));
end

for pairIdx = 1:numPairs
    for grpIdx = 1:2
        if grpIdx == 1; gl = 'Med'; else; gl = 'Con'; end
        for pIdx = 1:nP
            fname = fullfile(savedDataFolder, ...
                [pairedSubjectNameList{pairIdx,grpIdx} '_' protocols{pIdx} '_stateSpace.mat']);
            if ~exist(fname,'file'); continue; end
            r = load(fname);
            if isnan(K_std) && isfield(r,'K_std'); K_std = r.K_std; end
            data.k.(gl)(pairIdx,pIdx) = r.k;
            data.velocity.(gl)(pairIdx,pIdx) = r.metrics.velocity;
            data.logVolume.(gl)(pairIdx,pIdx) = r.metrics.logVolume;
            data.radius.(gl)(pairIdx,pIdx) = r.metrics.radius;
            data.TNN_d.(gl)(pairIdx,pIdx) = r.tnnD;
            if isfield(r.metrics,'PR'); data.PR.(gl)(pairIdx,pIdx) = r.metrics.PR; end
            if isfield(r.metrics,'totalVariance'); data.totalVariance.(gl)(pairIdx,pIdx) = r.metrics.totalVariance; end
        end
    end
end

medIdx = find(strcmp(protocols, medState), 1);

fprintf('\n## %s-Band Expanded Results (K_std = %d)\n\n', bandLabel, K_std);
fprintf('Protocols: %s\n\n', strjoin(protocols, ', '));

%% Between-group at meditation state
fprintf('### Between-Group at %s: Meditators vs Controls\n\n', medState);
fprintf('| Metric | Meditators (mean +/- SEM) | Controls (mean +/- SEM) | p-value |\n');
fprintf('|---|---|---|---|\n');
displayMetrics = {'k','Velocity','Log Volume','Radius','PR','TNN d'};
metricFields  = {'k','velocity','logVolume','radius','PR','TNN_d'};
for m = 1:length(metricFields)
    M = data.(metricFields{m}).Med(:,medIdx);
    C = data.(metricFields{m}).Con(:,medIdx);
    v = ~isnan(M) & ~isnan(C); n = sum(v);
    if n < 5; continue; end
    p = signrank(M(v), C(v));
    fprintf('| %s | %.2f +/- %.2f | %.2f +/- %.2f | %.4f%s |\n', ...
        displayMetrics{m}, ...
        nanmean(M), nanstd(M)/sqrt(n), nanmean(C), nanstd(C)/sqrt(n), p, sigStar(p));
end

%% ANCOVA at meditation state
fprintf('\n### ANCOVA at %s: LogVolume ~ Group + TotalVariance\n\n', medState);
fprintf('| Predictor | Coefficient | SE | p-value |\n');
fprintf('|---|---|---|---|\n');
volMed = data.logVolume.Med(:,medIdx); volCon = data.logVolume.Con(:,medIdx);
varMed = data.totalVariance.Med(:,medIdx); varCon = data.totalVariance.Con(:,medIdx);
v = ~isnan(volMed) & ~isnan(volCon) & ~isnan(varMed) & ~isnan(varCon);
n = sum(v);
if n >= 6
    volV = [volMed(v); volCon(v)];
    grpV = [ones(n,1); zeros(n,1)];
    varV = [varMed(v); varCon(v)];
    tbl = table(volV, grpV, varV, 'VariableNames', {'LogVolume','Group','TotalVariance'});
    lm = fitlm(tbl, 'LogVolume ~ Group + TotalVariance');
    ct = lm.Coefficients;
    fprintf('| Group (Med vs Con) | %+.4f | %.4f | %.4f%s |\n', ...
        ct{'Group','Estimate'}, ct{'Group','SE'}, ct{'Group','pValue'}, sigStar(ct{'Group','pValue'}));
    fprintf('| TotalVariance | %+.5f | %.5f | %.4g%s |\n', ...
        ct{'TotalVariance','Estimate'}, ct{'TotalVariance','SE'}, ct{'TotalVariance','pValue'}, sigStar(ct{'TotalVariance','pValue'}));
    fprintf('| (R^2 = %.3f, n = %d per group) | | | |\n', lm.Rsquared.Ordinary, n);
end

%% Within-group transitions
% All 3 pairwise transitions
transitions = nchoosek(1:nP, 2);

for t = 1:size(transitions,1)
    a = transitions(t,1); b = transitions(t,2);
    fprintf('\n### Within-Group: %s vs %s\n\n', protocols{a}, protocols{b});
    fprintf('| Metric | Group | %s (mean +/- SEM) | %s (mean +/- SEM) | p-value |\n', protocols{a}, protocols{b});
    fprintf('|---|---|---|---|---|\n');
    for m = 1:length(metricFields)
        for grpIdx = 1:2
            if grpIdx == 1; gl = 'Med'; gLabel = 'Meditators'; else; gl = 'Con'; gLabel = 'Controls'; end
            D = data.(metricFields{m}).(gl);
            v = ~isnan(D(:,a)) & ~isnan(D(:,b)); n = sum(v);
            if n < 5; continue; end
            p = signrank(D(v,a), D(v,b));
            fprintf('| %s | %s | %.2f +/- %.2f | %.2f +/- %.2f | %.4f%s |\n', ...
                displayMetrics{m}, gLabel, ...
                nanmean(D(:,a)), nanstd(D(:,a))/sqrt(n), ...
                nanmean(D(:,b)), nanstd(D(:,b))/sqrt(n), p, sigStar(p));
        end
    end
end

end

function s = sigStar(p)
    if p < 0.005; s = ' ★★★';
    elseif p < 0.01; s = ' ★★';
    elseif p < 0.05; s = ' ★';
    else; s = '';
    end
end
