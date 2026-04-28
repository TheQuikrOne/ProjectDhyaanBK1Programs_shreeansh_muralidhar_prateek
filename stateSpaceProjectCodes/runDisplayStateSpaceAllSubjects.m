% runDisplayStateSpaceAllSubjects
% Configuration and runner script for displaying state-space analysis results.
% Run this after runSaveStateSpaceData.m has completed.

%% ========================= Configuration ================================

% Protocols analyzed
protocolsToAnalyze = {'EO1', 'M1'};

% Path to saved results
savedDataFolder = fullfile(fileparts(mfilename('fullpath')), 'savedDataAlpha');

% Get paired subject list
pairedSubjectNameList = getPairedSubjectsBK1;

%% ========================= Display Settings =============================
displaySettings.parametricTest = 0;   % 0 = non-parametric (Wilcoxon/Mann-Whitney)
displaySettings.pairedDataFlag = 1;   % 1 = show paired lines in violin plots
displaySettings.medianFlag = 1;       % 1 = show median, 0 = show mean

%% ========================= Run Display ==================================
displayStateSpaceAllSubjects(pairedSubjectNameList, protocolsToAnalyze, ...
    savedDataFolder, displaySettings);

%% ========================= Auto-Save All Figures ========================
figs = findall(0, 'Type', 'figure');
fprintf('\nSaving %d figures to %s\n', length(figs), savedDataFolder);
for i = 1:length(figs)
    name = figs(i).Name;
    if isempty(name); name = ['figure_' num2str(figs(i).Number)]; end
    safeName = regexprep(name, '[^a-zA-Z0-9_]', '_');
    saveas(figs(i), fullfile(savedDataFolder, [safeName '.png']));
    saveas(figs(i), fullfile(savedDataFolder, [safeName '.fig']));
    fprintf('  Saved: %s.png + .fig\n', safeName);
end
fprintf('All figures saved.\n');
