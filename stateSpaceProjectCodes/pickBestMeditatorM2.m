function [subjectName, controlName, deltaVol, allDeltas] = pickBestMeditatorM2(dataFolder)
% pickBestMeditatorM2 - Pick the meditator whose gamma trajectory expands
% most dramatically from EO1 to M2 (largest LogVolume increase).
%
% Used to select the most visually compelling subject for the side-by-side
% animated 3D trajectory comparison.
%
% INPUT:
%   dataFolder - absolute path to savedDataGamma_AllStates folder
%
% OUTPUTS:
%   subjectName - meditator subject ID (e.g. '050UR')
%   controlName - matched control subject ID (kept for future side-by-side
%                 with control rather than EO1)
%   deltaVol    - LogVolume(M2) - LogVolume(EO1) for the picked meditator
%   allDeltas   - per-pair table (for inspection)

pairedSubjectNameList = getPairedSubjectsBK1;
numPairs = size(pairedSubjectNameList, 1);

deltas = NaN(numPairs, 1);
volEO1 = NaN(numPairs, 1);
volM2  = NaN(numPairs, 1);

for pairIdx = 1:numPairs
    medName = pairedSubjectNameList{pairIdx, 1};
    fEO1 = fullfile(dataFolder, [medName '_EO1_stateSpace.mat']);
    fM2  = fullfile(dataFolder, [medName '_M2_stateSpace.mat']);
    if ~exist(fEO1, 'file') || ~exist(fM2, 'file'); continue; end
    rEO1 = load(fEO1, 'metrics');
    rM2  = load(fM2,  'metrics');
    volEO1(pairIdx) = rEO1.metrics.logVolume;
    volM2(pairIdx)  = rM2.metrics.logVolume;
    deltas(pairIdx) = volM2(pairIdx) - volEO1(pairIdx);
end

[deltaVol, bestIdx] = max(deltas);
subjectName = pairedSubjectNameList{bestIdx, 1};
controlName = pairedSubjectNameList{bestIdx, 2};

allDeltas = table((1:numPairs)', pairedSubjectNameList(:,1), volEO1, volM2, deltas, ...
    'VariableNames', {'pairIdx', 'meditator', 'volEO1', 'volM2', 'deltaVol'});

fprintf('Picked meditator: %s (pair %d)  ΔLogVolume EO1→M2 = %+.2f\n', ...
    subjectName, bestIdx, deltaVol);
fprintf('  EO1 LogVolume = %.2f, M2 LogVolume = %.2f\n', volEO1(bestIdx), volM2(bestIdx));

end
