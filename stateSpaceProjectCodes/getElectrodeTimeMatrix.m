function [dataMatrix, goodElectrodeList, Fs_out, numGoodTrials] = ...
    getElectrodeTimeMatrix(folderSourceString, subjectName, gridType, expDate, ...
    protocolName, badEyeCondition, badTrialVersion, stRange, gammaRange, targetFs)
% getElectrodeTimeMatrix - Build gamma-filtered, downsampled, trial-concatenated
% electrode-time matrix for state-space analysis.
%
% INPUTS:
%   folderSourceString - root path (e.g., 'C:\AIML\patagonia')
%   subjectName        - subject ID (e.g., '003S')
%   gridType           - electrode grid type (default: 'EEG')
%   expDate            - experiment date string (e.g., '221021')
%   protocolName       - protocol name ('EO1' or 'M1')
%   badEyeCondition    - 'ep' (eye position) or 'wo' (without eye data)
%   badTrialVersion    - bad trial version string (e.g., 'v8')
%   stRange            - [startTime endTime] in seconds (e.g., [0.25 1.25])
%   gammaRange         - [lowFreq highFreq] in Hz (e.g., [30 80])
%   targetFs           - target sampling rate after downsampling (e.g., 250)
%
% OUTPUTS:
%   dataMatrix         - (N_good_elec x T_total) concatenated trajectory
%   goodElectrodeList  - list of electrode numbers used (1-indexed)
%   Fs_out             - actual output sampling rate after resampling
%   numGoodTrials      - number of good trials used
%
% NOTE: Trials are CONCATENATED, not averaged. Spontaneous gamma is not
% phase-locked to trial onset, so averaging would cancel the signal.

if ~exist('gridType','var');        gridType = 'EEG';           end
if ~exist('badEyeCondition','var'); badEyeCondition = 'ep';     end
if ~exist('badTrialVersion','var'); badTrialVersion = 'v8';     end
if ~exist('stRange','var');         stRange = [0.25 1.25];      end
if ~exist('gammaRange','var');      gammaRange = [30 80];       end
if ~exist('targetFs','var');        targetFs = 250;             end

electrodeList = 1:64;
cutoffNumTrials = 30;
cutoffNumElectrodes = 10;

% Build path to segmented data
folderSegment = fullfile(folderSourceString, 'data', 'segmentedData', ...
    subjectName, gridType, expDate, protocolName, 'segmentedData');

% Load timing information
timingFile = fullfile(folderSegment, 'LFP', 'lfpInfo.mat');
if ~exist(timingFile, 'file')
    warning('Timing file does not exist: %s', timingFile);
    dataMatrix = []; goodElectrodeList = []; Fs_out = targetFs; numGoodTrials = 0;
    return;
end

t = load(timingFile);
timeVals = t.timeVals;
Fs = round(1 / (timeVals(2) - timeVals(1)));

% Time indices for stimulus window
stPos = find(timeVals >= stRange(1) & timeVals < stRange(2));

% Get bad trials and bad electrodes
[badTrials, badElecs] = getBadTrialsAndElectrodes(subjectName, expDate, ...
    protocolName, folderSourceString, badEyeCondition, badTrialVersion);

% Determine good electrodes
allBadElecs = badElecs(:);
goodElectrodeList = setdiff(electrodeList, allBadElecs);

if length(goodElectrodeList) < cutoffNumElectrodes
    warning('Subject %s protocol %s: only %d good electrodes (< %d cutoff)', ...
        subjectName, protocolName, length(goodElectrodeList), cutoffNumElectrodes);
    dataMatrix = []; Fs_out = targetFs; numGoodTrials = 0;
    return;
end

% Load one electrode to determine trial count
e = load(fullfile(folderSegment, 'LFP', 'elec1.mat'));
totalTrials = size(e.analogData, 1);
goodTrials = setdiff(1:totalTrials, badTrials);
numGoodTrials = length(goodTrials);

if numGoodTrials < cutoffNumTrials
    warning('Subject %s protocol %s: only %d good trials (< %d cutoff)', ...
        subjectName, protocolName, numGoodTrials, cutoffNumTrials);
    dataMatrix = []; Fs_out = targetFs;
    return;
end

% Design bandpass filter for gamma isolation (30-80 Hz)
filterOrder = 4;
[b, a] = butter(filterOrder, gammaRange / (Fs/2), 'bandpass');

% Compute downsampling parameters
% resample(x, P, Q) resamples at P/Q * original rate
[P, Q] = rat(targetFs / Fs);
T_ds = length(resample(zeros(1, length(stPos)), P, Q)); % downsampled length per trial

% Pre-allocate output: electrodes x (numGoodTrials * T_ds)
T_total = numGoodTrials * T_ds;
numGoodElecs = length(goodElectrodeList);
dataMatrix = zeros(numGoodElecs, T_total);

Fs_out = targetFs;

% Process each good electrode
for elecIdx = 1:numGoodElecs
    elecNum = goodElectrodeList(elecIdx);
    e = load(fullfile(folderSegment, 'LFP', ['elec' num2str(elecNum) '.mat']));

    % Extract good trials, stimulus window
    goodTrialData = e.analogData(goodTrials, stPos); % (nGoodTrials x T)

    % Filter and downsample each trial
    concatenated = zeros(1, T_total);
    for trIdx = 1:numGoodTrials
        trialSignal = goodTrialData(trIdx, :); % 1 x T

        % Zero-phase bandpass filter (gamma: 30-80 Hz)
        filteredSignal = filtfilt(b, a, trialSignal);

        % Downsample (resample includes anti-aliasing filter)
        dsSignal = resample(filteredSignal, P, Q);

        % Handle edge case: resample output length may vary by 1 sample
        nSamples = min(length(dsSignal), T_ds);
        startIdx = (trIdx - 1) * T_ds + 1;
        concatenated(startIdx : startIdx + nSamples - 1) = dsSignal(1:nSamples);
    end

    dataMatrix(elecIdx, :) = concatenated;
end

end
