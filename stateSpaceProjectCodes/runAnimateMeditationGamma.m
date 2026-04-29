% runAnimateMeditationGamma
% Generates a side-by-side animated 3D trajectory of gamma-band PC1-PC3
% evolution for a single auto-picked meditator: EO1 baseline (left) vs M2
% meditation (right). Saves as MP4.

%% ====== Tunables ========================================================
STEP        = 10;       % temporal subsampling factor (1 = full 250Hz)
PCTILE_CLIP = 1.0;      % clip axis bounds to [PCTILE, 100-PCTILE] of data per axis (avoids outlier blow-up)
VIEW_ANGLE  = [45 25];  % azimuth, elevation
FRAME_RATE  = 60;       % output video FPS
QUALITY     = 90;       % VideoWriter Quality (1-100)
OUTPUT_NAME = 'gamma_meditation_evolution.mp4';

%% ====== Setup ===========================================================
addpath(genpath('C:\AIML\patagonia\nsp_project'));
thisDir = fileparts(mfilename('fullpath'));
dataFolder = fullfile(thisDir, 'savedDataGamma_AllStates');
outputPath = fullfile(thisDir, '..', OUTPUT_NAME);

%% ====== Pick meditator with the strongest EO1->M2 expansion =============
[subjectName, ~, deltaVol] = pickBestMeditatorM2(dataFolder);

%% ====== Load EO1 and M2 projected trajectories ==========================
rEO1 = load(fullfile(dataFolder, [subjectName '_EO1_stateSpace.mat']));
rM2  = load(fullfile(dataFolder, [subjectName '_M2_stateSpace.mat']));

projEO1 = rEO1.projected(:, 1:3);
projM2  = rM2.projected(:,  1:3);
nTrEO1  = rEO1.numGoodTrials;
nTrM2   = rM2.numGoodTrials;

% Subsample for speed (the comet doesn't need 250Hz to look smooth)
projEO1_ds = projEO1(1:STEP:end, :);
projM2_ds  = projM2(1:STEP:end,  :);
% Trial counts stay the same; per-trial sample count shrinks proportionally
% but boundary positions need to be recomputed from the new T
nTrEO1_ds = nTrEO1;
nTrM2_ds  = nTrM2;

fprintf('EO1: %d frames (%d trials)  |  M2: %d frames (%d trials)\n', ...
    size(projEO1_ds,1), nTrEO1_ds, size(projM2_ds,1), nTrM2_ds);

%% ====== Per-panel axis bounds via percentiles (outlier-safe) ============
% Per-panel autoscale gets blown up by a few outlier trials, making the
% bulk of the trajectory look small. Instead, use the [PCTILE_CLIP,
% 100-PCTILE_CLIP] percentile per axis so 98% of the trajectory fills the
% box (the rare outliers just clip outside, but the head will still be
% drawn there — MATLAB shows it fine even outside axis limits if 'Clipping'
% is off, but we'll just live with occasional clipping for the bulk view).
mgFrac = 0.05;
mn = prctile(projEO1_ds, PCTILE_CLIP,        1);
mx = prctile(projEO1_ds, 100 - PCTILE_CLIP,  1);
limsEO1 = [mn - mgFrac*(mx-mn); mx + mgFrac*(mx-mn)];
mn = prctile(projM2_ds,  PCTILE_CLIP,        1);
mx = prctile(projM2_ds,  100 - PCTILE_CLIP,  1);
limsM2  = [mn - mgFrac*(mx-mn); mx + mgFrac*(mx-mn)];

fprintf('EO1 axis ranges: PC1 [%.1f, %.1f]  PC2 [%.1f, %.1f]  PC3 [%.1f, %.1f]\n', ...
    limsEO1(1,1), limsEO1(2,1), limsEO1(1,2), limsEO1(2,2), limsEO1(1,3), limsEO1(2,3));
fprintf('M2  axis ranges: PC1 [%.1f, %.1f]  PC2 [%.1f, %.1f]  PC3 [%.1f, %.1f]\n', ...
    limsM2(1,1),  limsM2(2,1),  limsM2(1,2),  limsM2(2,2),  limsM2(1,3),  limsM2(2,3));

%% ====== Build figure ====================================================
fig = figure('Name', sprintf('Gamma trajectory evolution — %s', subjectName), ...
    'Position', [80 80 1500 720], 'Color', 'w');

axEO1 = subplot(1, 2, 1);
sEO1 = animateTrajectory3D(axEO1, projEO1_ds, nTrEO1_ds, VIEW_ANGLE, ...
    limsEO1, sprintf('%s  —  EO1 baseline    LogVol = %.2f', ...
    subjectName, rEO1.metrics.logVolume));

axM2  = subplot(1, 2, 2);
sM2  = animateTrajectory3D(axM2,  projM2_ds,  nTrM2_ds,  VIEW_ANGLE, ...
    limsM2, sprintf('%s  —  M2 meditation    LogVol = %.2f   (Δ = %+.2f)', ...
    subjectName, rM2.metrics.logVolume, deltaVol));

%% ====== Open VideoWriter (MP4 with AVI fallback) ========================
try
    vw = VideoWriter(outputPath, 'MPEG-4');
    vw.Quality = QUALITY;
catch
    outputPath = strrep(outputPath, '.mp4', '.avi');
    vw = VideoWriter(outputPath, 'Motion JPEG AVI');
    vw.Quality = QUALITY;
    fprintf('MP4 unavailable, falling back to AVI: %s\n', outputPath);
end
vw.FrameRate = FRAME_RATE;
open(vw);

%% ====== Animation loop ==================================================
nFrames = max(sEO1.numFrames, sM2.numFrames);
fprintf('Rendering %d frames -> %s\n', nFrames, outputPath);
tStart = tic;

for f = 1:nFrames
    if f <= sEO1.numFrames; animateTrajectory3D('step', sEO1, f); end
    if f <= sM2.numFrames;  animateTrajectory3D('step', sM2,  f); end

    if mod(f, 30) == 0
        sgtitle(fig, sprintf('Frame %d / %d   (%.0f%%)', f, nFrames, 100*f/nFrames));
    end

    drawnow limitrate;
    % getframe() is unreliable under -batch / no-display modes; print() works
    try
        frame = getframe(fig);
    catch
        cdata = print(fig, '-RGBImage', '-r0');
        frame = im2frame(cdata);
    end
    writeVideo(vw, frame);
end

close(vw);
fprintf('Done in %.1f s. Saved to: %s\n', toc(tStart), outputPath);
