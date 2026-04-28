function [d, mu_values] = computeTNNDimension(dataPoints)
% computeTNNDimension - Two-Nearest Neighbors (TNN) intrinsic dimension
% estimator (Facco et al., 2017).
%
% INPUTS:
%   dataPoints  - (T_total x D) matrix where D is the FULL dimensionality
%                 (electrodes). NOT truncated to k — TNN estimates the true
%                 non-linear intrinsic dimension which may differ from PCA's k.
%
% OUTPUTS:
%   d           - estimated intrinsic dimension (MLE)
%   mu_values   - ratio r2/r1 for each data point (for diagnostic plots)
%
% METHOD:
%   For each point, find distances to 1st and 2nd nearest neighbors.
%   The ratio mu = r2/r1 follows a power-law distribution with exponent d.
%   MLE: d = N / sum(log(mu_i))
%
% If D is very large, we first project with PCA retaining 99% of variance.
% This does NOT restrict to k dimensions — it preserves virtually all structure.
%
% Uses knnsearch (KD-tree) instead of pdist2 for memory efficiency.
% For very large T, subsamples to maxPoints to keep runtime tractable.

[T, D] = size(dataPoints);
maxPoints = 10000; % subsample if T exceeds this to avoid memory/time issues

% Optional variance-preserving dimensionality reduction
maxDimForFull = 500;
if D > maxDimForFull
    [~, scores, ~, ~, explained] = pca(dataPoints);
    cumVar = cumsum(explained);
    nComp = find(cumVar >= 99, 1);
    if isempty(nComp); nComp = length(explained); end
    dataPoints = scores(:, 1:nComp);
    fprintf('  TNN: Reduced from %d to %d dims (99%% variance)\n', D, nComp);
end

% Subsample if too many points (random without replacement)
if T > maxPoints
    subsampleIdx = randperm(T, maxPoints);
    dataPoints = dataPoints(subsampleIdx, :);
    T = maxPoints;
    fprintf('  TNN: Subsampled to %d points\n', T);
end

% Use knnsearch with KD-tree for memory-efficient nearest neighbor search
% Find 2 nearest neighbors (excluding self — knnsearch excludes self by default
% only if query == training, but here they are the same so we search for 3
% and skip the first which is self with distance 0)
[~, dists] = knnsearch(dataPoints, dataPoints, 'K', 3);
% dists columns: [self(=0), 1st-NN, 2nd-NN]
r1 = dists(:, 2); % 1st nearest neighbor distance
r2 = dists(:, 3); % 2nd nearest neighbor distance

% Compute mu = r2/r1
mu_values = r2 ./ r1;

% Remove degenerate cases (mu = 1 when points are equidistant, or r1 = 0)
validIdx = mu_values > 1 & isfinite(mu_values);
mu_valid = mu_values(validIdx);
N_valid = length(mu_valid);

if N_valid < 10
    warning('Too few valid mu values (%d) for reliable TNN estimation', N_valid);
    d = NaN;
    return;
end

% Maximum likelihood estimator for intrinsic dimension
% Under the TNN model: P(mu) = d * mu^(-(d+1)), so log(mu) ~ Exp(1/d)
% MLE: d_hat = N / sum(log(mu))
d = N_valid / sum(log(mu_valid));

end
