function metrics = computeStateSpaceMetrics(projected, K_std, eigenvalues, numTrials)
% computeStateSpaceMetrics - Compute dynamical markers in the standardized
% K_std-dimensional PCA space.
%
% INPUTS:
%   projected   - (T_total x K_std) matrix of PCA-projected data
%   K_std       - standardized embedding dimension (for volume formula)
%   eigenvalues - (optional) full PCA eigenvalues from computeParallelAnalysis,
%                 used for Participation Ratio. If omitted, PR is computed
%                 from the covariance of the projected subspace.
%   numTrials   - (optional) number of concatenated trials. When provided,
%                 velocity is computed only within-trial, ignoring the
%                 discontinuous jump between the last sample of trial N and
%                 the first sample of trial N+1.
%
% OUTPUTS:
%   metrics     - struct with fields:
%     .velocity     - mean within-trial Euclidean step size (trajectory speed)
%     .logVolume    - log hyper-ellipsoid volume from covariance eigenvalues
%     .radius       - mean distance of trajectory points from centroid
%     .PR           - Participation Ratio: (sum λ)² / sum(λ²)
%     .totalVariance - sum of all PCA eigenvalues (total signal power)

T = size(projected, 1);

%% 1. Trajectory Velocity (within-trial only)
% When trials are concatenated end-to-end the diff() between the last sample
% of trial N and the first of trial N+1 is a meaningless "teleportation"
% spike. We mask those boundary steps out before averaging.
diffs = diff(projected, 1, 1);            % (T-1) x K_std displacement vectors
stepSizes = sqrt(sum(diffs.^2, 2));       % (T-1) x 1 Euclidean step sizes

if nargin >= 4 && ~isempty(numTrials) && numTrials > 1
    % Samples per trial (integer, guaranteed by getElectrodeTimeMatrix)
    T_ds = T / numTrials;
    % Boundary positions: the diff at row i corresponds to the step from
    % time point i to i+1. Trial boundaries fall at i = T_ds, 2*T_ds, ...
    boundaryRows = (1 : numTrials - 1) * T_ds;   % (numTrials-1) x 1
    stepSizes(boundaryRows) = [];                  % remove boundary jumps
end

metrics.velocity = mean(stepSizes);

%% 2. State Space Volume (log hyper-ellipsoid)
% Volume of the K_std-dimensional hyper-ellipsoid enclosing the trajectory.
% V = π^(k/2) / Γ(k/2+1) · ∏ sqrt(λ_i)
% log(V) = k/2·log(π) - gammaln(k/2+1) + 0.5·Σ log(λ_i)
covK    = cov(projected);                 % K_std x K_std covariance matrix
eigVals = eig(covK);
eigVals = sort(eigVals, 'descend');

% Clamp floor: 1e-6 rather than eps (~2.2e-16).
% Using eps forces log(eps) ≈ -36, which catastrophically dominates the
% volume sum for subjects whose true dimensionality is slightly below K_std.
% 1e-6 is still negligibly small relative to any real EEG eigenvalue yet
% imposes only a mild penalty (log(1e-6) ≈ -14 before the 0.5 factor).
eigVals(eigVals < 1e-6) = 1e-6;

metrics.logVolume = K_std/2 * log(pi) - gammaln(K_std/2 + 1) + ...
    0.5 * sum(log(eigVals(1:K_std)));

%% 3. Excursion Radius
% Mean distance of trajectory points from the geometric centroid
centroid  = mean(projected, 1);                              % 1 x K_std
distances = sqrt(sum((projected - centroid).^2, 2));         % T x 1
metrics.radius = mean(distances);

%% 4. Participation Ratio
% PR = (Σλ)² / Σ(λ²): amplitude-independent measure of effective
% linear dimensionality. PR = 1 when variance is concentrated in one
% component; PR = N when variance is spread uniformly across N components.
if nargin >= 3 && ~isempty(eigenvalues)
    % Use full PCA eigenvalue spectrum for the most informative PR
    lam = eigenvalues(:);
else
    % Fall back to projected-subspace covariance eigenvalues
    lam = eigVals(1:K_std);
end
metrics.PR = sum(lam)^2 / sum(lam.^2);

%% 5. Total Variance (signal power proxy for ANCOVA covariate)
metrics.totalVariance = sum(lam);

end
