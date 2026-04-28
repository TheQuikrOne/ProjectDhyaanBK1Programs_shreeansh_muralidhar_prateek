function [k, eigenvalues, threshold95, pcCoeffs, pcScores] = ...
    computeParallelAnalysis(dataMatrix, nPermutations)
% computeParallelAnalysis - PCA with Horn's Parallel Analysis to determine
% the embedding dimension k.
%
% INPUTS:
%   dataMatrix     - (N_elec x T_total) electrode-time matrix
%   nPermutations  - number of permutations for parallel analysis (default: 1000)
%
% OUTPUTS:
%   k              - embedding dimension (number of significant components)
%   eigenvalues    - PCA eigenvalues (explained variance), sorted descending
%   threshold95    - 95th percentile threshold from random permutations
%   pcCoeffs       - PCA coefficient matrix (N_elec x N_elec), columns = PCs
%   pcScores       - PCA scores (T_total x N_elec)
%
% METHOD:
%   Horn's Parallel Analysis with Phase Randomization surrogates. For each
%   permutation, the FFT phase of every electrode's time series is randomized
%   (with conjugate symmetry enforced so the IFFT is real). This destroys
%   cross-electrode correlations while preserving the bandpassed power
%   spectrum of each channel — producing "colored noise" surrogates that
%   match the autocorrelation structure of the original filtered data.
%   Time-shuffling (randperm) is NOT used because it produces white-noise
%   surrogates that artificially lower the eigenvalue threshold and inflate k.

if ~exist('nPermutations', 'var'); nPermutations = 1000; end

[N_elec, T_total] = size(dataMatrix);

% Center data: subtract temporal mean from each electrode
centeredData = dataMatrix - mean(dataMatrix, 2);

% PCA on transposed data (samples = time points, features = electrodes)
[pcCoeffs, pcScores, eigenvalues] = pca(centeredData');
% pcCoeffs: N_elec x N_elec (columns are PC loading vectors)
% pcScores: T_total x N_elec (projections)
% eigenvalues: N_elec x 1 (explained variance per component)

numComponents = length(eigenvalues);

% Pre-compute FFTs of all centered electrodes (reused each permutation)
centeredFFT = fft(centeredData, [], 2);   % N_elec x T_total complex

% Indices for independent phase bins (excludes DC at index 1, and Nyquist
% if T_total is even, which must stay real to keep the IFFT real).
% For a signal of length N: floor((N-1)/2) independent phases suffice.
nFreePhases = floor((T_total - 1) / 2);
posFreqIdx  = 2 : nFreePhases + 1;                       % positive frequencies
negFreqIdx  = T_total : -1 : T_total - nFreePhases + 1;  % conjugate mirror

% Horn's Parallel Analysis: phase-randomization surrogate threshold
permEigenvalues = zeros(nPermutations, numComponents);

for perm = 1:nPermutations
    % Draw one set of random phases per electrode (N_elec x nFreePhases)
    randPhases = 2 * pi * rand(N_elec, nFreePhases);

    % Build phase rotation vector for every electrode simultaneously
    phaseVec = zeros(N_elec, T_total);
    phaseVec(:, posFreqIdx) =  randPhases;
    phaseVec(:, negFreqIdx) = -randPhases;  % conjugate symmetry → real IFFT

    % Rotate phases: preserves |X(f)| (power spectrum) but destroys
    % cross-electrode phase coherence, i.e. genuine functional connectivity
    surrogate = real(ifft(centeredFFT .* exp(1i * phaseVec), [], 2));

    [~, ~, permEig] = pca(surrogate');
    permEigenvalues(perm, 1:length(permEig)) = permEig';
end

% 95th percentile threshold at each component position
threshold95 = prctile(permEigenvalues, 95, 1)';

% Embedding dimension k: number of real eigenvalues exceeding threshold
k = sum(eigenvalues > threshold95);

% Ensure k >= 1
k = max(k, 1);

end
