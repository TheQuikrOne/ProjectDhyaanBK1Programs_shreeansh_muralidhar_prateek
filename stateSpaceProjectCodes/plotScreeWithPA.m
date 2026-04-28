function plotScreeWithPA(eigenvalues, threshold95, k, K_std, subjectName, protocolName, ax)
% plotScreeWithPA - Scree plot with Horn's Parallel Analysis threshold.
%
% INPUTS:
%   eigenvalues  - PCA eigenvalues (sorted descending)
%   threshold95  - 95th percentile threshold from parallel analysis
%   k            - individual embedding dimension for this subject
%   K_std        - global standardized dimension (shown as reference)
%   subjectName  - subject ID for title
%   protocolName - protocol name for title
%   ax           - axes handle (optional, creates new figure if empty)

if ~exist('ax', 'var') || isempty(ax)
    figure;
    ax = gca;
end

numToPlot = min(30, length(eigenvalues));

axes(ax); hold on;

% Plot eigenvalues
plot(1:numToPlot, eigenvalues(1:numToPlot), 'b-o', 'LineWidth', 1.5, ...
    'MarkerSize', 5, 'MarkerFaceColor', 'b');

% Plot parallel analysis threshold
plot(1:numToPlot, threshold95(1:numToPlot), 'r--', 'LineWidth', 1.5);

% Mark individual k
if k <= numToPlot
    xline(k + 0.5, 'g-', ['k=' num2str(k)], 'LineWidth', 1.2, ...
        'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'left');
end

% Mark K_std
if K_std <= numToPlot
    xline(K_std + 0.5, 'm--', ['K_{std}=' num2str(K_std)], 'LineWidth', 1.2, ...
        'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'right');
end

hold off;

xlabel('Component');
ylabel('Eigenvalue (Variance Explained)');
title(sprintf('%s - %s (k=%d)', subjectName, protocolName, k));
legend('Data', 'PA Threshold (95%)', 'Location', 'northeast');
xlim([0 numToPlot + 1]);
grid on;

end
