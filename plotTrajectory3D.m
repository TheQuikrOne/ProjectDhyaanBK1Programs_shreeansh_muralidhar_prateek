function plotTrajectory3D(projected, subjectName, protocolName, groupLabel, ax)
% plotTrajectory3D - 3D scatter plot of neural trajectory in PC1-PC3 space.
%
% INPUTS:
%   projected    - (T x K_std) matrix of PCA-projected data (uses first 3 PCs)
%   subjectName  - subject ID for title
%   protocolName - protocol name for title
%   groupLabel   - 'Meditator' or 'Control'
%   ax           - axes handle (optional)

if ~exist('ax', 'var') || isempty(ax)
    figure;
    ax = gca;
end

if size(projected, 2) < 3
    warning('Need at least 3 PCs for 3D plot, only %d available', size(projected, 2));
    return;
end

T = size(projected, 1);
pc1 = projected(:, 1);
pc2 = projected(:, 2);
pc3 = projected(:, 3);

axes(ax); hold on;

% Color by time progression
timeColor = linspace(0, 1, T)';
colormap(ax, parula);

% Plot trajectory as thin line
plot3(pc1, pc2, pc3, '-', 'Color', [0.7 0.7 0.7 0.3], 'LineWidth', 0.3);

% Scatter points colored by time
scatter3(pc1, pc2, pc3, 8, timeColor, 'filled', 'MarkerFaceAlpha', 0.5);

% Mark centroid
centroid = mean(projected(:, 1:3), 1);
scatter3(centroid(1), centroid(2), centroid(3), 150, 'r', 'p', ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

% Mark start and end
scatter3(pc1(1), pc2(1), pc3(1), 80, 'g', 's', 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1);
scatter3(pc1(end), pc2(end), pc3(end), 80, 'm', 'd', 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1);

hold off;

xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title(sprintf('%s - %s (%s)', subjectName, protocolName, groupLabel));

cb = colorbar;
cb.Label.String = 'Time (normalized)';

grid on;
view(30, 20); % default 3D viewing angle

end
