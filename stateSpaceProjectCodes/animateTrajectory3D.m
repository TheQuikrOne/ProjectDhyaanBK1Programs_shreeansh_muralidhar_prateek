function out = animateTrajectory3D(varargin)
% animateTrajectory3D - Two-form 3D trajectory animation primitive.
%
% Persistent growing line + always-visible bright head.
% Trial boundaries become one-frame gaps in the line (no teleportation
% segments) but the head stays visible at every real data point.
%
% INIT FORM:
%   state = animateTrajectory3D(ax, projected, numTrials, viewAngle, lims, titleStr)
%
% STEP FORM:
%   animateTrajectory3D('step', state, frameIdx)

if ischar(varargin{1}) && strcmp(varargin{1}, 'step')
    state = varargin{2};
    fIdx  = varargin{3};
    pt    = state.pcView(fIdx, :);

    % Trial boundary: insert a NaN gap so the line doesn't draw a
    % teleportation segment. We do NOT clear the line, so all prior
    % history stays visible.
    if ismember(fIdx, state.boundaryFrames)
        addpoints(state.line, NaN, NaN, NaN);
    end

    % Always extend the line with the real point and move the head there.
    % The head is never set to NaN, so it never vanishes.
    addpoints(state.line, pt(1), pt(2), pt(3));
    set(state.headPoint, 'XData', pt(1), 'YData', pt(2), 'ZData', pt(3));
    return;
end

%% --- INIT ---
ax        = varargin{1};
projected = varargin{2};
numTrials = varargin{3};
viewAngle = varargin{4};
lims      = varargin{5};
if nargin >= 6; titleStr = varargin{6}; else; titleStr = ''; end

pcView = projected(:, 1:3);
T = size(pcView, 1);

% Trial boundary frame indices (computed in subsampled coordinates)
if numTrials > 1
    T_ds = T / numTrials;
    boundaryFrames = round((1 : numTrials - 1) * T_ds) + 1;
    boundaryFrames(boundaryFrames < 1 | boundaryFrames > T) = [];
else
    boundaryFrames = [];
end

axes(ax); hold(ax, 'on');

% Persistent growing line — orange, moderate weight, no max-points cap
state.line = animatedline(ax, ...
    'Color', [1.0 0.55 0.0], ...
    'LineWidth', 1.3);

% Bright head marker — never NaN, always at current position
state.headPoint = scatter3(ax, NaN, NaN, NaN, 130, ...
    [1.0 0.15 0.15], 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.7);

% Persistent green start marker
sp = pcView(1, :);
scatter3(ax, sp(1), sp(2), sp(3), 90, [0.1 0.7 0.1], 's', 'filled', ...
    'MarkerEdgeColor', 'k');

% Axes setup
if ~isempty(lims)
    xlim(ax, lims(:,1)); ylim(ax, lims(:,2)); zlim(ax, lims(:,3));
end
view(ax, viewAngle);
axis(ax, 'vis3d');
grid(ax, 'on');
xlabel(ax, 'PC1'); ylabel(ax, 'PC2'); zlabel(ax, 'PC3');
if ~isempty(titleStr); title(ax, titleStr, 'FontSize', 12); end
set(ax, 'BoxStyle', 'full', 'Box', 'on');

state.pcView         = pcView;
state.boundaryFrames = boundaryFrames;
state.numFrames      = T;
state.ax             = ax;

out = state;
end
