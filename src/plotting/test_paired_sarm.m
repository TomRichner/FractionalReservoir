% Example data: N subjects measured at two conditions
N = 100;
subj = (1:N)';
yA = randn(N,1);
yB = yA + 0.5*randn(N,1) + 0.4;

x = [repmat("A",N,1); repmat("B",N,1)];
y = [yA; yB];
pairID = [subj; subj];

figure; ax = axes; hold(ax,'on');

% Use numeric x positions for robustness, keep labels
xCat  = categorical(x);
xCats = categories(xCat);
xPos  = double(xCat);  % 1..G

hs = swarmchart(ax, xPos, y, 36, 'filled');  % markersize=36 (points^2)
hs.XJitter = 'density';
hs.XJitterWidth = 1;  % reduce from default 0.5 to decrease spread
drawnow;  % ensure jitter is computed

% Extract jittered coordinates
xyz = struct(hs).XYZJittered;          % undocumented, but works
xj  = xyz(:,1);
yj  = xyz(:,2);

% Draw paired lines
u = unique(pairID);
hl = gobjects(numel(u),1);
for k = 1:numel(u)
    idx = (pairID == u(k)) & ~isnan(yj) & ~isnan(xj);
    if nnz(idx) < 2, continue; end
    [xs,ord] = sort(xj(idx));
    ys = yj(idx); ys = ys(ord);
    hl(k) = line(ax, xs, ys, 'LineWidth', 0.5, 'Color', [0 0 0]); % adjust styling as desired
end

% Put points on top of lines
uistack(hs, 'top');

% Restore category labels
ax.XTick = 1:numel(xCats);
ax.XTickLabel = xCats;
xlim(ax, [0.5, numel(xCats)+0.5]);
box(ax,'on');
