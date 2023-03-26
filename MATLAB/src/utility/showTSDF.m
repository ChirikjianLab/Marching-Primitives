function [] = showTSDF(sdf, grid)
sdf(or(sdf <= grid.disp_range(1), sdf >= grid.disp_range(2))) = nan;

h = slice(grid.x, grid.y, grid.z, sdf, [], [], grid.z);

set(h, 'EdgeColor','none', 'FaceColor','interp')
set(h, 'EdgeColor','none')
alpha(.2)
axis equal
end