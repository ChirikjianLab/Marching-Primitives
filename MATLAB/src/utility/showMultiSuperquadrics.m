function [] = showMultiSuperquadrics(x, c, arclength)

n = size(x, 1);
rng(5);
% color = rand(100, 3);

for i = 1 : n
%     c = color(i, :);
    showSuperquadrics(x(i, :), 'Color', c, ...
        'FaceAlpha', 1, 'ShowAxis', 1, ...
        'Light', false, 'Arclength', arclength);
    hold on
end
hold off
rng('shuffle')
end