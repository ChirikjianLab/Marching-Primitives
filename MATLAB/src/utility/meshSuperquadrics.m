function [mesh_merge] = meshSuperquadrics(x, varargin)

n = size(x, 1);

% with tapering or not
taper = false;
% color = 'r';
arclength = 0.02;

for k = 1 : size(varargin, 2)
    if strcmp(varargin{k}, 'Taper')
        taper = varargin{k + 1};
    end
%     if strcmp(varargin{k}, 'Color')
%         color = varargin{k + 1};
%     end   
    if strcmp(varargin{k}, 'Arclength')
        arclength= varargin{k + 1};
    end    
end

mesh = cell(n, 1);
v_num = zeros(n, 1);
f_num = zeros(n, 1);


for i = 1 : n
    [x_mesh, y_mesh, z_mesh] = singleMeshSuperquadrics(x(i, :), ...
        'Arclength', arclength, 'Taper', taper);
    [mesh{i}.f, mesh{i}.v] = mesh2tri(x_mesh, y_mesh, z_mesh,'b'); % f b x
    v_num(i) = size(mesh{i}.v, 1);
    f_num(i) = size(mesh{i}.f, 1);
    if i > 1
        mesh{i}.f = mesh{i}.f + sum(v_num(1 : i - 1));
    end
end
mesh_merge.v = zeros(sum(v_num), 3);
mesh_merge.f = zeros(sum(f_num), 3);
mesh_merge.v(1 : v_num(1), :) = mesh{1}.v;
mesh_merge.f(1 : f_num(1), :) = mesh{1}.f;
for i = 2 : n
    mesh_merge.v(...
        sum(v_num(1 : i - 1)) + 1 : sum(v_num(1 : i - 1)) + v_num(i), :) = ...
        mesh{i}.v;
    mesh_merge.f(...
        sum(f_num(1 : i - 1)) + 1 : sum(f_num(1 : i - 1)) + f_num(i), :) = ...
        mesh{i}.f;
end
end

