% generate trimesh representation for save
%% Loading file paths
clear
close all
addpath('src/utility/')
addpath('src/')
addpath('../mesh2tri/')

[file,path] = uigetfile({'*.csv'},...
    'Please select the .csv file... containing SDF of an object', ...
    '', 'MultiSelect','off');
if isequal(file, 0)
    error('No file selected.');
end

file_path = [path, file];
sdf = csvread(file_path)';
voxelGrid.size = ones(1, 3) * sdf(1);
voxelGrid.range = sdf(2 : 7);
sdf = sdf(8 : end);

voxelGrid.x = linspace(voxelGrid.range(1), voxelGrid.range(2), voxelGrid.size(1));
voxelGrid.y = linspace(voxelGrid.range(3), voxelGrid.range(4), voxelGrid.size(2));
voxelGrid.z = linspace(voxelGrid.range(5), voxelGrid.range(6), voxelGrid.size(3));
[x, y, z] = ndgrid(voxelGrid.x, voxelGrid.y, voxelGrid.z);
voxelGrid.points = reshape(cat(4, x, y, z), [], 3)';
voxelGrid.interval = (voxelGrid.range(2) - voxelGrid.range(1)) / (voxelGrid.size(1) -1);
voxelGrid.truncation = 1.1 * voxelGrid.interval; %1.1 - 1.5
voxelGrid.disp_range = [-inf, voxelGrid.truncation];
voxelGrid.visualizeArclength = 0.01 * sqrt(voxelGrid.range(2) - voxelGrid.range(1));
clearvars x y z

sdf = min(max(sdf, -voxelGrid.truncation), voxelGrid.truncation);


%% marching-primitives
tic
[x] = MPS(sdf, voxelGrid, 'minArea', 5);
toc

%% triangularization and compression
[mesh_original] = meshSuperquadrics(x, 'Arclength', voxelGrid.visualizeArclength);
% compression
mesh = reducepatch(mesh_original.f, mesh_original.v, 0.05); %0.05-2
stl = triangulation(mesh.faces, mesh.vertices);

ifsave = true;
[pathname,name,ext] = fileparts([path, file]);
if ifsave
    x_save = single(x);
    save(fullfile(pathname, [name,'_sq.mat']),'x_save')
    stlwrite(stl,fullfile(pathname, [name,'_sq.stl']), 'binary')
end

%% visualize

close all
view_vector = [151, -40];
light_vector = [190,10];
camera_roll = 50;
color = [145,163,176] ./255;

% rearrange sdf to 3D array for region connection checking
sdf3d_region = reshape(sdf, voxelGrid.size(1), voxelGrid.size(2), voxelGrid.size(3));
% rearrange sdf to 3D array for visualization
sdf3d = permute(sdf3d_region, [2, 1, 3]);


[mesh_gt.f, mesh_gt.v] = plyread(...
    fullfile(pathname, [name,'_watertight.ply']),'tri');

figure(1)
trisurf(mesh_gt.f,mesh_gt.v(:,1),mesh_gt.v(:,2),mesh_gt.v(:,3), ...
    'FaceColor', color, 'FaceAlpha', 1, 'EdgeColor', 'none')
axis equal
view(view_vector)
camroll(camera_roll)
light
lightangle(light_vector(1), light_vector(2))
material dull
axis(voxelGrid.range)
grid off
axis off
title('ground truth mesh from marching cube')

figure(2)
trisurf(mesh.faces,mesh.vertices(:,1),mesh.vertices(:,2),mesh.vertices(:,3), ...
    'FaceColor', color, 'FaceAlpha', 1, 'EdgeColor', 'none')
axis equal
view(view_vector)
camroll(camera_roll)
light
lightangle(light_vector(1), light_vector(2))
material dull
axis(voxelGrid.range)
grid off
axis off
title('superquadrics representation from marching primitives')

figure(3)
trisurf(mesh_gt.f,mesh_gt.v(:,1),mesh_gt.v(:,2),mesh_gt.v(:,3), ...
    'FaceColor', 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none')
hold on
trisurf(mesh.faces,mesh.vertices(:,1),mesh.vertices(:,2),mesh.vertices(:,3), ...
    'FaceColor', color, 'FaceAlpha', 1, 'EdgeColor', 'none')
view(view_vector)
camroll(camera_roll)
axis equal
light
lightangle(light_vector(1), light_vector(2))
material dull
axis(voxelGrid.range)
grid off
axis off
hold off
title('overlapping recovered representation with the ground truth')


