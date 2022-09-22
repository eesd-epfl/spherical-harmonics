clear all; close all; clc;
load('reconstruction_data.mat')

x = spharm_verts(:, 1);
y = spharm_verts(:, 2);
z = spharm_verts(:, 3);

[phi, theta, r] = cart2sph(spharm_verts(:,1), spharm_verts(:,2), spharm_verts(:,3));
theta = pi/2 - theta; neg_angles = find(phi<0); phi(neg_angles) = phi(neg_angles) + 2*pi;

fig1=figure(1);
fig1.Renderer='Painters'; % To force saving vectorized figures only in Matlab not Octave

tri = delaunay(phi,theta);

red_color = zeros(1, 255); green_color = zeros(1, 255); blue_color = zeros(1, 255);
red_color(1:127) = linspace(42,220,127); red_color(128:255) = linspace(220,174,127+1);
green_color(1:127) = linspace(63,220,127); green_color(128:255) = linspace(220,0,127+1);
blue_color(1:127) = linspace(181,220,127); blue_color(128:255) = linspace(220,22,127+1);

Paraviewmap = [red_color', green_color', blue_color']./255;

colormap(Paraviewmap);

% colormap ('default')
% trisurf(tri,phi,theta, r, 'LineWidth', 0.000001, 'LineStyle', 'none')
% trimesh(tri,phi,theta, r, 'LineWidth', 0.000001)
trisurf(tri,phi,theta, r, 'LineWidth', 0.000001)

xlim([0, 2*pi])
ylim([0, pi])
set(gca, 'xtick', [0.  pi 2*pi])
set(gca, 'ytick', [0.   0.5 * pi pi])
axis('equal')

% ##colorbar;
cb = colorbar; 
% set(cb,'position',[.15 .05 .05 .3])

grid('off')
view(0,90)

