%% The main funtion to implement the sphercial harmonic analysis
clc; clear; close all; tic
addpath('LinearSphericalConformalMap_stone')
addpath('LinearSphericalConformalMap_stone/mfile')
addpath('LinearSphericalConformalMap_stone/extension')
addpath('LinearSphericalConformalMap_stone/stlTools')
%% The inputs
parametrize_input = false;
maxDeg = 20 ; % the degree for spherical basis calculation
reconDeg = [3, 5, 10, 20]; % tuntitled2he degree used for spherical harmonic reconstrution
% reconDeg = 1:20;
% reconDeg = [50];
reconstruction_mesh_refinment = 2; % Possible inputs 1,2,3,4,5,6, inf? for icospheres
plot_figures = true;

input_mesh_name = 'Stone.stl';
% input_mesh_name = 'StanfordBunny.stl';

% input_stl_files = '/home/mahmoudshaqfa/3d-microstructure-generator/3d-walls-generator/Modular Programming/Saved Layers/STLMESH/';
input_stl_files = '/home/mahmoudshaqfa/3d-microstructure-generator/3d-walls-generator/Modular Programming/SphericalHarmonics/SH_finalized/SH_main/test_geometry';

export_video = false;
video_output = '/home/mahmoudshaqfa/3d-microstructure-generator/3d-walls-generator/Modular Programming/SphericalHarmonics/SH_finalized/SH_main/videos/';

reconstruction_output_mesh = '/home/mahmoudshaqfa/3d-microstructure-generator/3d-walls-generator/Modular Programming/SphericalHarmonics/SH_finalized/SH_main/reconstructed/';
%% Collect the files list and organise inputs/outputs
% Surface parametrization
input_mesh = strcat(input_stl_files, '/', input_mesh_name);
output_mesh = strcat(pwd, '/parametrized/', input_mesh_name);

if parametrize_input
    parametrize_surface(input_mesh, output_mesh)
end

% Read the original mesh
[v, f] = stlRead(input_mesh);
faces = f;
vertices = v;
% Record readings from the origin of the particle
center = mean(vertices,1);
vertices = vertices - center(ones(size(vertices,1),1),:);
[mesh_volume, surface_area] = stlVolume(vertices,faces);

% Read the parametrized input sphere
[para_v, para_f] = stlRead(output_mesh);
spharm_verts = para_v;
spharm_verts = spharm_verts - center(ones(size(spharm_verts,1),1),:);

% Plot the original input mesh - to compare
if plot_figures
    r_original = sqrt(vertices(:,1).^2 + vertices(:,2).^2 + vertices(:,3).^2);
    paraview_patch(vertices, faces, r_original);
    [A, V] = stlVolume(vertices, faces);
    title(sprintf('The original (input) mesh\n Area: %f, Volume: %f'...
        ,A, V ))
end
%% spherical expansion and reconstructon
clear sph_verts faces Z;
[sph_verts, faces] = icosphere(reconstruction_mesh_refinment);
reconDeg = min(reconDeg, maxDeg);
fvec = SH_expansion(vertices,maxDeg); % Spherical expansion

% Shape descriptors (2-norm) for frequency accumulates at a certain
% frequency degree..
Dl = zeros([1, maxDeg+1]);
for l = 1:maxDeg
    for m = -l:1:l
        Dl(1, l) = Dl(1, l) + (real(fvec(l^2 + l + m + 1)))^2 + (imag(fvec(l^2 + l + m + 1)))^2;
    end
    Dl(1, l) = sqrt(Dl(1, l)); % Export only the last object stone.stl
end
Dl(1, :) = Dl(1, :)/Dl(1, 1);
if plot_figures
   figure
   x_temp = 1:length(Dl);
   loglog(x_temp(3:end), Dl(3:end), 'LineWidth', 2)
   title('The shape descriptors')
   xlabel('Freqency degree (l)')
   ylabel('Normalized amplitude')
   grid on
end
%% Reconstruction starts here
idx = 0;
for deg = reconDeg
    idx = idx +1;
    spharm_verts=SH_reconstruction(sph_verts, deg, fvec);
    % Save the reconstruction results
    output_mesh = strcat(reconstruction_output_mesh, 'rec_',...
        num2str(deg), input_mesh_name);
    stlWrite(output_mesh, faces, spharm_verts, 'mode','ascii')
    if plot_figures
        r = sqrt(spharm_verts(:, 1).^2 + spharm_verts(:, 2).^2 + spharm_verts(:, 3).^2);
        paraview_patch(spharm_verts, faces, r);
        [A, V] = stlVolume(spharm_verts, faces);
        title(sprintf('Reconstruction Degree: %1.0f \n Area: %f, Volume: %f'...
        ,deg, A, V ))
        caxis([0.5 0.7]) % color axis limits
        if export_video == true
            saveas(gcf, fullfile(video_output, strcat(num2str(idx), '.png')))
        end
    end
end

if export_video == true
   close all 
end

fprintf('Elapsed time is %1.2f min.\n', toc/60)