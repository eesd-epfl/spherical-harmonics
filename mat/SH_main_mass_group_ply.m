%% The main funtion to implement the sphercial harmonic analysis
clc; clear; close all; tic
% addpath('LinearSphericalConformalMap_stone')
addpath('LinearSphericalConformalMap_stone/mfile')
addpath('LinearSphericalConformalMap_stone/extension')
addpath('LinearSphericalConformalMap_stone/stlTools')

%% The inputs
parametrize_input = false;
maxDeg = 15 ; % the degree for spherical basis calculation
reconDeg = 15; % tuntitled2he degree used for spherical harmonic reconstrution
reconstruction_mesh_refinment = 2; % Possible inputs 1,2,3,4,5,6 for icospheres
plot_figures = true;

input_stl_files = '/home/mahmoudshaqfa/Desktop/PLYS/';

% Parametrized meshes if it doesn't exist use the default directory
input_para_stl = "";
%% Collect the files list and organise inputs/outputs
fileList = dir(fullfile(input_stl_files, '*.ply'));
fileList = struct2cell(fileList);
n = size(fileList, 2);

Dl = zeros([n, maxDeg+1]);

for stl_idx=1:n
    % Surface parametrization
    input_mesh_name = fileList{1, stl_idx};
    input_mesh = strcat(input_stl_files, '/', input_mesh_name);
    reconstruction_output_mesh = strcat(input_stl_files, '/Output/', input_mesh_name);
    reconstruction_output_mesh = strcat(reconstruction_output_mesh(1:end-3), 'stl');
    
    if (input_para_stl == "")
        output_mesh = strcat(pwd, '/parametrized/', input_mesh_name);
    else
        output_mesh = strcat(input_para_stl, input_mesh_name);
    end
    
    disp(strcat('Analysing stone: ', input_mesh_name))
    
    % Read the original mesh
    pc_data = pcread(input_mesh);
    vertices = pc_data.Location;
        
    % Record readings from the origin of the particle
    center = mean(vertices,1);
    vertices = vertices - center(ones(size(vertices,1),1),:);
    %% spherical expansion and reconstructon
    clear sph_verts faces Z;
    [sph_verts, faces] = icosphere(reconstruction_mesh_refinment);
    reconDeg = min(reconDeg, maxDeg);
    
    fvec = SH_expansion(vertices,maxDeg); % Rotational spherical expansion
    spharm_verts=SH_reconstruction(sph_verts,reconDeg,fvec);
    if plot_figures
        paraview_patch(spharm_verts, faces);
    end
    % Save the reconstruction results
    stlWrite(reconstruction_output_mesh, faces, spharm_verts, 'mode','ascii')
end
disp("Done.")
fprintf('Elapsed time is %1.2f min.\n', toc/60)% Export only the last object stone.stl