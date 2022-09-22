function parametrize_surface(input_file, output_file)
% Demo for mapping a stone surface onto a sphere conformally
%
% Please refer to the following papers:
%
% spherical_conformal_map:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% mobius_area_correction_spherical:
% [2] G. P. T. Choi, Y. Leung-Liu, X. Gu, and L. M. Lui, 
%     "Parallelizable global conformal parameterization of simply-connected surfaces via partial welding."
%     SIAM Journal on Imaging Sciences, 2020.
close all; %clc;

if nargin == 0
    % Just change the name of the file here .. and run this script directly
    input_file = 'test_amir.stl';
    output_file = strcat(pwd, '/test_stone_output.stl');
elseif nargin == 1
    output_file = strcat(pwd, '/test_stone_output.stl');
end

addpath('LinearSphericalConformalMap_stone')
addpath('LinearSphericalConformalMap_stone/mfile')
addpath('LinearSphericalConformalMap_stone/extension')
addpath('LinearSphericalConformalMap_stone/stlTools')

%% Read the input file
disp (input_file)
[v,f] = stlRead(input_file);
% plot_mesh(v,f);

%% spherical conformal map [1]
% map = spherical_conformal_map(v,f);
% 
% plot_mesh(map,f); 
% 
% % evaluate the angle distortion
% angle_distortion(v,f,map);

%% Spherical conformal map [1] combined with the method in our latest work [2] 
% the Mobius area correction method further reduces the area distortion of the spherical parameterization while preserving the conformality
map = spherical_conformal_map(v,f);
map2 = mobius_area_correction_spherical(v,f,map);
stlWrite(output_file, f, map2, 'mode','ascii')

% plot_mesh(map2,f);
% evaluate the angle distortion
%angle_distortion(v,f,map2);
end