%% The main funtion to implement the sphercial harmonic analysis
% by Zhou Bo on 2014/06/25
clc; clear;
%% The input of particle Number and result dictionary
Pid = 'LBS01';   % input "LBS02" or "HDG01"
% Input the Data after spherical paremetrazation
inputname = [Pid, '_CALD_ini.mat'];
inputdir = fullfile(cd, 'Surf_input\Parameterize1Part_Origin');
load([inputdir '\' inputname]);

%[ini_verts,Volume,Surf]=Cal_Geometry(faces,vertices);
subplot(1,2,1); 
patch_lighta(vertices, faces);
%axis([-1.5,1.5,-1.5,1.5,-1.5,1.5]); 
break;
clear sph_verts faces Z; load('L4_icosa.mat');

NormRot_verts = Surface_reconstruction(Pid,faces,sph_verts,ini_verts);   % spherical reconstruction
subplot(1,2,2); 
patch_lightmesh(NormRot_verts, faces);
axis([-1.5,1.5,-1.5,1.5,-1.5,1.5]); 
break;