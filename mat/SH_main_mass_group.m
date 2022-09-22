%% The main funtion to implement the sphercial harmonic analysis
clc; clear; close all; tic
addpath('LinearSphericalConformalMap_stone')
addpath('LinearSphericalConformalMap_stone/mfile')
addpath('LinearSphericalConformalMap_stone/extension')
addpath('LinearSphericalConformalMap_stone/stlTools')

%% The inputs
parametrize_input = false;
maxDeg = 20 ; % the degree for spherical basis calculation
reconDeg = 20; % tuntitled2he degree used for spherical harmonic reconstrution
reconstruction_mesh_refinment = 4; % Possible inputs 1,2,3,4,5,6 for icospheres
plot_figures = false;
% input_stl_files = '/home/mahmoudshaqfa/3d-microstructure-generator/3d-walls-generator/Modular Programming/Saved Layers/STLMESH/';
% input_stl_files = '/home/mahmoudshaqfa/3d-microstructure-generator/3d-walls-generator/Modular Programming/SphericalHarmonics/SH_finalized/SH_main/Amir_stones_analysis/STLMESH';
% input_stl_files = '/home/mahmoudshaqfa/3d-microstructure-generator/3d-walls-generator/Modular Programming/SphericalHarmonics/SH_finalized/SH_main/Amir_stones_simplified_stones/STLMESH';
% input_stl_files = '/home/mahmoudshaqfa/3d-microstructure-generator/3d-walls-generator/Modular Programming/SphericalHarmonics/SH_finalized/SH_main/test_geometry';
input_stl_files = '/home/mahmoudshaqfa/3d-microstructure-generator/3d-walls-generator/Modular Programming/SphericalHarmonics/SH_finalized/SH_main/Amir_extra_simplified';

% Parametrized meshes if it doesn't exist use the default directory
% input_para_stl = "";
% input_para_stl = '/home/mahmoudshaqfa/3d-microstructure-generator/3d-walls-generator/Modular Programming/SphericalHarmonics/SH_finalized/SH_main/Amir_stones_analysis/parametrized/';
input_para_stl = '/home/mahmoudshaqfa/3d-microstructure-generator/3d-walls-generator/Modular Programming/SphericalHarmonics/SH_finalized/SH_main/Amir_stones_simplified_stones/parametrized/';

%% Collect the files list and organise inputs/outputs
fileList = dir(fullfile(input_stl_files, '*.stl'));
fileList = struct2cell(fileList);
n = size(fileList, 2);

Dl = zeros([n, maxDeg+1]);

for stl_idx=1:n
    % Surface parametrization
    input_mesh_name = fileList{1, stl_idx};
    input_mesh = strcat(input_stl_files, '/', input_mesh_name);
    if (input_para_stl == "")
        output_mesh = strcat(pwd, '/parametrized/', input_mesh_name);
    else
        output_mesh = strcat(input_para_stl, input_mesh_name);
    end
    reconstruction_output_mesh = strcat(pwd, '/', 'reconstructed/', 'rec', input_mesh_name);
    if parametrize_input
       parametrize_surface(input_mesh, output_mesh)
    end
    disp(strcat('Analysing stone: ', input_mesh_name))
    
    % Read the original mesh
    [vertices, faces] = stlRead(input_mesh);
    % Record readings from the origin of the particle
    center = mean(vertices,1);
    vertices = vertices - center(ones(size(vertices,1),1),:);
    [mesh_volume, surface_area] = stlVolume(vertices,faces);
    
    % Read the parametrized input sphere
%     [para_v, para_f] = stlRead(output_mesh);
%     spharm_verts = para_v;
%     center = mean(spharm_verts,1);
%     spharm_verts = spharm_verts - center(ones(size(spharm_verts,1),1),:);
    
    % Plot the original input mesh - to compare
    if plot_figures
        figure(1)
        subplot(1,2,1)
        paraview_patch(vertices, faces);
    end
    %% spherical expansion and reconstructon
    clear sph_verts faces Z;
    [sph_verts, faces] = icosphere(reconstruction_mesh_refinment);
    reconDeg = min(reconDeg, maxDeg);
    
    fvec = SH_expansion(vertices,maxDeg); % Rotational spherical expansion
    
    % Shape descriptors (2-norm) for frequency accumulates in a certain
    % frequency degree..
%     for l = 1:maxDeg
%         for m = -l:1:l
%             Dl(stl_idx, l) = Dl(stl_idx, l) + (real(fvec(l^2 + l + m + 1)))^2 + (imag(fvec(l^2 + l + m + 1)))^2;
%         end
%         Dl(stl_idx, l) = sqrt(Dl(stl_idx, l));
%     end
%     Dl(stl_idx, :) = Dl(stl_idx, :)/Dl(stl_idx, 1);
    
    spharm_verts=SH_reconstruction(sph_verts,reconDeg,fvec);
    if plot_figures
        subplot(1,2,2)
        paraview_patch(spharm_verts, faces);
    end
    % Save the reconstruction results
    stlWrite(reconstruction_output_mesh, faces, spharm_verts, 'mode','ascii')
    if plot_figures
        figure(2)
        subplot(1,2,1)
        bar(fvec)
        subplot(1,2,2)
        bar(Dl(stl_idx, :))
        hold on
        plot_figures(1:length(Dl(stl_idx, :)), Dl(stl_idx, :), 'k', 'LineWidth', 2)
        hold off
    end
    %% Spherical harmonics toplogy plot
    x = spharm_verts(:, 1);
    y = spharm_verts(:, 2);
    z = spharm_verts(:, 3);
    [phi, theta, r] = cart2sph(spharm_verts(:,1), spharm_verts(:,2), spharm_verts(:,3));
    theta = pi/2 - theta; neg_angles = find(phi<0); phi(neg_angles) = phi(neg_angles) + 2*pi;
    
    if plot_figures
        figure(3)
        % tri = delaunay(phi,theta);
        % triplot(tri,phi,theta)
        % trisurf(tri,phi,theta, sqrt(z.^2 + x.^2 + y.^2), 'LineWidth', 0.01)
        % view(0,90)
        % scatter3(x,y,z, 'filled')
        
        tri = delaunay(x,y,z);
        h = trisurf(tri, x, y, z);
        axis vis3d
        axis off
        l = light('Position',[-50 -15 29]);
        % set(gca,'CameraPosition',[208 -50 7687])
        % lighting phong
        % shading interp
        colorbar EastOutside
    end
end
disp("Done.")

% figure
figure('Renderer', 'painters', 'Position', [10 10 900 600])
subplot(1,2,1)
data_std = std(Dl, 1);
data_mean = mean(Dl, 1);
x = 1:length(data_mean);
loglog(x(3:end), data_mean(3:end), '-ks', 'MarkerSize',7,...
    'MarkerFaceColor',[1 .6 .6], 'LineWidth', 2)
if n > 1
    hold on
    loglog(x(3:end), data_mean(3:end) + data_std(3:end), '--b', 'LineWidth', 2)
    hold on
    loglog(x(3:end), data_mean(3:end) - data_std(3:end), '--r', 'LineWidth', 2)
end
hold on
% fitness = fit(x(3:end)', data_mean(3:end)','b*x^m');
% plot(fitness,'k')
m = -1.462;
b = 12.66;
r2 = 0.9624;
loglog(x(3:end), b * x(3:end).^m)

fit_text = sprintf('$D_{n} = %f~ n ^ {%f}$', b, m);
text(mean(x(3:end))/3, mean(data_mean(3:end))*2.5, fit_text, 'Interpreter',...
    'Latex', 'FontSize', 13, 'Color', 'k')

fit_text = sprintf('$R^{2} = %f $', r2);

text(mean(x(3:end))/3, mean(data_mean(3:end))*2.5, fit_text, 'Interpreter',...
    'Latex', 'FontSize', 13, 'Color', 'k')
fit_text2 = sprintf('$ Fractal~dimension = %f$', (6 + m) * 0.5);
text(mean(x(3:end))/2, mean(data_mean(3:end))*3.5, fit_text2, 'Interpreter',...
    'Latex', 'FontSize', 12, 'Color', 'k')
hold off
xlim([3, length(x)])
ylim([0,10])
xlabel('Frequency degree (n)')
ylabel('Normalized amplitude')
legend('Mean','Upper bound', 'Lower bound', 'Curve fitting')

subplot(1,2,2)
for figs =  1:size(Dl, 1)
    loglog(x, Dl(figs, :), 'Color', rand(1,3), 'LineWidth', 1.5)
    hold on
end
hold off
xlim([3, length(x)])
ylim([0,10])
xlabel('Frequency degree (n)')
ylabel('Normalized amplitude')
fprintf('Elapsed time is %1.2f min.\n', toc/60)% Export only the last object stone.stl
% save('reconstruction_data.mat', 'spharm_verts');
% save('coefficients.mat', 'Dl');