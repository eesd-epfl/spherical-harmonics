%% This function is used for calculating the spherical harmonic coefficiencs

function  fvec = SH_expansion(vertices,maxDeg,~)
vertnum = size(vertices,1);
max_d = maxDeg;

deg = max(1, floor(sqrt(vertnum)*1/2));
deg = min(deg, max_d);
fprintf('\nUse spharm up to %d degree (vec_len=%d).\n',deg,(deg+1)^2);
center = mean(vertices,1);
rads = zeros(length(vertices(:, 1)), 1);

for i = 1:size(vertices,1)
    vertices(i,:) = vertices(i,:)-center; % the polar radius from the center
    rads(i,:) = norm(vertices(i,:));
end
% calculate the spherical basis
Z = SH_basis(vertices, deg); %% vertices is the direct mapping

[x,y] = size(Z);
fprintf('\nLeast square for %d equations and %d unknowns\n',x,y);

% Least square fitting (Both give identical results)
fvec = Z\rads;         % The least square using the pesudo inverse (faster in Matlab)
% fvec = inv(Z'*Z)*Z'*rads; % The least square equation (slower)