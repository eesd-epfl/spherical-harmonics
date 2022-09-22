function spharm_verts = SH_reconstruction(sph_verts, recDeg, fvec)
[phi,theta] = cart2sph(sph_verts(:,1), sph_verts(:,2), sph_verts(:,3));
Z = SH_basis(sph_verts, recDeg); 
r = real(Z(:, 1: (recDeg+1)^2)*fvec(1: (recDeg+1)^2, :));
[rx, ry, rz] = sph2cart(phi, theta, r);
spharm_verts = [rx, ry, rz];