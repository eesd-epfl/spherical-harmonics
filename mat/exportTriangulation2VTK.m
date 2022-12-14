function exportTriangulation2VTK(file,XYZ,tri,k1,k2,km,kg)
%
% This function takes as input a 2D unrestricted triangulation and export
% it to an ASCII VTK file which can be oppened with the viewer Paraview.
%
% Input :
%           "file" is the name without extension of the file (string).
%           "XYZ" is the coordinate of the vertex of the triangulation (nx3 matrix).
%           "tri" is the list of triangles which contain indexes of XYZ (mx3 matrix).
%
% Sample example :
%
%   [X,Y,Z]=peaks(25);
%   X=reshape(X,[],1);
%   Y=reshape(Y,[],1);
%   Z=0.4*reshape(Z,[],1);
%   tri = delaunay(X,Y);
%   exportTriangulation2VTK('sampleExampleTri',[X Y Z],tri)
%
% Note : If the triangulation doesn't have Z component (a plane), put the
% third column of XYZ with all zeros. Paraview only deals with 3D object.
%
% David Gingras, January 2009


X=XYZ(:,1);
Y=XYZ(:,2);
Z=XYZ(:,3);
nbpoint=length(X);
if mod(nbpoint,3)==1
    X(nbpoint+1:nbpoint+2,1)=[0;0];
    Y(nbpoint+1:nbpoint+2,1)=[0;0];
    Z(nbpoint+1:nbpoint+2,1)=[0;0];
	k1(nbpoint+1:nbpoint+2,1)=[0;0];
	k2(nbpoint+1:nbpoint+2,1)=[0;0];
	km(nbpoint+1:nbpoint+2,1)=[0;0];
	kg(nbpoint+1:nbpoint+2,1)=[0;0];	
elseif mod(nbpoint,3)==2
    X(nbpoint+1,1)=0;
    Y(nbpoint+1,1)=0;
    Z(nbpoint+1,1)=0;
	k1(nbpoint+1,1)=0;
	k2(nbpoint+1,1)=0;
	km(nbpoint+1,1)=0;
	kg(nbpoint+1,1)=0;
end
nbpoint=length(X);

fid=fopen([file '.vtk'],'wt');
fprintf(fid,'# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\n');
fprintf(fid,'POINTS %d float\n',nbpoint);
fprintf(fid,'%3.7f %3.7f %3.7f %3.7f %3.7f %3.7f %3.7f %3.7f %3.7f\n',[X(1:3:end-2) Y(1:3:end-2) Z(1:3:end-2) X(2:3:end-1) Y(2:3:end-1) Z(2:3:end-1) X(3:3:end) Y(3:3:end) Z(3:3:end)]');
ntri=length(tri);
fprintf(fid,'POLYGONS %d %d\n',ntri,4*ntri);

fprintf(fid,'3 %d %d %d\n',(tri-ones(ntri,3))');

%append some scalar data
fprintf(fid, 'POINT_DATA   %d\n', nbpoint); %ASCII header
fprintf(fid, 'SCALARS k1 float\n'); %ASCII header
fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
fprintf (fid, '%3.7f \n', k1); %binary data

%append some scalar data
%fprintf(fid, 'POINT_DATA   %d\n', size(k2,1)); %ASCII header
fprintf(fid, 'SCALARS k2 float\n'); %ASCII header
fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
fprintf (fid, '%3.7f \n', k2); %binary data

%append some scalar data
%fprintf(fid, 'POINT_DATA   %d\n', size(km,1)); %ASCII header
fprintf(fid, 'SCALARS km float\n'); %ASCII header
fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
fprintf (fid, '%3.7f \n', km); %binary data

%append some scalar data
%fprintf(fid, 'POINT_DATA   %d\n', size(kg,1)); %ASCII header
fprintf(fid, 'SCALARS kg float\n'); %ASCII header
fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
fprintf (fid, '%3.7f \n', kg); %binary data

fclose(fid);
end