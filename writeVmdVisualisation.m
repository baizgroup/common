function [] = writeVmdVisualisation (fileName, points, radiusLine, radiusSphere, colorCode)
% This function writes coordinates as textfile to visualize them in VMD.
% The content of the textfile can then be copy&pasted into the VMD / extensions / TK console.
%
% input:
% - fileName: the file where the vmd script will be written to
% -points is a nx3 array containing the x,y,z coordinate of n points in
% space. At each point a sphere will be drawn and every 2 points are
% connected with a line (cylinder).
% - radiusLine: the radius of the line used to connect the points
% - radiusSphere: the radius of the sphere drawn at each point
% - colorCode: the color code for the drawn lines and spheres according to
% vmd (e.g. 0=blue, 1=red, 2=gray, 3=orange, ...)
%
%  B. Knapp 2014-01-15

if nargin < 3
    radiusLine = 0.1;
end
if nargin < 4
    radiusSphere = 0.9;
end
if nargin < 5
    colorCode = 1;
end

points = nanometer2angstroem(points); % we use nm while vmd uses angstroem
fid = fopen(fileName,'w');

if fid == -1
    error ('Could not open file "%s" for writing.', fileName);
end

fprintf(fid, 'mol load graphics %s\n', 'distances');
fprintf(fid, 'draw color %1.0f\n', colorCode);

resolutionCurve = 15;
filledCurve = 'yes';

for pointIter = 1:2:size(points,1)-1
    fprintf(fid, 'graphics top cylinder {%1.4f %1.4f %1.4f} {%1.4f %1.4f %1.4f} radius %1.1f resolution %1.0f filled %s\n',...
        points(pointIter,1), points(pointIter,2), points(pointIter,3), points(pointIter+1, 1), points(pointIter+1, 2), points(pointIter+1, 3),...
        radiusLine, resolutionCurve, filledCurve);
end

fprintf(fid, 'mol load graphics %s\n', 'points');
fprintf(fid, 'draw color %1.0f\n', colorCode);

for pointIter = 1:size(points,1)
    fprintf(fid, 'graphics top sphere {%1.4f %1.4f %1.4f} radius %1.1f resolution %1.0f\n',...
        points(pointIter,1), points(pointIter,2), points(pointIter,3), radiusSphere, resolutionCurve);
end

fprintf(fid, 'display projection Orthographic\n');
fprintf(fid, 'display depthcue off\n');
fprintf(fid, 'axes location Off\n');
fprintf(fid, 'display resetview\n\n');

fclose(fid);
