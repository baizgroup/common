function [d] = calcDist (a,b, doPlot)
% function calculates the distance between two points in 3d space and
% optionally plots them
%
% params:
% - a: the first point
% - b: the second point
% - d: the distance calculated
%
% example for usage:
% a = rand (3,1);
% b = rand (3,1);
% d = calcDist(a, b, true); title (sprintf ('dist=%2.3f', d));
%
%
% B.Knapp 2013-06-05

if nargin == 2
    doPlot = false;
end

d = sqrt((a(1)-b(1))^2 + (a(2)-b(2))^2 + (a(3)-b(3))^2);

if doPlot
    hold on;
    plot3(a(1), a(2), a(3), 'ro');
    plot3(b(1), b(2), b(3), 'ro');
    plot3([a(1),b(1)], [a(2),b(2)], [a(3),b(3)], 'r--');
    hold off;
end