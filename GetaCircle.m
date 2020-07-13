function [xyz] = GetaCircle(Center,normal,R,N)
% Given a center point, a normal vector and a radius, output the
% coordinates of the points and plot the circle
% CENTER is the center point as an array [x,y,z]
% NORMAL is the normal vector array as [nx,ny,nz]
% R is the radius value
% N is thr number of points that you want

theta2 = 0:2*pi/N:(2*pi)-(2*pi/N);
v=null(normal);
%points=repmat(Center',1,size(theta,2))+R*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
%plot3(points(1,:),points(2,:),points(3,:),'o','MarkerFaceColor', 'r'); hold on
xyz = repmat(Center',1,size(theta2,2))+R*(v(:,1)*cos(theta2)+v(:,2)*sin(theta2));
end