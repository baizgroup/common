function [newData] = angstroem2nanometer(myData) 
% This function converts angrom to namometer (e.g. make a pdb file compareable with an xtc file)
%
% example for usage: 
%
% angstroem2nanometer(3)  
% angstroem2nanometer([5 8 2])  
%
% B. Knapp 2013-10-10

%%THIS CHECK SLOWS DOWN THE CODE
if ~isnumeric(myData)
    error ('You tried to convert a non-number.');
end
if isnan(myData)
    warning('There are nan-values in your data');
end

newData=myData./10; % to nanometer is a division by 10
