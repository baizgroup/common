function [newData] = nanometer2angstroem (myData) 
% This function converts namometer to angstroem (e.g. make an xtc file compareable with a pdb file)
%
% example for usage: 
%
% nanometer2angstroem(0.3)  
% nanometer2angstroem([0.5 0.8 0.2])  
%
% B. Knapp 2013-10-10

if ~isnumeric(myData)
    error ('You tried to convert a non-number.');
end
if isnan(myData)
    warning('There are nan-values in your data');
end
   
newData=myData.*10; % to angststroem is a multiplication by 10
