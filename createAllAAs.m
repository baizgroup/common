function [AAsCell, AAsArr] = createAllAAs ()
% Returns a cell array as well an array containing the 20 amino acids in
% the 1 letter code. You might want to use createAAHashTable.m to convert
% them to 3 letter codes.
%
% example for call: [AAsCell, AAsArr] = createAllAAs;
%
% B. Knapp 2011-05-04 


AAsArr = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'];
AAsCell = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
