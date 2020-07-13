function [trj] = convertTrj(fnXtc)
% This function combines the output of parsePdb with parseXtc. It parses
% the atom information from a pdb file and adds the trajectory information
% from an xtc file.
%

xtc = parseXtc(fnXtc);
