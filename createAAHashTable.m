function [h] = createAAHashTable ()
% function creates an associative array for the amino acids to translate
% them from 1 to 3 letter code and vice versa
%
% params:
%       h: the hashtable containing the amino acids
%
%   example for call:
%       h = createAAHashTable;
%       h.get('A')
%       h.get('ALA')
%
%  B.Knapp 2011-05-04


h =  java.util.Hashtable;

h.put('ALA', 'A');
h.put('ARG', 'R');
h.put('ASN', 'N');
h.put('ASP', 'D');
h.put('CYS', 'C');
h.put('GLN', 'Q');
h.put('GLU', 'E');
h.put('GLY', 'G');
h.put('HIS', 'H');
h.put('ILE', 'I');
h.put('LEU', 'L');
h.put('LYS', 'K');
h.put('MET', 'M');
h.put('PHE', 'F');
h.put('PRO', 'P');
h.put('SER', 'S');
h.put('THR', 'T');
h.put('TRP', 'W');
h.put('TYR', 'Y');
h.put('VAL', 'V');

h.put('A', 'ALA');
h.put('R', 'ARG');
h.put('N', 'ASN');
h.put('D', 'ASP');
h.put('C', 'CYS');
h.put('Q', 'GLN');
h.put('E', 'GLU');
h.put('G', 'GLY');
h.put('H', 'HIS');
h.put('I', 'ILE');
h.put('L', 'LEU');
h.put('K', 'LYS');
h.put('M', 'MET');
h.put('F', 'PHE');
h.put('P', 'PRO');
h.put('S', 'SER');
h.put('T', 'THR');
h.put('W', 'TRP');
h.put('Y', 'TYR');
h.put('V', 'VAL');

