function [pdb] = constructorPdb ()
% this constructor creates an empty structure to store pdb data according
% to the official nomenclature
%
% COLUMNS        DATA  TYPE    FIELD        DEFINITION
% -------------------------------------------------------------------------------------
% 1 -  6         Record name   "ATOM  "
% 7 - 11         Integer       serial       Atom  serial number.
% 13 - 16        Atom          name         Atom name.
% 17             Character     altLoc       Alternate location indicator.
% 18 - 20        Residue name  resName      Residue name.
% 22             Character     chainID      Chain identifier.
% 23 - 26        Integer       resSeq       Residue sequence number.
% 27             AChar         iCode        Code for insertion of residues.
% 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
% 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
% 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
% 55 - 60        Real(6.2)     occupancy    Occupancy.
% 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
% 77 - 78        LString(2)    element      Element symbol, right-justified.
% 79 - 80        LString(2)    charge       Charge  on the atom.
% 
% http://www.wwpdb.org/documentation/format33/sect9.html
% 
%   B. Knapp 2011-05-04




% example for a pdb line in the flatfile:
% ATOM   2853  CA  HIS B 163     157.846  20.965 -32.973  1.00 70.57           C 


% array of the size [n*3*t]: where n is the number of atoms; 3 means x, y and z; t is the n*3 information for several timesteps
pdb.coords = [];

pdb.coordsUnit = 'nm';                  % needs to multiply original pdb coords by 10 to be in the same unit as xtc and trr. nanometer instead of angstrom
pdb.recordName = {};                    % e.g. "ATOM"
pdb.atomNumber = [];                    % e.g. "2853"
pdb.atomName   = {};                    % e.g. "CA"
pdb.alternateLocationIndicator = {};    % e.g. ???
pdb.residueName = {};                   % e.g. "HIS"
pdb.chainIdentifier = {};               % e.g. "B"
pdb.residueNumber = [];                 % e.g. "163"
pdb.insertionCode = {};                 % e.g. ???
pdb.occupancy = [];                     % e.g. "1.00"
pdb.temperatureFactor = [];             % e.g. "70.57"
pdb.element = {};                       % e.g. "C"
pdb.atomCharge = [];                    % e.g. ???
pdb.information = {};
pdb.name = '';                          % the file name or the pdb id where the data comes from (field for free use)

