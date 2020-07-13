function [pdb] = parsePdb (fileName, doPreallocation)
% function parses a *.pdb file according to the below mentioned
% official pdb-nomenclature.
% the results are stored as pdb struct which is created by the
% function "constructorPdb".
%
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
%
% examples for usage:
% pdb = parsePdb ('/data/housemartin/knapp/projects/dien_rmsf/testFiles/testsForHung/gro2mat_package/functions_andTests/pdbStructure.pdb');
% plot3(pdb.coords(:,1), pdb.coords(:,2), pdb.coords(:,3), 'rx');
% plotPdb(pdb, 'coloringMethod', 'chain', 'drawingMethod', 'lines', 'plotOnlyCA', true, 'plotHetAtm', false, 'turnLightOn', false)
%
%
% Bernhard Knapp 2011-08-02


if nargin == 1
    doPreallocation = true;
end

pdb = constructorPdb();
data = readFile(fileName);

% preallocate arrays
if doPreallocation
    nrOfAtomEntries = 0;
    for dataIter = 1:numel(data)
        line = data{dataIter};
        if ((length(line)>=4) && (strcmp(line(1:4), 'ATOM') || strcmp(line(1:6), 'HETATM')))
            nrOfAtomEntries=nrOfAtomEntries+1;
        end
    end
    pdb.coordsUnit = zeros(nrOfAtomEntries,3,1);
    pdb.recordName = cell(nrOfAtomEntries,1);
    pdb.atomNumber = zeros(nrOfAtomEntries,1);
    pdb.atomName   = cell(nrOfAtomEntries,1);
    pdb.alternateLocationIndicator = cell(nrOfAtomEntries,1);
    pdb.residueName = cell(nrOfAtomEntries,1);
    pdb.chainIdentifier = cell(nrOfAtomEntries,1);
    pdb.residueNumber = zeros(nrOfAtomEntries,1);
    pdb.insertionCode = cell(nrOfAtomEntries,1);
    pdb.occupancy = zeros(nrOfAtomEntries,1);
    pdb.temperatureFactor = zeros(nrOfAtomEntries,1);
    pdb.element = cell(nrOfAtomEntries,1);
    pdb.atomCharge = zeros(nrOfAtomEntries,1);
    pdb.information = {};
    pdb.coords = zeros(nrOfAtomEntries,3);
end

% do the parsing stuff
pos = 0;
for dataIter = 1:numel(data)
    line = data{dataIter};
    if ((length(line)>=4) && (strcmp(line(1:4), 'ATOM') || strcmp(line(1:6), 'HETATM')))  % atom line found
        pos = pos + 1;
        
        pdb.coords(pos, 1, 1)       = angstroem2nanometer(str2double(line(31:38))); % change unit to nanometer!
        pdb.coords(pos, 2, 1)       = angstroem2nanometer(str2double(line(39:46)));
        pdb.coords(pos, 3, 1)       = angstroem2nanometer(str2double(line(47:54)));
        pdb.coordsUnit              = 'nm';
        pdb.recordName{pos, 1}         = strtrim(line(1:6));
        pdb.atomNumber(pos, 1)         = str2double(line(7:11));
        pdb.atomName{pos, 1}           = strtrim(line(13:16));
        pdb.alternateLocationIndicator{pos, 1} = line(17);
        pdb.residueName{pos, 1}        = strtrim(line(18:20));
        pdb.chainIdentifier{pos, 1}    = line(22);
        pdb.residueNumber(pos, 1)      = str2double(line(23:26));
        pdb.insertionCode{pos, 1}      = line(27);
        pdb.occupancy(pos, 1)          = str2double(line(55:60));
        pdb.temperatureFactor(pos, 1)  = str2double(line(61:66));
        if length(line) >= 67
            pdb.element{pos, 1}            = strtrim(line(77:78));
            if length(line) > 79
                pdb.atomCharge(pos, 1)         = str2double(line(79:80));
            end
        else
            pdb.element{pos, 1}            = '';
            pdb.atomCharge(pos, 1)         = 0;
        end
    else
        pdb.information{end+1} = line;
    end
end

pdb.sourceFileName = fileName;
