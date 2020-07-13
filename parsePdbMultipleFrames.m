function [pdb] = parsePdbMultipleFrames (fileName, doPreallocation)
% function parses a multiple *.pdb file (for example created from an xtc
% file via the trjconv command of Gromacs) according to the below mentioned
% official pdb-nomenclature.
%
% Please note that the input file needs to contain lines determining
% when the frame is done: "ENDMDL" at the end of each frame
%
%
% The results are stored as pdb struct which is created by the
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
% pdb = parsePdbMultipleFrames ('/data/housemartin/knapp/projects/dien_rmsf/testFiles/testsForHung/4_peptide/peptide_model1.md.500ns.fitRotTrans.pbcCluster.100ps.xtc.pdb');
% pdb = parsePdbMultipleFrames ('/data/housemartin/knapp/projects/dien_rmsf/testFiles/testsForHung/4_peptide/peptide_model1.md.500ns.fitRotTrans.pbcCluster.1ns.xtc.pdb');
% pdb = parsePdbMultipleFrames ('/data/housemartin/knapp/projects/dien_rmsf/testFiles/testsForHung/4_peptide/peptide_model1.md.500ns.fitRotTrans.pbcCluster.xtc.pdb');
% figure; hold on;
% plot3(pdb.coords(:,1,1), pdb.coords(:,2,1), pdb.coords(:,3,1), 'bx');
% plot3(pdb.coords(:,1,5), pdb.coords(:,2,5), pdb.coords(:,3,5), 'rx');
% title('frame 1 and frame 5');
%
% B. Knapp 2014-01-14

if nargin == 1
    doPreallocation = true;
end

pdb = constructorPdb();
data = readFile(fileName);
secondFrameFound = false;
line = '';
%  check for the number of atoms and preallocate arrays
sizeFirstFrame = -1;
if doPreallocation
    nrOfAtomEntries = 0;
    dataIter = 1;
    while dataIter < numel(data) && not(secondFrameFound)
        line = data{dataIter};
        if length(line) >= 6
            secondFrameFound = strcmp(line(1:6), 'ENDMDL');
            if secondFrameFound
                sizeFirstFrame = dataIter;
            end
        end
        if ((length(line)>=4) && (strcmp(line(1:4), 'ATOM') || strcmp(line(1:6), 'HETATM')))
            nrOfAtomEntries=nrOfAtomEntries+1;
        end
        dataIter = dataIter + 1;
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
end

% do the actual atom parsing
atomPos = 0;
nrFramesParsed = 0;
for dataIter = 1:numel(data)
    line = data{dataIter};
    if ((length(line)>=4) && (strcmp(line(1:4), 'ATOM') || strcmp(line(1:6), 'HETATM')))  % atom line found
        atomPos = atomPos + 1;
        
        
        pdb.coords(atomPos, 1, nrFramesParsed+1)       = angstroem2nanometer(str2double(line(31:38))); % change unit to nanometer!
        pdb.coords(atomPos, 2, nrFramesParsed+1)       = angstroem2nanometer(str2double(line(39:46)));
        pdb.coords(atomPos, 3, nrFramesParsed+1)       = angstroem2nanometer(str2double(line(47:54)));
        
        if nrFramesParsed > 0 % we need the meta information only once
            pdb.coordsUnit              = 'nm';
            pdb.recordName{atomPos, 1}         = strtrim(line(1:6));
            pdb.atomNumber(atomPos, 1)         = str2double(line(7:11));
            pdb.atomName{atomPos, 1}           = strtrim(line(13:16));
            pdb.alternateLocationIndicator{atomPos, 1} = line(17);
            pdb.residueName{atomPos, 1}        = strtrim(line(18:20));
            pdb.chainIdentifier{atomPos, 1}    = line(22);
            pdb.residueNumber(atomPos, 1)      = str2double(line(23:26));
            pdb.insertionCode{atomPos, 1}      = line(27);
            pdb.occupancy(atomPos, 1)          = str2double(line(55:60));
            pdb.temperatureFactor(atomPos, 1)  = str2double(line(61:66));
            if length(line) >= 67
                pdb.element{atomPos, 1}            = strtrim(line(77:78));
                if length(line) > 79
                    pdb.atomCharge(atomPos, 1)         = str2double(line(79:80));
                end
            else
                pdb.element{atomPos, 1}            = '';
                pdb.atomCharge(atomPos, 1)         = 0;
            end
        end % if nrFramesParsed > 0
    else
        pdb.information{end+1} = line;
    end % length line and ATOM or HETATM
    
    if length(line) >= 6 && strcmp(line(1:6), 'ENDMDL'); % end of one frame found
        nrFramesParsed = nrFramesParsed + 1;
        atomPos = 0;
    end
end

pdb.sourceFileName = fileName;
