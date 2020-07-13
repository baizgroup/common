function [pdb] = writePdb (pdb, pdbFileName)
% function writes an object which matches "constructorPdb" to the
% speciefied filename in pdb-format.
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
%
%
% B. Knapp 2013-07-04

% that would haven been too easy...
% matlabPdb = ourPdb2matlabPdb(pdbData, false);
% pdbwrite(pdbFileName, matlabPdb);
%
%
% example for useage: 
% trjWildType = parseTrj('/run/media/knapp/Seagate Expansion Drive/1mi5_172_100ns_sims/wildtype/7Y.firstFrame.pdb', '/run/media/knapp/Seagate Expansion Drive/1mi5_172_100ns_sims/wildtype/7Y.final.md.xtc');
% trjWildType.coords(:,:, 4:end) = []; 
% writePdb(trjWildType, 'test.multipe.pdb');

if ~isfield(pdb, 'coords')
    error ('The first argument is not a valid pdb structure because the "coords" field is missing.');
end

[nrOfCoords, myDims, nrOfFrames] = size(pdb.coords);

% these values might not be set (e.g. in models)
if isfield(pdb, 'temperatureFactor')
    pdb.temperatureFactor(isnan(pdb.temperatureFactor)) = 0;
else
    pdb.temperatureFactor = [];
    pdb.temperatureFactor(1:nrOfCoords) = 0;
end
if isfield(pdb, 'atomCharge')
    pdb.atomCharge(isnan(pdb.atomCharge)) = 0;
else
    pdb.atomCharge = [];
    pdb.atomCharge(1:nrOfCoords) = 0;
end




fid = fopen(pdbFileName,'w');            % open the file
if fid ~= -1
    for frameIter = 1:nrOfFrames
        if nrOfFrames > 1
            fprintf(fid, 'MODEL %i\n', frameIter);
        end
        for atomIter = 1:nrOfCoords
            
            line = sprintf('%6s%5.0f %4s%1s%3s %1s%4.0f%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2d\n',...
                'ATOM  ', ...
                pdb.atomNumber(atomIter), ...
                pdb.atomName{atomIter}, ...
                ' ', ...
                pdb.residueName{atomIter}, ...
                pdb.chainIdentifier{atomIter}, ...
                pdb.residueNumber(atomIter), ...
                pdb.insertionCode{atomIter}, ...
                nanometer2angstroem(pdb.coords(atomIter, 1, frameIter)), ...
                nanometer2angstroem(pdb.coords(atomIter, 2, frameIter)), ...
                nanometer2angstroem(pdb.coords(atomIter, 3, frameIter)), ...
                pdb.occupancy(atomIter), ...
                pdb.temperatureFactor(atomIter), ...
                pdb.element{atomIter}, ...
                pdb.atomCharge(atomIter));
            fprintf(fid, '%s', line);
            
        end % atomIter
        fprintf(fid, 'TER\n');
        if nrOfFrames > 1
            fprintf(fid, 'ENDMDL\n');
        end
    end % frameIter
    fclose(fid);                     %# Close the file
else
    fprintf ('Invalid file name "%s"', pdbFileName);
end % fid




