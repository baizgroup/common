function [pdbOut] = selectAAs (pdbIn, chainIDs, aaStarts, aaEnds)
% This functions select specific residues from a molecules given in the
% matrix pdb format (created by constructorPdb.m).
%
% params:
% - pdbIn: a structure created by constructorPdb
% - chainIDs: chain of the complex, where the selected part of the molecule is
% located (e.g. chainID='A'). Also multiple chains can be given. Their
% number must match the starts and ends.
% - aaStarts: Number of amino acid, where the specific part starts (e.g.
% aaStart=46). Also multipe starts can be given. Their number must match
% the chains and ends.
% - aaEnds: Number of amino acid, where the specific part ends (e.g.
% aaEnd=77).  Also multipe ends can be given. Their number must match
% the chains and starts.
% - pdbOut: the same as pdbIn but only the selected part
%
%
% Example for usage:
% pdb = parsePdb ('pdbStructure.pdb');
%
% pdb_C = selectionPart(pdb, 'C', 1, 9);
% plotPdb(pdb_C, 'coloringMethod', 'residue', 'drawingMethod', 'o', 'plotOnlyCA', false, 'plotHetAtm', false);
% 
% pdb_CDR3interaction = selectionPart(pdb, 'ACDE', [1 1 93 95], [180 9 104 107]);
% plotPdb(pdb_CDR3interaction, 'coloringMethod', 'chain', 'drawingMethod', 'lines', 'plotOnlyCA', true, 'plotHetAtm', false);
%
% B. Knapp 2013-01-16

if length(chainIDs) ~= length(aaStarts) || length(aaEnds) ~= length(aaStarts)
    error ('The chains, starts, and ends must have the same number of members.');
end

[numberOfAtoms,nrOfDims,nrOfTimeSteps] = size(pdbIn.coords);

pdbOut = constructorPdb();
pdbOut.coords = zeros(1,nrOfDims,nrOfTimeSteps);

fillCtr = 1;
for atomIter = 1:1:numberOfAtoms % for all atoms
    for chainIter = 1:length(chainIDs) % do it for all parameterized chains ...
        currChain = chainIDs(chainIter);
        
        if (isequal(currChain, pdbIn.chainIdentifier{atomIter})) % only if the chain matches
            
            if (pdbIn.residueNumber(atomIter) >= aaStarts(chainIter)) % only if we are in the range of the wanted AAs
                if (pdbIn.residueNumber(atomIter) <= aaEnds(chainIter))
                    pdbOut.coords(fillCtr,:,:) = pdbIn.coords(atomIter,:,:);
                    pdbOut.coordsUnit = pdbIn.coordsUnit;
                    pdbOut.recordName{fillCtr} = pdbIn.recordName{atomIter};
                    pdbOut.atomNumber(fillCtr) = pdbIn.atomNumber(atomIter);
                    pdbOut.atomName{fillCtr} = pdbIn.atomName{atomIter};
                    pdbOut.alternateLocationIndicator{fillCtr} = pdbIn.alternateLocationIndicator{atomIter};
                    pdbOut.residueName{fillCtr} = pdbIn.residueName{atomIter};
                    pdbOut.chainIdentifier{fillCtr} = pdbIn.chainIdentifier{atomIter};
                    pdbOut.residueNumber(fillCtr) = pdbIn.residueNumber(atomIter);
                    if isfield (pdbIn, 'insertionCode') % field might not exists if it is computationally created data e.g. with MOSAICS
                        pdbOut.insertionCode{fillCtr} = pdbIn.insertionCode{atomIter};
                    end
                    if isfield (pdbIn, 'occupancy')
                        pdbOut.occupancy(fillCtr) = pdbIn.occupancy(atomIter);
                    end
                    if isfield (pdbIn, 'temperatureFactor')
                        pdbOut.temperatureFactor(fillCtr) = pdbIn.temperatureFactor(atomIter);
                    end
                    if isfield (pdbIn, 'element')
                        pdbOut.element{fillCtr} = pdbIn.element{atomIter};
                    end
                    if isfield (pdbIn, 'atomCharge')
                        pdbOut.atomCharge(fillCtr) = pdbIn.atomCharge(atomIter);
                    end
                    if isfield (pdbIn, 'information')
                        pdbOut.information = pdbIn.information;
                    end
                    fillCtr = fillCtr+1;
                end % end
            end % start
            
        end % chain
    end % for chain
end %for

