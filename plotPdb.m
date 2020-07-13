function [] = plotPdb(pdb, varargin)
% function plots an input created by "constructorPdb" and threedimensionally
% displays its content in matlab.
%
% Further parameters can be set via:
% - coloringMethod: 'chain', 'atom', 'residueSeq', 'residueName' or any
%   specific color given as array e.g. [1 0.5 0.11]
% - drawingMethod: 'lines', 'vdw', or any of matlabs plotting type (see
%   help plot)
% - turnLightOn: (true or false) mostly for the vdw representation. It looks
%   better but takes longer
% - plotOnlyCA: (true or false) plot only the C-alpha atoms
% - plotHetAtm: (true or false) plot also the heterogenous atoms
% - timeStep: the time step which should be plotted (default=1)
%
%   B. Knapp 2011-08-02
%   B. Knapp revised 2014-04-10
%
% example for usage:
%
%   pdb = parsePdb ('/data/housemartin/knapp/projects/dien_rmsf/testFiles/testsForHung/gro2mat_package/functions_andTests/pdbStructure.pdb');
%   plotPdb(pdb, 'coloringMethod', 'chain', 'drawingMethod', 'lines', 'plotOnlyCA', true, 'plotHetAtm', false, 'turnLightOn', false, 'timeStep', 1)
%   plotPdb(pdb, 'coloringMethod', 'atom', 'drawingMethod', 'dots', 'plotOnlyCA', false, 'plotHetAtm', false, 'turnLightOn', false, 'timeStep', 1)
%   plotPdb(pdb, 'coloringMethod', 'chain', 'drawingMethod', 'vdw', 'plotOnlyCA', false, 'plotHetAtm', false, 'turnLightOn', true, 'timeStep', 1)
%

hold on;
argStruct = struct('coloringMethod', 'chain', ...
    'drawingMethod', 'o', ...
    'turnLightOn', false, ...
    'plotOnlyCA',  false, ...
    'plotHetAtm', false, ...
    'timeStep', 1 ...
    );

argStruct = parseArgs(varargin, argStruct);

if argStruct.timeStep > size(pdb.coords,3)
    error ('Frame number %i can not be drawn because there are only %i frames.', argStruct.timeStep, size(pdb.coords,3));
end

if argStruct.turnLightOn
    light;
    light('Position',[-1 -1 -2]);
end



[x,y,z] = sphere(20);
[c, ~] = create20colors();
aaHash = createAAHashTable;
[~, arrAAs] = createAllAAs;
chainNr = 0;
oldChain = '';
plotList = {};

for atomIter = 1:size(pdb.coords, 1) % for all atoms in the pdb file
    
    if strcmp(pdb.recordName{atomIter}, 'ATOM') || (argStruct.plotHetAtm && strcmp(pdb.recordName{atomIter}, 'HETATM')) % enter if it is ATOM or HETATM (+we want to draw HETATMs)
        if strcmp(pdb.atomName{atomIter}, 'CA') || not(argStruct.plotOnlyCA) % enter if it is a CA or no CA and we want to draw all
            
            if (strcmp(pdb.recordName{atomIter}, 'HETATM')) % it is a heterogenious ATOM (they are always plotted as the same
                
                r = 0.03;
                surface('XData',pdb.coords(atomIter,1,argStruct.timeStep) + r*x, ...
                    'YData',pdb.coords(atomIter,2,argStruct.timeStep)     + r*y, ...
                    'ZData',pdb.coords(atomIter,3,argStruct.timeStep)     + r*z, ...
                    'FaceColor','m',...
                    'EdgeColor','none','FaceLighting','gouraud');
                
            elseif strcmp(pdb.recordName{atomIter}, 'ATOM') % it is an ATOM here all our settings apply
                
                % decide about the color
                if strcmp(argStruct.coloringMethod, 'atom') % coloring according to the atom type
                    switch(pdb.atomName{atomIter}(1))
                        case  'H', color = [0.0 0.0 0.0]; % white
                        case  'C', color = [0.8 0.8 0.8]; % gray
                        case  'O', color = [1.0 0.0 0.0]; % red
                        case  'N', color = [0.0 0.0 1.0]; % blue
                        case  'S', color = [1.0 1.0 0.0]; % yellow
                        otherwise, color = [1.0 1.0 1.0];
                    end
                elseif strcmp(argStruct.coloringMethod, 'chain') % coloring according to the chain
                    colorPos = unicode2native(pdb.chainIdentifier{atomIter})-unicode2native('A')+1;
                    color = c{colorPos};
                elseif strcmp(argStruct.coloringMethod, 'residueSeq') % coloring according to the residue sequence (1 to n)
                    colorPos = pdb.residueNumber(atomIter);
                    if colorPos > length(c)
                        colorPos = mod(colorPos, length(c))+1;
                    end
                    color = c{colorPos};
                elseif strcmp(argStruct.coloringMethod, 'residueName') || strcmp(argStruct.coloringMethod, 'residue') % coloring according to the residue name
                    currAA = aaHash.get(pdb.residueName{atomIter}); % transform 3 letter to one letter
                    aaNumber = strfind (arrAAs, currAA); % transform to number
                    color = c{aaNumber};
                elseif isnumeric (argStruct.coloringMethod) % a color code is given
                    color = argStruct.coloringMethod;
                else % invalid name for color
                    error ('Coloring type "%s" is not known.', argStruct.coloringMethod);
                end
                
                % decide about the type
                if strcmp(argStruct.drawingMethod, 'vdw')
                    switch(pdb.atomName{atomIter}(1)) % set size for vdw
                        case  'H',  r = 0.06;
                        case  'C',  r = 0.10;
                        case  'O',  r = 0.10;
                        case  'N',  r = 0.08;
                        otherwise,  r = 0.10;
                    end
                    surface('XData',pdb.coords(atomIter,1,argStruct.timeStep) + r*x, ...
                        'YData',pdb.coords(atomIter,2,argStruct.timeStep) + r*y, ...
                        'ZData',pdb.coords(atomIter,3,argStruct.timeStep) + r*z, ...
                        'FaceColor',color,...
                        'EdgeColor','none','FaceLighting','gouraud');
                elseif length(argStruct.drawingMethod) == 1 % one of matlabs 1 letter symbols
                    plot3(pdb.coords(atomIter,1,argStruct.timeStep), pdb.coords(atomIter,2,argStruct.timeStep), pdb.coords(atomIter,3,argStruct.timeStep), argStruct.drawingMethod, 'MarkerEdgeColor', color);
                elseif strcmp(argStruct.drawingMethod, 'lines') % a meaning full line i.e. chains are not connected
                    % in this case we do not directly plot but the save the coords
                    % for a continous plot later
                    if (~strcmp(pdb.chainIdentifier{atomIter}, oldChain)) % we found a new chain
                        oldChain = pdb.chainIdentifier{atomIter};
                        chainNr = chainNr+1;
                        plotList{chainNr} = [];
                    end
                    plotList{chainNr}(end+1, :) = pdb.coords(atomIter,:,argStruct.timeStep);
                else
                    error ('Unkown drawing type "%s".', argStruct.drawingMethod);
                end % type
                
            end % HEATM or ATOM
        end % CA or not CA
    end % is it a valid ATOM entry
end % for all atoms


if strcmp(argStruct.drawingMethod, 'lines')
    if not(strcmp(argStruct.coloringMethod, 'chain') || isnumeric (argStruct.coloringMethod))
        warning ('If you want a continuos line per chain it has to be the same color');
    end
    for chainIter = 1:length(plotList)
        if strcmp(argStruct.coloringMethod, 'chain')
            plot3(plotList{chainIter}(:,1), plotList{chainIter}(:,2), plotList{chainIter}(:,3), '-', 'color',c{chainIter});
        else
            plot3(plotList{chainIter}(:,1), plotList{chainIter}(:,2), plotList{chainIter}(:,3), '-', 'color', color);
        end
    end
end


hold off;



