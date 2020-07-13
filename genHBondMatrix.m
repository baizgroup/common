%GENHBONDMATRIX This function determines if there is a hydrogen bond between two water molecules.
% Version 0.1 2020-01-23
% Author: Me, Joe Shirley
% 
% Required Parameters:
%   ptmData      - Data that comes out of parseTrjMat and that has been
%                  futher processed by parseTopMat. Both are required.
% 
% References: None
% Notes: None
%
% Header finished.

function result = genHBondMatrix(ptmData)
    % We start by assuming no hydrogen bonding
    bondStatus = false;
    % A good defualt angle for hydrogen bonding    
    defaultAngle = 30; 
    % We assume the waters are close enough to hydrogen bond.
    assumeClose = -1; 
    
    % Lets parse our inputs
    p = inputParser;
    
    % Valid data from parseTrjMat is a structure and some more things
    errMsg = "Valid data must be from parseTrjMat and parseTopMat." + newline + ...
             "It must be a structure with a field called trjFrames and a field in trjFrames called atomName." + newline + ...
             "It must also have a field called topol and a field called bondMatrix in topol.";
    validData = @(x) assert(isstruct(x) && isfield(x, 'trjFrames') && isfield(x.trjFrames, 'atomName'), sprintf(errMsg));
    addRequired(p, 'ptmData', validData);
    
    parse(p, ptmData);

    % If the user inputs a bond length, we will perform the closeness test
    if (p.Results.bondLength ~= -1)
        % Calculate the distance between the oxygens
        componentDifference = p.Results.watCoords1(1,:) - p.Results.watCoords2(1,:);
        differenceMagnitude = norm(componentDifference);
        
        % If the distance is larger than the bondLength, we know they are
        % not hydrogen bonded, and can quit the function
        if (differenceMagnitude > p.Results.bondLength)
            result = false;
            return
        end
        % In any other case, the waters have passed the closeness test and
        % can now be treated as such
    end
    
    % Calculate the sum displacement between each hydrogen and both oxygens
    %% This needs to take into account PBC....... How can we do that without including the PBC function or having to pass the box to the function....
    sumDisplacement = zeros(4, 1);
    hydrogenCoords = [p.Results.watCoords1(2:end, :); p.Results.watCoords2(2:end, :)];
    for n = [1:4]
        sumDisplacement(n) = abs(norm(PBC(hydrogenCoords(n, :) - p.Results.watCoords1(1, :), boxSize)) + norm(PBC(hydrogenCoords(n, :) - p.Results.watCoords2(1, :), p.Results.boxSize)));
    end
    
    % Find the hydrogen that is closest (on average) to both oxygens
    [~, HofInterest] = min(sumDisplacement);

    % Determine the bond angle
    Ang = acosd((PBC(hydrogenCoords(HofInterest,:) - p.Results.watCoords1(1,:), p.Results.boxSize)*PBC(hydrogenCoords(HofInterest,:) - p.Results.watCoords2(1,:), p.Results.boxSize)')./(norm(PBC(hydrogenCoords(HofInterest,:) - p.Results.watCoords1(1,:), p.Results.boxSize)).*norm(PBC(hydrogenCoords(HofInterest,:) - p.Results.watCoords2(1,:), p.Results.boxSize))));

    % If the bond angle is acceptable, we have a hydrogen bond
    if real(Ang) < p.Results.angle || real(Ang) > (180-p.Results.angle)
        bondStatus = true;
    end
    
    result = bondStatus;
end








%PBC Returns a vector given a vector and a box for periodic boundary.
function vec = PBC(vec, boxSize)
    % calculate vec in periodic boundary conditions
    for n = 1:3
        Lx = boxSize(n);
        hLx = Lx/2;

        more =  vec(:,n) > hLx;
        less =  vec(:,n) <= -hLx;
    
        vec(more,n) = vec(more,n) - Lx;
        vec(less,n) = vec(less,n) + Lx;
    end
    
end













