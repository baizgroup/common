%%DISTANCEFINDER
% Version 1.0 2020-06-24
% Author: Joe Shirley
%
% 
% Required Parameters:
%   data         - These are  the X and Y values to fit. This should be a
%                  matrix which consistst of two column vectors.
%
%   func         - This is the function you are fitting to. This is a
%                  character vector. 'a*exp(-x/b) + c' is an example.
%
%   fitOpts      - These are fit options such as the fit algorithm and
%                  tolerance. They should be generated with 'fitoptions.'
%                  Startpoints will be determined from the upper and lower
%                  bounds, so don't set the 'Startpoint' parameter.
% 
% Optional Parameters:
%   numAvg       - This is should be an integer which is the number of
%                  times (loop iterations) you want to fit the data. By
%                  default this value is 100.
%
%   startPoint   - This specifies the algorithm for finding the start point
%                  for each iteration of the loop. By default it is a
%                  random number between the upper and lower bounds.
%                  Options include: 'rand' 'mean' 'low'... (see code)
% 
% Examples:
%   result = hyperfit(data, 'ax + b', fitOpts)
%       The data is fit to the linear equation 'ax + b' 100 times with
%       randomized starting values. A structure is returned that contains:
%       all fitting parameters, the average fit, a final fit starting at
%       the average, the goodness of the final fit, the coefficients
%       associated with the final fit, and the 95 % confidence interval on
%       those coefficients.
%
%   result = hyperfit(data, 'ax + b', fitOpts, 'numAvg', 10)
%       Same as the first, but we only average 10 fits now.
%
%   result = hyperfit(data, 'ax + b', fitOpts, 'startPoint', 'mean')
%       Same as the first, but instead of starting at a random point for
%       each fit, we start at the mean. This is useless.
%
%   result = hyperfit(data, 'ax + b', fitOpts, 'numAvg', 1000, 'startPoint', 'mean 50')
%       Same as the first, but we average 1000 fits and start each fit at
%       the mean +- randomly up to 25 % of the range. For example, if the
%       limits were [100 200], 'mean 50' would translate to 150 +- a random
%       number from 0 to 25. 50 indicates the total variation is 50 % of
%       the total range. 
%
% Header Finished

function result = distanceFinder(ref, sel, varargin)
%% Defaults
defaultBox = NaN;
defaultBonusCalcs = 'none';
defaultDownsample = 1;

%% Lets parse our inputs
p = inputParser;

% Valid coordinates are numbers in the form Nx3xF, where N is the number of
% atoms and F is the number of frames
errMsg = "Valid coordinate data (ref/sel) should be an Nx3xF numeric array, where N is the number of atoms and F is the number of frames.";
validCoord = @(x) assert(isnumeric(x) && (size(x,2) == 3), sprintf(errMsg));
addRequired(p, 'ref', validCoord);
addRequired(p, 'sel', validCoord);

% Valid box is numeric
errMsg = "Box needs to be numeric.";
validBox = @(x) assert(isnumeric(x), sprintf(errMsg));
addOptional(p, 'box', defaultBox, validBox);

% 
errMsg = "A valid input for 'numAvg' should be a scalar, positive integer on the range [1,inf).";
validDownsample = @(x) assert(isreal(x) && isscalar(x) && (x >= 1) && (rem(x,1) == 0), errMsg);
addParameter(p, 'downsample', defaultDownsample, validDownsample);

% Valid bonus calcs haven't been determined
errMsg = "A valid input for 'numAvg' should be a single, positive integer.";
validBonusCalcs = @(x) assert(isreal(x) && isscalar(x) && (sign(x) == 1) && (rem(x,1) == 0), errMsg);
addParameter(p, 'numAvg', defaultBonusCalcs, validBonusCalcs);

parse(p, ref, sel, varargin{:});

box = p.Results.box;

if isnan(box)
    warning(['You''ve failed to provide a box for the application of the periodic boundary condition.' newline ...
        'Calculations will be performed assuming there is no periodic boundary to apply...' newline ...
        'To avoid this message, add box dimensions or explicitly set them to [] for no periodic boundary.']);
    box = [];
end

ax1Coords = permute(p.Results.sel(:,:,1:p.Results.downsample:end), [1 4 2 3]);
ax2Coords = permute(p.Results.ref(:,:,1:p.Results.downsample:end), [4 1 2 3]);

box = permute(box(:,:,1:p.Results.downsample:end), [1 4 2 3]);

halfBox = box ./ 2;

numFrames = size(box, 4);
    
%% Let's do this part in a for loop so we don't run out of memory.
% It works outside of a loop if we remove the indexing, but more than a
% couple of frames will cause your computer to run out of RAM and be
% very sad. :-(   
for j = (1:5:numFrames-5)
    disp(['Frame: ' num2str(j)]);
    % The periodic boundary condition says vectors that are greater than
    % halfBox are actually box - vec. Vectors that are negative and
    % have a greater magnitude are box + vec.        
    workingDiffMatrix = ax1Coords(:,:,:,j:j+5) - ax2Coords(:,:,:,j:j+5);
    workingHalfBox = halfBox(:,:,:,j:j+5);
    workingBox = box(:,:,:,j:j+5);
    tooBig = workingDiffMatrix > workingHalfBox;
    shift = workingBox - workingDiffMatrix;
    workingDiffMatrix(tooBig) = shift(tooBig);
    clear tooBig; clear shift
    tooBig = (-workingDiffMatrix) > workingHalfBox;
    shift = workingBox + workingDiffMatrix;
    workingDiffMatrix(tooBig) = shift(tooBig);
    clear tooBig; clear shift

    %% Going from XYZ to r
    % We need to sqrt(sum(square all of these differences. First, square
    % all of the XYZ components. Then sum along the matrix dimension those
    % are stored in (3). Then square root the result. This allows us to
    % return to a 3D matrix. We only want the minimum distance, so then
    % we do a sum through the crowder atoms.
    diffSqMatrix = workingDiffMatrix .^ 2;
    clear workingDiffMatrix;
    diffSumSqMatrix = sum(diffSqMatrix,3);
    clear diffSqMatrix;
    minXYDistSq = min(diffSumSqMatrix,[],2);


    result.minDistMatrix(:,j:j+5) = minXYDistSq .^ (0.5);
end

disp(['It took ' num2str(toc, 2) ' seconds to analyze ' num2str(size(box,4)) ' frames.']);
end











