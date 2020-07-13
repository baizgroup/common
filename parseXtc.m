function [xtc, tElapsed] = parseXtc( filename, startFrame, endFrame)
% Mex-version of the gromacs xtc-parser: This function loads a given xtc
% file into matlab using the a modified gromacs C code via MEX.
%
% params:
% - filename: the xtc file to read
% - startFrame: the first frame to be read
% - endFrame: the last frame to be read
% - xtc: the resulting xtc trajectory
% - tElapsed: how long the parsing took
%
% H. Dien 2012 
% H. Dien last modified 2014-01-14

    xtc=[];
    
%
% check arguments
%
    if(exist('filename') == 0  || isnumeric(filename) == 1)
       disp('ERROR: First argument is not a string!');
       return;        
    end
    
    if(exist('startFrame') && isnumeric(startFrame) == 0)
       disp('ERROR: Second argument is not a number!');
       return;        
    end

    if(exist('endFrame') && isnumeric(endFrame) == 0)
       disp('ERROR: Third argument is not a number!');
       return;        
    end
    
    if(exist('startFrame') && exist('endFrame') == 0)
        % testing with xUnit test Framework -> use disp instead error()
        disp('ERROR: Please enter an endframe');
        return;
    end
    
    if(exist('startFrame') == 0)
        startFrame=1;
    end
    
    if(exist('endFrame') == 0)
        endFrame=0;
    else
       if(endFrame < 1)
            disp('ERROR: Endframe has to be higher than zero!'); 
            return;
       end 
    end
        
    if(startFrame < 1)
        disp('ERROR: Startframe has to be higher than zero!');
         return;
    end

    if(startFrame > endFrame && endFrame ~= 0) 
         disp('ERROR: Endframe has to be higher than startframe!');
          return;
    end
%
% end chk args
%
    tStart = tic; 
    
   [coords, steps, time, precision, box] = gromacsMex(startFrame, endFrame, filename);
  
    xtc = constructorPdb();
    xtc.coords = coords;
    xtc.trajectoryData.step = steps;
    xtc.trajectoryData.time = time;
    xtc.trajectoryData.precision = precision;
    xtc.trajectoryData.box = box;
   
    tElapsed = toc(tStart);
     fprintf('\n xtcread Done Time: %.2f \n', tElapsed);
     fprintf('\n xtc read: %.0f frames\n', size(coords,3));
end


