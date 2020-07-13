function [trj] = parseTrjMat(fnPdb, fnXtc, maxTol, startFrame, endFrame)
% This function combines the output of parsePdb with parseXtc. It parses
% the atom information from a pdb file and adds the trajectory information
% from an xtc file.
%
% params:
% - fnPdb: filename of the pdb file
% - fnXtc: filename of the xtc trajectory file
% - maxTol: maximum tolerance of the absoltue coordinate difference between
% the pdb file and the first frame of the xtc file. Note that this
% difference will not be exactly zero due to floating point operations and
% binary vs ascii saving of the data.
% - startFrame: the first frame to be read
% - endFrame: the last frame to be read
% - trj: the resulting trajectory
%
% examples for usage:
%
% trj7Y = parseTrj('/data/housemartin/not-backed-up/knapp/1mi5_172_100ns_sims/wildtype/7Y.firstFrame.pdb', '/data/housemartin/not-backed-up/knapp/1mi5_172_100ns_sims/wildtype/7Y.final.md.xtc');
% trj1A = parseTrj('/data/housemartin/not-backed-up/knapp/1mi5_172_100ns_sims/1A/1A.firstFrame.pdb', '/data/housemartin/not-backed-up/knapp/1mi5_172_100ns_sims/1A/1A.final.md.xtc');
% trj7A = parseTrj('/data/housemartin/not-backed-up/knapp/1mi5_172_100ns_sims/7A/7A.firstFrame.pdb', '/data/housemartin/not-backed-up/knapp/1mi5_172_100ns_sims/7A/7A.final.md.xtc');
%
% % crash due to wrong atom number:
% trjCrash = parseTrj('/data/housemartin/not-backed-up/knapp/1mi5_172_100ns_sims/1A/1A.firstFrame.pdb', '/data/housemartin/not-backed-up/knapp/1mi5_172_100ns_sims/7A/7A.final.md.xtc');
%
% % the number of atom coords matches but not the coords itself
% trjCrash = parseTrj('/data/housemartin/not-backed-up/knapp/1mi5_172_100ns_sims/3W/3W.firstFrame.pdb', '/data/housemartin/not-backed-up/knapp/1mi5_172_100ns_sims/5W/5W.final.md.xtc');
%
%
% B. Knapp 2013-09-05

%load matlab file (pre processed)
p=strtrim(fnXtc);
q=strsplit(p,'.');
if strcmp(q(end),'mat');
     load(fnXtc);
else
    %load xtc
    if nargin == 2
        maxTol = 0.001;
        xtc = parseXtc(fnXtc);
    elseif nargin == 3
        xtc = parseXtc(fnXtc);
    elseif nargin == 5
        xtc = parseXtc(fnXtc, startFrame, endFrame);
    else
        error ('invalid number of input arguments');
    end
end
    

pdb = parsePdb(fnPdb);
sPdb = size(pdb.coords, 1);
sXtc = size(xtc.coords(:,:,1),1);
if (sPdb ~= sXtc)
    error('The number of atoms does not match between the pdb of the file "%s" (n=%i) and the first frame of the xtc of "%s" (n=%i)). Are you sure that those two files match?', fnPdb, sPdb, fnXtc, sXtc);
end

trj = pdb;
trj.coords = xtc.coords;
trj.trajectoryData = xtc.trajectoryData;
trj.sourceFileName = fnPdb;
trj.sourceFileNameXtc = fnXtc;