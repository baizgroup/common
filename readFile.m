function [data, nrOfLines] = readFile (fileName, delimiter)
% [data, nrOfLines] = readFile (fileName, delimiter);
% function reads a file at once (much better time performance, but more
% memory neccesary), returns a cell array of data-lines
%
% B.Knapp 2007-11-02
% 
%   Created:            $LastChangedBy: bknapp $
%   Last modified:      $LastChangedDate: 2011-05-04 14:35:21 +0200 (Mi, 04 Mai 2011) $
%   Revision:           $LastChangedRevision: 6 $
%   Version:            $Id: readFile.m 10 2011-05-05 09:50:54Z bknapp $
%
% example for call:
%
% fileName = 'W:\bknapp\pdb-files\_experimentelle_Bindingsdaten\MHCI_Community_Set\hla_a0201\comparisonWithStructureBasedMethods\YMMGIEYGL_vs_GLFVLLAFL\peptX_24_10_2007-16_46_14.log.xml';
% [data, nrOfLines] = readFile (fileName, '\n');
% 

fid=fopen(fileName,'r');

if fid==-1
    error('Could not read file "%s".', fileName);
end
if nargin == 1
    delimiter = '\n';
end

% fprintf('reading data.\n');
fullText     = fread(fid,'char=>char')';
fclose(fid);

% fprintf('splitting data.\n');
data = strread(fullText,'%s','delimiter', delimiter); % use textscan instead?
clear('fulltext');


nrOfLines = numel(data);
%fprintf('done reading data.\n');

