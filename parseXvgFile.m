function [result] = parseXvgFile (fileName)
% function [] = parseXvgFile (fileName)
% function parses an n-dim *.xvg file (output from gromacs (g_rms, g_rmsf,
% g_sas, g_gyrate, g_energy, g_hbond) and decides internally which
% one it is.
%
% params:
%   fileName: the xvg file
%   result: struct consisting of the values (nrOfValues x dimension) and
%   several additional information
%
%
%
%
% created: B.Knapp 2007-02-06
% revised: B.Knapp 2014-01-14
%
% example for usage:
%
% myPath = '/data/housemartin/knapp/projects/dien_rmsf/testFiles/testsForHung/gro2mat_package/functions_andTests/';
% test.rmsd = parseXvgFile([myPath, 'rmsd.xvg']);
% test.rmsf = parseXvgFile([myPath, 'rmsf.xvg']);
% test.sasa = parseXvgFile([myPath, 'sasa.xvg']);
% test.enery = parseXvgFile([myPath, 'energy.xvg']);
% test.hbond = parseXvgFile([myPath, 'hbond.xvg']);
% test.gyration = parseXvgFile([myPath, 'gyration.xvg']);

[data, ~] = readFile (fileName);

typeFound = false;
type = 'undefined';
line = '#begin';
lineIter = 0;
while ischar(line) && (strcmp(line(1), '#')) && ~typeFound % parse for the type of xvg-file
    lineIter = lineIter+1;
    line = data{lineIter};
    % e.g. g_rms is part of G R O M A C S
    pat = '(?<type>\w+)\sis\spart\sof\sG\sR\sO\sM\sA\sC\sS';
    m = regexp(line, pat, 'names');
    if length(m) > 0
        type = m.type;
        typeFound = true;
    end
end
fprintf('starting evaluating "%s" for "%s"-type, depending on the filesize this may take some time.\n', fileName, type);

if strcmp (type(end-3:end), '_mpi') % remove the mpi appendix
   type = type(1:end-4); 
end

switch type % decide which type must be parsed
    case 'g_rms'
        result = parseRSMD (data, lineIter);
    case 'g_rmsf'
        result = parseRSMF (data, lineIter);
    case 'g_sas'
        result = parseSAS (data, lineIter);
    case 'g_energy'
        result = parseEnergy (data, lineIter);
    case 'g_hbond'
        result = parseHBond (data, lineIter);
    case 'g_gyrate'
        result = parseGyrate (data, lineIter);
    case 'g_dist'
        result = parseDist (data, lineIter);
    otherwise
        error ('unknown type: "%s"', type);
end

result.type = type;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result] = parseRSMD (data, lineIter) % root mean square deviation
legends = {};
vals = [];
line = data{lineIter};
while ischar(line) && lineIter < numel(data)
    lineIter = lineIter+1;
    line = data{lineIter};
    if (isempty(line) || strcmp(line(1), '#')) % unimportant header
        %ignore
    elseif (strcmp(line(1), '@')) % important header
        % e.g. @    title "RMSD"
        pat = '@\s+title\s+"(?<title>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.title = m.title;
        end
        % e.g. @ subtitle "of chA_&_Backbone after lsq fit to Protein"
        pat = 'subtitle\s+"(?<subtitle>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.subtitle = m.subtitle;
        end
        % parse for the legends
        % eg: @ s0 legend "chA"
        pat = 'legend\s+".+"';
        m = regexp(line, pat, 'match');
        if (length(m)>0)
            quotes = find (m{1}=='"');
            legends{length(legends)+1} = m{1}(quotes(1)+1:quotes(2)-1);
        end
        %parse for the labels of the axis
        % e.g. @    xaxis  label "Time (ps)"
        pat = 'xaxis\s+label\s+"(?<label>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.labelX = m.label;
        end
        % e.g. @    yaxis  label "RMSD (nm)"
        pat = 'yaxis\s+label\s+"(?<label>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.labelY = m.label;
        end
    else % parse for the values
        if line ~= -1
            % eg:   0.0000   0.0485   0.0482   0.0418
            pos = size(vals,1)+1;
            tmp = textscan(line, '%f');
            vals(pos, :) = tmp{1};
        end
    end
end
result.legends = legends;
if isempty (result.legends)
    result.legends{1} = 'n/a';
end
result.vals = vals;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result] = parseRSMF (data, lineIter) % root mean square fluctuation
vals = [];
line = data{lineIter};
while ischar(line) && lineIter < numel(data)
    lineIter = lineIter+1;
    line = data{lineIter};
    if (strcmp(line(1), '#')) % unimportant header
        %ignore
    elseif (strcmp(line(1), '@')) % important header
        % e.g. @    title "RMS fluctuation"
        pat = 'title\s+"(?<title>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.title = m.title;
        end
        %parse for the labels of the axis
        % e.g. @    xaxis  label "Atom"
        pat = 'xaxis\s+label\s+"(?<label>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.labelX = m.label;
        end
        % e.g. @    yaxis  label "(nm)"
        pat = 'yaxis\s+label\s+"(?<label>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.labelY = m.label;
        end
    else % parse for the values
        if line ~= -1
            % eg:  182   0.0546 (the RMSF has always just 2 columns)
            pat = '(?<id>\d+)\s+(?<val>\d+\.\d+)';
            m = regexp(line, pat, 'names');
            pos = size(vals,1)+1;
            vals(pos, 1) = str2double(m.id);
            vals(pos, 2) = str2double(m.val);
        end
    end
end

result.legends{1} = 'RMSF'; % no legends available for RMSF data
result.vals = vals;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result] = parseSAS (data, lineIter) % Solvent Accessible Surface
legends = {};
vals = [];
line = data{lineIter};
while ischar(line) && lineIter < numel(data)
    lineIter = lineIter+1;
    line = data{lineIter};
    if (strcmp(line(1), '#')) % unimportant header
        %ignore
    elseif (strcmp(line(1), '@')) % important header
        % e.g. @    title "Solvent Accessible Surface"
        pat = 'title\s+"(?<title>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.title = m.title;
        end
        % parse for the legends
        % eg: @ s0 legend "chA"
        pat = 'legend\s+".+"';
        m = regexp(line, pat, 'match');
        if (length(m)>0)
            quotes = find (m{1}=='"');
            legends{length(legends)+1} = m{1}(quotes(1)+1:quotes(2)-1);
        end
        %parse for the labels of the axis
        % e.g. @    xaxis  label "Time (ps)"
        pat = 'xaxis\s+label\s+"(?<label>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.labelX = m.label;
        end
        % e.g. @    yaxis  label "Area (nm\S2\N)"
        pat = 'yaxis\s+label\s+"(?<label>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.labelY = m.label;
        end
    else % parse for the values
        if line ~= -1
            % eg:       100           0   0.0235718   0.0235718   -0.147936
            % (the sasa has always 4 colums (time, hydrophobic, hydrophilic, and totla)))
            pat = '(?<vals>-?\d+\.?\d*)';
            m = regexp(line, pat, 'names');
            pos = size(vals,1)+1;
            for ii=1:size(m,2)
                vals(pos, ii) = str2double(m(ii).vals);
            end
            
        end
    end
end
result.legends = legends;
if isempty (result.legends)
    result.legends{1} = 'n/a';
end
result.vals = vals;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result] = parseEnergy (data, lineIter) % energy
legends = {};
vals = [];
line = data{lineIter};
while ischar(line) && lineIter < numel(data)
    lineIter = lineIter+1;
    line = data{lineIter};
    if (strcmp(line(1), '#')) % unimportant header
        %ignore
    elseif (strcmp(line(1), '@')) % important header
         % e.g. @    title "Gromacs Energies"
        pat = 'title\s+"(?<title>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.title = m.title;
        end
        % parse for the legends
        % eg: @ s0 legend "G96Angle"
        pat = 'legend\s+".+"';
        m = regexp(line, pat, 'match');
        if (~isempty(m))
            quotes = find (m{1}=='"');
            legends{length(legends)+1} = m{1}(quotes(1)+1:quotes(2)-1);
        end
        %parse for the labels of the axis
        % e.g. @    xaxis  label "Time (ps)"
        pat = 'xaxis\s+label\s+"(?<label>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.labelX = m.label;
        end
        % e.g. @    yaxis  label "(kJ/mol), (K), (bar)"
        pat = 'yaxis\s+label\s+"(?<label>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.labelY = m.label;
        end
    else % parse for the values
        if line ~= -1
            % eg:  182   0.0546
            pos = size(vals,1)+1;
            tmp = textscan(line, '%f');
            vals(pos, :) = tmp{1};
        end
    end
end
result.legends = legends;
if isempty (result.legends)
    result.legends{1} = 'n/a';
end
result.vals = vals;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result] = parseHBond (data, lineIter) % hydrogen bond
legends = {};
vals = [];
line = data{lineIter};
while ischar(line) && lineIter < numel(data)
    lineIter = lineIter+1;
    line = data{lineIter};
    if (strcmp(line(1), '#')) % unimportant header
        %ignore
    elseif (strcmp(line(1), '@')) % important header
         % e.g. @    title "Hydrogen Bonds"
        pat = 'title\s+"(?<title>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.title = m.title;
        end
        %parse for the labels of the axis
        % e.g. @    xaxis  label "Time (ps)"
        pat = 'xaxis\s+label\s+"(?<label>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.labelX = m.label;
        end
        % e.g. @    yaxis  label "Number"
        pat = 'yaxis\s+label\s+"(?<label>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.labelY = m.label;
        end
        pat = 'legend\s+".+"';
        m = regexp(line, pat, 'match');
        if (~isempty(m))
            quotes = find (m{1}=='"');
            legends{length(legends)+1} = m{1}(quotes(1)+1:quotes(2)-1);
        end
    else % parse for the values
        if line ~= -1
            % eg:  182   0.0546
            pos = size(vals,1)+1;
            tmp = textscan(line, '%f');
            vals(pos, :) = tmp{1};

        end
    end
end
result.legends = legends;
if isempty (result.legends)
    result.legends{1} = 'n/a';
end
result.vals = vals;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result] = parseGyrate (data, lineIter) % radius of gyration parsing
legends = {};
vals = [];
line = data{lineIter};
while ischar(line) && lineIter < numel(data)
    lineIter = lineIter+1;
    line = data{lineIter};
    if (isempty(line) || strcmp(line(1), '#')) % unimportant header
        %ignore
    elseif (strcmp(line(1), '@')) % important header
         % e.g. @    title "Radius of Gyration"
        pat = 'title\s+"(?<title>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.title = m.title;
        end
        % parse for the legends
        % eg: @ s0 legend "Rg"
        pat = 'legend\s+".+"';
        m = regexp(line, pat, 'match');
        if (length(m)>0)
            quotes = find (m{1}=='"');
            legends{length(legends)+1} = m{1}(quotes(1)+1:quotes(2)-1);
        end
        %parse for the labels of the axis
        % e.g. @    xaxis  label "Time (ps)"
        pat = 'xaxis\s+label\s+"(?<label>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.labelX = m.label;
        end
        % e.g. @    yaxis  label "Rg (nm)"
        pat = 'yaxis\s+label\s+"(?<label>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.labelY = m.label;
        end
    else % parse for the values
        if line ~= -1
            % eg:   0.0000   0.0485   0.0482   0.0418 (always 4 values)
            pos = size(vals,1)+1;
            tmp = textscan(line, '%f');
            vals(pos, :) = tmp{1};
        end
    end
end
result.legends = legends;
if isempty (result.legends)
    result.legends{1} = 'n/a';
end
result.vals = vals;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result] = parseDist (data, lineIter) % radius of gyration parsing
legends = {};
vals = [];
line = data{lineIter};
while ischar(line) && lineIter < numel(data)
    lineIter = lineIter+1;
    line = data{lineIter};
    if (isempty(line) || strcmp(line(1), '#')) % unimportant header
        %ignore
    elseif (strcmp(line(1), '@')) % important header
         % e.g. @    title "Distance"
        pat = 'title\s+"(?<title>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.title = m.title;
        end
        % parse for the legends
        % eg: 
        % @ s0 legend "|d|"
        % @ s1 legend "d\sx\N"
        % @ s2 legend "d\sy\N"
        % @ s3 legend "d\sz\N"
        pat = 'legend\s+".+"';
        m = regexp(line, pat, 'match');
        if (length(m)>0)
            quotes = find (m{1}=='"');
            legends{length(legends)+1} = m{1}(quotes(1)+1:quotes(2)-1);
        end
        %parse for the labels of the axis
        % e.g. @    xaxis  label "Time (ps)"
        pat = 'xaxis\s+label\s+"(?<label>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.labelX = m.label;
        end
        % e.g. @    yaxis  label "Distance (nm)"
        pat = 'yaxis\s+label\s+"(?<label>.+)"';
        m = regexp(line, pat, 'names');
        if length(m) > 0
            result.labelY = m.label;
        end
    else % parse for the values
        if line ~= -1
            % eg:  0.0000000    2.1721919    0.0792522   -2.0870516    0.5969520
            pos = size(vals,1)+1;
            tmp = textscan(line, '%f');
            vals(pos, :) = tmp{1};
        end
    end
end
result.legends = legends;
if isempty (result.legends)
    result.legends{1} = 'n/a';
end
result.vals = vals;


