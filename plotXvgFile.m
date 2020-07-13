function [] = plotXvgFile (r, nXSubplot, nYSubplot, subplotID, columsToPlot, fontSize)
% function takes the the output of parseXvgFile and plots it
% params
%   r: the result from parseXvgFile
%   nXSubplot: nr of subplotts in x direction
%   nYSubplot: nor of subplots in y direction
%   subplotID: where to place the current plot
%   maxSize: maxium size of the data taken, e.g. 2 means 2 columns from the
%   xvg file
%   fontSize: the size of the font for text and axis (default 16)
%
% created: B.Knapp 2007-02-01
% revised: B.Knapp 2014-01-15
%
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
% figure;
% plotXvgFile(testNew.rmsd);
% figure;
% plotXvgFile(testNew.rmsd, 2,3,1, 4);
% plotXvgFile(testNew.rmsf, 2,3,2, 2);
% plotXvgFile(testNew.sasa, 2,3,3, 4);
% plotXvgFile(testNew.energy, 2,3,4, [2 5 6]);
% plotXvgFile(testNew.hbond, 2,3,5, 2);
% plotXvgFile(testNew.gyration, 2,3,6, 2);

if nargin == 4
    columsToPlot = 2:size(r.vals,2); % all but not the time
end
if nargin == 1
    columsToPlot = 2:size(r.vals,2); % all but not the time
    nXSubplot = 1;
    nYSubplot = 1;
    subplotID = 1;
end

if nargin < 6
    fontSize = 16;
end

if sum(columsToPlot==1) > 0
    warning ('Column 1 contains only the time scale (or atom id). Are you sure that you want to plot this column on the y Axis?');
end


% the underscore will cause a different formating therefore it is replaced
% with a whitespace
for legendIter = 1:length(r.legends)
    tmp = r.legends{legendIter};
    tmp(tmp=='_') = ' ';
    r.legends{legendIter} = tmp;
end

colors = create20colors;
subplot(nYSubplot, nXSubplot, subplotID);
hold on;

switch r.type % decide which type must be plotted
    case 'g_rms'
        plotRMSD (r, colors, columsToPlot, fontSize);
    case 'g_rmsf'
        plotRMSF (r, colors, columsToPlot, fontSize);
    case 'g_sas'
        plotSAS (r, colors, columsToPlot, fontSize);
    case 'g_energy'
        plotEnergy (r, colors, columsToPlot, fontSize);
    case 'g_hbond'
        plotHBond (r, colors, columsToPlot, fontSize);
    case 'g_gyrate'
        plotGyr (r, colors, columsToPlot, fontSize);
    case 'g_dist'
        plotDist (r, colors, columsToPlot, fontSize);
    otherwise
        error ('unknown type: "%s"', type);
end

hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plotGyr (r, colors, columsToPlot, fontSize)
plotData (r, colors, columsToPlot, fontSize)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plotRMSD (r, colors, columsToPlot, fontSize)
plotData (r, colors, columsToPlot, fontSize);
r.subtitle = strrep (r.subtitle, '_', ' ');
idx = strfind(r.subtitle, ' after lsq');
if idx > 1
    myTitle = [r.title, r.subtitle(idx:end)];
    title (myTitle);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plotRMSF (r, colors, columsToPlot, fontSize)
r.labelY = 'nm';
tmp = r.vals(:,1); % backup for chain finding
r.vals(:,1) = 1:size(r.vals(:,1)); % give it a consecutive ordering
plotData (r, colors, columsToPlot, fontSize);
chainIdxs = [];
for aaIter = 2:size(tmp, 1)
    if tmp(aaIter,1) < tmp(aaIter-1,1)
        chainIdxs(end+1) = aaIter;
    end
end
maxVal = max(r.vals(:,2));
for newChainIter = 1:length (chainIdxs)
   plot ([chainIdxs(newChainIter) chainIdxs(newChainIter)],[0 maxVal],'r--'); 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plotSAS (r, colors, columsToPlot, fontSize)
r.labelY = 'Area (nm2)';
plotData (r, colors, columsToPlot, fontSize);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plotEnergy (r, colors, columsToPlot, fontSize)
plotData (r, colors, columsToPlot, fontSize);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plotHBond (r, colors, columsToPlot, fontSize)
plotData (r, colors, columsToPlot, fontSize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plotDist (r, colors, columsToPlot, fontSize)
for legendIter = 1:length(r.legends)
    tmp = r.legends{legendIter};
    tmp = strrep (tmp, '\s', ' in ');
    tmp = strrep (tmp, '\N', '');
    r.legends{legendIter} = tmp;
end
plotData (r, colors, columsToPlot, fontSize);




function [] = plotData (r, colors, columsToPlot, fontSize)
myLegends = {};
for colIter = 1:length(columsToPlot)
    col = columsToPlot(colIter);
    plot(r.vals(:,1), r.vals(:, col),  'color', colors{colIter});
    myLegends{end+1} = r.legends{col-1};
end

xlabel(r.labelX, 'FontSize', fontSize);
ylabel(r.labelY, 'FontSize', fontSize);
title(r.title, 'FontSize', fontSize);
legend(myLegends, 'FontSize', fontSize);
set(gca,'FontSize', fontSize);
