% This file contains some simple test cases for the gro2mat package. You
% can adapt them for your needs:


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% XVG TESTING %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GROMACS commands to create some (not very meaningful) trajectory analysis
% g_make_ndx -f test_pdbStructure.pdb -o test_index.ndx
% g_rms -f test_trajectoryCoords.xtc -s test_pdbStructure.pdb -o test_rmsd.xvg -n test_index.ndx -ng 4
% g_rmsf -f test_trajectoryCoords.xtc -s test_pdbStructure.pdb -o test_rmsf.xvg -n test_index.ndx
% g_sas -f test_trajectoryCoords.xtc -s test_pdbStructure.pdb -o test_sasa.xvg -n test_index.ndx
% g_energy -f test_energy.edr -s test_pdbStructure.pdb -o test_energy.xvg
% g_hbond -f test_trajectoryCoords.xtc -s test_structure.tpr -num test_hbond.xvg -n test_index.ndx
% g_gyrate -f test_trajectoryCoords.xtc -s test_pdbStructure.pdb -o test_gyration.xvg -n test_index.ndx


% do the new parsing
testNew.rmsd = parseXvgFile('./testFiles/test_rmsd.xvg');
testNew.rmsf = parseXvgFile('./testFiles/test_rmsf.xvg');
testNew.sasa = parseXvgFile('./testFiles/test_sasa.xvg');
testNew.energy = parseXvgFile('./testFiles/test_energy.xvg');
testNew.hbond = parseXvgFile('./testFiles/test_hbond.xvg');
testNew.gyration = parseXvgFile('./testFiles/test_gyration.xvg');

% % create a new reference set (execute this only if you are sure the values
% % are correct!)
% test.rmsd = parseXvgFile('test_rmsd.xvg');
% test.rmsf = parseXvgFile('test_rmsf.xvg');
% test.sasa = parseXvgFile('test_sasa.xvg');
% test.energy = parseXvgFile('test_energy.xvg');
% test.hbond = parseXvgFile('test_hbond.xvg');
% test.gyration = parseXvgFile('test_gyration.xvg');
% save test_referenceData.mat

% check the reference vs the new results
load ./testFiles/test_referenceData.mat
fprintf('comparison between the reference data and the actual version:\n');
fprintf ('rmsd : %1.0f\n', isequal(test.rmsd, testNew.rmsd));
fprintf ('rmsf : %1.0f\n', isequal(test.rmsf, testNew.rmsf));
fprintf ('sasa : %1.0f\n', isequal(test.sasa, testNew.sasa));
fprintf ('energy : %1.0f\n', isequal(test.energy, testNew.energy));
fprintf ('hbond : %1.0f\n', isequal(test.hbond, testNew.hbond));
fprintf ('gyration : %1.0f\n', isequal(test.gyration, testNew.gyration));


% plot the data
figure; plotXvgFile(testNew.rmsd);
figure; plotXvgFile(testNew.rmsf);
figure; plotXvgFile(testNew.sasa);
figure; plotXvgFile(testNew.energy);
figure; plotXvgFile(testNew.hbond);
figure; plotXvgFile(testNew.gyration);

figure;
plotXvgFile(testNew.rmsd, 2,3,1, 4);
plotXvgFile(testNew.rmsf, 2,3,2, 2);
plotXvgFile(testNew.sasa, 2,3,3, 4);
plotXvgFile(testNew.energy, 2,3,4, [2 5 6]);
plotXvgFile(testNew.hbond, 2,3,5, 2);
plotXvgFile(testNew.gyration, 2,3,6, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% PDB TESTING %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
pdb = parsePdb ('./testFiles/test_pdbStructure.pdb');
pdb_C = selectAAs(pdb, 'C', 1, 9);
plotPdb(pdb_C, 'coloringMethod', 'chain', 'drawingMethod', 'x', 'plotOnlyCA', false, 'plotHetAtm', false);

figure;
pdb_C_Ahead = selectAAs(pdb, 'AC', [1 1], [180 9]);
plotPdb(pdb_C_Ahead, 'coloringMethod', 'chain', 'drawingMethod', 'lines', 'plotOnlyCA', true, 'plotHetAtm', false);

figure;
pdb_CDR3interaction = selectAAs(pdb, 'ACDE', [1 1 93 95], [180 9 104 107]);
plotPdb(pdb_CDR3interaction, 'coloringMethod', 'chain', 'drawingMethod', 'lines', 'plotOnlyCA', true, 'plotHetAtm', false);

figure;
plotPdb(pdb_CDR3interaction);
figure;
plotPdb(pdb_CDR3interaction, 'coloringMethod', 'chain', 'drawingMethod', 'lines', 'plotOnlyCA', true, 'plotHetAtm', false, 'turnLightOn', false, 'timeStep', 1);
figure;
plotPdb(pdb_CDR3interaction, 'coloringMethod', 'chain', 'drawingMethod', 'o', 'plotOnlyCA', true, 'plotHetAtm', false, 'turnLightOn', false, 'timeStep', 1);
figure;
plotPdb(pdb_CDR3interaction, 'coloringMethod', 'chain', 'drawingMethod', 'vdw', 'plotOnlyCA', false, 'plotHetAtm', false, 'turnLightOn', false, 'timeStep', 1);
figure;
plotPdb(pdb_CDR3interaction, 'coloringMethod', 'chain', 'drawingMethod', 'vdw', 'plotOnlyCA', false, 'plotHetAtm', false, 'turnLightOn', true, 'timeStep', 1);
figure;
plotPdb(pdb_C, 'coloringMethod', 'residue', 'drawingMethod', 'vdw', 'plotOnlyCA', false, 'plotHetAtm', false, 'turnLightOn', true, 'timeStep', 1);
figure;
plotPdb(pdb_CDR3interaction, 'coloringMethod', 'chain', 'drawingMethod', 'lines', 'plotOnlyCA', false, 'plotHetAtm', false, 'turnLightOn', false, 'timeStep', 1);
figure;
plotPdb(pdb_C, 'coloringMethod', [1 0 0], 'drawingMethod', 'x', 'plotOnlyCA', false, 'plotHetAtm', false, 'turnLightOn', false, 'timeStep', 1);
figure;
plotPdb(pdb_C, 'coloringMethod', 'residueSeq', 'drawingMethod', 'vdw', 'plotOnlyCA', false, 'plotHetAtm', false, 'turnLightOn', true, 'timeStep', 1);
figure;
plotPdb(pdb_C, 'coloringMethod', 'residueName', 'drawingMethod', 'vdw', 'plotOnlyCA', false, 'plotHetAtm', false, 'turnLightOn', true, 'timeStep', 1);
figure;
plotPdb(pdb_C, 'coloringMethod', 'atom', 'drawingMethod', 'd', 'plotOnlyCA', false, 'plotHetAtm', false, 'turnLightOn', false, 'timeStep', 1);
figure; hold on;
plotPdb(pdb_C, 'coloringMethod', 'residue', 'drawingMethod', 'x', 'plotOnlyCA', false, 'plotHetAtm', false, 'turnLightOn', true, 'timeStep', 1);
plotPdb(pdb_C, 'coloringMethod', [0 0 1], 'drawingMethod', 'lines', 'plotOnlyCA', true, 'plotHetAtm', false, 'turnLightOn', true, 'timeStep', 1);
plotPdb(pdb_C, 'coloringMethod', 'residue', 'drawingMethod', 'o', 'plotOnlyCA', true, 'plotHetAtm', false, 'turnLightOn', true, 'timeStep', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% XTC testing %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trjFrames_1_10 = parseTrj('./testFiles/test_pdbStructure.pdb', './testFiles/test_trajectoryCoords.xtc', 0.001, 1, 10); % pdb matches with the first frame
trjFrames_10_20 = parseTrj('./testFiles/test_pdbStructure.pdb', './testFiles/test_trajectoryCoords.xtc', 5, 10, 19); % pdb does not match with the 10th frame
trjFrames_all = parseTrj('./testFiles/test_pdbStructure.pdb', './testFiles/test_trajectoryCoords.xtc');

figure; hold on;
plotPdb(trjFrames_1_10, 'coloringMethod', [1 0 0], 'drawingMethod', 'o', 'plotOnlyCA', false, 'timeStep', 10); % the 10th frame of 1-10 is equivalent to the 1st frame of 10-19
plotPdb(trjFrames_10_20, 'coloringMethod', [1 0 0], 'drawingMethod', 'x', 'plotOnlyCA', false, 'timeStep', 1);


trjFrames_1_10_C = selectAAs(trjFrames_1_10, 'C', 1, 9);
figure; hold on;
plotPdb(trjFrames_1_10_C, 'coloringMethod', [0.8 0 0], 'drawingMethod', 'vdw', 'plotOnlyCA', false, 'timeStep', 1, 'turnLightOn', true);
plotPdb(trjFrames_1_10_C, 'coloringMethod', [0 0.9 0], 'drawingMethod', 'vdw', 'plotOnlyCA', false, 'timeStep', 5, 'turnLightOn', true);
plotPdb(trjFrames_1_10_C, 'coloringMethod', [0 0 0.7], 'drawingMethod', 'vdw', 'plotOnlyCA', false, 'timeStep', 10, 'turnLightOn', true);

figure; hold on;
for timeIter = 1:100:size(trjFrames_all.coords, 3)
    plotPdb(trjFrames_all, 'coloringMethod', 'chain', 'drawingMethod', 'lines', 'plotOnlyCA', true, 'timeStep', timeIter);
end

figure; hold on;
trjFrames_all_chainC = selectAAs(trjFrames_all, 'C', 1, 9);
for timeIter = 1:5:size(trjFrames_all_chainC.coords, 3)
    plotPdb(trjFrames_all_chainC, 'coloringMethod', 'chain', 'drawingMethod', 'lines', 'plotOnlyCA', true, 'timeStep', timeIter);
end

figure; hold on;
trjFrames_all_chainC = selectAAs(trjFrames_all, 'C', 1, 9);
plotPdb(trjFrames_all, 'coloringMethod', 'chain', 'drawingMethod', 'lines', 'plotOnlyCA', true, 'timeStep', 1);
for timeIter = 1:5:size(trjFrames_all_chainC.coords, 3)
    plotPdb(trjFrames_all_chainC, 'coloringMethod', [1 0 0], 'drawingMethod', 'lines', 'plotOnlyCA', true, 'timeStep', timeIter);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% remaining AA tools testing %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[AAsCell, AAsArr] = createAllAAs;
h = createAAHashTable;
threeLetter = h.get(AAsArr(1));
oneLetter = h.get(threeLetter);

fprintf ('aa hash : %1.0f\n', isAA(oneLetter) && isequal(AAsArr(1), oneLetter));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% USE CASE TESTING %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% read the pdb and xtc file
trj = parseTrj('./testFiles/test_pdbStructure.pdb', './testFiles/test_trajectoryCoords.xtc');
% plot the first frame in 3d to check if everything is ok
figure;
plotPdb(trj, 'coloringMethod', 'chain', 'drawingMethod', 'lines', 'plotOnlyCA', true, 'plotHetAtm', false, 'turnLightOn', false, 'timeStep', 1);

% extract the relevant points
h1StartAllFrames = selectAAs(trj, 'A', 57, 61);
h1StartAllFrames = mean(h1StartAllFrames.coords, 1);
h1EndAllFrames = selectAAs(trj, 'A', 81, 85);
h1EndAllFrames = mean (h1EndAllFrames.coords, 1);
h2StartAllFrames = selectAAs(trj, 'A', 138, 142);
h2StartAllFrames = mean(h2StartAllFrames.coords, 1);
h2EndAllFrames = selectAAs(trj, 'A', 165, 170);
h2EndAllFrames = mean (h2EndAllFrames.coords, 1);

% calculate the distances for all time frames
nrOfFrames = size(trj.coords,3);
for frameIter = 1:nrOfFrames
    distsAllFramesH1sH2e(frameIter) = calcDist (h1StartAllFrames(1, :, frameIter), h2EndAllFrames  (1, :, frameIter)); %
    distsAllFramesH1eH2s(frameIter) = calcDist (h1EndAllFrames  (1, :, frameIter), h2StartAllFrames(1, :, frameIter));
end

% plot distance with moving average
figure; hold on;
plot (distsAllFramesH1sH2e, 'r-','LineWidth',1);
plot (mvgAvg(distsAllFramesH1sH2e, 5), 'r-','LineWidth',2);
xlabel('time steps'); ylabel('distance (nm)');

% plot the autocorrelation of the distance
figure; hold on;
plot (autocorr(distsAllFramesH1sH2e, 800), 'r-','LineWidth',1);
plot ([0, 801], [0, 0], 'k-'); % add a line at zero
xlabel('time steps'); ylabel('autocorrelation');

% plot the correlation between the distances
figure; hold on;
plot(distsAllFramesH1sH2e, distsAllFramesH1eH2s, 'kx');
lsline;  % add a least-squares fit line to the plot
cc = corr(distsAllFramesH1sH2e', distsAllFramesH1eH2s');
myTitle = sprintf ('cc=%1.2f', cc);
title(myTitle);
xlabel('H1s/H2e');
ylabel('H1e/H2s');

points = [h1StartAllFrames(:,:,1); h2EndAllFrames(:,:,1); h1EndAllFrames(:,:,1); h2StartAllFrames(:,:,1)];
writeVmdVisualisation ('vmdTclTkFile.txt', points, 0.2, 1.0, 1); % file resulting from here can be used in VMD. e.g. vmd ./testFiles/test_pdbStructure.pdb => extensions => TK console => paste the file content