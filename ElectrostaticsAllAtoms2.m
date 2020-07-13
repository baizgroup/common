function [field pos] = ElectrostaticsAllAtoms2(data,in)

%this function also takes two variables
% EFieldMap(data,in)

%indices(osc, atom)
%https://en.wikipedia.org/wiki/Atomic_units
Const.FourPiEpsilonNaught = 4*pi*8.8541878E-12;%Epsilon naught in (F/m)
Const.AUConvField = 1./5.142E11; % Conversion between SI and Atomic Units (Eh/e ao);
Const.AUConvPot = 1./2.72114E1; %This is the correct conversion
Const.AUConvGrad = 1/9.7173E21; %Unit conversion from Wolfram Alpha

%build atom index
SoluteAtoms = strcmp(data.trjFrames.residueName,in.SoluteMolName);
SoluteAtomsAll = data.trjFrames.atomName(SoluteAtoms);
%This is a general method for calculating the number of solute molecules -
%cb may 2019
FirstSoluteAtom = min(data.trjFrames.atomNumber(SoluteAtoms));
NumSoluteMols = max(data.trjFrames.residueNumber(SoluteAtoms)) - min(data.trjFrames.residueNumber(SoluteAtoms)) +1;
NumAtomsSolute = sum(SoluteAtoms)/NumSoluteMols;

for Osc = 1:size(in.atomName,1)
    for Atom = 1:size(in.atomName,2)
        %build the index for current atom
        
        k = strcmp(SoluteAtomsAll,in.atomName{Osc,Atom});
        [OscInd(Osc,Atom,:),~] = find(k);
    end
end

OscInd = OscInd+FirstSoluteAtom-1;

%exclude vector
ExcludeVecAll = zeros([size(in.ExcludeAtoms,1) size(data.trjFrames.residueName)]);
    k2 = strcmp(data.trjFrames.residueName,in.SoluteMolName);
for Osc = 1:size(in.ExcludeAtoms,1)
    for Atom=1:size(in.ExcludeAtoms,2)
    k1 = strcmp(data.trjFrames.atomName,in.ExcludeAtoms{Osc,Atom});
    ExcludeVecAll(Osc,:) = ExcludeVecAll(Osc,:)+(k1' & k2');
    end
end


tic
%these are in box coordinates
field = PotFieldGradAllAtoms(data,OscInd, ExcludeVecAll, Const, in.CalcGrad, NumSoluteMols, FirstSoluteAtom, NumAtomsSolute);
toc

field.index.OscInd = OscInd;
%atompos has the same naming convention as AtomName
pos.AtomPos = data.trjFrames.coords(OscInd, :,:);
pos.box = data.trjFrames.trajectoryData.box;

