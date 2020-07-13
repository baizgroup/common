function [field pos] = ElectrostaticsAllAtoms(data,in)

%this function also takes two variables
% EFieldMap(data,in)

%indices(osc, atom)
%https://en.wikipedia.org/wiki/Atomic_units
Const.FourPiEpsilonNaught = 4*pi*8.85E-12;%Epsilon naught in (F/m)
Const.AUConvField = 1./5.142E11; % Conversion between SI and Atomic Units (Eh/e ao);
Const.AUConvPot = 1./2.72114E1; %This is the correct conversion
Const.AUConvGrad = 1/9.7173E21; %Unit conversion from Wolfram Alpha
SoluteAtoms = data.trjFrames.atomName(1:data.topol.NumAtomsinMol(in.MolGroup)); %THISLINE
SoluteAtomsAll = repmat(SoluteAtoms,[data.topol.NumMolsInBox(in.MolGroup) 1]);
%build index

for Osc = 1:size(in.atomName,1)
    for Atom = 1:size(in.atomName,2)
        %build the index for current atom
        
        k = strcmp(SoluteAtomsAll,in.atomName{Osc,Atom});
        [OscInd(Osc,Atom,:),~] = find(k);
    end
end

%field.index.OscInd = OscInd;


%exclude vector
ExcludeVecAll = zeros([size(in.ExcludeAtoms,1) data.topol.NumMolsInBox(in.MolGroup)*data.topol.NumAtomsinMol(in.MolGroup)]);
 for Osc = 1:size(in.ExcludeAtoms,1)
     for Atom=1:size(in.ExcludeAtoms,2)
        k = strcmp(SoluteAtomsAll,in.ExcludeAtoms{Osc,Atom});
        ExcludeVecAll(Osc,:) = ExcludeVecAll(Osc,:)+k';
     end
 end
 
ExcludeVecAll = [ExcludeVecAll zeros([size(in.ExcludeAtoms,1) sum(data.topol.NumMolsInBox.*data.topol.NumAtomsinMol)-size(ExcludeVecAll,2)])];

tic
%these are in box coordinates

field = PotFieldGradAllAtoms(data,OscInd, ExcludeVecAll, Const, in.CalcGrad,in.excludeCutOff);
toc

field.index.OscInd = OscInd;
%atompos has the same naming convention as AtomName
pos.AtomPos = data.trjFrames.coords(OscInd, :,:);
pos.box = data.trjFrames.trajectoryData.box;

