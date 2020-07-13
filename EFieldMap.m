function output = EFieldMap(data,in, Const)

%this function also takes two variables
% EFieldMap(data,in)

%indices(osc, atom)
Const.FourPiEpsilonNaught = 4*pi*8.85E-12;%Epsilon naught in (F/m)
Const.AUConv = 1./5.142E11; % Conversion between SI and Atomic Units (Eh/e ao);
%Const.MapCoeff = 2371.4;

SoluteAtoms = data.trjFrames.atomName(1:data.topol.NumAtomsinMol(1));
    

%initialize output matrix
output.FreqShift = zeros(data.topol.NumMolsInBox(1),size(in.atomName,1),size(data.trjFrames.coords,3));
output.COcentArray = zeros(3, data.topol.NumMolsInBox(1),size(in.atomName,1),size(data.trjFrames.coords,3));
output.COArray = zeros(3, data.topol.NumMolsInBox(1),size(in.atomName,1),size(data.trjFrames.coords,3));

%build index
for Osc = 1:size(in.atomName,1)
    for Atom = 1:size(in.atomName,2)
        %build the index for current atom
        k = strcmp(SoluteAtoms,in.atomName{Osc,Atom});
        [OscInd(Osc,Atom,:),~] = find(k);
    end
end

%First build index
for Vec = 1:size(in.atomVecName,1) %each vector
    for Osc = 1:size(in.atomVecName,2)
        %split string
        StrSpl = strsplit(in.atomVecName{Vec,Osc},'-');
        [IndexVec(Vec,Osc,1,:),~] = find(strcmp(SoluteAtoms,char(StrSpl(1))));
        [IndexVec(Vec,Osc,2,:),~] = find(strcmp(SoluteAtoms,char(StrSpl(2))));
    end
end


if isfield(in,'ExcludeAtomRadius')
ExcludeVec = in.ExcludeAtomRadius;
disp(['Excluding atoms within: ' num2str(ExcludeVec) ' nm'])
else
    %Exclude Atom Index;
    

for Osc = 1:size(in.ExcludeAtoms,1) %each vector
    for Atom = 1:size(in.ExcludeAtoms,2)
        [IndexVecExclAtoms(Osc,Atom,:),~] = find(strcmp(SoluteAtoms,in.ExcludeAtoms{Osc,Atom}));
    end
end
%Create the Atom Exclude Vector

ExcludeVec = zeros(size(in.ExcludeAtoms,1),size(data.trjFrames.atomCharge,1));

for Osc = 1:size(in.ExcludeAtoms,1)
tempVec = IndexVecExclAtoms(Osc,:,:);
ExcludeVec(size(in.ExcludeAtoms,1),tempVec(:)) = 1;
end

end
%This function calculates the electric field (x,y,z) at
%all the atom positions. OscInd(Osc, Atom)
%EFieldAtom(:,Atom,Osc,Mol,frame)
EFieldAtom = EfieldAllAtoms(data,OscInd, ExcludeVec, Const);

EFieldContribution = zeros(size(in.atomVecName,1),1);

%Index atoms and build vector index
for frame = 1:size(data.trjFrames.coords,3);
    box = data.trjFrames.trajectoryData.box(frame,:);
    %box(3) = 10000; %make pseudo 2D PBC
    
    for Mol  = 1:data.topol.NumMolsInBox(1);
        %disp(['Mol:' num2str(Mol)]);
        for Osc = 1:size(in.atomVecName,2)
            %disp(['Osc:' num2str(Osc)]);
            for Vec = 1:size(in.atomVecName,1) %each vector
                %disp(['Vec:' num2str(Vec)]);
                
                CurrAtomPos1 = data.trjFrames.coords(IndexVec(Vec,Osc,1,Mol),:,frame);
                CurrAtomPos2 = data.trjFrames.coords(IndexVec(Vec,Osc,2,Mol),:,frame);
                Bond = PBC(CurrAtomPos2-CurrAtomPos1,box);
                Bond = Bond./sqrt(sum(Bond.^2));
                
                EFieldContribution(Vec) = Const.MapCoeff(Vec)*EFieldAtom(:,Vec,Osc,Mol,frame)'*Bond';
            end
            
            %Sum all map contributions to a specific oscillator
            output.FreqShift(Mol,Osc,frame) = sum(EFieldContribution);
        end
    end
    
end


%Build the CO atom index FOR SOLUTE ONLY
   for Osc = 1:size(in.COatomName,1);
        %build the index for current atom
        k = strcmp(SoluteAtoms,char(in.COatomName{Osc,1}));
        [OscIndAtom1(Osc,:),~] = find(k);
        %build the index for current atom
        k = strcmp(SoluteAtoms,char(in.COatomName{Osc,2}));
        [OscIndAtom2(Osc,:),~] = find(k);
   end

   %repackage the variables for output
for frame = 1:size(data.trjFrames.coords,3)
    box = data.trjFrames.trajectoryData.box(frame,:);
    for Osc = 1:size(in.COatomName,1)
        for Mol  = 1:data.topol.NumMolsInBox(1)
            
            CurrAtomPos0 = data.trjFrames.coords(OscIndAtom1(Osc,Mol),:,frame);
            CurrAtomPos1 = data.trjFrames.coords(OscIndAtom2(Osc,Mol),:,frame);
            CObond = PBC(CurrAtomPos1-CurrAtomPos0,box);
            CObondNorm = CObond./sqrt(sum(CObond.^2));
            output.COArray(:,Mol,Osc,frame) = CObondNorm;
            output.COcentArray(:,Mol,Osc,frame) = CurrAtomPos1;
            
        end
    end
end


