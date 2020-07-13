function EFieldAtom = EfieldFunction(AllAtomPos, CurrAtomPos, ThisCharges, ExcludeVec, box, ResidueVec,Const,Osc)
%Inputs: AllAtomPos CurrAtomPos Cutoff

CurrAtomPosMat = repmat(CurrAtomPos, [size(AllAtomPos,1) 1]);

%calculate the distances
xyzDiffVec = PBC(CurrAtomPosMat - AllAtomPos, box); %this is a super slow line. Should work to speed it up
rDiffDist = sqrt(sum(xyzDiffVec.^2,2));
xyzDiffVecNorm = xyzDiffVec./repmat(rDiffDist, [1 3]);
%exclude atoms within cutoff and same molecule

CutoffAtoms = ExcludeVec(Osc,:)' & ResidueVec;

%Electric field in units of Hartree/(borh electron charge)
EField = Const.AUConvField*(1./Const.FourPiEpsilonNaught)*repmat(1.60217657e-19.*ThisCharges./((1E-9*rDiffDist).^2),[1 3]).*(xyzDiffVecNorm);

%Kill the contribution for atoms within cutoff
EField(CutoffAtoms,:) = 0;


%EField(ResidueVec,:) = 0;

EFieldAtom = sum(EField);