function EPotAtom = EPotFunction(AllAtomPos, CurrAtomPos, ThisCharges, ExcludeVec, box, ResidueVec,Const,Osc)
%Inputs: AllAtomPos CurrAtomPos Cutoff

CurrAtomPosMat = repmat(CurrAtomPos, [size(AllAtomPos,1) 1]);

%calculate the distances
xyzDiffVec = PBC(CurrAtomPosMat - AllAtomPos, box);
rDiffDist = sqrt(sum(xyzDiffVec.^2,2));
%exclude atoms within cutoff and same molecule

CutoffAtoms = ExcludeVec(Osc,:)' & ResidueVec;

%Electric field in units of Hartree/(borh electron charge)
EPot = Const.AUConvPot*(1./Const.FourPiEpsilonNaught)*(1.60217657e-19.*ThisCharges./((1E-9*rDiffDist)));

%Kill the contribution for atoms within cutoff
%EPot(isinf(EPot) | isnan(EPot)) = 0;
EPot(CutoffAtoms,:) = 0;
EPotAtom = sum(EPot);