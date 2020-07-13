function EGradAtom = EGradFunction(AllAtomPos, CurrAtomPos, ThisCharges, ExcludeVec, box, ResidueVec,Const,Osc)
%Inputs: AllAtomPos CurrAtomPos Cutoff

CurrAtomPosMat = repmat(CurrAtomPos, [size(AllAtomPos,1) 1]);

%calculate the distances
xyzDiffVec = PBC(CurrAtomPosMat - AllAtomPos, box);
rDiffDist = sqrt(sum(xyzDiffVec.^2,2));
%exclude atoms within cutoff and same molecule

if length(ExcludeVec)==1
CutoffAtoms = rDiffDist<ExcludeVec;
else
CutoffAtoms = ExcludeVec(Osc,:)' & ResidueVec;
end

%Electric field in units of Hartree/(borh electron charge)
% we included the conversion fron nm to meter (1E27)
cv =Const.AUConvGrad*(1./Const.FourPiEpsilonNaught)*1.60217657e-19*1E27;
x = xyzDiffVec(:,1);
y = xyzDiffVec(:,2);
z = xyzDiffVec(:,3);
xyzc = (x.^2 + y.^2 + z.^2).^(5/2);

EGrad(:,1) = (cv*ThisCharges').*(2*x.^2 - y.^2 - z.^2)./xyzc;
EGrad(:,2) = (cv*ThisCharges').*(3*x.*y)./xyzc;
EGrad(:,3) = (cv*ThisCharges').*(3*x.*z)./xyzc;
EGrad(:,4) = (cv*ThisCharges').*(3*x.*y)./xyzc;
EGrad(:,5) = (cv*ThisCharges').*(-x.^2 + 2*y.^2 - z.^2)./xyzc;
EGrad(:,6) = (cv*ThisCharges').*(3*z.*y)./xyzc;
EGrad(:,7) = (cv*ThisCharges').*(3*x.*z)./xyzc;
EGrad(:,8) = (cv*ThisCharges').*(3*y.*z)./xyzc;
EGrad(:,9) = (cv*ThisCharges').*(-x.^2 + -y.^2 + 2*z.^2)./xyzc;
%EGRAD 
% dxdx
% dxdy
% dxdz
% dydx
% dydy
% dydz
% dzdx
% dzdy
% dzdz

%Kill the contribution for atoms within cutoff
EGrad(CutoffAtoms,:) = 0;
EGradAtom = sum(EGrad);