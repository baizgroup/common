%units: coupling->cm-1 ;  dipole->Debye/(angstrom amu^(1/2))
% length -> angstrom
% amide-I transition dipole in NMA 3.466

function TDC = doTDC(geomVects)

dip.CO = 3.466*geomVects.CO;

%next calculate Rij(i,j,dim)
Rij = zeros([size(geomVects.COcent,1) size(geomVects.COcent,1) 3]);
for n=1:3;
    Vi=repmat(geomVects.COcent(:,n), [1 size(geomVects.COcent,1)]);
    Vj=repmat(geomVects.COcent(:,n)',[size(geomVects.COcent,1) 1]);
    
    Rij(:,:,n) = Vi-Vj;
end

%distance between units
Distij = sqrt(sum(Rij.^2,3));
TDC1 = 0*Distij;
TDC2 = 0*Distij;

temp(:,1,:) = dip.CO;
COall=repmat(temp,1,size(temp,1),1);
COallTr = permute(COall, [2 1 3]);
dipCOs = (COall.*COallTr);

TDC1f = sum(dipCOs,3)./(Distij.^3);
TDC2f = (sum(COall.*Rij,3).*sum(COallTr.*Rij,3))./(Distij.^5);

%old code(slow but correct)
% for n1=1:size(geomVects.COcent,1)
%     for n2=n1:size(geomVects.COcent,1)
%         %first term
%         TDC1(n1,n2) = (dip.CO(n1,:)*dip.CO(n2,:)') ./ ((Distij(n1,n2)).^3);
%         %second term
%         TDC2(n1,n2) = (dip.CO(n1,:)*squeeze(Rij(n1,n2,:)))*(dip.CO(n2,:)*squeeze(Rij(n1,n2,:)))...
%             ./((Distij(n1,n2)).^5);
%     end
% end
%this is the final TDC coupling in cm-1
TDC = (8.48E5/1650) * (9.79E-2 * (TDC1f - 3*TDC2f));

TDC(1:(size(TDC,1)+1):end)=0;%replace the diagonal elems with zeros



