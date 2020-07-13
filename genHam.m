function Ham = genHam(FreqShift, Dipoles, CouplingMats,in)



if size(CouplingMats,1)==1 %Single oscillator
    Ham.Ham = FreqShift+in.SiteEner;
    Ham.Evects = squeeze(sum(Dipoles.^2,1));
    Ham.Evals = FreqShift+in.SiteEner;
    Ham.TrDipoles = squeeze(sum(Dipoles.^2,1));
else
    Ham.Ham = 0*CouplingMats;
    Ham.Evects = 0*CouplingMats;
    Ham.Evals = zeros(size(CouplingMats,2),size(CouplingMats,3));
    Ham.TrDipoles = zeros(3,size(CouplingMats,2),size(CouplingMats,3));
    
    
    disp('genHam');
    c=progress('init'); 
    for frame=1:size(CouplingMats,3)
        
    c=progress(c,frame/size(CouplingMats,3)); 
        
        diagElems = FreqShift(:,:,frame);
        diagFreqs = in.SiteEner+reshape(diagElems, 1,[] );
        Ham.Ham(:,:,frame) =  diag(diagFreqs)+CouplingMats(:,:,frame);
        [Ham.Evects(:,:,frame) temp] = eig(Ham.Ham(:,:,frame));
        
        Ham.Evals(:,frame) = diag(temp);
        
        evects = squeeze(Ham.Evects(:,:,frame));
        moldipoles = squeeze(Dipoles(:,frame))';
        
        for n=1:size(CouplingMats,2)
            Ham.TrDipoles(:,n,frame) = [sum(moldipoles(:,1).*evects(:,n)) sum(moldipoles(:,2).*evects(:,n)) sum(moldipoles(:,3).*evects(:,n))];
        end
        
    end
end