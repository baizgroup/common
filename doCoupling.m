function [CouplingMats Dipoles] = doCoupling(in,data)

%%Calculate CO Dipoles
data.Dipoles = zeros(3,size(in.atomVecName,2),size(data.field.index.OscInd,3),size(data.trjFrames.coords,3));
data.DipoleCenters = zeros(3,size(in.atomVecName,2),size(data.field.index.OscInd,3),size(data.trjFrames.coords,3));


disp('Coupling');
c=progress('init');

for frame=1:size(data.trjFrames.coords,3)
    
    c=progress(c,frame/size(data.trjFrames.coords,3));
    
    for Mol =1:size(data.field.index.OscInd,3)
        for Osc = 1:size(in.atomVecName,2)
            
            AtomPos1 = 10*data.trjFrames.coords(data.field.index.OscInd(Osc,1,Mol),:,frame);
            AtomPos2 = 10*data.trjFrames.coords(data.field.index.OscInd(Osc,2,Mol),:,frame);
            AtomPos3 = 10*data.trjFrames.coords(data.field.index.OscInd(Osc,3,Mol),:,frame);
            Bond1 = (AtomPos1-AtomPos2)./sqrt(sum((AtomPos1-AtomPos2).^2));
            Bond2 = (AtomPos3-AtomPos2)./sqrt(sum((AtomPos3-AtomPos2).^2));
            
            Mu = 0.94*(Bond1) + 0.23*(Bond2);
            Mu = Mu./sqrt(sum(Mu.^2));
            %disp('debug mode');
            data.Dipoles(:,Osc,Mol,frame) = Mu;
            data.DipoleCenters(:,Osc,Mol,frame) = AtomPos2;
            
        end
    end
end


%%arrange oscillators by molecule. First molecule, then oscillator
%%Lipid1, Osc1, Osc2, Lipid2, Osc1, Osc2;
geomVects.AllCO = zeros(3,size(in.atomVecName,2)*size(data.field.index.OscInd,3),size(data.trjFrames.coords,3));
geomVects.AllCOCent = zeros(3,size(in.atomVecName,2)*size(data.field.index.OscInd,3),size(data.trjFrames.coords,3));
CouplingMats = zeros(size(in.atomVecName,2)*size(data.field.index.OscInd,3),size(in.atomVecName,2)*size(data.field.index.OscInd,3),size(data.trjFrames.coords,3));

disp('doTDC');
c=progress('init');

if size(CouplingMats,1)>1
    for frame=1:size(data.trjFrames.coords,3)
        
        c=progress(c,frame/size(data.trjFrames.coords,3));
        
        p=1;
        for Mol =1:size(data.field.index.OscInd,3)
            for Osc = 1:size(in.atomVecName,2)
                
                geomVects.AllCO(:,p,frame) = data.Dipoles(:,1+mod(p,size(in.atomVecName,2)),Mol,frame);
                geomVects.AllCOCent(:,p,frame) = data.DipoleCenters(:,1+mod(p,size(in.atomVecName,2)),Mol,frame);
                
                p=p+1;
            end
        end
        
        geomVects.CO = squeeze(geomVects.AllCO(:,:,frame))';
        geomVects.COcent = squeeze(geomVects.AllCOCent(:,:,frame))';
        
        %if there is only one oscillator, do not do TDC
        if size(size(geomVects.COcent,1)~=1)
            CouplingMats(:,:,frame) = doTDC(geomVects);
        end
        
    end
    
    Dipoles = geomVects.AllCO;
else
    Dipoles(:,1,:) = squeeze(data.Dipoles);
end




