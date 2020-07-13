function IR = gen1DSpec(Ham,in)
IRLines = zeros(size(in.FreqAx));


    disp('gen1DSpec');
    c=progress('init');
    
for frame=1:size(Ham.Ham,3)
    
    c=progress(c,frame/size(Ham.Ham,3)); 
    
    if size(Ham.TrDipoles,2)>1
        Dip = squeeze(Ham.TrDipoles(:,:,frame))';
    else
        Dip = Ham.TrDipoles(frame);
    end
    Eigenvals = Ham.Evals(:,frame);
    
    dipStrength = sum(Dip.^2,2);
    IRLines = assign1DGrid(Eigenvals,dipStrength,IRLines,in.FreqAx);
end

IR =  conv(IRLines,voight(in.FreqAx,in.FreqAx(floor(end/2)),3,0),'same');