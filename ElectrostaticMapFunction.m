function FreqShift = ElectrostaticMapFunction(field,in)

disp('Calculating Frequencies');
%Osc, Mol, frame
FreqShift = zeros(size(field.EField,3),size(field.EField,4), size(field.EField,5));

for frame = 1:size(field.EField,5)
    for Mol  = 1:size(field.EField,4);
        %disp(['Mol:' num2str(Mol)]);
        for Osc = 1:size(field.EField,3)
            ThisShift = 0;
            for Atom = 1:size(field.EField,3)
                
                EsAtom = [field.EPotential(1,Atom,Osc,Mol,frame) field.EField(:,Atom,Osc,Mol,frame)' field.EGrad(:,Atom,Osc,Mol,frame)'];
                ThisShift = ThisShift + dot(squeeze(in.MapCoeff(Osc, Atom,:)),EsAtom);
            end
            FreqShift(Osc,Mol,frame) = ThisShift;
            
        end
    end
end
