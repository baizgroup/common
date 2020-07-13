function FreqShift = ElectrostaticMapFunctionFast(field,in)

disp('Calculating Frequencies');
%Osc, Mol, frame
%FreqShift = zeros(size(field.EField,3),size(field.EField,4), size(field.EField,5));

%Reshape and replicatte map for simple element-by-element multiplication
PerMap = permute(in.MapCoeff,[3 2 1]);
MapTraj = repmat(PerMap, 1,1,1, size(field.EField,4), size(field.EField,5));

disp(['Map before application = [' num2str(MapTraj(1,:,1,7,1)) ']']);

%concatenate pot field and gradient
FieldTraj = cat(1,field.EPotential,field.EField);

disp(['Potential before application = [' num2str(FieldTraj(1,:,1,7,1)) ']']);

shift2display = FieldTraj(1,:,1,7,1).*MapTraj(1,:,1,7,1);

disp(['Shifts (cm^{-1}) = [' num2str(shift2display) ']']);

disp(['Atom Shifts (cm^{-1}) = [' num2str(sum(shift2display(1:3))) ']']);
disp(['Virtual Shifts (cm^{-1}) = [' num2str(sum(shift2display(4:end))) ']']);

%sum over xyz dimensions and atoms
freqTrajs = (sum(sum(MapTraj.*FieldTraj,1),2));

disp(['Sum of shifts (cm^{-1}) = ' num2str(freqTrajs(1,1,1,7,1))]);


%remove first two singleton dimensions
%FreqShift(Osc,Mol,frames);
FreqShift = permute(freqTrajs, [3 4 5 1 2]);

%%THIS WAS THE OLD VERSION
% %% field (:,Atom,Osc,Mol,frame)

% for frame = 1:size(field.EField,5)
%     for Mol  = 1:size(field.EField,4);
%         %disp(['Mol:' num2str(Mol)]);
%         for Osc = 1%1:size(field.EField,3)
%             ThisShift = 0;
%             for Atom = 1:size(field.EField,2)
%                 
%                 EsAtom = [field.EPotential(1,Atom,Osc,Mol,frame) field.EField(:,Atom,Osc,Mol,frame)' field.EGrad(:,Atom,Osc,Mol,frame)'];
%                 ThisShift = ThisShift + dot(squeeze(in.MapCoeff(Osc, Atom,:)),EsAtom);
%             end
%             FreqShift(Osc,Mol,frame) = ThisShift;
%             
%         end
%     end
% end