function FreqShift = ElectrostaticMapFunctionFast_FSquared(field,in)

%Reshape and replicatte map for simple element-by-element multiplication
PerMapF = permute(in.MapCoeff(:,:,1:3),[3 2 1]);
MapTrajF = repmat(PerMapF, [1,1,1, size(field.EField,4), size(field.EField,5)]);

PerMapSq = permute(in.MapCoeff(:,:,4:6),[3 2 1]);
MapTrajSq = repmat(PerMapSq, [1,1,1, size(field.EField,4), size(field.EField,5)]);



%sum over xyz dimensions and atoms
freqTrajsF =  (sum(sum(MapTrajF.*field.EField,1),2));
freqTrajsSq = (sum(sum(MapTrajSq.*(field.EField.^2),1),2));

FreqShift = permute(freqTrajsF + freqTrajsSq, [3 4 5 1 2]);

