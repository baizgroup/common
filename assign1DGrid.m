
function outputmat = assign1DGrid(freq1,vals,InputMat,Freqax)

%first index
outsideRange =  freq1<Freqax(1) | freq1>Freqax(end);

freq1(outsideRange) = [];
vals(outsideRange) = [];

FF = repmat(freq1,size(InputMat));
VV = repmat(Freqax,size(vals));

[junk Inds] = min(abs(FF-VV)');

outputmat = InputMat;
if ~isempty(Inds)
outputmat(Inds) = outputmat(Inds) + vals';
end