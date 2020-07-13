function [mData] = mvgAvg (data, meanRange)
% functions calculates the moving average for a given dataset
%
%   params:
%   - data: the data to calculate the mvgAvg on
%   - meanRange: the range to the left and to the right from the current
%   value. e.g. beeing at value number 15 with range 5 values means that the mean is calculated between the values number 10 and 20
%   - mData: the moving average data   
%
% example for usage: 
%   a = rand (100,1); plot (a); hold on ; plot(mvgAvg(a, 5), 'r-');
%
% B. Knapp 2013-03-22

if meanRange > length(data)/2
   error ('mean range of %i is too big for data size %i', meanRange, length(a)); 
end

mData = zeros(size(data));

sum = 0;
for i=1:meanRange % sum up the right side for the first element
    sum = sum + data(i);
end
xIter = 1;
while (xIter <= meanRange + 1) % move as long as we are hitting the left boundary, no left-side removing necessary
   sum = sum + data(xIter+meanRange); 
   mData(xIter) = sum / (xIter+meanRange);
   xIter=xIter+1;
end
while (xIter < length(data)-meanRange+1) % now we are somewhere in the middle
   sum = sum-data(xIter-meanRange-1)+data(xIter+meanRange);
   mData(xIter) = sum / (2*meanRange+1);
   xIter=xIter+1;
end
while (xIter <= length(data)) % only removing on the left side while on the right no further addition is possible
    sum = sum - data(xIter-meanRange-1);
    mData(xIter) = sum / ((length(data)-xIter)+meanRange + 1);
    xIter=xIter+1;
end




