function [peak,dist] = peakdetect(A, B, binsize)
%% PEAKDETECT help
%
% Determines the peak value of a distribution.
%
% INPUTS:
% A: x values of distribution.
% B: distribution values.
% binsize: Number of points to be included in a bin.
% OUTPUTS:
% peak: center of bin with the maximum distribution value.
% dist: highest average distribution value for a bin.
%
% EXAMPLE:
% A = (x values go here)
% B = (distribution values go here)
% binsize = 10
% [peak, dist] = peakdetect(A,B,binsize);
%
%% PEAKDETECT
%The number of points to extend forward or backward
%Step size of the distribution
dx = A(2)-A(1);
%Initialize integral array
integral = zeros(length(B)-binsize+1,1);
%Calculate the integrals
for i = 1:length(B)- binsize + 1
%Take the sum of the x values multiplied by y values in the bin
integral(i) = sum(dx×B(i:i+binsize-1));
end
%Find the index of the center of the bin with the max integral
j = find(integral==max(integral))+round(binsize/2);
%Find the value of the max integral
peak = A(j);
%Return the average dist around the max distribution value
dist = mean(B(j-round(binsize/2):j+round(binsize/2)));
end