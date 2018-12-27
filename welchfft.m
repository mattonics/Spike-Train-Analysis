function [f, spec] = welchfft(Spikes, dt, numBins, pOverlap, Padding)
%% WELCHFFT
%
% Calculates power spectrum by averaging Fast Fourier Transforms of
% overlapping window divisions
%
% INPUTS:
% Spikes (s): Times of spikes
% dt: Sampling rate for new sampled train
% numBins: Number of window divisions
% pOverlap: Percent overlap between windows
% Padding (optional): Number of zeros to zero-pad signal (increases frequency
% resolution)
% NOTE: Can change window type manually in the code (SEE LINE 69).
% OUTPUTS:
% f: Column vector of frequencies
% spec: Column vector of normalized spectrum values
%
% EXAMPLE values used in this paper:
% Spikes = (Spike times go here);
% dt = 0.001;
% numBins = 15;
% pOverlap = 50;
% Padding = 50;
% [f, spec] = welchfft(Spikes, dt, numBins, pOverlap, Padding);
% plot(f,spec);
%
% Reference: Welch (1967) The use of fast Fourier transform for the
% estimation of power spectra: a method based on time averaging over
% short, modified periodograms. IEEE Trans Audio Electroacoust 15:70–73.
%%
if nargin == 4
Padding = 0;
end
nTime = (Spikes(1):dt:Spikes(end)+dt); %New uniformly sampled times
nTime(2,:) = 0; %Initialize spike train
%Find equivalent times and set value at that time to 1
nTime(2,ismember(round(nTime(1,:).*(1/dt)).*dt, round(Spikes.*(1/dt)).*dt)) …
= 1;
%Spike train equals binary train (1,0)
spikeTrain = nTime(2,:);
N=length(spikeTrain); %Length of spike train
fs = 1/dt; %Sampling rate of new spike train
L = floor(N/(numBins-pOverlap/100×numBins+pOverlap/100)); %Number of points in
%length of window
Overlap_N = floor(L×pOverlap/100); %Number of points in each overlapping section
%Set endpoints of the windows
windowEndpoints = [linspace(1,1+(L-Overlap_N)×(numBins-1),numBins); …
L:L-Overlap_N:length(spikeTrain)+(L-Overlap_N)/2];
%Fix error from floor due to percent input
spikeTrain = [spikeTrain zeros(1,windowEndpoints(2,end)-length(spikeTrain))];
%Initalize windows
windows = zeros(numBins, 1);
%Compute fft for all windows
for i = 1:numBins
%Take current window from the spike train
currTrain = [spikeTrain(windowEndpoints(1,i):windowEndpoints(2,i)) …
zeros(1,Padding)];
%Can change window here by replacing line 71 with
%options below:
windowed = hanning(length(currTrain))'.*(currTrain);
%windowed = blackman(length(currTrain))'.*(currTrain);
%windowed = flattopwin(length(currTrain))'.*(currTrain);
%windowed = hamming(length(currTrain))'.*(currTrain);
%Take zero-mean, fft, and magnitude-squared
windows(i,1:length(windowed)) = abs(fft(windowed-mean(windowed))).^2;
%Normalize to the mean
windows(i,:) = windows(i,:)./mean(windows(i,:));
end
%Frequencies up to the nyquist
f = linspace(0,fs/2,floor(length(windows)/2))';
%Average results of windows of corresponding spectra values
spec = mean(windows(:,1:floor(length(windows)/2)))';
end