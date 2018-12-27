function [BinCenters, Dist] = ISIFreq(Spikes, Steps)
%% ISIFREQ help
%
% Returns a histogram of 1/ISI frequencies.
%
% INPUTS:
% Spikes: Times of spikes.
% Steps (1/units of Spikes): Steps (boundaries) for the histogram
% frequency count.
% OUTPUTS:
% BinCenters: Returns centers of the bins used for the histogram.
% Dist: Probability distribution of the ISI Frequencies
%
% EXAMPLE values used in this paper
% Spikes = (Spike times go here);
% Steps = 0:0.005:10;
% [BinCenters, Dist] = ISIFreq(Spikes, Steps);
% figure
% bar(BinCenters, Dist);
% title 'ISI Frequency Distribution'
% xlabel 'Frequencies'
% ylabel 'Distribution'
%
%% ISIFREQ
%Offset spikes by 1 and subtract for ISI's;
ISIs = Spikes(2:end) - Spikes(1:end-1);
%Calculate ISI Frequencies
ISIFreqs = 1./ISIs;
%Locates the bin centers of the steps
BinCenters = [(Steps(1:end-1) + Steps(2:end))./2 Steps(end)];
%Run histogram count
Counts = histc(ISIFreqs,Steps);
%Convert to probability distribution
Dist = Counts./sum(Counts);
end