function NLISITrain = NLISI(Spikes, p)
%%% NORMALIZED LOG ISI DISTRIBUTION
ISIs = Spikes(2:end)-Spikes(1:end-1); %Offset spikes by 1 and subtract 
%for ISI's
LogISIs = log(ISIs); %Take the log10 of ISI's
 
N = length(LogISIs); %N = number of ISIs
Q = max([20,floor(0.2*N)]); %Set window length as max of 20 and 20% of N
 
NLISITrain = zeros(1,length(LogISIs))'; %Initialize Normalized Log ISI Train
 
%Central Location of first 2*Q+1 ISIs
CentralLoc1 = ComputeCL(LogISIs(1:2*Q+1), p); 
%Subtract Central Location (Normalize)
NLISITrain(1:Q) = LogISIs(1:Q) - CentralLoc1; 
%Central Location of last 2*Q+1 ISIs
CentralLoc2 = ComputeCL(LogISIs(N-2*Q:N), p);
%Subtract Central Location (Normalize)
NLISITrain(N-Q+1:end) = LogISIs(N-Q+1:end) - CentralLoc2;
 
%For the middle portion
for i = Q+1:N-Q
   %Compute central location for portion of index +/- Q and subtract
   NLISITrain(i) = LogISIs(i) - ComputeCL(LogISIs(i-Q:i+Q), p);
end
 
%Get statistics of the NLISI train
% med = median(NLISITrain);
% pool_MAD = mad(NLISITrain);
% CentralDistBounds = [med - pool_MAD*2.58 med + pool_MAD*2.58];
 
% %Plot the NLISI Distribution
% figure
% hold on
% %Run a smoothing pdf kernel.
% NLISIpdf = pdf(fitdist(NLISITrain,'Kernel'), Steps);
% NLISIpdf = NLISIpdf./sum(NLISIpdf);
% plot(Steps, NLISIpdf,'g')
% %Plot treshold lines
% plot([CentralDistBounds(1) CentralDistBounds(1)], [0 max(NLISIpdf)], '--r')
% plot([CentralDistBounds(2) CentralDistBounds(2)], [0 max(NLISIpdf)], '--b')
%  
% xlabel 'Normalized Log ISIs'
% ylabel 'Probability'
% title 'Normalized Log ISI Distribution'
end