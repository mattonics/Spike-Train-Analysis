function [CentralLocation] = ComputeCL(LogISIs, p)
%%    COMPUTECL help
%
%    Subroutine required for Matlab code for RGS burst and pause detection
%   (Table 4). This subroutine computes the central location given an ISI train
%     using robust measures
%    of the median absolute difference (MAD), median, and central set.
%
%    INPUTS:
%        ISIs (sec): Lengths of ISIs in seconds.
%        Steps: Bin edges for histogram count (histc) of the ISIs.
%        p: Bottom and top p% used as outliers to calculate central
%        location; keep p in range [0.05, 0.30] (default 0.05).
%        NOTE: RGSDetect inputs log scale ISIs and Steps.
%    OUTPUTS:
%        CentralLocation: Central location of the ISI distribution
%
%    REFERENCE: Ko D, Wilson CJ, Lobb CJ, Paladini CA (2012) Detection of bursts and
%      pauses in spike trains. J Neurosci Methods 211:145-158
%
%%  COMPUTECL
 %
    B = sort(LogISIs);
    burstquantid = ceil(length(B)*p);
    pausequantid = floor(length(B)*(1-p));
    E_Center = (B(burstquantid) + B(pausequantid))/2;
    smad = mad(LogISIs-E_center, 1);
    lx = LogISIs(abs(LogISIs-E_Center)<=norminv(0.95)*smad);
    %Caclulates E-Center as average of 2 thresholds
    
    %Calculates C1 set using MAD
    C1Set = B((B > E_Center + norminv(p)*smad) & (B < E_Center + norminv(1-p)*smad));
    
    %Calculates median of C1 set
    C1 = median(C1Set);
    %Calculates central set using MAD
    CentralSet = B((B > C1 - 1.64*smad) & (B < C1 + 1.64*smad));
    %Calculates median of central set
    CentralLocation = median(CentralSet);
end
