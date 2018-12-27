function [Bursts, Pauses] = RGSAlgo(NLISIs, Spikes, N_min, alpha)
%% BURST AND PAUSE STRING DETECTION
%Get index and ISI lengths of all ISIs that satisfy burst threshold
NLISITrain = NLISIs;
Burst_Thresh = CentralDistBounds(1);
BurstINDXS = 1:length(NLISITrain);
BurstINDXS(NLISITrain >= Burst_Thresh) = []; %Delete all indexes greater than
%the burst threshold
if ~isempty(BurstINDXS)
    %Matrix of all potential burst ISIs and their indexes
    BurstsM = [NLISITrain(NLISITrain < Burst_Thresh)';BurstINDXS];
    [~,c] = size(BurstsM);
    Burst_Seed = mat2cell(BurstsM,2,ones(1,c,1));
else
    Burst_Seed = {};
    Bursts = Burst_Seed;
end
 
%Loop through each potential burst ISI (Burst Seed)
for i = 1:length(Burst_Seed)
    %Go forwards and backwards from the current burst until both conditions
    %are unsatisfied
    forwards = 1;
    backwards = 1;
    while forwards || backwards
        currBurst = cell2mat(Burst_Seed(i)); 
        %Go forwards 1 ISI
        if forwards
            %Set current ISI as end of the current burst
            currSpike = currBurst(:,end);
            if currSpike(2) ~= length(NLISITrain)
                %q is number of spikes
                [~,q] = size(currBurst);
                %P1 is probability burst will occur assuming Gaussian
                %distribution with mean, mu*q, and std, sqrt(q)*sigma
                P1 = normcdf(sum(currBurst(1,:)), mu*q, sqrt(q).*sigma);
                testBurst = [currBurst [NLISITrain(currSpike(2)+1);...
                    currSpike(2)+1]];
                %P2 is the same probability with the next ISI added to the
                %burst
                P2 = normcdf(sum(testBurst(1,:)), mu*(q+1), ...
                    sqrt(q+1).*sigma);
                %If the next ISI increased the probability of the burst
                %occurring
                if P2 >= P1
                    %Stop going forwards
                    forwards = 0;
                else
                    %Otherwise, set the current burst seed to the tested
                    %burst
                    Burst_Seed{i} = testBurst;
                end
            else
                %Stop going forwards if at the end of the ISI train
                forwards = 0;
            end
        end
        
        currBurst = cell2mat(Burst_Seed(i)); 
        %Go backwards 1 ISI
        if backwards
            %Set current ISI as end of the current burst
            currSpike = currBurst(:,1);
            if currSpike(2) ~= 1
                %q is number of spikes
                [~,q] = size(currBurst);
                %P1 is probability burst will occur assuming Gaussian
                %distribution with mean, mu*q, and std, sqrt(q)*sigma
                P1 = normcdf(sum(currBurst(1,:)), mu*q, sqrt(q).*sigma);
                testBurst = [[NLISITrain(currSpike(2)-1);currSpike(2)-1] ...
                    currBurst];
                %P2 is the same probability with the next ISI added to the
                %burst
                P2 = normcdf(sum(testBurst(1,:)), mu*(q+1), ...
                    sqrt(q+1).*sigma);
                %If the next ISI increased the probability of the burst
                %occurring
                if P2 >= P1
                    %Stop going backwards
                    backwards = 0;
                else
                    %Otherwise, set the current burst seed to the tested
                    %burst
                    Burst_Seed{i} = testBurst;
                end
            else
                %Stop going backwards if at the end of the ISI train
                backwards = 0;
            end
        end
    end
end
 
Bursts.BurstingSpikes = [];
Bursts.Windows = [];
Bursts.NumSpikes = [];
Bursts.IBF = [];

if ~isempty(Burst_Seed)
    %Initialize BurstInfo
    BurstInfo = zeros(length(Burst_Seed),3);
    %Get start index of each burst
    BurstInfo(:,1) = cellfun(@(x) x(2,1),Burst_Seed);
    %Get end index of each burst
    BurstInfo(:,2) = cellfun(@(x) x(2,end),Burst_Seed);
    %Get P-value of each burst (probability of occurence assuming Gaussian
    %distribution)
    BurstInfo(:,3) = cellfun(@(x) normcdf(sum(x(1,:)), mu*length(x), ...
        sqrt(length(x)).*sigma),Burst_Seed);
    %Filter out bursts less than minimum number of spikes specified by N_min
    BurstInfo(BurstInfo(:,2)-BurstInfo(:,1)+2 < N_min,:) = [];
    
    %Filter out overlapping bursts
    no_overlap = 0;
    i=1;
    if ~isempty(BurstInfo)
        [r,~] = size(BurstInfo);
        if r ~= 1
            while ~no_overlap
                %If the indexes of the burst ISIs don't intersect
                if isempty(intersect(BurstInfo(i,1):BurstInfo(i,2),...
                    BurstInfo(i+1,1):BurstInfo(i+1,2)))
                    %move to the next burst
                    i = i+1;
                else
                    %If they intersect, choose the burst with the lower P
                    %value
                    if BurstInfo(i,3) <= BurstInfo(i+1,3)
                        BurstInfo(i+1,:) = [];
                    else
                        BurstInfo(i,:) = [];
                    end
                end
                %When the end is reached, stop
                [r,~] = size(BurstInfo);
                if i == r
                    no_overlap = 1;
                end
            end
        end
        
        %Bonferroni correction for false positives
        KB = length(find(BurstInfo(:,3) < alpha));
        BurstInfo(BurstInfo(:,3)*KB >= alpha,:) = [];

        %r is the number of rows or the number of bursts
        [r,~] = size(BurstInfo);
        %for each burst, append the burst spikes
        for i = 1:r
            Bursts.BurstingSpikes = [Bursts.BurstingSpikes;...
                Spikes(BurstInfo(i,1):BurstInfo(i,2)+1)];
        end
        
        %Use the indexes in burst info to find the burst windows
        Bursts.Windows = [Spikes(BurstInfo(:,1)) Spikes(BurstInfo(:,2)+1)];
        %Use the indexes to find the number of spikes in each burst
        Bursts.NumSpikes = BurstInfo(:,2) - BurstInfo(:,1) + 2;
        %Use the number of spikes and windows to calculate the IBF
        Bursts.IBF = Bursts.NumSpikes./(Bursts.Windows(:,2) - ...
        Bursts.Windows(:,1));
    end
end
 
%Get index and ISI lengths of all NLISIs that satisfy pause threshold
Pause_Thresh = CentralDistBounds(2);
PauseINDXS = 1:length(NLISITrain);
%Delete all indexes less than the pause threshold
PauseINDXS(NLISITrain <= Pause_Thresh) = [];
 
if ~isempty(PauseINDXS)
    %Matrix of all potential pause string NLISIs and their indexes
    PausesM = [NLISITrain(NLISITrain > Pause_Thresh)';PauseINDXS];
    [~,c] = size(PausesM);
    Pause_Seed = mat2cell(PausesM,2,ones(1,c,1));
else
    Pause_Seed = {};
    Pauses = [];
end
%Loop through each potential pause string NLISI (Pause Seed)
for i = 1:length(Pause_Seed)
    %Go forwards and backwards from the current pause string until both conditions
    %are unsatisfied
    forwards = 1;
    backwards = 1;
    while forwards || backwards
        currPause = cell2mat(Pause_Seed(i));
        %Go forwards 1 ISI
        if forwards
            %Set current ISI as end of the current pause string
            currPauseind = currPause(:,end);
            if currPauseind(2) ~= length(NLISITrain)
                [~,q] = size(currPause);
                %P1 is probability pause string will occur assuming Gaussian
                %distribution with mean, mu*q, and std, sqrt(q)*sigma
                P1 = (1-normcdf(sum(currPause(1,:)), mu*q, sqrt(q).*sigma));
                testPause = [currPause [NLISITrain(currPauseind(2)+1);...
                    currPauseind(2)+1]];
                %P2 is the same probability with the next ISI added to the
                %pause string
                P2 = (1-normcdf(sum(testPause(1,:)), mu*(q+1), ...
                    sqrt(q+1).*sigma));
                %If the next ISI increased the probability of the pause
                %string occurring
                if P2 >= P1
                    %Stop going forwards
                    forwards = 0;
                else
                    %Otherwise, set the current pause seed to the tested
                    %pause string
                    Pause_Seed{i} = testPause;
                end
            else
                forwards = 0;
            end
        end
        currPause = cell2mat(Pause_Seed(i));
        %Go backwards 1 ISI
        if backwards
            currPauseind = currPause(:,1);
            if currPauseind(2) ~= 1
                [~,q] = size(currPause);
                %P1 is probability pause string will occur assuming Gaussian
                %distribution with mean, mu*q, and std, sqrt(q)*sigma
                P1 = (1-normcdf(sum(currPause(1,:)), mu*q, sqrt(q).*sigma));
                testPause = [[NLISITrain(currPauseind(2)-1);...
                    currPauseind(2)-1] currPause];
                %P2 is the same probability with the next ISI added to the
                %pause string
                P2 = (1-normcdf(sum(currPause(1,:)), mu*(q+1), ...
                    sqrt(q+1).*sigma));
                %If the next ISI increased the probability of the pause
                %string occurring
                if P2 >= P1
                    %Stop going forwards
                    backwards = 0;
                else
                    %Otherwise, set the current pause seed to the tested
                    %pause string
                    Pause_Seed{i} = testPause;
                end
            else
                backwards = 0;
            end
        end
    end
end
 
Pauses.Windows = [];
Pauses.NumSpikes = [];
Pauses.IPF = [];
Pauses.PausingSpikes = [];
Pauses.AllSpikes = [];
Pauses.AllLengths = [];

if ~isempty(Pause_Seed)
    %Initialize PauseInfo variable
    PauseInfo = zeros(length(Pause_Seed),3);
    %Starting indexes of pause strings
    PauseInfo(:,1) = cellfun(@(x) x(2,1),Pause_Seed);
    %Ending indexes of pause strings
    PauseInfo(:,2) = cellfun(@(x) x(2,end),Pause_Seed);
    %P-value of the pause strings
    PauseInfo(:,3) = cellfun(@(x) (1-normcdf(sum(x(1,:)), mu*length(x), ...
        sqrt(length(x)).*sigma)),Pause_Seed);
    %Minimum number of spikes filter
    PauseInfo(PauseInfo(:,2)-PauseInfo(:,1) + 2 < N_min,:) = [];
    
    %Filter out overlaps
    no_overlap = 0;
    i=1;
    if ~isempty(PauseInfo)
        %r is number of current pause strings
        [r,~] = size(PauseInfo);
        if r ~= 1
            while ~no_overlap
                %If the indexes of the burst ISIs don't intersect
                if isempty(intersect(PauseInfo(i,1):PauseInfo(i,2),...
                    PauseInfo(i+1,1):PauseInfo(i+1,2)))
                    %Move to next pause string
                    i = i+1;
                else
                    %Choose the pause string with lower P-value
                    if PauseInfo(i,3) <= PauseInfo(i+1,3)
                        PauseInfo(i+1,:) = [];
                    else
                        PauseInfo(i,:) = [];
                    end
                end
                %End if the last pause string is reached
                [r,~] = size(PauseInfo);
                if i == r
                    no_overlap = 1;
                end
            end
        end
        %Bonferroni correction
        KP = length(find(PauseInfo(:,3) < alpha));
        PauseInfo(PauseInfo(:,3)*KP >= alpha,:) = [];

        %Use indexes to find start and end times
        Pauses.Windows = [Spikes(PauseInfo(:,1)) Spikes(PauseInfo(:,2)+1)];
        %Use indexes to determine number of spikes
        Pauses.NumSpikes = PauseInfo(:,2) - PauseInfo(:,1) + 2;
        %Use windows and numspikes to calculate IPF
        Pauses.IPF = Pauses.NumSpikes./(Pauses.Windows(:,2) - ...
        Pauses.Windows(:,1));

        %Add pausing spikes from each pause string to the pausingspikes element
        [r,~] = size(PauseInfo);
        for i = 1:r
            Pauses.PausingSpikes = [Pauses.PausingSpikes;...
                Spikes(PauseInfo(i,1):PauseInfo(i,2)+1)];
        end

        %Get all pauses using the pause indexes
        Pauses.AllSpikes = Spikes(PauseINDXS);
        %Get the lengths of all the pauses that satisfy the threshold
        Pauses.AllLengths = ISIs(PauseINDXS);
    end
end