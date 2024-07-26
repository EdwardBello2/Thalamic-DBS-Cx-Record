function [spikeburst] = burstanalysis(all_sp)

% function [spikeburst] = burstanalysis(all_sp)
%
% Modified from the function 'findexcitation' written by Gary Russo as part of the SpikeAnalyze program.  
% The fuction was a modification of the algorithm originally published by Legendy and Salcman 1985
% and later modified by Hanes et al 1995.
%
% The probability that a spike is part of a burst is defined using a possion surprise index that
% determines whether or not the occurance of a spike would be expected given a poisson distribution.
% The algorithm finds the time where the response is highest and returns its begin time and end time.  However,
% the spike rate may still be significantly above baseline before the onset time as when there is some
% sort of delay-period activity.  Thus a beginning of activation time (boa_time) is returned, where
% the spike train before a response first goes below sig_level.
%
% Input:   base_rate - baseline spike rate of the neuron (spikes/ms)
%          sp_time  - a list of spike times
%          MaxNumberExtraSpikes - the maximum number of spikes allowed in the search process where the surprise index does not increase
%          MaxExtraTime - maximum time allowed that does not raise the poisson surprise index
%          MinBurstSpikes - the minimum number of spikes in a burst that is considered significant
%          sig_level - the significance level of the burst
%          verbose_flag (optional) - flag to make function print zillions of messages as to exactly
%          what it is doing.  Mostly for debugging (0 = off, 1 = on).  Default is off.
%
% Output:  bor_time - beginning of response relative to start_time
%          boa_time - beginning of activation relative to start_time
%          eor_time - end of response relative to start_time
%          P_response  -  poisson probability of chance occurence
%          status - 1 if a complete burst was found, 0 otherwise
%          
%
% Note: when computing surprise index, spikes are counted by including the first spike and excluding the last in a burst
%   (Legendy and Salcman 1985)
%
% The following is from Tom Lane at Mathworks:
% The definition of poisscdf is just p = sum(poisspdf(0:x,m))
% If we did the sum up to x=Infinity we would get 1.0, so loosely speaking
%    1-p = sum(poisspdf(x+1:Inf,m))
%
% In practice you just have to figure out what upper bound to use in place of
% Inf.  For example, using your distribution x=33 gives the value 1-p=0
% because p is smaller than eps.  You could sum up values from 34 to 50 (say)
% to avoid the truncation caused by the subtraction, and you would get
%   sum(poisspdf(34:50,5))
%   ans =
%        1.54867851062915e-017
%
% It's hard for me to give detailed advice about your algorithm without
% knowing much about it.  Perhaps you could use a rule such as
%
%     p = 1 - poisscdf(x,m)
%     if (p<0.01)
%        p = sum(poisspdf(x+1:xmax, m))
%     end
%
% The trick is to decide what to use in place of xmax.  If execution time is
% not much of an issue, you could compute the probabilities one by one, and
% sum them up until the ratio of the next probability to the sum so far is
% very small.  Or you could estimate the number of terms you need using the
% fact that the standard deviation of the distribution is sqrt(m).

%Burst searching parameters 
%*****CHANGE AS NEEDED*****
BurstMaxSearchSpikes = 3;
BurstMinSpikes = 3;
BurstMaxGapTime = 10;
BurstSigLevel = 0.05;
BurstMaxTime = 150;    %Use to limit bursts to short intervals
BurstMinRate = 100;    %Use to limit bursts to short intervals

%Calculate the mean rate to use for burst detection
all_sp = all_sp*1e3; % all spikes should be in ms!
s = all_sp; status = 1; bob = []; eob = []; all_burst_rates = []; all_burst_sp_number = [];
burst_detection_mean_rate = (length(s)-1)/s(end); %Mean rate of entire data stream (not counting the first spike), in spikes/ms

%Parse through the spike data looking for bursts
h_wait = waitbar(0,'Finding bursts.  Please wait...');
while ~isempty(s) && status
    
    waitbar((length(all_sp)-length(s))/length(all_sp));
    
    %Call the findbursts function to sequentially find spikes in putative bursts
    [bob_time,eob_time, boa_time, P, status] = findbursts(burst_detection_mean_rate, s, ...
        BurstMaxSearchSpikes, BurstMaxGapTime, BurstMinSpikes, BurstSigLevel); %Mean rate in seconds, but function wants spikes/ms

    if status

        %Remove burst spikes if the duration is longer than BurstMaxTime
        if eob_time-bob_time > BurstMaxTime
            s(s <= eob_time) = [];
            continue
        end

        %Remove burst spikes if the burst rate is less than BurstMinRate
        current_sp_in_burst = length(find(all_sp >= bob_time & all_sp <= eob_time)); %Include spikes at beginning AND end of burst
        b_rate = current_sp_in_burst/(eob_time-bob_time)*1000; %Spikes/seconds within the burst
        if b_rate < BurstMinRate
            s(s <= eob_time) = [];
            continue
        end
        
        bob = [bob; bob_time];
        eob = [eob; eob_time];
        all_burst_rates = [all_burst_rates; b_rate];
        all_burst_sp_number = [all_burst_sp_number; current_sp_in_burst]; %Number of spikes in each burst
        s(s <= eob_time) = []; %remove burst spikes and continue finding bursts
        
    end
    
end

bob_times = bob;
eob_times = eob;
close(h_wait);

% Define output parameters
sp = all_sp(all_sp>0);
record_duration = all_sp(end)-all_sp(1);
bursts_to_keep = bob_times>0 & eob_times>0;
burst_start = bob_times(bursts_to_keep);
burst_end = eob_times(bursts_to_keep);
burst_rates = all_burst_rates(bursts_to_keep);
burst_sp_number = all_burst_sp_number(bursts_to_keep);
if ~isempty(burst_start)
    burst_duration = mean(burst_end - burst_start);
    burst_rate = mean(burst_rates);
    burst_freq = length(burst_rates)*1000 * 60 / record_duration; %Bursts per minute.  Not corrected for removed data because individual bursts can span removed sections
    percent_spikes_in_bursts = sum(burst_sp_number) / length(sp) * 100;
else
    burst_duration = 0;
    burst_rate = 0;
    burst_freq = 0;
    percent_spikes_in_bursts = 0;
end

% Reassign the parameters into an output structure
spikeburst.rate = burst_detection_mean_rate*1e3;
spikeburst.burst_duration = burst_duration;
spikeburst.burst_rate = burst_rate;
spikeburst.burst_freq = burst_freq;
spikeburst.percent_bursts = percent_spikes_in_bursts;
spikeburst.bob_times = bob_times;
spikeburst.eob_times = eob_times;

% end of function




%**************************************************************************
%*********************** Find bursts function  ****************************
%**************************************************************************
function [bor_time, eor_time, boa_time, P_response, status] = findbursts(base_rate, sp_time, MaxNumberExtraSpikes,...
   MaxExtraTime, MinBurstSpikes, sig_level)

    %Start out assuming the worst.  Will return these if a burst was not found.
    status = 0;
    bor_time = NaN;
    eor_time = NaN;
    boa_time = NaN;
    P_response = NaN;

    sp_time = sort(sp_time);

    if isempty(sp_time)
        disp('Warning: No data to search.');
        return
    end

    %Add .1 to all spike times that are duplicates (when analyzing collapsed data) so no isi is zero
    dup_count = 0;
    while any(~diff(sp_time)) %as long as there are duplicates.  May need to do more than once if triplets, etc
       dups = find(~diff(sp_time)) + 1; %indexes to duplicate spike times NOT including first one
       dup_count = dup_count + length(dups);
       sp_time(dups) = sp_time(dups - 1) + .1; %Make it .1 plus previous
       sp_time = sort(sp_time);
    end
    if dup_count > 0
       disp('%1.0f simultaneous spikes',dup_count);
    end

    if isempty(sp_time) || length(sp_time) < MinBurstSpikes
       disp('Warning: Insufficient spikes for finding neural events.');
       return
    end
    isi = diff(sp_time);  %interspike intervals

    ExtraSpikes = 0;
    ExtraTime = 0;

    %Find putative bors by computing vector of 1st two isi averages and find elements greater than avg_rate
    running_average_isi = (isi(1:end-1) + isi(2:end))./2;       %average of 2 consecutive isi's, in ms
    running_average_rates = 1./running_average_isi;             %average rate of each 3 consecutive spikes across entire spike train in spikes/ms
    putative_bor = find(running_average_rates > base_rate);     %indexes to spikes in sp_time that are at the beginning of 2 consecutive ISIs that are bor_isi_factor times above the average baseline rate

    if isempty(putative_bor)
        disp('No putative bursts found');
        return
    end
    putative_bor_times = sp_time(putative_bor); %Times of putative BORs.  Used to find the next one in the while loop if the first isn't significant 
    current_putative_bor = putative_bor(1);

    while current_putative_bor ~= putative_bor(end)
       pdf_flag = 0; %start examining the burst using cdf distribution

       %Compute SI's of current spike train
       current_sp_time = sp_time(current_putative_bor:end);     %Spike times from search point on
       t = (current_sp_time(2:end) - current_sp_time(1));       %time intervals of 2nd spike on relative to the first
       num_spikes = 1:length(current_sp_time) - 1;              %number of spikes in each interval t, not counting the first
       num_spikes = num_spikes';

       %initializations
       bor_index = current_putative_bor; %First spike of current_sp_time
       eor_index = bor_index + 1;        %Third spike of current_sp_time
       boa_index = bor_index;            %Index to the beginning of activation

       [SI_max,P,pdf_flag] = burstsurprise(t(1),num_spikes(1),base_rate,pdf_flag); %First SI of the first time interval (second spike)
       for i = 2:length(t) %Start comparing SI on the 2nd time interval (3rd spike)
          [SI,P,pdf_flag] = burstsurprise(t(i),num_spikes(i),base_rate,pdf_flag);
          if boa_index == bor_index && P <= sig_level
             boa_index = boa_index + i + 1;
          end

          if SI >= SI_max
             %If poisscdf saturates and gives this outrageously huge SI, count it as an
             %increasing SI, but don't use this ridiculous number as the max SI
             if SI < 744
                SI_max = SI;
             else
                SI_max = SI_max + 1;
             end
             eor_index = bor_index + i;  %index of current eor
             ExtraSpikes = 0;            %reset extra spike cound
             ExtraTime = 0;              %reset extra time found
          else
             ExtraSpikes = ExtraSpikes + 1; %add to extra spike count
             ExtraTime = ExtraTime + current_sp_time(i + 1) - current_sp_time(i); %add to extra time count current search spike time minus current eor_index time 
          end

          %if the extra spike count is exceeded, the eor is found
          if ExtraSpikes >= MaxNumberExtraSpikes
             break
          end

          %if the extra time count is exceeded, the eor is found
          if ExtraTime >= MaxExtraTime
             break
          end
          
       end

       %Now compute SI of burst each time after taking a spike off the beginning and find max   
       ExtraSpikes = 0;
       ExtraTime = 0;
       test_bors = (bor_index:eor_index - 1)'; %indexes to potential new bors
       t = (sp_time(eor_index) - sp_time(test_bors)); %Time of the bursts with those potential new bors

       %Below is a fix for a special case when collapsing trials.  If two spikes are at the same ms, then add .1 to
       %one of them.  But if this .1 is at the eob, then this last ISI will have a very large SI, and the
       %algorithm will remove all spikes from the beginning to this point in an effort to maximize the SI.
       %Then of course the burst is only two spikes and will be rejected every time.  So round up to 1 ms to
       %remove this bogusly short time interval.
       t(end) = ceil(t(end));

       num_spikes = eor_index - test_bors;
       pdf_flag = 0; %try again using cdf distribution
       [SI,P,pdf_flag] = burstsurprise(t(1),num_spikes(1),base_rate,pdf_flag); %SI of current burst

       remove_spikes = 0;
       for i = 2:length(t)
          [SI_bor,P,pdf_flag] = burstsurprise(t(i), num_spikes(i), base_rate, pdf_flag);
          if SI_bor >= SI
             SI = SI_bor; %Make SI the maximum possible for the burst
             remove_spikes = i - 1; % t = 2 is best, remove one spike
          else
             % do nothing
          end
       end

       bor_index = bor_index + remove_spikes;  %New bor index that maximizes SI
       eor_time = sp_time(eor_index);
       bor_time = sp_time(bor_index);
       if bor_index < boa_index %In some cases, the optimize response start is earlier than the first found increase above sig_level.  In those cases, boa = bor
          boa_index = bor_index;
       end
       boa_time = sp_time(boa_index);

       P_response = 1 - poisscdf(eor_index - bor_index, base_rate * (eor_time  - bor_time));

       %If a found response is not significant, the while loop iterates again, but starting at the spike after the nonsignificant eor
       if P_response > sig_level || eor_index - bor_index + 1 < MinBurstSpikes

          %if a response is not significant, find putative_bor greater than the end of the last burst found
          next_bor = find(putative_bor_times > sp_time(eor_index)); %Indexes to all the putative BORs with times after the end of the last response attempt

          if ~isempty(next_bor)
             current_putative_bor = putative_bor(next_bor(1)); %The first one.  Will start searching here on next iteration of the while loop
          else
             current_putative_bor = putative_bor(end);  %This will disable the while loop
          end 

       else
          status = 1;
          break %Found acceptable response
       end

       %If the while loop ran out, indicating that it got to the end of the isi list
       bor_time = NaN;
       eor_time = NaN;
       boa_time = NaN;
       P_response = NaN;
       status = 0;
    end

% end of function


%**************************************************************************
%******************** Surprise index function  ****************************
%**************************************************************************
function [SI,P,pdf_flag] = burstsurprise(t, num_spike, base_rate, pdf_flag)

    lambda = base_rate * t; %lambda, the expected number of spikes based on the average rate (base_rate) in the given interval (t)

    if pdf_flag
       P = poisspdf(num_spike,lambda);
    else
       P = 1 - poisscdf(num_spike - 1,lambda);  %Corrected by Phil Hahn
    end

    %If a huge burst with very little baseline, poissoncdf saturates to zero.
    %When P is zero, will make it the smallest number possible in matlab.  But this can be
    %a problem when evaluating a burst's bor.  If the burst is huge, then it
    %is saturated as you remove spikes, even if the spikes should not be part of the burst!

    if ~P
       P = poisspdf(num_spike,lambda);
       pdf_flag = 1;
    end

    P(~P) = exp(-745); %If any zeros left, will get an error when taking the log, so just make it the smallest number possible and hope for the best.
    SI = -log(P);

% end of function