%Written by Seth Konig on 2/8/15
%Script attemtps to detect if Cortex lags after collecting large amounts of data
%using duration of ITI and reward period. The reward period has eye data 
% being collected during it while the ITI period does not. FYI This does not
% account for lag that the system cannot measure though, just lag in the
% amount of time it takes to compute something and the processor can still
% count. 

% % list of files for RED
% init = 'RR';
% cortexFiles = {'150123.3','150126.3','150127.3','150128.3','150129.3',...
%     '150130.3','150202.3','150204.3','150205.3'};
% itmFiles = {'411','412','413','415','425','426','428','429','430'};


% list of files for Tobii
init = 'TO';
cortexFiles = {'1503 127.1','150128.2','150129.1'};
itmFiles = {'416','417','418'};

iti_start_code = 15; %start of the inter trial interval
iti_end_code = 16;%end of the inter trial interval
reward_code = 3; %cortex delivered a reward pulse, usually 4-5 of them
eye_on_code = 100;%cortex starts collecting data
cross_on_code = 35;%cortex displays fixation cross 

iti_dur = NaN(length(itmFiles),1800);
reward_dur = NaN(length(itmFiles),1800);
eye_start_to_cross_on_dur = NaN(length(itmFiles),1800);

for i=1:length(itmFiles);
    datfil = strcat(init,cortexFiles{i});
    if strcmp(init,'IW')==1 || strcmp(init,'iw')==1
        datfil=['R:\Buffalo Lab\Cortex Data\Irwin\' datfil];
    elseif strcmp(init,'MP')==1 || strcmp(init,'mp')==1
        datfil=['R:\Buffalo Lab\Cortex Data\Peepers\' datfil];
    elseif strcmp(init,'WR')==1 || strcmp(init,'wr')==1
        datfil=['R:\Buffalo Lab\Cortex Data\Wilbur\' datfil];
    elseif strcmp(init,'TT')==1 || strcmp(init,'tt')==1
        datfil=['R:\Buffalo Lab\Cortex Data\Timmy\' datfil];
    elseif strcmp(init,'JN')==1 || strcmp(init,'jn')==1
        datfil=['R:\Buffalo Lab\Cortex Data\Guiseppe\' datfil];
    elseif strcmp(init,'TD')==1 || strcmp(init,'td')==1
        datfil=['R:\Buffalo Lab\Cortex Data\Theodore\' datfil];
    elseif strcmp(init,'PW')==1 || strcmp(init,'pw')==1  
        datfil=['R:\Buffalo Lab\Cortex Data\Vivian\' datfil];
    elseif strcmp(init,'RR')==1 || strcmp(init,'rr')==1  
        datfil=['R:\Buffalo Lab\Cortex Data\Red\' datfil];
    elseif strcmp(init,'TO')==1 || strcmp(init,'to')==1  
        datfil=['R:\Buffalo Lab\Cortex Data\Tobii\' datfil];
    end
    [time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata(datfil);
    
    for trial = 1:size(event_arr,2);
        iti_dur(i,trial) = time_arr(event_arr(:,trial) == iti_end_code,trial)...
            - time_arr(event_arr(:,trial) == iti_start_code,trial); 
        eye_start_to_cross_on_dur(i,trial) = time_arr(event_arr(:,trial) == cross_on_code,trial)...
            - time_arr(event_arr(:,trial) == eye_on_code,trial);
        
        rewards = find(event_arr(:,trial) == reward_code); %find if trial had a reward
        if ~isempty(rewards)
        reward_dur(i,trial)= time_arr(rewards(end),trial) - time_arr(rewards(1),trial); 
        end
    end
    
end

%% For fake data collected without a monkey, but added pupil data collection 


[time_arr,event_arr,eog_arr,epp_arr, header,trialcount] = get_ALLdata('PUPILRED.1');


iti_start_code = 15; %start of the inter trial interval
iti_end_code = 16;%end of the inter trial interval
reward_code = 3; %cortex delivered a reward pulse, usually 4-5 of them
eye_on_code = 100;%cortex starts collecting data
cross_on_code = 35;%cortex displays fixation cross 
stimulus_on_code = 23;%turn TL display on
stimulus_off_code = 24;%turn TL display off

iti_dur = NaN(1,760);
stimulus_dur = NaN(1,760);
eye_start_to_cross_on_dur =  NaN(1,760);

for trial = 1:size(event_arr,2);
    iti_dur(trial) = time_arr(event_arr(:,trial) == iti_end_code,trial)...
        - time_arr(event_arr(:,trial) == iti_start_code,trial);
    stimulus_dur(trial) = time_arr(event_arr(:,trial) == stimulus_off_code,trial)...
        - time_arr(event_arr(:,trial) == stimulus_on_code,trial);
    eye_start_to_cross_on_dur(trial) = time_arr(event_arr(:,trial) == cross_on_code,trial)...
        - time_arr(event_arr(:,trial) == eye_on_code,trial);
end
