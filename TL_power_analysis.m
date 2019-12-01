% Code written by Seth Konig 2/15/15. Based mainly on Analyze_All_TLs_Data.m
% but code does power analysis at the end. Code removes trials with reaction
% times less than 150 ms because these are likley false positives. Trials
% with reaction times > 5000 ms are also removed because the monkey is
% probably not paying attention. Code looks at the 1st max_desired_blocks blocks of data
% for individual monkeys and across all monkeys. Below is a list of 100
% potentail sessions for analysis. Sessions with less than max_desired_blocks blocks are
% removed (~ 4-5 out of the 100). Also, sessions in which more than 10% of
% trials (96/(max_desired_blocks*24 = 960) are removed from the power analysis as well.
%
% [1] Import task data and detect fixations and saccades.
% [2] Analyze reaction times and number of saccades to find target.
% [3] Compute power analysis for individual monkeys.
% [4] Compute power anlaysis across monkeys.
%     power analysis assumes novel contexts are null distribution (P0) and
%     repeat contexts are the alternative distribution.

%%
%---[1] Import task data and detect fixations and saccades---%

%data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\TL_relationalmemory\Eyedat\';
%data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\TL_relationalmemory\Eyedat\Power Analysis\';%for sets ran power analysis on
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\TL_relationalmemory\Eyedat\'; %where to find preprocessed data

%cortex data files for TLs
% TL_files = {
%     'MP101216.2','MP110204.2','MP110215.2','MP110217.4','MP110218.2',...
%     'MP110224.2','MP110228.2','MP110301.2','MP110303.3','MP110308.2',...
%     'MP110314.2','MP110405.2','MP110425.2','MP110426.2','MP110427.2',...
%     'MP110428.2','MP110502.2','MP110509.2','MP110510.2','MP110511.2',...
%     'MP110513.2','MP110516.2','MP110518.2',...
%     'PW130410.2','PW130411.2','PW130412.2','PW130415.3','PW130416.2',...
%     'PW130417.3','PW130418.3','PW130422.3','PW130424.3','PW130425.3',...
%     'PW130426.3','PW130429.3','PW130501.3','PW130502.4','PW130506.3',...
%     'PW130507.3','PW130508.4','PW130510.3','PW130513.3','PW130515.3',...
%     'PW130516.3','PW130520.3','PW130521.3','PW130522.3',...
%     'JN110908.1','JN110909.1','JN110913.1','JN110914.1','JN110919.1',...
%     'JN110921.2','JN110926.1','JN110928.1','JN111005.2','JN111010.2',...
%     'JN111012.2','JN111019.2','JN111024.2','JN111026.2','JN111027.2',...
%     'JN111031.2','JN111101.2','JN111107.2','JN111108.2','JN111109.2',...
%     'JN111111.2','JN111116.2','JN111121.2',...
%     'IW110816.1','IW110819.1','IW110823.1','IW110829.4','IW110830.1',...
%     'IW110906.1','IW110907.1','IW110909.1','IW110912.1','IW110921.1',...
%     'IW110926.1','IW110928.1','IW111010.1','IW111012.1','IW111019.1',...
%     'IW111021.1','IW111026.1','IW111107.1','IW111110.1','IW111208.1',...
%     'IW120124.1','IW120203.1','IW120221.1','IW120223.1','IW120320.1'
%    };
% 
% %which set e.g. tl425.itm
% item_num = [
%     425 437 442 444 445, ... %MP's sets
%     449 451 452 454 455, ... %MP's sets
%     459 470 481 482 483, ... %MP's sets
%     484 485 489 490 491, ... %MP's sets
%     493 494 496,...          %MP's sets
%     421 422 423 424 425, ... %PW's sets
%     426 428 430 432 433, ... %PW's sets
%     434 435 437 438 440, ... %PW's sets
%     441 442 444 445 447, ... %PW's sets
%     448 449 451 450,...      %PW's sets
%     429 430 431 432 434, ... %JN's sets
%     436 437 438 440 442, ... %JN's sets
%     443 446 447 448 449, ... %JN's sets
%     450 451 454 455 456, ... %JN's sets
%     458 459 460,...          %JN's sets
%     421 424 426 432 433, ... %IW's sets
%     437 438 440 441 449, ... %IW's sets
%     452 453 460 462 467, ... %IW's sets
%     469 472 479 482 498, ... %IW's sets
%     422 431 434 435 443, ... %IW's sets
%     ];   

TL_files =      {'MF160920.2','MF160921.2','MF160922.2','MF160923.2','MF160926.2',...
                 'MF160927.2','MF161003.2','MF161004.2','MF161006.2','MF161007.2',...
                 'MF161011.2','MF161013.2','MF161014.2','MF161017.2','MF161018.2',...
                 'MF161019.2'};
item_num = [425 426 428 429 430 ...
            431 432 433 434 435 ...
            436 437 438 439 440 ...
            441 ] ;


%---PW sessions after had run through all 400s Ts and Ls once---%
% TL_files = {'PW131203.2','PW131204.2','PW131206.3','PW131209.3',...
%     'PW131210.4','PW140520.4','PW140521.4','PW140522.4',...
%     'PW140528.3','PW140529.4','PW140605.4','PW140610.5',...
%     'PW140611.5','PW140613.5','PW140616.5','PW140617.5',...
%     'PW140703.5','PW140707.3','PW140710.4','PW140714.3'};
% 
% item_num = [421 422 424 425,...
%     428 430 431 432,...
%     434 435 437 439,...
%     440 442 443 444,...
%     450 451 454 456];


% going to use a generic 1 for all, calibration doesn't really matter for
% for this anlaysis, but the code that imports data requires it. Generic
% file so probably doesnt have correct gain, which doesn't matter.
% clrchng_files = {'RR150210.1'};

% for file = 1:9%length(TL_files)
%     getTLsData(TL_files{file},clrchng_files{1},item_num(file))
%     close all
% end

%%
%---[2] Analyze reaction times and number of saccades to find target---%

max_desired_blocks = 40;%which blocks 1-max_desired_blocks to analyze

repeat_reactiontimes = NaN(length(TL_files),max_desired_blocks);%avearge reaction times by block
novel_reactiontimes  = NaN(length(TL_files),max_desired_blocks);%avearge reaction times by block
repeat_num_saccades  = NaN(length(TL_files),max_desired_blocks);%average number of saccades to find T
novel_num_saccades   = NaN(length(TL_files),max_desired_blocks);%average number of saccades to find T
num_trials_removed = zeros(2,length(TL_files));%number of trials with rts < 150 ms or > 5000 ms that are removed
which_monkey = NaN(1,length(TL_files));%which monkey is this data for
for  file = 1:length(TL_files)
    %load preprocessed Ts and Ls file
%     load([data_dir TL_files{file}(1:2) '\' TL_files{file}(1:end-2) 'set'...
%         num2str(item_num(file)) 'eyedat.mat']);
  load([data_dir TL_files{file}(1:2) '\' TL_files{file}(1:end-2) 'set'...
        num2str(item_num(file)) 'eyedat.mat']);
    
    %---quickly extract block number without using a for-loop---%
    block = struct2cell(per);
    block = reshape(block,size(block,1),size(block,3));
    block = cell2mat(block(4,:));
    
    maxblock = max(block);
    if maxblock < max_desired_blocks
        disp([TL_files{file}(1:end-2) 'set' num2str(item_num(file)) 'eyedat.mat' ...
            ' Did not complete at least ' num2str(max_desired_blocks) ' blocks ' ...
            'Skipping Session!']);%should have but don't know for sure
        continue%skip analysis on this session
    end
    maxblock(maxblock > max_desired_blocks) = max_desired_blocks; %if did more than max_desired_blocks blocks just look at the 1st max_desired_blocks blocks
    
    if strcmpi(TL_files{file}(1:2),'IW')
        which_monkey(file) = 1;
    elseif strcmpi(TL_files{file}(1:2),'JN')
        which_monkey(file) = 2;
    elseif strcmpi(TL_files{file}(1:2),'MP')
        which_monkey(file) = 3;
    elseif strcmpi(TL_files{file}(1:2),'PW')
        which_monkey(file) = 4;
    end
    
    for blk = 1:maxblock
        if sum(block == blk) > 20 %so completed nearly all trials in this block, kind of arbitrary
            rts = NaN(2,12); %reaction times
            num_sacs = NaN(2,12);%number of saccades to find T
            row_index = [0 0]; %row 1 repeat, row 2 novel
            
            block_trials = find(block == blk);%find which trials go with this block
            for bt = 1:length(block_trials)
                trialtype =  per(block_trials(bt)).trialtype;
                events = per(block_trials(bt)).allval;
                times = per(block_trials(bt)).alltim;
                if per(block_trials(bt)).blk ~= blk
                    error('Blocks organized wrong')
                end
                
                if trialtype == 1 %repeated context
                    row_index(1) = row_index(1)+1;
                else %novel context
                    row_index(2) = row_index(2)+1;
                end
                
                trial_start = times(find(events == 15));
                eyedata_start = times(find(events == 100));
                fix_on_cross = times(find(events == 8))-trial_start-eyedata_start;
                if ~isempty(find(events == 3))%so sucessful trial
                    fix_on_T = fix_on_cross(end); %last fixation is on T
                else
                    fix_on_T = times(find(events == 24))-trial_start-eyedata_start;
                    %not really fixated but image turns off
                end
                fix_on_cross = fix_on_cross(1);%first fixation is on crosshair to start trial
                stimulus_on =  times(find(events == 23))-trial_start-eyedata_start;
                stimulus_off = times(find(events == 24))-trial_start-eyedata_start;
                
                saccadetimes = fixationstats{block_trials(bt)}.saccadetimes;
                if ~isempty(saccadetimes)
                    saccadetimes(:,saccadetimes(1,:) < stimulus_on | saccadetimes(1,:) > fix_on_T) = [];
                end
                %find saccades that started after the stimulus turned on
                %but before the monkey fixated the T
                
                rts(trialtype,row_index(trialtype)) = fix_on_T-stimulus_on;
                num_sacs(trialtype,row_index(trialtype)) = size(saccadetimes,2);
                
                %throw out unrealistic reaction times. Reaction times < 150
                %mn should not be possible/are false positives while
                %reaction times > 5000 ms probably indicate a lack of
                %engagemetn in that trial. 13 distractors + 1 target at a
                %saccade rate of 4 Hz should only take a max of 3500 ms to
                %fixate all items.
                if rts(trialtype,row_index(trialtype)) < 150 || rts(trialtype,row_index(trialtype)) > 5000
                    rts(trialtype,row_index(trialtype)) = NaN;
                    num_sacs(trialtype,row_index(trialtype)) = NaN;
                    num_trials_removed(trialtype,file) = ...
                        num_trials_removed(trialtype,file)+1;
                end
            end
            if ~all(row_index == 12)
                disp(['Error missing a few trials in block ' num2str(blk)])
            end
            
            repeat_reactiontimes(file,blk) = nanmean(rts(1,:));%avearge reaction times by block
            novel_reactiontimes(file,blk)  = nanmean(rts(2,:));%avearge reaction times by block
            repeat_num_saccades(file,blk)  = nanmean(num_sacs(1,:));%average number of saccades to find T
            novel_num_saccades(file,blk)   = nanmean(num_sacs(2,:));%average number of saccades to find T
            
        end
    end
end

%---Remove sessions in which we had to throw out more than 10% oftrials---%
poor_quality_sessions = find(sum(num_trials_removed) > 0.1*24*max_desired_blocks);
%# is 96 for 40 blocks,114 for 60 blocks

for pqs = 1:length(poor_quality_sessions)
    disp(['Removing ' TL_files{poor_quality_sessions(pqs)}(1:end-2) 'set' ...
        num2str(item_num(poor_quality_sessions(pqs)))...
        'eyedat.mat from analyis. Had to throw out more than 10% of trails!'])
end

repeat_reactiontimes(poor_quality_sessions,:) = NaN;
novel_reactiontimes(poor_quality_sessions,:) = NaN;
repeat_num_saccades(poor_quality_sessions,:) = NaN;
novel_num_saccades(poor_quality_sessions,:) = NaN;
which_monkey(poor_quality_sessions) = NaN;

%%
%---[3] Compute power analysis for individual monkeys---%

monkey_novel_rt_means = NaN(max_desired_blocks,4);% mean reaction times by monkey for novel contexts
monkey_novel_rt_stds =  NaN(max_desired_blocks,4);% std reaction times by monkey for novel contexts
monkey_repeat_rt_means = NaN(max_desired_blocks,4);% mean reaction times by monkey for repeat contexts
monkey_repeat_rt_stds =  NaN(max_desired_blocks,4);% std reaction times by monkey for repeat contexts
monkey_novel_sac_means = NaN(max_desired_blocks,4);% mean # sacs by monkey for novel contexts
monkey_novel_sac_stds =  NaN(max_desired_blocks,4);% std # sacs by monkey for novel contexts
monkey_repeat_sac_means = NaN(max_desired_blocks,4);% mean # sacs by monkey for repeat contexts
monkey_repeat_sac_stds =  NaN(max_desired_blocks,4);% std # sacs by monkey for repeat contexts
monkey_n = NaN(1,4);%number of sessions by monkey
monkey_inits = {'IW','JN','MP','PW'};

for monk = 1:4
    their_sessions = find(which_monkey == monk);
    monkey_n(monk) = length(their_sessions);
    
    monkey_novel_rt_means(:,monk) = mean(novel_reactiontimes(their_sessions,:));
    monkey_novel_rt_stds(:,monk) = std(novel_reactiontimes(their_sessions,:));
    monkey_repeat_rt_means(:,monk) = mean(repeat_reactiontimes(their_sessions,:));
    monkey_repeat_rt_stds(:,monk) = std(repeat_reactiontimes(their_sessions,:));
    
    monkey_novel_sac_means(:,monk) = mean(novel_num_saccades(their_sessions,:));
    monkey_novel_sac_stds(:,monk) = std(novel_num_saccades(their_sessions,:));
    monkey_repeat_sac_means(:,monk) = mean(repeat_num_saccades(their_sessions,:));
    monkey_repeat_sac_stds(:,monk) = std(repeat_num_saccades(their_sessions,:));
end

for monk = 1:4
    figure
    subplot(2,2,1)
    hold on
    errorbar(monkey_novel_rt_means(:,monk), monkey_novel_rt_stds(:,monk)...
        ./sqrt(monkey_n(monk)),'b')
    errorbar(monkey_repeat_rt_means(:,monk), monkey_repeat_rt_stds(:,monk)...
        ./sqrt(monkey_n(monk)),'r')
    hold off
    xlabel('Block #')
    ylabel('Reaction time (ms)')
    legend('Novel','Repeat')
    yl = ylim;
    ylim([0 yl(2)]);
    xlim([0 41])
    
    subplot(2,2,2)
    hold on
    errorbar(monkey_novel_sac_means(:,monk), monkey_novel_sac_stds(:,monk)...
        ./sqrt(monkey_n(monk)),'b')
    errorbar(monkey_repeat_sac_means(:,monk), monkey_repeat_sac_stds(:,monk)...
        ./sqrt(monkey_n(monk)),'r')
    hold off
    xlabel('Block #')
    ylabel('# of Saccades to Find the T')
    legend('Novel','Repeat')
    yl = ylim;
    ylim([0 yl(2)]);
    xlim([0 41])
    
    subplot(2,2,3)
    rt_n_80 = NaN(1,max_desired_blocks);%number of sessions needed for power 0.8
    rt_n_90 = NaN(1,max_desired_blocks);%number of sessions needed for power 0.9
    for b = 1:max_desired_blocks
        rt_n_80(b) = sampsizepwr('t',[monkey_novel_rt_means(b,monk) monkey_novel_rt_stds(b,monk)],...
            monkey_repeat_rt_means(b,monk),0.8);%power analysis for reaction times with power of 0.8
        rt_n_90(b) = sampsizepwr('t',[monkey_novel_rt_means(b,monk) monkey_novel_rt_stds(b,monk)],...
            monkey_repeat_rt_means(b,monk),0.9);%power analysis for reaction times with power of 0.9
    end
    hold on
    plot(rt_n_80)
    plot(rt_n_90,'r')
    xlabel('Block #')
    ylabel('# of Sessions Needed')
    legend('Power 0.8','Power 0.9')
    title(['Recommend n_{0.8} = ' num2str(mean(rt_n_80(end-9:end)))...
        ' &  n_{0.9} = ' num2str(mean(rt_n_90(end-9:end)))]);%recommend average of last 10 blocks interested in
    ylim([0 25])
    
    subplot(2,2,4)
    sac_n_80 = NaN(1,max_desired_blocks);%number of sessions needed for power 0.8
    sac_n_90 = NaN(1,max_desired_blocks);%number of sessions needed for power 0.9
    for b = 1:max_desired_blocks
        sac_n_80(b) = sampsizepwr('t',[monkey_novel_sac_means(b,monk) monkey_novel_sac_stds(b,monk)],...
            monkey_repeat_sac_means(b,monk),0.8);%power analysis for reaction times with power of 0.8
        sac_n_90(b) = sampsizepwr('t',[monkey_novel_sac_means(b,monk) monkey_novel_sac_stds(b,monk)],...
            monkey_repeat_sac_means(b,monk),0.9);%power analysis for reaction times with power of 0.9
    end
    hold on
    plot(sac_n_80)
    plot(sac_n_90,'r')
    xlabel('Block #')
    ylabel('# of Sessions Needed')
    title(['Recommend n_{0.8} = ' num2str(mean(sac_n_80(end-9:end)))...
        ' &  n_{0.9} = ' num2str(mean(sac_n_90(end-9:end)))]);%recommend average of last 10 blocks interested in
    legend('Power 0.8','Power 0.9')
    ylim([0 25])
    
    subtitle(['Monkey: ' monkey_inits{monk} ' n = ' num2str(monkey_n(monk))]);
end
%%
%---[4] Compute power anlaysis across monkeys---%

all_novel_rt_means = nanmean(novel_reactiontimes);% mean reaction times across all monkeys for all novel contexts
all_novel_rt_stds =  nanstd(novel_reactiontimes);% std reaction times across all monkeys for all novel contexts
all_repeat_rt_means = nanmean(repeat_reactiontimes);% mean reaction times across all monkeys for all repeat contexts
all_repeat_rt_stds =  nanstd(repeat_reactiontimes);% std reaction times across all monkeys for all repeat contexts
all_novel_sac_means = nanmean(novel_num_saccades);% mean # sacs across all monkeys for all novel contexts
all_novel_sac_stds =  nanstd(novel_num_saccades);% std # sacs across all monkeys for all novel contexts
all_repeat_sac_means = nanmean(repeat_num_saccades);% mean # sacs across all monkeys for all repeat contexts
all_repeat_sac_stds =  nanstd(repeat_num_saccades);% std # sacs across all monkeys for all repeat contexts
all_n = sum(~isnan(which_monkey)); %how many sessions across all monkeys

figure
subplot(2,2,1)
hold on
errorbar(all_novel_rt_means,all_novel_rt_stds./sqrt(all_n),'b')
errorbar(all_repeat_rt_means,all_repeat_rt_stds./sqrt(all_n),'r')
hold off
xlabel('Block #')
ylabel('Reaction time (ms)')
legend('Novel','Repeat')
yl = ylim;
ylim([0 yl(2)]);
xlim([0 41])

subplot(2,2,2)
hold on
errorbar(all_novel_sac_means,all_novel_sac_stds./sqrt(all_n),'b')
errorbar(all_repeat_sac_means,all_repeat_sac_stds./sqrt(all_n),'r')
hold off
xlabel('Block #')
ylabel('# of Saccades to Find the T')
legend('Novel','Repeat')
yl = ylim;
ylim([0 yl(2)]);
xlim([0 41])

subplot(2,2,3)
sac_n_80 = NaN(1,max_desired_blocks);%number of sessions needed for power 0.8
sac_n_90 = NaN(1,max_desired_blocks);%number of sessions needed for power 0.9
for b = 1:max_desired_blocks
    sac_n_80(b) = sampsizepwr('t',[all_novel_rt_means(b) all_novel_rt_stds(b)],...
        all_repeat_rt_means(b),0.8);%power analysis for reaction times with power of 0.8
    sac_n_90(b) = sampsizepwr('t',[all_novel_rt_means(b) all_novel_rt_stds(b)],...
        all_repeat_rt_means(b),0.9);%power analysis for reaction times with power of 0.9
end
hold on
plot(sac_n_80,'b')
plot(sac_n_90,'r')
hold off
xlabel('Block #')
ylabel('# of Sessions Needed')
title(['Recommend n_{0.8} = ' num2str(mean(sac_n_80(end-9:end)))...
    ' &  n_{0.9} = ' num2str(mean(sac_n_90(end-9:end)))]);%recommend average of last 10 blocks interested in
legend('Power 0.8','Power 0.9')
ylim([0 25])


subplot(2,2,4)
sac_n_80 = NaN(1,max_desired_blocks);%number of sessions needed for power 0.8
sac_n_90 = NaN(1,max_desired_blocks);%number of sessions needed for power 0.9
for b = 1:max_desired_blocks
    sac_n_80(b) = sampsizepwr('t',[all_novel_sac_means(b) all_novel_sac_stds(b)],...
        all_repeat_sac_means(b),0.8);%power analysis for reaction times with power of 0.8
    sac_n_90(b) = sampsizepwr('t',[all_novel_sac_means(b) all_novel_sac_stds(b)],...
        all_repeat_sac_means(b),0.9);%power analysis for reaction times with power of 0.9
end
hold on
plot(sac_n_80,'b')
plot(sac_n_90,'r')
hold off
xlabel('Block #')
ylabel('# of Sessions Needed')
title(['Recommend n_{0.8} = ' num2str(mean(sac_n_80(end-9:end)))...
    ' &  n_{0.9} = ' num2str(mean(sac_n_90(end-9:end)))]);%recommend average of last 10 blocks interested in
legend('Power 0.8','Power 0.9')
ylim([0 25])

subtitle(['Data Across monkeys, n_{sessions} = ' num2str(all_n)])
%% Quick statistical analysis
pval = [];
for block = 1:40
    [~,p] = ttest2(novel_reactiontimes(:,block),repeat_reactiontimes(:,block),'tail','both');
    pval(block) = p;
end