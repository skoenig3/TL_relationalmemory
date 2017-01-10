% code written to analyze pre-vs-post analysis
% Modified from Analyze_All_TLs_Data on August 25, 2016 by Seth Konig

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\TL_relationalmemory\Eyedat\';


%---Vivian Pre-Lesion Data---%
pre_files = {'PW150518.2','PW150519.2','PW150520.2','PW150521.2','PW150522.2',...
             'PW150526.2','PW150527.2','PW150528.2','PW150529.2','PW150601.2',...
             'PW150602.2','PW150603.2','PW150604.2','PW150605.2','PW150608.2'};
%didn't do any reaclimation may need to do use original pass through TLs or
%throw out first 5 sets :( :(
pre_item_num = [500 501 502 503 504,...
    505 506 507 508 509,....
    510 511 512 513 514];

post_files = {'PW160705.2','PW160706.2','PW160707.2','PW160708.3','PW160713.2',...
              'PW160720.2','PW160722.2','PW160727.2','PW160728.2','PW160804.2',...
              'PW160805.2','PW160811.2','PW160818.2','PW160819.2','PW160823.2'};
post_item_num =       [520 521 522 523 524 ...
    529 531 534 535 538 ...
    539 542 547 548 549];


max_desired_blocks = 40;%which blocks 1-max_desired_blocks to analyze


%%%---Pre Lesion Stuff---%%%

pre_novel_saccade_rate =  NaN(length(pre_files),max_desired_blocks);%avearge saccade rate by block
pre_repeat_saccade_rate =  NaN(length(pre_files),max_desired_blocks);%avearge saccade rate by block
pre_repeat_reactiontimes = NaN(length(pre_files),max_desired_blocks);%avearge reaction times by block
pre_novel_reactiontimes  = NaN(length(pre_files),max_desired_blocks);%avearge reaction times by block
pre_novel_num_saccades= NaN(length(pre_files),max_desired_blocks);%average number of saccades to find T
pre_repeat_num_saccades= NaN(length(pre_files),max_desired_blocks);%average number of saccades to find T

for  file = 1:length(pre_files)
    %load preprocessed Ts and Ls file
    load([data_dir pre_files{file}(1:2) '\' pre_files{file}(1:end-2) 'set'...
        num2str(pre_item_num(file)) 'eyedat.mat']);
    
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
    
    for blk = 1:maxblock
        if sum(block == blk) > 20 %so completed nearly all trials in this block
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
                end
            end
            if ~all(row_index == 12)
                disp(['Error missing a few trials in block ' num2str(blk)])
            end
            
            pre_repeat_reactiontimes(file,blk) = nanmean(rts(1,:));%avearge reaction times by block
            pre_novel_reactiontimes(file,blk)  = nanmean(rts(2,:));%avearge reaction times by block
            pre_repeat_num_saccades(file,blk)  = nanmean(num_sacs(1,:));%average number of saccades to find T
            pre_novel_num_saccades(file,blk)   = nanmean(num_sacs(2,:));%average number of saccades to find T
            pre_novel_saccade_rate(file,blk)   = nanmean(1000*num_sacs(1,:)./rts(1,:));%avearge saccade rate
            pre_repeat_saccade_rate(file,blk)   = nanmean(1000*num_sacs(2,:)./rts(2,:));%avearge saccade rate
            
        end
    end
end

% Get Percent change from block 1 for last 10 blocks
pre_last_block_rt_means_repeat = 100*mean(pre_repeat_reactiontimes(:,30:40),2)./mean(pre_novel_reactiontimes(:,1));
pre_last_block_rt_means_novel = 100*mean(pre_novel_reactiontimes(:,30:40),2)./mean(pre_novel_reactiontimes(:,1));

pre_last_block_sac_means_repeat = 100*mean(pre_repeat_num_saccades(:,30:40),2)...
    ./mean(pre_novel_num_saccades(:,1));
pre_last_block_sac_means_novel = 100*mean(pre_novel_num_saccades(:,30:40),2)...
    ./mean(pre_novel_num_saccades(:,1));

% Get Saccade Rate by Block and Last 10 blocks

pre_last_block_sr_means_repeat = mean(pre_repeat_saccade_rate(:,30:40),2);
pre_last_block_sr_means_novel = mean(pre_novel_saccade_rate(:,30:40),2);


%%%---Post Lesion Stuff---%%%
post_novel_saccade_rate =  NaN(length(post_files),max_desired_blocks);%avearge saccade rate by block
post_repeat_saccade_rate =  NaN(length(post_files),max_desired_blocks);%avearge saccade rate by block
post_repeat_reactiontimes = NaN(length(post_files),max_desired_blocks);%avearge reaction times by block
post_novel_reactiontimes  = NaN(length(post_files),max_desired_blocks);%avearge reaction times by block
post_novel_num_saccades= NaN(length(post_files),max_desired_blocks);%average number of saccades to find T
post_repeat_num_saccades= NaN(length(post_files),max_desired_blocks);%average number of saccades to find T

for  file = 1:length(post_files)
    %load preprocessed Ts and Ls file
    load([data_dir post_files{file}(1:2) '\' post_files{file}(1:end-2) 'set'...
        num2str(post_item_num(file)) 'eyedat.mat']);
    
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
    
    for blk = 1:maxblock
        if sum(block == blk) > 20 %so completed nearly all trials in this block
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
                end
            end
            if ~all(row_index == 12)
                disp(['Error missing a few trials in block ' num2str(blk)])
            end
            
            post_repeat_reactiontimes(file,blk) = nanmean(rts(1,:));%avearge reaction times by block
            post_novel_reactiontimes(file,blk)  = nanmean(rts(2,:));%avearge reaction times by block
            post_repeat_num_saccades(file,blk)  = nanmean(num_sacs(1,:));%average number of saccades to find T
            post_novel_num_saccades(file,blk)   = nanmean(num_sacs(2,:));%average number of saccades to find T
            post_novel_saccade_rate(file,blk)   = nanmean(1000*num_sacs(1,:)./rts(1,:));%avearge saccade rate
            post_repeat_saccade_rate(file,blk)   = nanmean(1000*num_sacs(2,:)./rts(2,:));%avearge saccade rate
            
        end
    end
end

% Get Percent change from block 1 for last 10 blocks
post_last_block_rt_means_repeat = 100*mean(post_repeat_reactiontimes(:,30:40),2)./mean(post_novel_reactiontimes(:,1));
post_last_block_rt_means_novel = 100*mean(post_novel_reactiontimes(:,30:40),2)./mean(post_novel_reactiontimes(:,1));

post_last_block_sac_means_repeat = 100*mean(post_repeat_num_saccades(:,30:40),2)...
    ./mean(post_novel_num_saccades(:,1));
post_last_block_sac_means_novel = 100*mean(post_novel_num_saccades(:,30:40),2)...
    ./mean(post_novel_num_saccades(:,1));

% Get Saccade Rate by Block and Last 10 blocks

post_last_block_sr_means_repeat = mean(post_repeat_saccade_rate(:,30:40),2);
post_last_block_sr_means_novel = mean(post_novel_saccade_rate(:,30:40),2);
%% Plot average reaction times by block

figure
hold on
errorbar(1:40,nanmean(pre_repeat_reactiontimes(:,1:40)),nanstd(pre_repeat_reactiontimes(:,1:40))...
    ./sqrt(sum(~isnan(pre_repeat_reactiontimes(:,1:40)))),'r');
errorbar(nanmean(pre_novel_reactiontimes(:,1:40)),nanstd(pre_novel_reactiontimes(:,1:40))...
    ./sqrt(sum(~isnan(pre_novel_reactiontimes(:,1:40)))),'b');
errorbar(1:40,nanmean(post_repeat_reactiontimes(:,1:40)),nanstd(post_repeat_reactiontimes(:,1:40))...
    ./sqrt(sum(~isnan(post_repeat_reactiontimes(:,1:40)))),'k--');
errorbar(nanmean(post_novel_reactiontimes(:,1:40)),nanstd(post_novel_reactiontimes(:,1:40))...
    ./sqrt(sum(~isnan(post_novel_reactiontimes(:,1:40)))),'g--');
hold off
xlim([0 41])
xlabel('Block #')
ylabel('Reaction time (ms)')
legend('Pre-Repeated','Pre-Novel','Post-Repeated','Post-Novel','Location','NorthEast')
xlim([0 41])
%% Plot average # of saccades to find T by block

figure
hold on
errorbar(1:40,nanmean(pre_repeat_num_saccades(:,1:40)),nanstd(pre_repeat_num_saccades(:,1:40))...
    ./sqrt(sum(~isnan(pre_repeat_num_saccades(:,1:40)))),'r');
errorbar(nanmean(pre_novel_num_saccades(:,1:40)),nanstd(pre_novel_num_saccades(:,1:40))...
    ./sqrt(sum(~isnan(pre_novel_num_saccades(:,1:40)))),'b');
errorbar(1:40,nanmean(post_repeat_num_saccades(:,1:40)),nanstd(post_repeat_num_saccades(:,1:40))...
    ./sqrt(sum(~isnan(post_repeat_num_saccades(:,1:40)))),'k--');
errorbar(nanmean(post_novel_num_saccades(:,1:40)),nanstd(post_novel_num_saccades(:,1:40))...
    ./sqrt(sum(~isnan(post_novel_num_saccades(:,1:40)))),'g--');
hold off
xlabel('Block #')
ylabel('# of Saccades')
legend('Pre-Repeated','Pre-Novel','Post-Repeated','Post-Novel','Location','NorthEast')
yl = ylim;
ylim([0 yl(2)]);
xlim([0 41])
%% Plot average saccade rate by block
figure
hold on
errorbar(1:40,nanmean(pre_repeat_saccade_rate(:,1:40)),nanstd(pre_repeat_saccade_rate(:,1:40))...
    ./sqrt(sum(~isnan(pre_repeat_saccade_rate(:,1:40)))),'r');
errorbar(nanmean(pre_novel_saccade_rate(:,1:40)),nanstd(pre_novel_saccade_rate(:,1:40))...
    ./sqrt(sum(~isnan(pre_novel_saccade_rate(:,1:40)))),'b');
errorbar(1:40,nanmean(post_repeat_saccade_rate(:,1:40)),nanstd(post_repeat_saccade_rate(:,1:40))...
    ./sqrt(sum(~isnan(post_repeat_saccade_rate(:,1:40)))),'k--');
errorbar(nanmean(post_novel_saccade_rate(:,1:40)),nanstd(post_novel_saccade_rate(:,1:40))...
    ./sqrt(sum(~isnan(post_novel_saccade_rate(:,1:40)))),'g--');
hold off
xlabel('Block #')
ylabel('Saccade Rate (Hz)')
yl = ylim;
ylim([3 yl(2)]);
xlim([0 41])
legend('Pre-Repeated','Pre-Novel','Post-Repeated','Post-Novel','Location','SouthEast')
%% Compare Reaction times for the last 10 blocks
means = [mean(pre_last_block_rt_means_novel)  mean(post_last_block_rt_means_novel);...
    mean(pre_last_block_rt_means_repeat) mean(post_last_block_rt_means_repeat)];

npre = sqrt(length(pre_last_block_rt_means_repeat));
npost = sqrt(length(post_last_block_rt_means_repeat));
sems = [std(pre_last_block_rt_means_novel)/npre std(post_last_block_rt_means_novel)/npost;...
    std(pre_last_block_rt_means_repeat)/npre std(post_last_block_rt_means_repeat)/npost];

[~,p_rt_nov] = ttest(pre_last_block_rt_means_novel,post_last_block_rt_means_novel,'tail','both');
[~,p_rt_rep] = ttest(pre_last_block_rt_means_repeat,post_last_block_rt_means_repeat,'tail','both');

figure
hold on
bar(means)
errorb(means,sems)
if p_rt_nov < 0.05
   plot(1,100,'*k')
end
if p_rt_rep < 0.05
   plot(2,100,'*k')
end
hold off
set(gca,'Xtick',1:2)
set(gca,'XtickLabel',{'Novel','Repeat'})
xlabel('Context')
ylabel('% of RT of novel block 1')
legend('Pre','Post')

%% Compare Number of Saccades to find T for last 10 blocks
means = [mean(pre_last_block_sac_means_novel)  mean(post_last_block_sac_means_novel);...
    mean(pre_last_block_sac_means_repeat) mean(post_last_block_sac_means_repeat)];

npre = sqrt(length(pre_last_block_sac_means_repeat));
npost = sqrt(length(post_last_block_sac_means_repeat));
sems = [std(pre_last_block_sac_means_novel)/npre std(post_last_block_sac_means_novel)/npost;...
    std(pre_last_block_sac_means_repeat)/npre std(post_last_block_sac_means_repeat)/npost];

[~,p_sac_nov] = ttest(pre_last_block_sac_means_novel,post_last_block_sac_means_novel,'tail','both');
[~,p_sac_rep] = ttest(pre_last_block_sac_means_repeat,post_last_block_sac_means_repeat,'tail','both');

figure
hold on
bar(means)
errorb(means,sems)
if p_sac_nov < 0.05
   plot(1,100,'*k')
end
if p_sac_rep < 0.05
   plot(2,100,'*k')
end
hold off
set(gca,'Xtick',1:2)
set(gca,'XtickLabel',{'Novel','Repeat'})
xlabel('Context')
ylabel('% of Saccades of novel block 1')
legend('Pre','Post')


