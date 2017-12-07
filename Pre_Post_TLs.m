% code written to analyze pre-vs-post analysis
% Modified from Analyze_All_TLs_Data on August 25, 2016 by Seth Konig
% runs basic analysis using ANOVA. Different statistics may be more
% informative.

clar
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\TL_relationalmemory\Eyedat\'; %where to find preprocessed data


%---Vivian Pre-Lesion Data---%
% pre_files = {'PW150518.2','PW150519.2','PW150520.2','PW150521.2','PW150522.2',...
%              'PW150526.2','PW150527.2','PW150528.2','PW150529.2','PW150601.2',...
%              'PW150602.2','PW150603.2','PW150604.2','PW150605.2','PW150608.2'};
% %didn't do any reaclimation may need to do use original pass through TLs or
% %throw out first 5 sets :( :(
% pre_item_num = [500 501 502 503 504,...
%     505 506 507 508 509,....
%     510 511 512 513 514];
% 
% post_files = {'PW160705.2','PW160706.2','PW160707.2','PW160708.3','PW160713.2',...
%               'PW160720.2','PW160722.2','PW160727.2','PW160728.2','PW160804.2',...
%               'PW160805.2','PW160811.2','PW160818.2','PW160819.2','PW160823.2'};
% post_item_num =       [520 521 522 523 524 ...
%     529 531 534 535 538 ...
%     539 542 547 548 549];

%---Red---%
% pre_files =      {'RR150129.3','RR150130.3','RR150202.3','RR150205.3','RR150206.3',...
%                  'RR150209.2','RR150210.2','RR150211.2','RR150213.2','RR150217.2',...
%                  'RR150219.2','RR150220.2','RR150223.2','RR150224.2','RR150225.2'};
% pre_item_num = [425 426 428 430 431,...
%             432 433 434 435 436,...
%             438 439 440 441 442];
% 
% post_files =      {'RR161012.2','RR161013.3','RR161014.2','RR161017.2','RR161018.2',...
%                  'RR161019.2','RR161020.2','RR161024.2','RR161025.2','RR161026.2',...
%                  'RR161027.2','RR161028.2','RR161101.2','RR161103.2'};
% post_item_num =      [489 491 492 493 494 ...
%                  495 496 498 499 500 ...
%                  501 502 504 505];


%---Tobii---%

pre_files =      {'TO150211.2','TO150218.2','TO150219.2','TO150223.2','TO150224.2',...
                'TO150225.2','TO150226.2','TO150302.2','TO150303.2','TO150305.2',...
                'TO150306.2','TO150309.2','TO150311.2','TO150312.2','TO150313.2'};

pre_item_num = [475 478 479 480 481,...
            482 483 485 486 488,...
            489 490 492 493 494];


post_files =      {'TO170510.2','TO170511.2','TO170512.2','TO170515.2','TO170516.2',...
                'TO170517.2','TO170518.2','TO170519.2','TO170523.2','TO170524.2',...
                'TO170526.2','TO170530.2','TO170601.2','TO170605.2','TO170606.2'};
post_item_num = [425 426 428 429 430,...
            431 432 433 436 437,...
            439 440 442 443 444];



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
                    fix_on_T = fix_on_T(end);
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
%% Plot Raw Data and Do Stats

% Plot average reaction times by block
figure
subplot(2,3,[1 2])
hold all
errorbar(1:40,nanmean(pre_repeat_reactiontimes(:,1:40)),nanstd(pre_repeat_reactiontimes(:,1:40))...
    ./sqrt(sum(~isnan(pre_repeat_reactiontimes(:,1:40)))));
errorbar(nanmean(pre_novel_reactiontimes(:,1:40)),nanstd(pre_novel_reactiontimes(:,1:40))...
    ./sqrt(sum(~isnan(pre_novel_reactiontimes(:,1:40)))));
errorbar(1:40,nanmean(post_repeat_reactiontimes(:,1:40)),nanstd(post_repeat_reactiontimes(:,1:40))...
    ./sqrt(sum(~isnan(post_repeat_reactiontimes(:,1:40)))));
errorbar(nanmean(post_novel_reactiontimes(:,1:40)),nanstd(post_novel_reactiontimes(:,1:40))...
    ./sqrt(sum(~isnan(post_novel_reactiontimes(:,1:40)))));
hold off
xlim([0 41])
xlabel('Block #')
ylabel('Reaction time (ms)')
legend('Pre-Repeated','Pre-Novel','Post-Repeated','Post-Novel','Location','SouthWest')
xlim([0 41])


% Plot average # of saccades to find T by block
subplot(2,3,[4 5])
hold all
errorbar(1:40,nanmean(pre_repeat_num_saccades(:,1:40)),nanstd(pre_repeat_num_saccades(:,1:40))...
    ./sqrt(sum(~isnan(pre_repeat_num_saccades(:,1:40)))));
errorbar(nanmean(pre_novel_num_saccades(:,1:40)),nanstd(pre_novel_num_saccades(:,1:40))...
    ./sqrt(sum(~isnan(pre_novel_num_saccades(:,1:40)))));
errorbar(1:40,nanmean(post_repeat_num_saccades(:,1:40)),nanstd(post_repeat_num_saccades(:,1:40))...
    ./sqrt(sum(~isnan(post_repeat_num_saccades(:,1:40)))));
errorbar(nanmean(post_novel_num_saccades(:,1:40)),nanstd(post_novel_num_saccades(:,1:40))...
    ./sqrt(sum(~isnan(post_novel_num_saccades(:,1:40)))));
hold off
xlabel('Block #')
ylabel('# of Saccades')
%legend('Pre-Repeated','Pre-Novel','Post-Repeated','Post-Novel','Location','NorthEast')
yl = ylim;
%ylim([yl yl(2)]);
xlim([0 41])


%---Run 2-way ANOVA on Raw Data for the last 10 blocks---%
raw_pre_last_block_rt_means_novel = mean(pre_novel_reactiontimes(:,30:40),2);
raw_pre_last_block_rt_means_repeat = mean(pre_repeat_reactiontimes(:,30:40),2);
raw_post_last_block_rt_means_novel = mean(post_novel_reactiontimes(:,30:40),2);
raw_post_last_block_rt_means_repeat = mean(post_repeat_reactiontimes(:,30:40),2);

raw_pre_last_block_sac_means_novel = mean(pre_novel_num_saccades(:,30:40),2);
raw_pre_last_block_sac_means_repeat = mean(pre_repeat_num_saccades(:,30:40),2);
raw_post_last_block_sac_means_novel = mean(post_novel_num_saccades(:,30:40),2);
raw_post_last_block_sac_means_repeat = mean(post_repeat_num_saccades(:,30:40),2);

npre = length(raw_pre_last_block_rt_means_novel);
npost = length(raw_post_last_block_rt_means_novel);

rt_means = [mean(raw_pre_last_block_rt_means_novel)  mean(raw_post_last_block_rt_means_novel);...
    mean(raw_pre_last_block_rt_means_repeat) mean(raw_post_last_block_rt_means_repeat)];
rt_sems = [std(raw_pre_last_block_rt_means_novel)/sqrt(npre) std(raw_post_last_block_rt_means_novel)/sqrt(npost);...
    std(raw_pre_last_block_rt_means_repeat)/sqrt(npre) std(raw_post_last_block_rt_means_repeat)/sqrt(npost)];


sac_means = [mean(raw_pre_last_block_sac_means_novel)  mean(raw_post_last_block_sac_means_novel);...
    mean(raw_pre_last_block_sac_means_repeat) mean(raw_post_last_block_sac_means_repeat)];
sac_sems = [std(raw_pre_last_block_sac_means_novel)/sqrt(npre) std(raw_post_last_block_sac_means_novel)/sqrt(npost);...
    std(raw_pre_last_block_sac_means_repeat)/sqrt(npre) std(raw_post_last_block_sac_means_repeat)/sqrt(npost)];


[P_ANOVA_rt,T,STATS,TERMS]=anovan([raw_pre_last_block_rt_means_novel; raw_pre_last_block_rt_means_repeat; ...
    raw_post_last_block_rt_means_novel; raw_post_last_block_rt_means_repeat],...
    {[ones(npre,1); ones(npre,1); 2*ones(npost,1); 2*ones(npost,1)],...
    [ones(npre,1); 2*ones(npre,1); ones(npost,1); 2*ones(npost,1)]},...
    'model','interaction','varnames',{'Pre/Post','Nov/Repeat'});


[P_ANOVA_sac,T,STATS,TERMS]=anovan([raw_pre_last_block_sac_means_novel; raw_pre_last_block_sac_means_repeat; ...
    raw_post_last_block_sac_means_novel; raw_post_last_block_sac_means_repeat],...
    {[ones(npre,1); ones(npre,1); 2*ones(npost,1); 2*ones(npost,1)],...
    [ones(npre,1); 2*ones(npre,1); ones(npost,1); 2*ones(npost,1)]},...
    'model','interaction','varnames',{'Pre/Post','Nov/Repeat'});


%---Run T-test---%
[~,p_rt_nov_raw] = ttest2(raw_pre_last_block_rt_means_novel,raw_post_last_block_rt_means_novel);
[~,p_rt_rep_raw] = ttest2(raw_pre_last_block_rt_means_repeat,raw_post_last_block_rt_means_repeat);
[~,p_sac_nov_raw] = ttest2(raw_pre_last_block_sac_means_novel,raw_post_last_block_sac_means_novel);
[~,p_sac_rep_raw] = ttest2(raw_pre_last_block_sac_means_repeat,raw_post_last_block_sac_means_repeat);


subplot(2,3,3)
errorb(rt_means,rt_sems)
hold on
if p_rt_nov_raw < 0.05
   plot(1,max(rt_means(1,:)+rt_sems(1,:))*1.05,'*k')
end
if p_rt_rep_raw < 0.05
   plot(2,max(rt_means(2,:)+rt_sems(2,:))*1.05,'*k')
end
hold off
set(gca,'Xtick',1:2)
set(gca,'XtickLabel',{'Novel','Repeat'})
xlabel('Context')
ylabel('Reaction Time (ms)')
legend('Pre','Post')
title(['ANOVA: p_{lesion} = ' num2str(P_ANOVA_rt(1),2) ', p_{novelty} = ' num2str(P_ANOVA_rt(2),3) ...
    ', p_{interact} = ' num2str(P_ANOVA_rt(3),3)])
box off

subplot(2,3,6)
errorb(sac_means,sac_sems)
hold on
if p_sac_nov_raw < 0.05
   plot(1,max(sac_means(1,:)+sac_sems(1,:))*1.05,'*k')
end
if p_sac_rep_raw < 0.05
   plot(2,max(sac_means(2,:)+sac_sems(2,:))*1.05,'*k')
end
hold off
set(gca,'Xtick',1:2)
set(gca,'XtickLabel',{'Novel','Repeat'})
xlabel('Context')
ylabel('# of Saccades')
title(['ANOVA: p_{lesion} = ' num2str(P_ANOVA_sac(1),2) ', p_{novelty} = ' num2str(P_ANOVA_sac(2),3) ...
    ', p_{interact} = ' num2str(P_ANOVA_sac(3),3)])
box off

subtitle(['Raw Data TLs: ' pre_files{1}(1:2)])

%% Plot Data and Do Stats
%---Run 2-way ANOVA on Normalized Data for the last 10 blocks---%
%normalize by average of all novel trials
% norm_pre_novel_reactiontimes = pre_novel_reactiontimes./mean(pre_novel_reactiontimes(:));
% norm_pre_repeat_reactiontimes = pre_repeat_reactiontimes./mean(pre_novel_reactiontimes(:));
% norm_post_novel_reactiontimes = post_novel_reactiontimes./mean(post_novel_reactiontimes(:));
% norm_post_repeat_reactiontimes = post_repeat_reactiontimes./mean(post_novel_reactiontimes(:));
% 
% norm_pre_novel_num_saccades = pre_novel_num_saccades./mean(pre_novel_num_saccades(:));
% norm_pre_repeat_num_saccades = pre_repeat_num_saccades./mean(pre_novel_num_saccades(:));
% norm_post_novel_num_saccades = post_novel_num_saccades./mean(post_novel_num_saccades(:));
% norm_post_repeat_num_saccades = post_repeat_num_saccades./mean(post_novel_num_saccades(:));

%normalize to mean reaction time of first 10 blocks on novel trials only
norm_pre_novel_reactiontimes = pre_novel_reactiontimes./mean(mean(pre_novel_reactiontimes(:,1:10)));
norm_pre_repeat_reactiontimes = pre_repeat_reactiontimes./mean(mean(pre_novel_reactiontimes(:,1:10)));
norm_post_novel_reactiontimes = post_novel_reactiontimes./mean(mean(post_novel_reactiontimes(:,1:10)));
norm_post_repeat_reactiontimes = post_repeat_reactiontimes./mean(mean(post_novel_reactiontimes(:,1:10)));

norm_pre_novel_num_saccades = pre_novel_num_saccades./mean(mean(pre_novel_num_saccades(:,1:10)));
norm_pre_repeat_num_saccades = pre_repeat_num_saccades./mean(mean(pre_novel_num_saccades(:,1:10)));
norm_post_novel_num_saccades = post_novel_num_saccades./mean(mean(post_novel_num_saccades(:,1:10)));
norm_post_repeat_num_saccades = post_repeat_num_saccades./mean(mean(post_novel_num_saccades(:,1:10)));

Normalized_pre_last_block_rt_means_norm_novel = mean(norm_pre_novel_reactiontimes(:,30:40),2);
Normalized_pre_last_block_rt_means_norm_repeat = mean(norm_pre_repeat_reactiontimes(:,30:40),2);
Normalized_post_last_block_rt_means_norm_novel = mean(norm_post_novel_reactiontimes(:,30:40),2);
Normalized_post_last_block_rt_means_norm_repeat = mean(norm_post_repeat_reactiontimes(:,30:40),2);

Normalized_pre_last_block_sac_means_norm_novel = mean(norm_pre_novel_num_saccades(:,30:40),2);
Normalized_pre_last_block_sac_means_norm_repeat = mean(norm_pre_repeat_num_saccades(:,30:40),2);
Normalized_post_last_block_sac_means_norm_novel = mean(norm_post_novel_num_saccades(:,30:40),2);
Normalized_post_last_block_sac_means_norm_repeat = mean(norm_post_repeat_num_saccades(:,30:40),2);

npre = length(Normalized_pre_last_block_rt_means_norm_novel);
npost = length(Normalized_post_last_block_rt_means_norm_novel);

rt_means_norm = [mean(Normalized_pre_last_block_rt_means_norm_novel)  mean(Normalized_post_last_block_rt_means_norm_novel);...
    mean(Normalized_pre_last_block_rt_means_norm_repeat) mean(Normalized_post_last_block_rt_means_norm_repeat)];
rt_sems_norm = [std(Normalized_pre_last_block_rt_means_norm_novel)/sqrt(npre) std(Normalized_post_last_block_rt_means_norm_novel)/sqrt(npost);...
    std(Normalized_pre_last_block_rt_means_norm_repeat)/sqrt(npre) std(Normalized_post_last_block_rt_means_norm_repeat)/sqrt(npost)];


sac_means_norm = [mean(Normalized_pre_last_block_sac_means_norm_novel)  mean(Normalized_post_last_block_sac_means_norm_novel);...
    mean(Normalized_pre_last_block_sac_means_norm_repeat) mean(Normalized_post_last_block_sac_means_norm_repeat)];
sac_sems_norm = [std(Normalized_pre_last_block_sac_means_norm_novel)/sqrt(npre) std(Normalized_post_last_block_sac_means_norm_novel)/sqrt(npost);...
    std(Normalized_pre_last_block_sac_means_norm_repeat)/sqrt(npre) std(Normalized_post_last_block_sac_means_norm_repeat)/sqrt(npost)];


%---Run T-test---%
[~,p_rt_nov_Normalized] = ttest2(Normalized_pre_last_block_rt_means_norm_novel,Normalized_post_last_block_rt_means_norm_novel);
[~,p_rt_rep_Normalized] = ttest2(Normalized_pre_last_block_rt_means_norm_repeat,Normalized_post_last_block_rt_means_norm_repeat);
[~,p_sac_nov_Normalized] = ttest2(Normalized_pre_last_block_sac_means_norm_novel,Normalized_post_last_block_sac_means_norm_novel);
[~,p_sac_rep_Normalized] = ttest2(Normalized_pre_last_block_sac_means_norm_repeat,Normalized_post_last_block_sac_means_norm_repeat);


% Plot average reaction times by block
figure
subplot(2,3,[1 2])
hold all
errorbar(1:40,nanmean(norm_pre_repeat_reactiontimes(:,1:40)),nanstd(norm_pre_repeat_reactiontimes(:,1:40))...
    ./sqrt(sum(~isnan(norm_pre_repeat_reactiontimes(:,1:40)))));
errorbar(nanmean(norm_pre_novel_reactiontimes(:,1:40)),nanstd(norm_pre_novel_reactiontimes(:,1:40))...
    ./sqrt(sum(~isnan(norm_pre_novel_reactiontimes(:,1:40)))));
errorbar(1:40,nanmean(norm_post_repeat_reactiontimes(:,1:40)),nanstd(norm_post_repeat_reactiontimes(:,1:40))...
    ./sqrt(sum(~isnan(norm_post_repeat_reactiontimes(:,1:40)))));
errorbar(nanmean(norm_post_novel_reactiontimes(:,1:40)),nanstd(norm_post_novel_reactiontimes(:,1:40))...
    ./sqrt(sum(~isnan(norm_post_novel_reactiontimes(:,1:40)))));
hold off
xlim([0 41])
xlabel('Block #')
ylabel('Proportion of #  RT')
legend('Pre-Repeated','Pre-Novel','Post-Repeated','Post-Novel','Location','SouthWest')
xlim([0 41])


% Plot average # of saccades to find T by block
subplot(2,3,[4 5])
hold all
errorbar(1:40,nanmean(norm_pre_repeat_num_saccades(:,1:40)),nanstd(norm_pre_repeat_num_saccades(:,1:40))...
    ./sqrt(sum(~isnan(norm_pre_repeat_num_saccades(:,1:40)))));
errorbar(nanmean(norm_pre_novel_num_saccades(:,1:40)),nanstd(norm_pre_novel_num_saccades(:,1:40))...
    ./sqrt(sum(~isnan(norm_pre_novel_num_saccades(:,1:40)))));
errorbar(1:40,nanmean(norm_post_repeat_num_saccades(:,1:40)),nanstd(norm_post_repeat_num_saccades(:,1:40))...
    ./sqrt(sum(~isnan(norm_post_repeat_num_saccades(:,1:40)))));
errorbar(nanmean(norm_post_novel_num_saccades(:,1:40)),nanstd(norm_post_novel_num_saccades(:,1:40))...
    ./sqrt(sum(~isnan(norm_post_novel_num_saccades(:,1:40)))));
hold off
xlabel('Block #')
ylabel('Proportion of # of Saccades')
%legend('Pre-Repeated','Pre-Novel','Post-Repeated','Post-Novel','Location','NorthEast')
yl = ylim;
%ylim([yl yl(2)]);
xlim([0 41])


subplot(2,3,3)
errorb(rt_means_norm,rt_sems_norm)
hold on
if p_rt_nov_Normalized < 0.05
   plot(1,max(rt_means_norm(1,:)+rt_sems_norm(1,:))*1.05,'*k')
end
if p_rt_rep_Normalized < 0.05
   plot(2,max(rt_means_norm(2,:)+rt_sems_norm(2,:))*1.05,'*k')
end
hold off
set(gca,'Xtick',1:2)
set(gca,'XtickLabel',{'Novel','Repeat'})
xlabel('Context')
ylabel('Reaction Time (ms)')
legend('Pre','Post')
title(['p_{t-test repeat pre vs post} = ' num2str(p_rt_rep_Normalized,3)])
box off

subplot(2,3,6)
errorb(sac_means_norm,sac_sems_norm)
hold on
if p_sac_nov_Normalized < 0.05
   plot(1,max(sac_means_norm(1,:)+sac_sems_norm(1,:))*1.05,'*k')
end
if p_sac_rep_Normalized < 0.05
   plot(2,max(sac_means_norm(2,:)+sac_sems_norm(2,:))*1.05,'*k')
end
hold off
set(gca,'Xtick',1:2)
set(gca,'XtickLabel',{'Novel','Repeat'})
xlabel('Context')
ylabel('# of Saccades')
title(['p_{t-test repeat pre vs post} = ' num2str(p_sac_rep_Normalized,3)])
box off

subtitle(['Normalized Data TLs: ' pre_files{1}(1:2)])




%% Plot average saccade rate by block
figure
hold all
errorbar(1:40,nanmean(pre_repeat_saccade_rate(:,1:40)),nanstd(pre_repeat_saccade_rate(:,1:40))...
    ./sqrt(sum(~isnan(pre_repeat_saccade_rate(:,1:40)))));
errorbar(nanmean(pre_novel_saccade_rate(:,1:40)),nanstd(pre_novel_saccade_rate(:,1:40))...
    ./sqrt(sum(~isnan(pre_novel_saccade_rate(:,1:40)))));
errorbar(1:40,nanmean(post_repeat_saccade_rate(:,1:40)),nanstd(post_repeat_saccade_rate(:,1:40))...
    ./sqrt(sum(~isnan(post_repeat_saccade_rate(:,1:40)))));
errorbar(nanmean(post_novel_saccade_rate(:,1:40)),nanstd(post_novel_saccade_rate(:,1:40))...
    ./sqrt(sum(~isnan(post_novel_saccade_rate(:,1:40)))));
hold off
xlabel('Block #')
ylabel('Saccade Rate (Hz)')
yl = ylim;
ylim([3 6]);
xlim([0 41])
legend('Pre-Repeated','Pre-Novel','Post-Repeated','Post-Novel','Location','SouthEast')
subtitle(['Saccade Rate: ' pre_files{1}(1:2)])



