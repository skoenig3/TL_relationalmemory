% code written to analyze all the Ts and Ls data
% written by Seth Konig 2/9/15
%
% [1] Import task data and detect fixations and saccades.
% [2] Analyze reaction times and number of saccades to find target.

%%
%---[1] Import task data and detect fixations and saccades---%

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\TL_relationalmemory\Eyedat\';

% TL_files =      {'RR150129.3','RR150130.3','RR150202.3','RR150205.3','RR150206.3',...
%                  'RR150209.2','RR150210.2','RR150211.2','RR150213.2','RR150217.2',...
%                  'RR150219.2','RR150220.2','RR150223.2','RR150224.2','RR150225.2'};
% clrchng_files = {'RR150129.1','RR150130.1','RR150202.1','RR150205.1','RR150206.1',...
%                  'RR150209.1','RR150210.1','RR150211.1','RR150213.1','RR150217.1',...
%                  'RR150219.1','RR150220.1','RR150223.1','RR150224.1','RR150225.1'};
% item_num = [425 426 428 430 431,...
%             432 433 434 435 436,...
%             438 439 440 441 442];

% TL_files =      {'TO150211.2','TO150218.2','TO150219.2','TO150223.2','TO150224.2',...
%                 'TO150225.2','TO150226.2','TO150302.2','TO150303.2','TO150305.2',...
%                 'TO150306.2','TO150309.2','TO150311.2','TO150312.2','TO150313.2'};
% 
% clrchng_files = {'TO150211.1','TO150218.1','TO150219.1','TO150223.1','TO150224.1',...
%                 'TO150225.1','TO150226.1','TO150302.1','TO150303.1','TO150305.1',...
%                 'TO150306.1','TO150309.1','TO150311.1','TO150312.1','TO150313.1'};
% 
% item_num = [475 478 479 480 481,...
%             482 483 485 486 488,...
%             489 490 492 493 494];

TL_files =       {'PW150518.2','PW150519.2','PW150520.2','PW150521.2','PW150522.2',...
                  'PW150526.2','PW150527.2','PW150528.2','PW150529.2','PW150601.2',...
                  'PW150602.2','PW150603.2','PW150604.2','PW150605.2','PW150608.2'};
    
clrchng_files =  {'PW150518.3','PW150519.1','PW150520.1','PW150521.1','PW150522.1',...
                  'PW150526.1','PW150527.1','PW150528.1','PW150529.1','PW150601.1',...
                  'PW150602.1','PW150603.1','PW150604.1','PW150605.1','PW150608.1'};

item_num = [500 501 502 503 504,...
            505 506 507 508 509,....
            510 511 512 513 514];

% for file = 11:length(TL_files)
%     getTLsData(TL_files{file},clrchng_files{file},item_num(file))
%     close all
% end
%%
%---[2] Analyze reaction times and number of saccades to find target---%

max_desired_blocks = 40;%which blocks 1-max_desired_blocks to analyze

repeat_reactiontimes = NaN(length(TL_files),max_desired_blocks);%avearge reaction times by block
novel_reactiontimes  = NaN(length(TL_files),max_desired_blocks);%avearge reaction times by block
repeat_num_saccades  = NaN(length(TL_files),max_desired_blocks);%average number of saccades to find T
novel_num_saccades   = NaN(length(TL_files),max_desired_blocks);%average number of saccades to find T
novel_pupil_data     = cell(length(TL_files),max_desired_blocks);%average pupil data
repeat_pupil_data    = cell(length(TL_files),max_desired_blocks);%average pupil data
novel_fix_duration =  cell(length(TL_files),max_desired_blocks);%avearge fixation duration
repeat_fix_duration = cell(length(TL_files),max_desired_blocks);%avearge fixation duration
num_trials_removed = zeros(2,length(TL_files));%number of trials with rts < 150 ms or > 5000 ms that are removed
for  file = 1:length(TL_files)
    %load preprocessed Ts and Ls file
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
    
    for blk = 1:maxblock
        if sum(block == blk) > 20 %so completed nearly all trials in this block
            rts = NaN(2,12); %reaction times
            num_sacs = NaN(2,12);%number of saccades to find T
            fixationdurations{1} = NaN(12,20);%fixation durations
            fixationdurations{2} = NaN(12,20);%fixation durations
            if ~isempty(pupildata)%horizontal pupil diameter
                pupil{1} = NaN(12,360);
                pupil{2} = NaN(12,360);
            end
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
                
                %get pupil data 500 ms before fixation on crosshair and 500 ms after fix_on_T
                if ~isempty(pupildata)%horizontal pupil diameter
                    pupil{trialtype}(row_index(trialtype),:) = ...
                        pupildata{block_trials(bt)}(ceil((stimulus_on-eyedata_start)/5)-100:...
                        ceil((stimulus_on-eyedata_start)/5)+259);
                end
                
                fixationtimes = fixationstats{block_trials(bt)}.fixationtimes;
                fixationtimes(:,fixationtimes(1,:) < stimulus_on | fixationtimes(1,:) >= fix_on_T) = [];
                %usally cortex is faster than Cluster Fix at stating fixation onset
                %find fixations that started after the stimulus turned on
                %but before the monkey fixated the T
                
                fixationdurations{trialtype}(row_index(trialtype),1:size(fixationtimes,2)) = ...
                    diff(fixationtimes)+1;
                
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
                    pupil{trialtype}(row_index(trialtype),:) = NaN;
                end
            end
            if ~all(row_index == 12)
                disp(['Error missing a few trials in block ' num2str(blk)])
            end
            
            repeat_reactiontimes(file,blk) = nanmean(rts(1,:));%avearge reaction times by block
            novel_reactiontimes(file,blk)  = nanmean(rts(2,:));%avearge reaction times by block
            repeat_num_saccades(file,blk)  = nanmean(num_sacs(1,:));%average number of saccades to find T
            novel_num_saccades(file,blk)   = nanmean(num_sacs(2,:));%average number of saccades to find T
            
            novel_fix_duration{file,blk}   = nanmean(fixationdurations{1}(:,1:5));
            repeat_fix_duration{file,blk}   = nanmean(fixationdurations{2}(:,1:5));
            if ~isempty(pupildata)
                %normalize pupil data by eye values 0-25 ms before stimulus
                %onset on a trial by trial basis
                novel_pupil_data{file,blk}  = nanmean(pupil{1}./...
                    abs(nanmean(pupil{1}(:,94:99)')'*ones(1,size(pupil{1},2))));%average pupil data
                repeat_pupil_data{file,blk} = nanmean(pupil{2}./...
                    abs(nanmean(pupil{2}(:,94:99)')'*ones(1,size(pupil{2},2))));%average pupil data
            end
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

%% Plot Pupil data averaged over last 10 blocks, blocks 30-40
all_rep = [];
all_nov = [];
for file = 1:length(TL_files)
    for blk = 30:40
        all_rep = [all_rep; repeat_pupil_data{file,blk}];
        all_nov = [all_nov; novel_pupil_data{file,blk}];
    end
end
figure
hold on
plot(nanmean(all_rep),'r')
plot(nanmean(all_nov),'b')
hold off
xlabel('Time from Fixation on Crosshair')
ylabel('Pupil Diameter (a.u.)')

%% plot fixation durations averaged over block 30 through 40
all_rep = [];
all_nov = [];
for  file = 1:length(TL_files)
    for blk = 30:40
        all_rep = [all_rep; repeat_fix_duration{file,blk}];
        all_nov = [all_nov; novel_fix_duration{file,blk}];
    end
end
figure
subplot(1,2,1)
hold on
errorbar(nanmean(all_rep),nanstd(all_rep)./sqrt(sum(~isnan(all_rep))),'r')
errorbar(nanmean(all_nov),nanstd(all_nov)./sqrt(sum(~isnan(all_nov))),'b')
hold off
xlabel('Fixation Number')
ylabel('Fixation Duration (ms)')
%% Computer statistical analysis by block for reaction times and number of saccade to target
rt_p = NaN(1,40); %p-value for reaciton times
sac_p = NaN(1,40);% p-value for number of saccades
for block = 1:40
    [~,p] = ttest2(novel_reactiontimes(:,block),repeat_reactiontimes(:,block),'tail','both');
    rt_p(block) = p;
    [~,p] = ttest2(novel_num_saccades(:,block),repeat_num_saccades(:,block),'tail','both');
    sac_p(block) = p;
end
%find which blocks have p-values less than ...
% rt_p05 = find(rt_p < 0.05 & rt_p >= 0.01);
% rt_p01 = find(rt_p < 0.01 & rt_p >= 0.001);
% rt_p001 = find(rt_p < 0.001);
% sac_p05 = find(sac_p < 0.05 & sac_p >= 0.01);
% sac_p01 = find(sac_p < 0.01 & sac_p >= 0.001);
% sac_p001 = find(sac_p < 0.001); 
rt_p05 = find(rt_p < 0.05) ;
sac_p05 = find(sac_p < 0.05);

%% Plot average reaction time and # of saccades to find T by block averaged by session

figure
subplot(1,2,1)
hold on
p(1) = errorbar(1:40,nanmean(repeat_reactiontimes(:,1:40)),nanstd(repeat_reactiontimes(:,1:40))...
    ./sqrt(sum(~isnan(repeat_reactiontimes(:,1:40)))),'r');
p(2) = errorbar(nanmean(novel_reactiontimes(:,1:40)),nanstd(novel_reactiontimes(:,1:40))...
    ./sqrt(sum(~isnan(novel_reactiontimes(:,1:40)))),'b');
p(3) = plot(rt_p05,1500*ones(1,length(rt_p05)),'k*');
% p(4) = plot(rt_p01,1500*ones(1,length(rt_p01)),'r*');
% p(5) = plot(rt_p001,1500*ones(1,length(rt_p001)),'g*');
hold off
xlabel('Block #')
ylabel('Reaction time (ms)')
legend('Repeated','Novel','p<0.05','p < 0.01','p < 0.001','Location','SouthEast')
yl = ylim;
ylim([0 yl(2)]);
xlim([0 41])

subplot(1,2,2)
hold on
p(1) = errorbar(1:40,nanmean(repeat_num_saccades(:,1:40)),nanstd(repeat_num_saccades(:,1:40))...
    ./sqrt(sum(~isnan(repeat_num_saccades(:,1:40)))),'r');
p(2) = errorbar(nanmean(novel_num_saccades(:,1:40)),nanstd(novel_num_saccades(:,1:40))...
    ./sqrt(sum(~isnan(novel_num_saccades(:,1:40)))),'b');
p(3) = plot(sac_p05,6*ones(1,length(sac_p05)),'k*');
% p(4) = plot(sac_p01,6*ones(1,length(sac_p01)),'r*');
% p(5) = plot(sac_p001,6*ones(1,length(sac_p001)),'g*');
hold off
xlabel('Block #')
ylabel('# of Saccades')
% legend('Repeated','Novel','p<0.05','p < 0.01','p < 0.001','Location','SouthEast')
legend('Repeated','Novel','p<0.05','Location','SouthEast')
yl = ylim;
ylim([0 yl(2)]);
xlim([0 41])
%% Get Percent change from block 1 for last 10 blocks 
last_block_rt_means_repeat = 100*mean(repeat_reactiontimes(:,30:40),2)./mean(novel_reactiontimes(:,1));
last_block_rt_means_novel = 100*mean(novel_reactiontimes(:,30:40),2)./mean(novel_reactiontimes(:,1));

last_block_sac_means_repeat = 100*mean(repeat_num_saccades(:,30:40),2)...
    ./mean(novel_num_saccades(:,1));
last_block_sac_means_novel = 100*mean(novel_num_saccades(:,30:40),2)...
    ./mean(novel_num_saccades(:,1));

[~,p_rt] = ttest(last_block_rt_means_repeat,last_block_rt_means_novel,'tail','both');
[~,p_sac] = ttest(last_block_sac_means_repeat,last_block_sac_means_novel,'tail','both');

figure
subplot(1,2,1)
hold on
bar([mean(last_block_rt_means_novel) mean(last_block_rt_means_repeat)])
errorb(1,mean(last_block_rt_means_novel),std(last_block_rt_means_novel)...
    ./sqrt(length(last_block_rt_means_novel)))
errorb(2,mean(last_block_rt_means_repeat),std(last_block_rt_means_repeat)...
    ./sqrt(length(last_block_rt_means_repeat)))
if p_rt < 0.05
   plot(1.5, 100,'*k')
end
hold off
set(gca,'Xtick',1:2)
set(gca,'XtickLabel',{'Novel','Repeat'})
xlabel('Context')
ylabel('% of RT of novel block 1')

subplot(1,2,2)
hold on
bar([mean(last_block_sac_means_novel) mean(last_block_sac_means_repeat)])
errorb(1,mean(last_block_sac_means_novel),std(last_block_sac_means_novel)...
    ./sqrt(length(last_block_sac_means_novel)))
errorb(2,mean(last_block_sac_means_repeat),std(last_block_sac_means_repeat)...
    ./sqrt(length(last_block_sac_means_repeat)))
if p_sac < 0.05
   plot(1.5, 100,'*k')
end
hold off
set(gca,'Xtick',1:2)
set(gca,'XtickLabel',{'Novel','Repeat'})
xlabel('Context')
ylabel('% of Saccades of novel block 1')