%Draft written by Seth Konig 5/1/18 modified from Analyzed_All_TLs_Data.m
%to look at #of context learned

%%
%---[1] Import task data and detect fixations and saccades---%

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\TL_relationalmemory\Eyedat\'; %where to find preprocessed data


TL_files =      {'RR150129.3','RR150130.3','RR150202.3','RR150205.3','RR150206.3',...
                 'RR150209.2','RR150210.2','RR150211.2','RR150213.2','RR150217.2',...
                 'RR150219.2','RR150220.2','RR150223.2','RR150224.2','RR150225.2',...
                'TO150211.2','TO150218.2','TO150219.2','TO150223.2','TO150224.2',...
                'TO150225.2','TO150226.2','TO150302.2','TO150303.2','TO150305.2',...
                'TO150306.2','TO150309.2','TO150311.2','TO150312.2','TO150313.2',...
                'PW150518.2','PW150519.2','PW150520.2','PW150521.2','PW150522.2',...
                'PW150526.2','PW150527.2','PW150528.2','PW150529.2','PW150601.2',...
                'PW150602.2','PW150603.2','PW150604.2','PW150605.2','PW150608.2',...
                'MF160920.2','MF160921.2','MF160922.2','MF160923.2','MF160926.2',...
                'MF160927.2','MF161003.2','MF161004.2','MF161006.2','MF161007.2',...
                'MF161011.2','MF161013.2','MF161014.2','MF161017.2','MF161018.2',...
                'MF161019.2'};
             
item_num = [425 426 428 430 431,...
            432 433 434 435 436,...
            438 439 440 441 442,...
            475 478 479 480 481,...
            482 483 485 486 488,...
            489 490 492 493 494,...
            500 501 502 503 504,...
            505 506 507 508 509,...
            510 511 512 513 514,...
            425 426 428 429 430,...
            431 432 433 434 435,...
            436 437 438 439 440,...
            441 ] ;

%%
CNDFile = 'Z:\eblab\Cortex Programs\Contextual\TL2.cnd';
cndfil=[];
h=fopen(CNDFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(cndfil)
        if length(tline)>size(cndfil,2)
            tline=tline(1:size(cndfil,2));
        end
    end
    tline = [tline ones(1,(size(cndfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        cndfil=[cndfil; tline];
    else
        break
    end
end
fclose(h);

cnd2itm = [];
for line = 1:size(cndfil)-1;
    str = textscan(cndfil(line+1,:),'%d');
    cnd2itm(line) = str{1}(2);
end

%%
%---[2] Analyze reaction times and number of saccades to find target---%

max_desired_blocks = 40;%which blocks 1-max_desired_blocks to analyze

all_repeat_reaction_times = cell(1,length(TL_files));
all_repeat_num_saccdaes = cell(1,length(TL_files));
repeat_reactiontimes = NaN(length(TL_files),max_desired_blocks);%avearge reaction times by block
novel_reactiontimes  = NaN(length(TL_files),max_desired_blocks);%avearge reaction times by block
repeat_num_saccades  = NaN(length(TL_files),max_desired_blocks);%average number of saccades to find T
novel_num_saccades   = NaN(length(TL_files),max_desired_blocks);%average number of saccades to find T
for  file = 1:length(TL_files)
    
    disp(['Loading file # ' num2str(file) '....'])
    %load preprocessed Ts and Ls file
    load([data_dir TL_files{file}(1:2) '\' TL_files{file}(1:end-2) 'set'...
        num2str(item_num(file)) 'eyedat.mat']);
    
    %---quickly extract block number without using a for-loop---%
    block = struct2cell(per);
    block = reshape(block,size(block,1),size(block,3));
    block = cell2mat(block(4,:));
    
    all_repeat_reaction_times{file} = NaN(12,40);
    all_repeat_reaction_times{file} = NaN(12,40);
    
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
            row_index = [0 0]; %row 1 repeat, row 2 novel
            
            block_trials = find(block == blk);%find which trials go with this block
            for bt = 1:length(block_trials)
                trialcnd = per(block_trials(bt)).cnd;
                context_num = cnd2itm(trialcnd);
                
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
                fix_on_T = fix_on_T(end);%in case theres's more than 1 code which happened in a few of Red's sets, extra 24 for something else
                fix_on_cross = fix_on_cross(1);%first fixation is on crosshair to start trial
                stimulus_on =  times(find(events == 23))-trial_start-eyedata_start;
                stimulus_off = times(find(events == 24))-trial_start-eyedata_start;
                
                saccadetimes = fixationstats{block_trials(bt)}.saccadetimes;
                if ~isempty(saccadetimes)
                    saccadetimes(:,saccadetimes(1,:) < stimulus_on | saccadetimes(1,:) > fix_on_T) = [];
                end
                %find saccades that started after the stimulus turned on
                %but before the monkey fixated the T
                
                
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
                end
                
                if context_num <= 12
                    all_repeat_reaction_times{file}(context_num,blk) = fix_on_T-stimulus_on;
                    all_repeat_num_saccdaes{file}(context_num,blk) = size(saccadetimes,2);
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

%% Plot average reaction time and # of saccades to find T by block averaged by session

figure
subplot(1,2,1)
hold on
p(1) = errorbar(1:40,nanmean(repeat_reactiontimes(:,1:40)),nanstd(repeat_reactiontimes(:,1:40))...
    ./sqrt(sum(~isnan(repeat_reactiontimes(:,1:40)))),'r');
p(2) = errorbar(nanmean(novel_reactiontimes(:,1:40)),nanstd(novel_reactiontimes(:,1:40))...
    ./sqrt(sum(~isnan(novel_reactiontimes(:,1:40)))),'b');
hold off
xlabel('Block #')
ylabel('Reaction time (ms)')
legend('Repeated','Novel','Location','SouthEast')
yl = ylim;
ylim([500 yl(2)]);
xlim([0 41])

subplot(1,2,2)
hold on
p(1) = errorbar(1:40,nanmean(repeat_num_saccades(:,1:40)),nanstd(repeat_num_saccades(:,1:40))...
    ./sqrt(sum(~isnan(repeat_num_saccades(:,1:40)))),'r');
p(2) = errorbar(nanmean(novel_num_saccades(:,1:40)),nanstd(novel_num_saccades(:,1:40))...
    ./sqrt(sum(~isnan(novel_num_saccades(:,1:40)))),'b');
hold off
xlabel('Block #')
ylabel('# of Saccades')
legend('Repeated','Novel','Location','SouthEast')
yl = ylim;
ylim([1 yl(2)]);
xlim([0 41])

%% Double check means are approximately the same
repeat_mean_rt = NaN(length(TL_files),40);
for  file = 1:length(TL_files)
    repeat_mean_rt(file,:) = mean(all_repeat_reaction_times{file});
end

subplot(1,2,1)
hold on
plot(1:40,mean(repeat_mean_rt),'k')
hold off

repeat_mean_sac = NaN(length(TL_files),40);
for  file = 1:length(TL_files)
    repeat_mean_sac(file,:) = mean(all_repeat_num_saccdaes{file});
end

subplot(1,2,2)
hold on
plot(1:40,mean(repeat_mean_sac),'k')
hold off

%%
count = 0;
ratio_rt = NaN(1,12*length(TL_files));
ratio_sac = NaN(1,12*length(TL_files));
for file = 1:length(TL_files)
    for c = 1:12
        count = count+1;
        ratio_rt(count) =  mean(all_repeat_reaction_times{file}(c,31:40))/mean(all_repeat_reaction_times{file}(c,1:10));
        ratio_sac(count) =  mean(all_repeat_num_saccdaes{file}(c,31:40))/mean(all_repeat_reaction_times{file}(c,1:10));
    end
end
