% Compare PW's data for the 1st time she went through Ts and Ls to the 2nd
% time she went through Ts and Ls. Vivian completed all 400 level Ts and Ls
% and then started on them again. This analysis is to determine if she her
% behavior had changed due to practice, memory for the stimuli, etc.
% between the 1st and 2nd time she completed 400 level Ts and Ls. Not all
% the sets being analyzed are the same but are as close to the same as
% possible. 

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\TL_relationalmemory\Eyedat\';

load([data_dir 'PW_1st_pass.mat']);
novel_reactiontimes1 = novel_reactiontimes(:,1:40);
repeat_reactiontimes1 = repeat_reactiontimes(:,1:40); 

load([data_dir 'PW_2nd_pass.mat']);
novel_reactiontimes2 = novel_reactiontimes(:,1:40);
repeat_reactiontimes2 = repeat_reactiontimes(:,1:40); 

load([data_dir 'PW_3rd_pass.mat']);
novel_reactiontimes3 = novel_reactiontimes(:,1:40);
repeat_reactiontimes3 = repeat_reactiontimes(:,1:40); 

clear novel_reactiontimes repeat_reactiontimes

figure
hold on
errorbar(nanmean(novel_reactiontimes1),nanstd(novel_reactiontimes1)./sqrt(sum(~isnan(novel_reactiontimes1))),'b');
errorbar(nanmean(novel_reactiontimes2),nanstd(novel_reactiontimes2)./sqrt(sum(~isnan(novel_reactiontimes2))),'g');
errorbar(nanmean(novel_reactiontimes3),nanstd(novel_reactiontimes3)./sqrt(sum(~isnan(novel_reactiontimes3))),'k');
errorbar(nanmean(repeat_reactiontimes1),nanstd(repeat_reactiontimes1)./sqrt(sum(~isnan(repeat_reactiontimes1))),'r');
errorbar(nanmean(repeat_reactiontimes2),nanstd(repeat_reactiontimes2)./sqrt(sum(~isnan(repeat_reactiontimes2))),'m');
errorbar(nanmean(repeat_reactiontimes3),nanstd(repeat_reactiontimes3)./sqrt(sum(~isnan(repeat_reactiontimes3))),'c');
hold off
legend('Novel 1st pass','Novel 2nd pass','Novel 3rd pass','Repeat 1st pass','Repeat 2nd pass','Repeat 3rd pass')
xlabel('Block #')
ylabel('Reaction Time (ms)')

change1 = novel_reactiontimes1-repeat_reactiontimes1;
change2 = novel_reactiontimes2-repeat_reactiontimes2; 
change3 = novel_reactiontimes3-repeat_reactiontimes3;
change2(5,:) = [];%all NaNs

figure
hold on
errorbar(nanmean(change1),nanstd(change1)./sqrt(sum(~isnan(change1))),'b');
errorbar(nanmean(change2),nanstd(change2)./sqrt(sum(~isnan(change2))),'r');
errorbar(nanmean(change3),nanstd(change3)./sqrt(sum(~isnan(change3))),'g');
hold off
xlabel('Block #')
ylabel('Change in Reaction Time (ms)')
legend('1st Pass','2nd Pass','3rd Pass')

%%
all_blocks1 = [];%which block the data came from 
all_sessions1 = [];%which sessions within the data, the data came from
first_vs_second1 = [];%1st or second time running through the data

for row = 1:size(change1,1)
    for col = 1:size(change1,2)
        all_blocks1(row,col) = col;
        all_sessions1(row,col) = row;
        first_vs_second1(row,col) = 1;
    end
end

all_blocks2= [];%which block the data came from 
all_sessions2 = [];%which sessions within the data, the data came from
first_vs_second2 = [];%1st or second time running through the data

for row = 1:size(change2,1)
    for col = 1:size(change2,2)
        all_blocks2(row,col) = col;
        all_sessions2(row,col) = row;
        first_vs_second2(row,col) = 2;
    end
end

all_blocks3= [];%which block the data came from 
all_sessions3 = [];%which sessions within the data, the data came from
first_vs_second3 = [];%1st or second time running through the data

for row = 1:size(change3,1)
    for col = 1:size(change3,2)
        all_blocks3(row,col) = col;
        all_sessions3(row,col) = row;
        first_vs_second3(row,col) = 3;
    end
end
%%
all_blocks = [all_blocks1(1:end) all_blocks2(1:end) all_blocks3(1:end)];
all_sessions = [all_sessions1(1:end) all_sessions2(1:end) all_sessions3(1:end)];
first_vs_second = [first_vs_second1(1:end) first_vs_second2(1:end) first_vs_second3(1:end) ];
change = [change1(1:end) change2(1:end) change3(1:end)]; 

[P,T,STATS,TERMS]=anovan(change',[first_vs_second' all_blocks' ]);

[P,T,STATS,TERMS]=anovan(change',[first_vs_second' all_blocks'],'model','interaction',...
    'varnames',{'1p2p','block'}); 

[P,T,STATS,TERMS]=anovan(change',[all_sessions' first_vs_second' all_blocks'],'model','interaction',...
    'varnames',{'sess','1p2p','block'}); 
%%

last_block_rt_means_repeat1 = 100*mean(repeat_reactiontimes1(:,30:40),2)./mean(novel_reactiontimes1(:,1));
last_block_rt_means_novel1 = 100*mean(novel_reactiontimes1(:,30:40),2)./mean(novel_reactiontimes1(:,1));

last_block_rt_means_repeat2 = 100*nanmean(repeat_reactiontimes2(:,30:40),2)./nanmean(novel_reactiontimes2(:,1));
last_block_rt_means_novel2 = 100*nanmean(novel_reactiontimes2(:,30:40),2)./nanmean(novel_reactiontimes2(:,1));

last_block_rt_means_repeat3 = 100*mean(repeat_reactiontimes3(:,30:40),2)./mean(novel_reactiontimes3(:,1));
last_block_rt_means_novel3 = 100*mean(novel_reactiontimes3(:,30:40),2)./mean(novel_reactiontimes3(:,1));



figure
hold on
bar([mean(last_block_rt_means_novel1) mean(last_block_rt_means_repeat1),...
    nanmean(last_block_rt_means_novel2) nanmean(last_block_rt_means_repeat2),...
    mean(last_block_rt_means_novel3) mean(last_block_rt_means_repeat3)])
errorb(1,mean(last_block_rt_means_novel1),std(last_block_rt_means_novel1)...
    ./sqrt(length(last_block_rt_means_novel1)))
errorb(2,mean(last_block_rt_means_repeat1),std(last_block_rt_means_repeat1)...
    ./sqrt(length(last_block_rt_means_repeat1)))

errorb(3,nanmean(last_block_rt_means_novel2),nanstd(last_block_rt_means_novel2)...
    ./sqrt(length(last_block_rt_means_novel2)))
errorb(4,nanmean(last_block_rt_means_repeat2),nanstd(last_block_rt_means_repeat2)...
    ./sqrt(length(last_block_rt_means_repeat2)))
errorb(5,mean(last_block_rt_means_novel3),std(last_block_rt_means_novel3)...
    ./sqrt(length(last_block_rt_means_novel3)))
errorb(6,mean(last_block_rt_means_repeat3),std(last_block_rt_means_repeat3)...
    ./sqrt(length(last_block_rt_means_repeat3)))
hold off
set(gca,'Xtick',1:6)
set(gca,'XtickLabel',{'Novel1','Repeat1','Novel2','Repeat2','Novel3','Repeat3'})
xlabel('Context/Pass')
ylabel('% of RT of novel block 1')
%%
[~,p_nov12] = ttest2(last_block_rt_means_novel1,last_block_rt_means_novel2,'tail','both')
[~,p_nov13] = ttest2(last_block_rt_means_novel1,last_block_rt_means_novel3,'tail','both')
[~,p_nov23] = ttest2(last_block_rt_means_novel2,last_block_rt_means_novel3,'tail','both')
%%
[~,p_rep12] = ttest2(last_block_rt_means_repeat1,last_block_rt_means_repeat2,'tail','both')
[~,p_rep13] = ttest2(last_block_rt_means_repeat1,last_block_rt_means_repeat3,'tail','both')
[~,p_rep23] = ttest2(last_block_rt_means_repeat2,last_block_rt_means_repeat3,'tail','both')