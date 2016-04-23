% Power analysis for TL's sets. written by Seth Konig 1/19/2015
% much of the code (the guts) is borrow from TLmodeyedat_avgs.m
% Scrapt 2/15/2015 based on code from network which had bugs. Removed bugs
% but started to use Cluster Fix and my own code which likely works better
% anyway.

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\TL_relationalmemory\Eyedat\';
monkeys = {'RR'};

rts = cell(2,length(monkeys)); %reaction time for finding the T
numsacs = cell(2,length(monkeys)); %# of saccades before finding the T
% row 1 novel trials, row 2 repeat trials

for m = 1:length(monkeys)
    cd([data_dir monkeys{m}]);%directory containing the individual monkeys' data
    files = what; %the files in this directory
    files = files.mat;%the .mat files in this directory
    
    %create array of NaNs to store future data
    for novrep = 1:2
        rts{novrep,m} = NaN(length(files),40);
        numsacs{novrep,m} = NaN(length(files),40); 
    end
    
    for f = 1:length(files);
        load(files{f}); %import the single days data
        
        trltyp1=cfg.trl((cfg.trl(:,6)==1),:); %repeated configurations
        trltyp2=cfg.trl((cfg.trl(:,6)==2),:); %novel configurations
        
        blkarr=unique(trltyp1(:,5));
        numblks=length(blkarr);
        if numblks < 40 %if didn't finish 40 blocks go to next data set
            continue
        else %will only analize 1st 40 blocks
            numblks = 40;
        end
        
        % get reaction times and number of saccdes by block 
        for k=1:numblks  
            % for repeated contexts
            sac1=trltyp1(trltyp1(:,5)==blkarr(k),7);
            rcttim1=trltyp1(trltyp1(:,5)==blkarr(k),2)-trltyp1(trltyp1(:,5)==blkarr(k),1);
            sac1=sac1(rcttim1<8000);
            rcttim1=rcttim1(rcttim1<8000);
            %remove trials in which took way to long may not be paying attention (was here before)

            %  novel contexts
            sac2=trltyp2(trltyp2(:,5)==blkarr(k),7);
            rcttim2=trltyp2(trltyp2(:,5)==blkarr(k),2)-trltyp2(trltyp2(:,5)==blkarr(k),1);
            sac2=sac2(rcttim2<8000);
            rcttim2=rcttim2(rcttim2<8000);
            %remove trials in which took way to long may not be paying attention (was here before)
            
            %store averages per block
            rts{1,m}(f,k) = mean(rcttim1);
            rts{2,m}(f,k) = mean(rcttim2); 
            numsacs{1,m}(f,k) = mean(sac1);
            numsacs{2,m}(f,k) = mean(sac2);
        end
    end
end
%% plot data by block number by monkey

for m = 1:length(monkeys)
   figure
   subplot(1,2,1)
   hold on
   errorbar(1:40,nanmedian(rts{1,m}),nanstd(rts{1,m})...
       ./sqrt(sum(~isnan(rts{1,m}))),'r')
     errorbar(nanmedian(rts{2,m}),nanstd(rts{2,m})...
       ./sqrt(sum(~isnan(rts{2,m}))),'b')
   hold off
   xlabel('Block #')
   ylabel('Reaction time (ms)')
   legend('Repeated','Novel')
   yl = ylim;
   ylim([0 yl(2)]);
   xlim([0 41])
   
   subplot(1,2,2)
   hold on
   errorbar(1:40,nanmedian(numsacs{1,m}),nanstd(numsacs{1,m})...
       ./sqrt(sum(~isnan(numsacs{1,m}))),'r')
   errorbar(nanmedian(numsacs{2,m}),nanstd(numsacs{2,m})...
       ./sqrt(sum(~isnan(numsacs{2,m}))),'b')
   hold off
   xlabel('Block #')
   ylabel('# of Saccades')
   legend('Repeated','Novel')
   yl = ylim;
   ylim([0 yl(2)]);
   xlim([0 41])
   
   subtitle(['Data for ' monkeys{m}])
end

%% plot data by block number across all monkeys
% combine data across all monkeys to get general population sets
allrts = cell(1,2);%combined reaction times across all monkeys
allnumsacs = cell(1,2);%combined number of saccades to T across all monkeys 
for m = 1:length(monkeys)
    for novrep = 1:2;
       allrts{novrep} = [allrts{novrep}; rts{novrep,m}];
       allnumsacs{novrep} = [allnumsacs{novrep}; numsacs{novrep,m}];
    end
end

figure
subplot(1,2,1)
hold on
errorbar(1:40,nanmedian(allrts{1}),nanstd(allrts{1})...
    ./sqrt(sum(~isnan(allrts{1}))),'r')
errorbar(nanmedian(allrts{2}),nanstd(allrts{2})...
    ./sqrt(sum(~isnan(allrts{2}))),'b')
hold off
xlabel('Block #')
ylabel('Reaction time (ms)')
legend('Repeated','Novel')
yl = ylim;
ylim([0 yl(2)]);
xlim([0 41])

subplot(1,2,2)
hold on
errorbar(1:40,nanmedian(allnumsacs{1}),nanstd(allnumsacs{1})...
    ./sqrt(sum(~isnan(allnumsacs{1}))),'r')
errorbar(nanmedian(allnumsacs{2}),nanstd(allnumsacs{2})...
    ./sqrt(sum(~isnan(allnumsacs{2}))),'b')
hold off
xlabel('Block #')
ylabel('# of Saccades')
legend('Repeated','Novel')
yl = ylim;
ylim([0 yl(2)]);
xlim([0 41])

subtitle('Combined Data Across All monkeys')
% %% calculate the number of samples needed using power analysis on individual monkeys
% % note the variability (i.e. std) is greater for novel contexts than
% % repeated contexts. I am using the std from the novel contexts for the
% % measure of the null distribution.
% rt_num_samps = NaN(length(monkeys),40);
% sacs_num_samps = NaN(length(monkeys),40);
% 
% for m = 1:length(monkeys)
%     for blk = 1:40
%         rt_num_samps(m,blk) = sampsizepwr('t',[nanmedian(rts{2,m}(:,blk)),...
%             nanstd(rts{2,m}(:,blk))],nanmedian(rts{1,m}(:,blk)),0.9);
%         sacs_num_samps(m,blk) = sampsizepwr('t',[nanmedian(numsacs{2,m}(:,blk)),...
%             nanstd(numsacs{2,m}(:,blk))],nanmedian(numsacs{1,m}(:,blk)),0.9);
%     end
% end
% %% calculate the number of samples needed using power analysis on individual monkeys
% % note the variability (i.e. std) is greater for novel contexts than
% % repeated contexts. I am using the std from the novel contexts for the
% % measure of the null distribution.
% rt_num_samps = NaN(length(monkeys),40);
% sacs_num_samps = NaN(length(monkeys),40);
% 
% for m = 1:length(monkeys)
%     for blk = 1:40
%         rt_num_samps(m,blk) = sampsizepwr('t',[nanmedian(rts{2,m}(:,blk)),...
%             nanstd(rts{2,m}(:,blk))],nanmedian(rts{1,m}(:,blk)),0.8);
%         sacs_num_samps(m,blk) = sampsizepwr('t',[nanmedian(numsacs{2,m}(:,blk)),...
%             nanstd(numsacs{2,m}(:,blk))],nanmedian(numsacs{1,m}(:,blk)),0.8);
%     end
% end
% %% calculate the number of samples needed using power analysis on
% % population data (i.e. combined across all monkeys)
% all_rt_num_samps = NaN(1,40);
% all_sacs_num_samps = NaN(1,40);
% 
% for blk = 20:40;
%     all_rt_num_samps(blk) = sampsizepwr('t',[nanmedian(allrts{2}(:,blk)),...
%         nanstd(allrts{2}(:,blk))],nanmedian(allrts{1}(:,blk)),0.8);
%     all_sacs_num_samps(blk) = sampsizepwr('t',[nanmedian(allnumsacs{2}(:,blk)),...
%         nanstd(allnumsacs{2}(:,blk))],nanmedian(allnumsacs{1}(:,blk)),0.8);
% end
