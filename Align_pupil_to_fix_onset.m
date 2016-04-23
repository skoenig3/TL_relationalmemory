datafile = 'RR150212.2';
init = datafile(1:2);
if strcmp(init,'IW')==1 || strcmp(init,'iw')==1
    datafile=['R:\Buffalo Lab\Cortex Data\Irwin\' datafile];
elseif strcmp(init,'MP')==1 || strcmp(init,'mp')==1
    datafile=['R:\Buffalo Lab\Cortex Data\Peepers\' datafile];
elseif strcmp(init,'WR')==1 || strcmp(init,'wr')==1
    datafile=['R:\Buffalo Lab\Cortex Data\Wilbur\' datafile];
elseif strcmp(init,'TT')==1 || strcmp(init,'tt')==1
    datafile=['R:\Buffalo Lab\Cortex Data\Timmy\' datafile];
elseif strcmp(init,'JN')==1 || strcmp(init,'jn')==1
    datafile=['R:\Buffalo Lab\Cortex Data\Guiseppe\' datafile];
elseif strcmp(init,'TD')==1 || strcmp(init,'td')==1
    datafile=['R:\Buffalo Lab\Cortex Data\Theodore\' datafile];
elseif strcmp(init,'PW')==1 || strcmp(init,'pw')==1
    datafile=['R:\Buffalo Lab\Cortex Data\Vivian\' datafile];
elseif strcmp(init,'RR')==1 || strcmp(init,'rr')==1
    datafile=['R:\Buffalo Lab\Cortex Data\Red\' datafile];
elseif strcmp(init,'TO')==1 || strcmp(init,'to')==1
    datafile=['R:\Buffalo Lab\Cortex Data\Tobii\' datafile];
end

[time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata(datafile);
pupil_diametery = NaN(size(epp_arr,2),1800);
pupil_diameterx = NaN(size(epp_arr,2),1800);
for trial = 1:size(epp_arr,2);
    if ~isempty(find(event_arr(:,trial) == 3)) %so trial was rewarded
        trial_start = time_arr(find(event_arr(:,trial) == 15),trial);
        epp_start = time_arr(find(event_arr(:,trial) == 100),trial);
        fix_on_cross = time_arr(find(event_arr(:,trial) == 8),trial);
        fix_on_cross = fix_on_cross(1);
        trial_end = time_arr(find(event_arr(:,trial) == 24),trial)-trial_start-epp_start;
        start_index = fix_on_cross-trial_start-epp_start;
        if trial_end-fix_on_cross > 1500
            trial_end = start_index+1499; 
        end
        epp_x = epp_arr(1:2:end,trial);
        epp_y = epp_arr(2:2:end,trial);

        epp_ind = round((start_index-300:trial_end)/5);
        epp_x = epp_x(epp_ind);
        epp_y = epp_y(epp_ind);
        pupil_diametery(trial,1:length(epp_y)) = epp_y;
        pupil_diameterx(trial,1:length(epp_x)) = epp_x;
    end
end

figure
plot(nanmean(pupil_diametery(:,1:1800)))
title(datafile(end-9:end))