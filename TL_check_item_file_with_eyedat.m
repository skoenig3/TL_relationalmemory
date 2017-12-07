% Written by Seth Konig 2/17/2015
% code used to check if item file location line up with eye data in case
% item files were confused or there was a typo in the log notes.

% code requires that file has already been preprocessed by getTLsData.m

TL_file = 'RR150212.2';
item_num = 435;%item # in log notes
real_item_num = 434;%suspected item #, can be same as above
cndfil = 'TL2.cnd';

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\TL_relationalmemory\Eyedat\';%directory where preprocessed data is stored

load([data_dir TL_file(1:2) '\' TL_file(1:end-2) 'set'...
    num2str(item_num) 'eyedat.mat']);

tloc = xlsread(['R:\Buffalo Lab\eblab\Tloc\Tloc'  num2str(real_item_num) '.xls']); %load TLs locations

%---Code reads in locations of T---%
fid=fopen(['R:\Buffalo Lab\eblab\Cortex Programs\Contextual\' cndfil],'rt');
tline=fgets(fid);
cndmrk=regexp(tline,'COND#'):regexp(tline,'COND#')+4;
tst0mrk=regexp(tline,'TEST0'):regexp(tline,'TEST0')+4;
tst1mrk=regexp(tline,'TEST1'):regexp(tline,'TEST1')+4;

cndind=[];
itmind=[];
tlinenew=0;
while tlinenew~=-1
    tlinenew=fgets(fid);
    if tlinenew~=-1
        cndind=[cndind str2double(tlinenew(cndmrk))];
        itmind=[itmind str2double(tlinenew(tst0mrk))];
    end
end

xpos=NaN(1,length(per));
ypos=NaN(1,length(per));
for trial = 1:24%length(per)
    cndnum=per(trial).cnd;
    itmnum=itmind(find(cndind==cndnum));
    xpos(trial)=tloc(itmnum,3);
    ypos(trial)=tloc(itmnum,4);
end
%%

for trial = 1:24
    XY = fixationstats{trial}.XY;
    events = per(trial).allval;
    times = per(trial).alltim;
    trial_start = times(find(events == 15));
    eyedata_start = times(find(events == 100));
    fix_on_cross = times(find(events == 8))-trial_start-eyedata_start;
    fix_on_T = fix_on_cross(end);
    fix_on_cross = fix_on_cross(1);
    stimulus_off = times(find(events == 24))-trial_start-eyedata_start;
    eyeind = [fix_on_cross stimulus_off]; 
    figure
    hold on
    plot(xpos(trial),ypos(trial),'k+','markersize',10)
    plot(XY(1,eyeind(1):eyeind(2)),XY(2,eyeind(1):eyeind(2)))
    plot(XY(1,fix_on_T),XY(2,fix_on_T),'r*','markersize',5)
    hold off
    xlim([-15 15])
    ylim([-12 12])
    pause
    close 
end
