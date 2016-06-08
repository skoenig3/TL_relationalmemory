function getTLsData(TL_cortexfile,clrchng_cortexfile,itemnum)
% written 2/9/15 by Seth Konig
% Imports a color change calibration file and cortex file from the TLs
% task into Matlab. Then the function parses the task data, detects fixations
% and saccades using Cluster Fix.
%
% Inputs:
%   1) TL_cortexfile: cortex file for the Ts & Ls task
%   2) clrchng_cortexfile: cortex file for the clrchng calibration task
%   3) itemnum: Item number for the Ts & Ls task e.g. 415
%
% Outputs:
%   1) Matlab file with the saved eye data including fixations and saccade
%   times as well as pupil data (if it exists). Save Directory on line 341!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calibration Stuff---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

samprate = 5;%number of ms between samples for ISCAN i.e. 200 Hz

%---Import Color Change Data for Calibration---%
init = clrchng_cortexfile(1:2); %monkey initials
if strcmpi(init,'IW')==1
    clrchng_cortexfile=['R:\Buffalo Lab\Cortex Data\Irwin\' clrchng_cortexfile];
elseif strcmpi(init,'MP')==1
    clrchng_cortexfile=['R:\Buffalo Lab\Cortex Data\Peepers\' clrchng_cortexfile];
elseif strcmpi(init,'WR')==1
    clrchng_cortexfile=['R:\Buffalo Lab\Cortex Data\Wilbur\' clrchng_cortexfile];
elseif strcmpi(init,'TT')==1 
    clrchng_cortexfile=['R:\Buffalo Lab\Cortex Data\Timmy\' clrchng_cortexfile];
elseif strcmpi(init,'JN')==1
    clrchng_cortexfile=['R:\Buffalo Lab\Cortex Data\Guiseppe\' clrchng_cortexfile];
elseif strcmpi(init,'TD')==1 
    clrchng_cortexfile=['R:\Buffalo Lab\Cortex Data\Theodore\' clrchng_cortexfile];
elseif strcmpi(init,'PW')==1
    clrchng_cortexfile=['R:\Buffalo Lab\Cortex Data\Vivian\' clrchng_cortexfile];
elseif strcmpi(init,'RR')==1
    clrchng_cortexfile=['R:\Buffalo Lab\Cortex Data\Red\' clrchng_cortexfile];
elseif strcmpi(init,'TO')==1
    clrchng_cortexfile=['R:\Buffalo Lab\Cortex Data\Tobii\' clrchng_cortexfile];
elseif strcmpi(init,'MF') == 1
    clrchng_cortexfile=['R:\Buffalo Lab\Cortex Data\Manfred\' clrchng_cortexfile];
end

ITMFile = 'R:\Buffalo Lab\eblab\Cortex Programs\ClrChng\cch25.itm';
CNDFile = 'R:\Buffalo Lab\eblab\Cortex Programs\ClrChng\cch25.cnd';
% this is different becasue the spacing is different and I don't have
% a new item file on the network for the new spacing
ind_spacex = [-6,-3,0,3,6]; %whats on the network
ind_spacey = [-6,-3,0,3,6];%whats on the network

if strcmpi(init,'RR') && str2double(clrchng_cortexfile(end-7:end-2)) < 150301
    % red's computer had smaller spacing matching network code
    spacex = [-6,-3,0,3,6]; %whats on the network
    spacey = [-6,-3,0,3,6];%whats on the network
else %the rest had spacing good for full screen
    spacex = [-12,-6,0,6,12];%what actually gets displayed
    spacey = [-8,-4,0,4,8];%what actually gets displayed
end

[time_arr,event_arr,eog_arr,~,~,~]  = ...
    get_ALLdata(clrchng_cortexfile);

%---read in item file---%
itmfil=[];
h =fopen(ITMFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(itmfil)
        if length(tline)>size(itmfil,2)
            tline=tline(1:size(itmfil,2));
        end
    end
    tline = [tline ones(1,(size(itmfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        itmfil=[itmfil; tline];
    else
        break
    end
end
fclose(h);

%---Read in condition file---%
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

itmlist = zeros(size(cndfil,1)-1,1);
for i = 2:size(cndfil,1);
    str = textscan(cndfil(i,:),'%d');
    itmlist(i-1) = str{1}(end);
end

numrpt = size(event_arr,2);
valrptcnt = 0;
clear per clrchgind
for rptlop = 1:numrpt
    if itmlist(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop)-1000) <= 189
        if size(find(event_arr(:,rptlop) == 200)) ~=0
            perbegind = find(event_arr(:,rptlop) == 24);%was originally 23, changed this and begtimdum line below to optimize
            perendind = find(event_arr(:,rptlop) == 24);
            if length( perbegind) > 1
                perbegind = perbegind(2);
                perendind = perendind(2);
            end
            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
            blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
            begtimdum = time_arr(perbegind,rptlop)-100;
            endtimdum = time_arr(perendind,rptlop);
            if endtimdum > begtimdum
                valrptcnt = valrptcnt + 1;
                clrchgind(valrptcnt)=rptlop;
                per(valrptcnt).begsmpind = begtimdum;
                per(valrptcnt).endsmpind = endtimdum;
                per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                per(valrptcnt).blk = event_arr(blknumind,rptlop);
                per(valrptcnt).allval = event_arr(:,rptlop);
                per(valrptcnt).alltim = time_arr(:,rptlop);
                per(valrptcnt).event = rptlop;
            end
        end
    end
end

clear cnd
numrpt = size(per,2);
cnd = zeros(1,numrpt);
for rptlop = 1:numrpt
    cnd(rptlop)=per(rptlop).cnd;
end

% Create structures x and y of the corresponding average eye data for each trial
% instance (l) of each condition (k)

x = cell(length(spacex),length(spacey));%For Calibration with Eye tracking data with cp2tform
y = cell(length(spacex),length(spacey));
control = NaN(length(cnd),2);
clr = ['rgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmk'];
figure
hold on
for k = 1:length(cnd)
    C = textscan(itmfil(itmlist(cnd(k)-1000)+5,:),'%d');
    control(k,:) = C{1}(9:10)';
    
    xi = find(C{1}(9) == ind_spacex);
    yi = find(C{1}(10) == ind_spacey);
    eyeind = floor(((per(k).begsmpind-1000)/samprate)*2):(floor((per(k).endsmpind-1000)/samprate))*2;
    evenind = eyeind(logical(~rem(eyeind,2)));
    oddind =  eyeind(logical(rem(eyeind,2)));
    x{xi,yi} = [x{xi,yi} mean(eog_arr(oddind,per(k).event))];
    y{xi,yi} = [y{xi,yi} mean(eog_arr(evenind,per(k).event))];
    plot(mean(eog_arr(oddind,per(k).event)),mean(eog_arr(evenind,per(k).event)),[clr(xi*yi) '+'])
end
title(['Calibration transformation for ' clrchng_cortexfile(end-9:end)])


%Test for errors%
count = zeros(length(spacey),length(spacex));
for xi = 1:length(spacex);
    for yi = 1:length(spacey);
        count(yi,xi) = sum(control(:,1) == ind_spacex(xi) & control(:,2) == ind_spacey(yi));
    end
end
if any(count < 5);
    disp('Calibration trial analysis incomplete or error')
    disp('Check number of calibration pionts or task not finished')
end

clear meanx meany
for k=1:numel(x)
    xss = x{k};
    low = mean(xss)-std(xss);
    high = mean(xss)+std(xss);
    xss(xss < low) = [];
    xss(xss > high) = [];
    meanx(k)=median(xss);
end
for k=1:numel(y)
    yss = y{k};
    low = mean(yss)-std(yss);
    high = mean(yss)+std(yss);
    yss(yss < low) = [];
    yss(yss > high) = [];
    meany(k)=median(y{k});
end

controlx = [];
controly = [];
for i = 1:length(spacex);
    for ii = 1:length(spacey);
        controly = [controly spacey(i)];
        controlx = [controlx spacex(ii)];
    end
end


tform = cp2tform([controlx' controly'], [meanx' meany'],'affine');
tform.forward_fcn = tform.inverse_fcn;

newx = [];
newy = [];
figure
hold on
for i = 1:length(controlx);
    plot(controlx(i),controly(i),'r+')
    [x,y] = tformfwd(tform,meanx(i),meany(i));
    plot(x,y,'*b')
    newx(i) = x;
    newy(i) = y;
end
if iscell(clrchng_cortexfile)
    title(['Calibration transformation for ' clrchng_cortexfile{1}(end-9:end)])
else
    title(['Calibration transformation for ' clrchng_cortexfile(end-9:end)])
end
xlim([-17.5 17.5])
ylim([-12.5 12.5])


%%%%%%%%%%%%%%%%%%%%%
%%%---TLs Stuff---%%%
%%%%%%%%%%%%%%%%%%%%%

%---Import Ts and Ls data---%
init = TL_cortexfile(1:2); %monkey initials
if strcmpi(init,'IW')==1 || strcmpi(init,'iw')==1
    TL_cortexfile=['R:\Buffalo Lab\Cortex Data\Irwin\' TL_cortexfile];
elseif strcmpi(init,'MP')==1 || strcmpi(init,'mp')==1
    TL_cortexfile=['R:\Buffalo Lab\Cortex Data\Peepers\' TL_cortexfile];
elseif strcmpi(init,'WR')==1 || strcmpi(init,'wr')==1
    TL_cortexfile=['R:\Buffalo Lab\Cortex Data\Wilbur\' TL_cortexfile];
elseif strcmpi(init,'TT')==1 || strcmpi(init,'tt')==1
    TL_cortexfile=['R:\Buffalo Lab\Cortex Data\Timmy\' TL_cortexfile];
elseif strcmpi(init,'JN')==1 || strcmpi(init,'jn')==1
    TL_cortexfile=['R:\Buffalo Lab\Cortex Data\Guiseppe\' TL_cortexfile];
elseif strcmpi(init,'TD')==1 || strcmpi(init,'td')==1
    TL_cortexfile=['R:\Buffalo Lab\Cortex Data\Theodore\' TL_cortexfile];
elseif strcmpi(init,'PW')==1 || strcmpi(init,'pw')==1
    TL_cortexfile=['R:\Buffalo Lab\Cortex Data\Vivian\' TL_cortexfile];
elseif strcmpi(init,'RR')==1 || strcmpi(init,'rr')==1
    TL_cortexfile=['R:\Buffalo Lab\Cortex Data\Red\' TL_cortexfile];
elseif strcmpi(init,'TO')==1 || strcmpi(init,'to')==1
    TL_cortexfile=['R:\Buffalo Lab\Cortex Data\Tobii\' TL_cortexfile];
elseif strcmpi(init,'MF') == 1
    TL_cortexfile=['R:\Buffalo Lab\Cortex Data\Manfred\' TL_cortexfile];
end
[time_arr,event_arr,eog_arr,epp_arr,~,~]  = get_ALLdata(TL_cortexfile);

%tloc = xlsread(['R:\Buffalo Lab\eblab\Tloc\Tloc'  num2str(itemnum) '.xls']);
%more lines are necessary for reading file information column 3 xpos and
%column 4 ypos. Row by column number


numrpt = size(event_arr,2);
new_eog_arr = [];
if ~isempty(epp_arr);%if there is pupil data
    new_epp_arr = [];
end
valrptcnt = 0;
for rptlop = 1:numrpt
    if size(find(event_arr(:,rptlop) == 284)) ~=0 %fine if trial finished not necessaril rewarded
        perbegind = find(event_arr(:,rptlop) == 100); % eye data starts at 100; might have predictive looking
        perendind = find(event_arr(:,rptlop) == 101); % eye data stops collecting after rewards so can stop here
        cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=4999);
        blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
        typnumind = find(event_arr(:,rptlop) == 1 | event_arr(:,rptlop) == 2);
        %1 for repeated contexts, 2 for novel contexts
        begtimdum = time_arr(perbegind,rptlop);
        endtimdum = time_arr(perendind(end),rptlop);
        if endtimdum > begtimdum
            valrptcnt = valrptcnt + 1;
            per(valrptcnt).begsmpind = begtimdum;
            per(valrptcnt).endsmpind = endtimdum(1);
            per(valrptcnt).begpos = 1;
            per(valrptcnt).cnd = event_arr(cndnumind,rptlop)-1000; %condition number
            per(valrptcnt).blk = event_arr(blknumind,rptlop)-500; %block number
            per(valrptcnt).allval = event_arr(:,rptlop); %cortex codes
            per(valrptcnt).alltim = time_arr(:,rptlop); %cortex timestamps for those codes
            per(valrptcnt).event = rptlop; %trial number (not necessairly successful trial number)
            per(valrptcnt).trialtype = event_arr(typnumind,rptlop); %novel (2)/repeat (1)
            new_eog_arr=cat(2,new_eog_arr,eog_arr(:,rptlop));
            if ~isempty(epp_arr);
                new_epp_arr = cat(2,new_epp_arr,epp_arr(:,rptlop));
            end
        end
    end
end

%---get eye data for only when fixation cross or picture is displayed---%
eyedat = cell(1,length(per));
if ~isempty(epp_arr);
    pupildata = cell(1,length(per)); %horizontal pupil diameter
end
cnd=[];
teststart = [];
for trlop=1:size(per,2)
    trleog=new_eog_arr(~isnan(new_eog_arr(:,trlop)),trlop); % eog for this trial
    horeog=trleog(1:2:size(trleog,1)); % horizontal eye dat
    vrteog=trleog(2:2:size(trleog,1)); % vertical eye dat
    
    picstart=1*samprate;
    picend=per(trlop).endsmpind-per(trlop).begsmpind;%added in case sometimes get weird
    %indexing artifacts that are off by 1 index due to cortex having a
    %clock speed with a 1 ms resoultion and the eye data collected at a 5 ms
    %resoultion
    
    eyedat{trlop}(1,:) = horeog(round(picstart/samprate):floor(picend/samprate));
    eyedat{trlop}(2,:) = vrteog(round(picstart/samprate):floor(picend/samprate));
    
    if ~isempty(epp_arr);
        trlepp=new_epp_arr(~isnan(new_epp_arr(:,trlop)),trlop); % eog for this trial
        pupildata{trlop} = trlepp(2:2:size(trlepp,1));%odd indexes contains nothing but noise
    end
end

%---Recalibrate and automatically scale eye data---%
for eye = 1:length(eyedat)
    x = eyedat{eye}(1,:);
    y = eyedat{eye}(2,:);
    [x,y] = tformfwd(tform,x,y);
    eyedat{eye} = [x;y];
end

%---Detect Fixations and Saccades with Cluster Fix---%
%re-saves eye data upsampled at 1 ms resultion. I don't usually plot this
%but can. The start of fixations/saccades is also at 1 ms resolution.
fixationstats = ClusterFixation_Short(eyedat);


%1 directory per monkey
data_dir = ['C:\Users\seth.koenig\Documents\MATLAB\TL_relationalmemory\Eyedat\' init '\'];
if ~isdir(data_dir)
    mkdir(data_dir); %if directory does not exist make one
end
filename = [TL_cortexfile(end-9:end-2)  'set' num2str(itemnum) 'eyedat.mat'];

if isempty(epp_arr)
    pupildata = [];
end
save([data_dir filename],'fixationstats','pupildata','per');
disp(['Generated:' filename])
end
