dataDir = './Data_PA2';
data = fullfile(dataDir,'ilham-act-new.mat');
load (data);
alpha = 0.25;
% IPPG1 = filter(alpha, [1 alpha-1], IPPG_1);
% IPPG2 = filter(alpha, [1 alpha-1], IPPG_2);
% IPPG3 = filter(alpha, [1 alpha-1], IPPG_3);
% IPPG4 = filter(alpha, [1 alpha-1], IPPG_4);
Peakwidth = 0.09;
minHeight = 0.0004;
minval = 0; % For finding peaks with positive value
minDist = 0.3; % For tresholding up to 180 bpm-> 3 bps = 1/3
filter1 = 0.0001;
filter2 = 0.0001;
filter3 = 0.0001;
filter4 = 0.0005;
%% For EVM
% filter1 = 0.02;
% filter2 = 0.015;
% filter3 = 0.005;
% filter4 = 0.02;
%%
% figure(1)
% [pks_1,locs_1]=findpeaks(IPPG_1,t,'MinPeakProminence',filter1);findpeaks(IPPG_1,t,'MinPeakProminence',filter1);title('HR Green Channel');
% figure(2)
% [pks_2,locs_2]=findpeaks(IPPG_2,t,'MinPeakProminence',filter2);findpeaks(IPPG_2,t,'MinPeakProminence',filter2);title('HR Green-Red');
% figure(3)
% [pks_3,locs_3]=findpeaks(IPPG_3,t,'MinPeakProminence',filter3);findpeaks(IPPG_3,t,'MinPeakProminence',filter3);title('HR CHROM');
% figure(4)
% [pks_4,locs_4]=findpeaks(IPPG_4,t,'MinPeakProminence',filter4);findpeaks(IPPG_4,t,'MinPeakProminence',filter4);title('HR POS');
%% Normal
% figure(1)
% [pks_1,locs_1]=findpeaks(iPPG_1,t,'MinPeakProminence',filter1);findpeaks(iPPG_1,t,'MinPeakProminence',filter1);title('HR Green Channel');
% figure(2)
% [pks_2,locs_2]=findpeaks(iPPG_2,t,'MinPeakProminence',filter2);findpeaks(iPPG_2,t,'MinPeakProminence',filter2);title('HR Green-Red');
% figure(3)
% [pks_3,locs_3]=findpeaks(iPPG_3,t,'MinPeakProminence',filter3);findpeaks(iPPG_3,t,'MinPeakProminence',filter3);title('HR CHROM');
% figure(4)
% [pks_4,locs_4]=findpeaks(iPPG_4,t,'MinPeakProminence',filter4);findpeaks(iPPG_4,t,'MinPeakProminence',filter4);title('HR POS');
%% Normal + MA(5)
startVal = min(t);
endVal = max(t);
t_interp=0.01:0.01:endVal;
PPG_1 = interp1(t, IPPG_1, t_interp, 'spline');
PPG_2 = interp1(t, IPPG_2, t_interp, 'spline');
PPG_3 = interp1(t, IPPG_3, t_interp, 'spline');
PPG_4 = interp1(t, IPPG_4, t_interp, 'spline');
PPG_1 = movmean(PPG_1,20);
PPG_2 = movmean(PPG_2,20);
PPG_3 = movmean(PPG_3,20);
PPG_4 = movmean(PPG_4,20);
%alpha = 0.28;
% PPG_1 = filter(alpha, [1 alpha-1], PPG_1);
% PPG_2 = filter(alpha, [1 alpha-1], PPG_2);
% PPG_3 = filter(alpha, [1 alpha-1], PPG_3);
% PPG_4 = filter(alpha, [1 alpha-1], PPG_4);
%%
figure(1)
axes1 = axes(figure(1));
[pks_1,locs_1]=findpeaks(PPG_1,t_interp,'MinPeakProminence',filter1,'MinPeakWidth',Peakwidth,'MinPeakDistance',minDist,'MinPeakHeight',minHeight);findpeaks(PPG_1,t_interp,'MinPeakProminence',filter1,'MinPeakWidth',Peakwidth,'MinPeakDistance',minDist,'MinPeakHeight',minHeight);title('HR Green Channel');
xlim([0 60]);
xlabel('Time (second)');
ylabel('Amplitude');
set(axes1,'FontSize',24);
set(axes1,'NextPlot','replaceall');
%%
figure(2)
axes2 = axes(figure(2));
[pks_2,locs_2]=findpeaks(PPG_2,t_interp,'MinPeakProminence',filter2,'MinPeakWidth',Peakwidth,'MinPeakDistance',minDist,'MinPeakHeight',minHeight);findpeaks(PPG_2,t_interp,'MinPeakProminence',filter2,'MinPeakWidth',Peakwidth,'MinPeakDistance',minDist,'MinPeakHeight',minHeight);title('HR Green-Red');
xlim([0 60]);
xlabel('Time (second)');
ylabel('Amplitude');
set(axes2,'FontSize',24);
set(axes2,'NextPlot','replaceall');
%%
figure(3)
axes3 = axes(figure(3));
[pks_3,locs_3]=findpeaks(PPG_3,t_interp,'MinPeakProminence',filter3,'MinPeakWidth',Peakwidth,'MinPeakDistance',minDist,'MinPeakHeight',minHeight);findpeaks(PPG_3,t_interp,'MinPeakProminence',filter3,'MinPeakWidth',Peakwidth,'MinPeakDistance',minDist,'MinPeakHeight',minHeight);title('HR CHROM');
xlim([0 60]);
xlabel('Time (second)');
ylabel('Amplitude');
set(axes3,'FontSize',24);
set(axes3,'NextPlot','replaceall');
%%
figure(4)
axes4 = axes(figure(4));
[pks_4,locs_4]=findpeaks(PPG_4,t_interp,'MinPeakProminence',filter4,'MinPeakWidth',Peakwidth,'MinPeakDistance',minDist,'MinPeakHeight',minHeight);findpeaks(PPG_4,t_interp,'MinPeakProminence',filter4,'MinPeakWidth',Peakwidth,'MinPeakDistance',minDist,'MinPeakHeight',minHeight);title('HR POS');
xlim([0 60]);
xlabel('Time (second)');
ylabel('Amplitude');
set(axes4,'FontSize',24);
set(axes4,'NextPlot','replaceall');
%%
figure(5)
axes5 = axes('Parent',figure(5));
plot1 = plot(t,rawColorSignal,'LineWidth',3,'Parent',axes5);
xlim([0 60]);
set(plot1(1),'DisplayName','Red Channel','LineWidth',3,'Color',[1 0 0]);
set(plot1(2),'DisplayName','Green Channel','LineStyle','--','Color',[0 1 0]);
set(plot1(3),'DisplayName','Blue Channel',...
    'MarkerIndices',1:50:3000,...
    'MarkerSize',8,...
    'Marker','o',...
    'LineStyle',':',...
    'Color',[0 0 1]);
ylabel('Pixel Intensity');
xlabel('Time (second)');
title('Raw Color Channel');
box(axes5,'on');
set(axes5,'FontSize',24);
set(axes5,'NextPlot','replaceall');
% Create legend
legend1 = legend(axes5,'show');
set(legend1,...
    'Position',[0.278884223030409 0.888951198203724 0.451069945744035 0.0952380972816831],...
    'Orientation','horizontal');
%% 
figure(6)
axes6 = axes('Parent',figure(6));
plot(t,iPPG_3);
xlim([0 60]);
xlabel('Time (second)');
ylabel('Amplitude');
set(axes6,'FontSize',24);
set(axes6,'NextPlot','replaceall');
%%
%% Using detrending-normalization
% figure(1)
% [pks_1,locs_1]=findpeaks(IPPG_1,t,'MinPeakProminence',0.0008);findpeaks(IPPG_1,t,'MinPeakProminence',0.0007);title('HR Green Channel');
% figure(2)
% [pks_2,locs_2]=findpeaks(IPPG_2,t,'MinPeakProminence',0.0006);findpeaks(IPPG_2,t,'MinPeakProminence',0.0007);title('HR Green-Red');
% figure(3)
% [pks_3,locs_3]=findpeaks(IPPG_3,t,'MinPeakProminence',0.0005);findpeaks(IPPG_3,t,'MinPeakProminence',0.0007);title('HR CHROM');
% figure(4)
% [pks_4,locs_4]=findpeaks(IPPG_4,t,'MinPeakProminence',0.0007);findpeaks(IPPG_4,t,'MinPeakProminence',0.0007);title('HR POS');

%% Using exponential MA
% figure(1)
% [pks_1,locs_1]=findpeaks(IPPG1,t,'MinPeakProminence',filter1);findpeaks(IPPG1,t,'MinPeakProminence',filter1);title('HR Green Channel');
% figure(2)
% [pks_2,locs_2]=findpeaks(IPPG2,t,'MinPeakProminence',filter2);findpeaks(IPPG2,t,'MinPeakProminence',filter2);title('HR Green-Red');
% figure(3)
% [pks_3,locs_3]=findpeaks(IPPG3,t,'MinPeakProminence',filter3);findpeaks(IPPG3,t,'MinPeakProminence',filter3);title('HR CHROM');
% figure(4)
% [pks_4,locs_4]=findpeaks(IPPG4,t,'MinPeakProminence',filter4);findpeaks(IPPG4,t,'MinPeakProminence',filter4);title('HR POS');
%% Calculate HR
total=numel(pks_1);
for i=1:total
    if i==2
        first=locs_1(i);
    end
    if i==total-1
        last=locs_1(i);
    end
end
length=last-first;
sample=((total-3)/length);
HR_1=sample*60;
%% Calculate HR
total=numel(pks_2);
for i=1:total
    if i==2
        first=locs_2(i);
    end
    if i==total-1
        last=locs_2(i);
    end
end
length=last-first;
sample=((total-3)/length);
HR_2=sample*60;
%% Calculate HR
total=numel(pks_3);
for i=1:total
    if i==2
        first=locs_3(i);
    end
    if i==total-1
        last=locs_3(i);
    end
end
length=last-first;
sample=((total-3)/length);
HR_3=sample*60;
%% Calculate HR
total=numel(pks_4);
for i=1:total
    if i==2
        first=locs_4(i);
    end
    if i==total-1
        last=locs_4(i);
    end
end
length=last-first;
sample=((total-3)/length);
HR_4=sample*60;

