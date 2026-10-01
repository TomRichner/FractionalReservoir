function [events,spatial] = ied_event_details(data,det)
% Additional waveform measures; leaves detector-v1 decisions unchanged.
t=data.t(:)'; fs=1/median(diff(t));
x=movmean(double(data.x),max(1,round(.025*fs)),2); groups=data.groups;
thin=x(:,1:max(1,round(fs/10)):end)'; valid=std(thin,0,1)>1e-8;
C=corrcoef(thin(:,valid)); labels=groups(valid); pairs=~eye(size(C));
same=labels==labels';
spatial=struct('within_corr',mean(C(pairs & same)), ...
    'between_corr',mean(C(pairs & ~same)));
names={'Time','DominantGroup','DetectedGroups','Recruitment','RecruitmentWidth', ...
    'WaveformWidth','WaveformAmplitude','CensoredWidth','PreLocalGrowth','PreExpansion'};
values=NaN(numel(det.events),numel(names));
for e=1:numel(det.events)
    tm=det.events(e).time; window=abs(t-tm)<=.10;
    active=any(det.active(:,window),2); counts=accumarray(groups,active);
    eligible=det.events(e).groups;
    [~,winner]=max(counts(eligible)); dominant=eligible(winner);
    trace=mean(x(groups==dominant,:),1);
    baseline=mean(det.baseline(groups==dominant));
    peak_window=find(abs(t-tm)<=.10); [height,k]=max(trace(peak_window));
    peak=peak_window(k); amplitude=height-baseline;
    left=peak; right=peak; half=baseline+amplitude/2;
    lower=find(t>=tm-.75,1); upper=find(t<=tm+.75,1,'last');
    while left>lower && trace(left)>half, left=left-1; end
    while right<upper && trace(right)>half, right=right+1; end
    censored=left==lower || right==upper || amplitude<=0;
    width=t(right)-t(left); if censored, width=NaN; end
    pre=det.growth.t>=tm-.50 & det.growth.t<tm-.10;
    pre_local=NaN; pre_q=NaN;
    if any(pre)
        pre_local=mean(det.growth.local(pre)); pre_q=mean(det.growth.expansion(pre));
    end
    values(e,:)=[tm dominant numel(det.events(e).groups) det.events(e).recruitment ...
        det.events(e).width width amplitude censored pre_local pre_q];
end
events=array2table(values,'VariableNames',names);
end
