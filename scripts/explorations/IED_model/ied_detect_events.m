function out = ied_detect_events(t,x,r,groups,lya)
% Candidate population events, not clinically validated IED classifications.
% Fixed detector v1. Thresholds are held constant across this pilot.
t=t(:)'; fs=1/median(diff(t));
xs=movmean(double(x),max(1,round(.025*fs)),2);
baseline=median(xs,2);
robust_sd=1.4826*median(abs(xs-baseline),2);
threshold=baseline+max(3*robust_sd,.12);
active=xs>threshold;
g=unique(groups(:))'; recruitment=zeros(numel(g),numel(t));
local_t=[]; local_group=[]; local_width=[];
long_count=0;
for j=1:numel(g)
    recruitment(j,:)=mean(active(groups==g(j),:),1);
    if max(recruitment(j,:))<=.20, continue; end
    [~,times,width]=findpeaks(recruitment(j,:),t, ...
        'MinPeakHeight',.20,'MinPeakProminence',.10,'MinPeakDistance',.15);
    long_count=long_count+nnz(width>1);
    keep=width>=.025 & width<=1;
    % Recruitment can have a flat maximum; findpeaks returns its first sample.
    % Locate the dendritic peak within that excursion instead of its onset.
    group_x=mean(xs(groups==g(j),:),1);
    for peak=1:numel(times)
        indices=find(abs(t-times(peak))<=max(width(peak),.025));
        [~,maximum]=max(group_x(indices));
        times(peak)=t(indices(maximum));
    end
    local_t=[local_t times(keep)]; %#ok<AGROW>
    local_group=[local_group repmat(g(j),1,nnz(keep))]; %#ok<AGROW>
    local_width=[local_width width(keep)]; %#ok<AGROW>
end
[local_t,order]=sort(local_t); local_group=local_group(order); local_width=local_width(order);
events=struct('time',{},'width',{},'recruitment',{},'groups',{},'peak_x',{});
j=1;
while j<=numel(local_t)
    k=j;
    while k<numel(local_t) && local_t(k+1)-local_t(j)<=.12, k=k+1; end
    tm=median(local_t(j:k)); window=abs(t-tm)<=.10;
    events(end+1)=struct('time',tm,'width',median(local_width(j:k)), ...
        'recruitment',mean(any(active(:,window),2)), ...
        'groups',unique(local_group(j:k)), ...
        'peak_x',max(mean(xs(:,window)-baseline,1))); %#ok<AGROW>
    j=k+1;
end
metrics=struct('event_count',numel(events),'events_per_min',numel(events)*60/(t(end)-t(1)), ...
    'local_event_count',numel(local_t),'long_group_excursions',long_count, ...
    'median_width',NaN,'median_recruitment',NaN,'global_fraction',NaN, ...
    'iei_cv',NaN,'mean_rate',mean(r,'all'),'saturation_fraction',mean(r>.95,'all'), ...
    'mean_active_fraction',mean(active,'all'),'mean_corr',NaN, ...
    'lambda1',NaN,'local_positive_fraction',NaN,'event_local_positive_fraction',NaN, ...
    'event_mean_local',NaN,'quiet_mean_local',NaN,'event_mean_expansion',NaN, ...
    'quiet_mean_expansion',NaN);
if ~isempty(events)
    metrics.median_width=median([events.width]);
    metrics.median_recruitment=median([events.recruitment]);
    metrics.global_fraction=mean([events.recruitment]>.60);
    if numel(events)>2
        intervals=diff([events.time]); metrics.iei_cv=std(intervals)/mean(intervals);
    end
end
thin=xs(:,1:max(1,round(fs/10)):end)';
valid=std(thin,0,1)>1e-8; C=corrcoef(thin(:,valid));
if size(C,1)>1, metrics.mean_corr=mean(C(~eye(size(C)))); end
growth=struct('t',[],'local',[],'accumulated',[],'expansion',[]);
if ~isempty(lya) && isfield(lya,'local_LE_spectrum_t')
    growth.t=lya.t_lya(:)+.025;
    growth.local=lya.local_LE_spectrum_t(:,1);
    growth.accumulated=lya.finite_LE_spectrum_t(:,1);
    growth.expansion=sum(max(lya.local_LE_spectrum_t,0),2)/log(2);
    keep=growth.t>=t(1) & growth.t<=t(end);
    event_mask=false(size(growth.t));
    for e=1:numel(events)
        event_mask=event_mask | (growth.t>=events(e).time-.10 & growth.t<=events(e).time+.20);
    end
    event_mask=event_mask & keep; quiet=keep & ~event_mask;
    metrics.lambda1=lya.LE_spectrum(1);
    metrics.local_positive_fraction=mean(growth.local(keep)>0);
    if any(event_mask)
        metrics.event_local_positive_fraction=mean(growth.local(event_mask)>0);
        metrics.event_mean_local=mean(growth.local(event_mask));
        metrics.event_mean_expansion=mean(growth.expansion(event_mask));
    end
    metrics.quiet_mean_local=mean(growth.local(quiet));
    metrics.quiet_mean_expansion=mean(growth.expansion(quiet));
end
out=struct('version','candidate_detector_v1','threshold',threshold, ...
    'baseline',baseline,'robust_sd',robust_sd,'active',active, ...
    'group_recruitment',recruitment,'events',events,'metrics',metrics,'growth',growth);
end
