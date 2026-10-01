function summary = summarize_ied_pilot(ids,best_id)
% Saved-data summary only: no trajectories or exponents are rerun.
arguments
    ids (1,:) double = 1:15
    best_id (1,1) double = 13
end
setup_paths(); root=fileparts(which('setup_paths'));
names={'Run','K','Events','RatePerMin','Recruitment','RecruitmentWidth', ...
    'WaveformWidth','WaveformAmplitude','MeanRate','Saturation','Lambda1', ...
    'WithinCorr','BetweenCorr','EventExpansion','QuietExpansion','SingleGroupFraction'};
values=NaN(numel(ids),numel(names));
for j=1:numel(ids)
    id=ids(j); D=load(fullfile(root,'data/IED_model/20260930',sprintf('run%02d',id),'run.mat'));
    [events,spatial]=ied_event_details(D.data,D.det); m=D.det.metrics;
    values(j,:)=[id D.cfg.K m.event_count m.events_per_min m.median_recruitment ...
        m.median_width median(events.WaveformWidth,'omitnan') ...
        median(events.WaveformAmplitude,'omitnan') m.mean_rate m.saturation_fraction ...
        m.lambda1 spatial.within_corr spatial.between_corr m.event_mean_expansion ...
        m.quiet_mean_expansion mean(events.DetectedGroups==1)];
    if id==best_id
        writetable(events,fullfile(root,'docs/IED_model',sprintf('events_run%02d.csv',id)));
        best=D;
    end
    clear D
end
summary=array2table(values,'VariableNames',names);
writetable(summary,fullfile(root,'docs/IED_model/summary_2026_09_30.csv'));
figdir=fullfile(root,'figs/IED_model/20260930/summary');
if ~isfolder(figdir), mkdir(figdir); end
fig=figure('Visible','off','Color','w','Position',[40 40 1400 900]);
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
quantities={'RatePerMin','Recruitment','WaveformWidth','MeanRate','Lambda1','WithinCorr'};
labels={'Candidates/min','Network Recruitment','Waveform FWHM (s)', ...
    'Mean Normalized Rate','Accumulated Growth (s^{-1})','Within/Between Correlation'};
for k=1:6
    ax=nexttile(tl); plot(ax,ids,summary.(quantities{k}),'-o','LineWidth',1.4); hold(ax,'on');
    if k==6
        plot(ax,ids,summary.BetweenCorr,'-s','LineWidth',1.4);
        legend(ax,{'Within Groups','Between Groups'},'Location','best');
    end
    if k==5, yline(ax,0,'--','Color',[0 .6 0]); end
    ylabel(ax,labels{k}); xlabel(ax,'Run'); xticks(ax,ids); box(ax,'off');
end
title(tl,'Sequential Single-Realization IED-like Pilot');
exportgraphics(fig,fullfile(figdir,'run_comparison.png'),'Resolution',140); close(fig);
plot_event_alignment(best,fullfile(figdir,'candidate_event_alignment.png'));
plot_candidate_detail(best,fullfile(figdir,'candidate_detail.png'));
fid=fopen(fullfile(root,'docs/IED_model/progress_2026_09_30.md'),'a');
guard=onCleanup(@()fclose(fid));
fprintf(fid,'\n## Final saved-data comparison\n\n');
fprintf(fid,'Waveform FWHM is measured on the dominant group mean x above its median baseline; it differs from detector recruitment width. Widths that do not return below half-height within ±0.75 s are censored. Within/between correlations use the same 10-Hz-subsampled x signals and fixed group labels. Random-run labels are diagnostic partitions, not anatomical modules.\n\n');
fprintf(fid,'| Run | K | Events/min | Recruitment | Recruitment width (s) | Waveform FWHM (s) | Mean rate | Leading growth (/s) | Within corr | Between corr |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for j=1:height(summary)
    fprintf(fid,'| %02d | %d | %.2f | %.3f | %.3f | %.3f | %.3f | %+.3f | %.3f | %.3f |\n', ...
        summary.Run(j),summary.K(j),summary.RatePerMin(j),summary.Recruitment(j), ...
        summary.RecruitmentWidth(j),summary.WaveformWidth(j),summary.MeanRate(j), ...
        summary.Lambda1(j),summary.WithinCorr(j),summary.BetweenCorr(j));
end
fprintf(fid,'\n![Run comparison](../../figs/IED_model/20260930/summary/run_comparison.png)\n\n');
fprintf(fid,'Candidate configuration selected for review: run%02d. [Event table](events_run%02d.csv).\n\n',best_id,best_id);
fprintf(fid,'![Candidate detail](../../figs/IED_model/20260930/summary/candidate_detail.png)\n\n');
fprintf(fid,'![Event-aligned dynamics](../../figs/IED_model/20260930/summary/candidate_event_alignment.png)\n\n');
fprintf(fid,'Alignment plots summarize all candidate peaks, not a selected successful event. Means/medians are descriptive; overlapping event windows and a single realization prevent causal or population-level inference. The top-K positive-growth sum remains a basis-dependent truncated diagnostic.\n');
disp(summary);
end

function plot_event_alignment(D,file)
[events,~]=ied_event_details(D.data,D.det); lag=-.75:.01:.75; n=height(events);
X=NaN(n,numel(lag)); L=X; Q=X;
for e=1:n
    group=events.DominantGroup(e);
    signal=mean(double(D.data.x(D.data.groups==group,:)),1)-mean(D.det.baseline(D.data.groups==group));
    X(e,:)=interp1(D.data.t,signal,events.Time(e)+lag,'linear',NaN);
    L(e,:)=interp1(D.det.growth.t,D.det.growth.local,events.Time(e)+lag,'linear',NaN);
    Q(e,:)=interp1(D.det.growth.t,D.det.growth.expansion,events.Time(e)+lag,'linear',NaN);
end
fig=figure('Visible','off','Color','w','Position',[40 40 1000 900]);
tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
series={X,L,Q}; labels={'Dominant-Group \Delta x (a.u.)','Local Leading Growth (s^{-1})', ...
    sprintf('Positive Top-%d Growth (bits/s)',D.cfg.K)};
for k=1:3
    ax=nexttile(tl); hold(ax,'on'); plot(ax,lag,series{k}','Color',[.8 .8 .8]);
    a=plot(ax,lag,mean(series{k},1,'omitnan'),'k','LineWidth',2);
    b=plot(ax,lag,median(series{k},1,'omitnan'),'Color',[.6 .15 .7],'LineWidth',1.5);
    xline(ax,0,'--','Color',[0 .6 0]); ylabel(ax,labels{k});
    if k==1, legend(ax,[a b],{'Mean','Median'},'Location','best'); end
end
xlabel(ax,'Time from Candidate Peak (s)');
title(tl,sprintf('Run %02d | %d Event Peaks | Individual Curves and Descriptive Summaries',D.cfg.id,n));
exportgraphics(fig,file,'Resolution',140); close(fig);
end

function plot_candidate_detail(D,file)
[events,~]=ied_event_details(D.data,D.det);
[~,e]=max(events.WaveformAmplitude); center=events.Time(e); group=events.DominantGroup(e);
t=D.data.t; window=abs(t-center)<=1;
fig=figure('Visible','off','Color','w','Position',[40 40 1200 1000]);
tl=tiledlayout(fig,4,1,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); hold(ax,'on');
for g=1:D.cfg.groups, plot(ax,t(window),mean(D.data.x(D.data.groups==g,window),1),'LineWidth',1.3); end
ylabel(ax,'Group Mean x (a.u.)'); legend(ax,compose('Group %d',1:D.cfg.groups),'Location','eastoutside');
ax=nexttile(tl); ix=find(D.data.groups==group);
imagesc(ax,t(window),1:numel(ix),double(D.data.x(ix,window))); axis(ax,'xy'); colorbar(ax);
ylabel(ax,sprintf('Group %d Neuron Index (n=%d)',group,numel(ix)));
ax=nexttile(tl); hold(ax,'on');
plot(ax,t(window),mean(D.det.active(:,window),1),'k','LineWidth',1.5);
plot(ax,t(window),mean(D.det.active(ix,window),1),'Color',[.6 .2 .7],'LineWidth',1.5);
ylabel(ax,'Recruitment Fraction'); legend(ax,{'Whole Network','Event Group'},'Location','eastoutside'); ylim(ax,[0 1]);
ax=nexttile(tl); hold(ax,'on'); time=D.det.growth.t; take=abs(time-center)<=1;
yyaxis(ax,'left'); a=plot(ax,time(take),D.det.growth.local(take),'Color',[.5 .5 .5]);
ylabel(ax,'Local Leading Growth (s^{-1})'); yline(ax,0,'--','Color',[0 .6 0],'HandleVisibility','off');
yyaxis(ax,'right'); b=plot(ax,time(take),D.det.growth.expansion(take),'Color',[.6 .2 .7],'LineWidth',1.4);
ylabel(ax,sprintf('Positive Top-%d Growth (bits/s)',D.cfg.K)); legend(ax,[a b],{'Local Growth','Positive-Growth Sum'},'Location','eastoutside');
xlabel(ax,'Time (s)'); title(tl,sprintf('Run %02d | Representative Local Event at %.3f s | Waveform FWHM %.3f s',D.cfg.id,center,events.WaveformWidth(e)));
exportgraphics(fig,file,'Resolution',150); close(fig);
end
