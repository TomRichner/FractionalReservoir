function out = fig_grouped_main(group, cfg)
% FIG_GROUPED_MAIN Compose the agreed main groups without running simulations.
% Groups 3--6 use saved numerical data. Legacy illustrative panels in groups
% 1, 2 and 8 use explicitly supplied archived PNGs; never silently simulate.
% Group 2 prefers the revised saved representative_dynamics stage.
arguments
    group (1,1) double
    cfg.run_dir (1,:) char = ''
    cfg.preset_name (1,:) char = ''
    cfg.source_fig_root (1,:) char = ''
    cfg.intro_native_root (1,:) char = ''
    cfg.pytorch_file (1,:) char = ''
    cfg.human_psd_file (1,:) char = ''
    cfg.out_dir (1,:) char = ''
    cfg.save (1,1) logical = true
    cfg.visible (1,1) logical = false
    cfg.verbose (1,:) char = 'minimal'
end
setup_paths();
root = fileparts(which('setup_paths'));
if ~isempty(cfg.source_fig_root) && ~is_absolute(cfg.source_fig_root)
    cfg.source_fig_root = fullfile(root,cfg.source_fig_root);
end
outdir = default_out_dir(cfg.out_dir,mfilename('fullpath'));
fig = figure('Visible',onoff(cfg.visible),'Color','w','Position',[40 40 1300 950]);
sources = {}; notes = {};
switch group
    case 1
        intro_root=cfg.intro_native_root;
        if isempty(intro_root), intro_root=fullfile(cfg.source_fig_root,'fig_introductory_concepts'); end
        a=fig_main1_concepts(cfg.preset_name,intro_root);
        native(a.figs,[0 0 1 1]); close(a.figs);
        sources=a.source; notes=a.notes;
    case 2
        a=fig_main2_dynamics(cfg.run_dir,cfg.source_fig_root);
        native(a.figs,[0 0 1 1]); close(a.figs);
        sources=a.source; notes=a.notes;
    case 3
        a=fig_main3_stability(cfg.run_dir);
        native(a.figs,[0 0 1 1]); close(a.figs);
        sources=a.source; notes=a.notes;
    case 4
        a=fig_sensitivity_medians('run_dir',cfg.run_dir,'preset_name',cfg.preset_name,'save',false,'visible',false);
        native(a.figs(1),[0 .34 1 .64]); close(a.figs);
        a=fig_lle_vs_rate('run_dir',cfg.run_dir,'save',false,'visible',false);
        native(a.figs(1),[0 .025 1 .29]); close(a.figs);
        sources{end+1}=cfg.run_dir;
        notes{end+1}='Medians and IQR are unchanged. Fixed sensitivity LLE display limits clip more negative/positive values; full distributions remain supplementary.';
    case 5
        a=fig_transient_amplification('run_dir',cfg.run_dir,'save',false,'visible',false);
        sources{end+1}=a.source;
        axm=findall(a.figs(2),'Type','axes');
        if isscalar(axm), axm.XTickLabel={'No adaptation','Single-timescale','Multiple-timescale'}; end
        native(a.figs(2),[0 .52 1 .46]); close(a.figs);
        f=fullfile(cfg.run_dir,'transient_gain','transient_gain_data.mat'); D=load(f); sources{end+1}=f;
        pan=panel([0 .06 1 .43]); tl=tiledlayout(pan,1,3,'TileSpacing','compact','Padding','compact');
        keys={'G_worst','G_noise','G_ei_diff','G_ei_sum','G_lyap'};
        labels={'Worst case','Noise RMS','Excite E / inhibit I','Excite E / excite I','Lyapunov direction'};
        colors=[0 0 0;0 .45 .7;.85 .33 .1;0 .6 .5;.55 .2 .65]; styles={'-','--','-.',':','--'}; marks={'none','none','o','s','^'};
        ia=find(strcmp(D.settings.variants,'active')); assert(isscalar(ia));
        gain_axes=gobjects(1,numel(D.results));
        for i=1:numel(D.results)
            ax=nexttile(tl); gain_axes(i)=ax; hold(ax,'on'); r=D.results(i); smp=r.samples(strcmp({r.samples.kind},'regular')); t=smp(1).t;
            hh=gobjects(1,5);
            for k=1:5
                Y=cell2mat(arrayfun(@(s)s.(keys{k})(ia,:),smp(:),'UniformOutput',false));
                hh(k)=plot(ax,t,median(Y,1),'Color',colors(k,:),'LineStyle',styles{k},'LineWidth',1.8,'Marker',marks{k},'MarkerIndices',unique(round(linspace(2,numel(t),8))),'MarkerSize',4);
            end
            set(ax,'YScale','log'); yline(ax,1,':','HandleVisibility','off'); title(ax,r.title,'FontWeight','normal'); xlabel(ax,'time after perturbation (s)'); box(ax,'off');
            if i==1, ylabel(ax,'active dendritic gain'); end
            if i==2, lg=legend(ax,hh,labels,'NumColumns',3,'Box','off','FontSize',10); lg.Layout.Tile='south'; end
        end
        linkaxes(gain_axes,'y');
        notes{end+1}=sprintf('Active propagator only, all five directions, regular-state medians. Actual saved horizon %g s; no extrapolation. Margin and gain stages do not sample matched states.',D.settings.horizon_s);
        notes{end+1}='Four panels in one row; 14-point fonts and axes linewidth 1.0. Condition labels/titles use manuscript colors. Margin display limits [0,100] assume the unspecified upper bound; ticks are 0,25,50,75. Gain panels retain shared logarithmic limits.';
    case 6
        a=fig_sfa_EOC_allStd('run_dir',cfg.run_dir,'preset_name',cfg.preset_name,'save',false,'visible',false);
        % Keep both axes linear: setting log on an imagesc axis would distort its bins.
        for j=1:numel(a.figs)
            axs=findall(a.figs(j),'Type','axes'); set(axs,'XScale','linear');
            native(a.figs(j),[(j-1)*.5 .43 .5 .55]);
        end
        close(a.figs);
        a=fig_memory_capacity('run_dir',cfg.run_dir,'save',false,'visible',false);
        ma=findall(a.figs(1),'Type','axes');
        for jj=1:numel(ma)
            if any(contains(string(ma(jj).XTickLabel),'Adaptation'))
                ma(jj).XTickLabel={'No adaptation','1TS','MTS'}; ma(jj).XTickLabelRotation=0;
            end
        end
        native(a.figs(1),[0 .10 1 .31]); close(a.figs); sources{end+1}=cfg.run_dir;
        notes{end+1}='Slowest-timescale axes are both linear. LLE distribution display clipping and MC trial/statistics remain those of the source run.';
        notes{end+1}='14-point fonts and axes linewidth 1.0. Row-major groups: (A) timescale LLE, (B) leading-vector fractions, (C) cumulative memory capacity and reconstruction, (D) memory horizon. All five scientific axes are retained; B-D use sparse y ticks.';
    case 7
        external(cfg.pytorch_file,[.03 .10 .94 .80],'PyTorch replacement experiment pending');
        banner('Provisional external learning figure');
        notes{end+1}='The intended three-condition learning experiment is deferred. Any provided old image is provisional and does not demonstrate faster MTS learning.';
    case 8
        raster('fig_stim_engages_adaptation/bursting_psd.png',[0 .10 .51 .80]);
        external(cfg.human_psd_file,[.52 .10 .47 .80],sprintf('Human SOZ PSD\nunavailable on\nthis computer'));
        banner('Model: uniform DC input                 Human: clinical 2-Hz stimulation');
        notes{end+1}='Different systems and protocols. Human PSD is external and is never substituted from a different local model run.';
    otherwise
        error('fig_grouped_main:BadGroup','Unknown main group.');
end
% Flatten layout containers into ordinary native axes for portable export.
flat=figure('Visible','off','Color','w','Position',fig.Position);
if ismember(group,[1 2 3])
    copy_axes(fig,flat,[.025 .025 .95 .95]);
else
    copy_axes(fig,flat,[0 0 1 1]);
end
close(fig); fig=flat;
if ismember(group,[1 2 3])
    set(findall(fig,'-property','FontSize'),'FontSize',14);
    if group==2
        tt=findall(fig,'Type','axes','-regexp','Tag','dynamics_r1_c[123]');
        for k=1:numel(tt), tt(k).Title.FontSize=20; end
    elseif group==3
        tt=findall(fig,'Type','axes','-regexp','Tag','stability_r1_c[123]');
        for k=1:numel(tt), tt(k).Title.FontSize=16; end
    end
else
    set(findall(fig,'Type','axes'),'FontSize',10);
end
if group==4, style_main4_grouped(fig); end
if group==5, style_main5_grouped(fig); end
if group==6, style_main6_grouped(fig); end
tag=sprintf('Fig_Main%d_Grouped',group);
out=struct('figs',fig,'files',{{}},'source',{sources});
if cfg.save
    if ~isfolder(outdir), mkdir(outdir); end
    % PNG plus editable MATLAB figure: local-rate traces make SVG needlessly large.
    drawnow;
    if ismember(group,[1 2 3 4 5 6])
        exportgraphics(fig,fullfile(outdir,[tag '.png']),'Resolution',200,'Padding',20);
    else
        exportgraphics(fig,fullfile(outdir,[tag '.png']),'Resolution',200);
    end
    savefig(fig,fullfile(outdir,[tag '.fig']));
    out.files=existing_outputs(outdir,tag);
    fid=fopen(fullfile(outdir,[tag '_provenance.md']),'w'); guard=onCleanup(@()fclose(fid));
    fprintf(fid,'# Main group %d\n\nRun: `%s`\n\n',group,cfg.run_dir);
    for j=1:numel(sources), fprintf(fid,'- Source: `%s`\n',string(sources{j})); end
    fprintf(fid,'\n'); for j=1:numel(notes), fprintf(fid,'%s\n\n',notes{j}); end
end
    function pan=panel(pos)
        pan=uipanel(fig,'Units','normalized','Position',pos,'BorderType','none','BackgroundColor','w');
    end
    function native(src,pos)
        drawnow; copy_axes(src,fig,pos);
    end
    function raster(rel,pos)
        assert(~isempty(cfg.source_fig_root),'fig_grouped_main:NoArchive','Supply source_fig_root for archived illustrative panels.');
        f=fullfile(cfg.source_fig_root,rel); assert(isfile(f),'fig_grouped_main:MissingArchive','Missing explicitly selected archive %s',f);
        ax=axes(panel(pos),'Position',[0 0 1 1]); imshow(imread(f),'Parent',ax); sources{end+1}=f;
        notes{end+1}=['Archived raster panel: ' rel];
    end
    function external(f,pos,msg)
        if ~isempty(f) && isfile(f)
            ax=axes(panel(pos),'Position',[0 0 1 1]); imshow(imread(f),'Parent',ax); sources{end+1}=f;
        else
            ax=axes(panel(pos)); axis(ax,'off'); text(ax,.5,.5,msg,'Units','normalized','HorizontalAlignment','center','FontSize',15); notes{end+1}=msg;
        end
    end
    function banner(msg)
        ax=axes(fig,'Position',[.02 .955 .96 .04]); axis(ax,'off'); text(ax,.5,.5,msg,'Units','normalized','HorizontalAlignment','center','FontSize',13,'FontWeight','bold');
    end
end
function tf=is_absolute(p)
tf=startsWith(p,'/') || ~isempty(regexp(p,'^[A-Za-z]:','once'));
end
function v=onoff(b)
if b, v='on'; else, v='off'; end
end
function copy_axes(src,dst,pos)
% Copy data graphics, legends and colorbars together so associations survive.
drawnow;
h=findall(src,'-isa','matlab.graphics.axis.Axes','-or','Type','legend','-or','Type','colorbar');
if isempty(h), return; end
rect=zeros(numel(h),4);
for k=1:numel(h), rect(k,:)=getpixelposition(h(k),true); end
sz=getpixelposition(src); rect=rect./[sz(3:4) sz(3:4)];
new=copyobj(h,dst);
for k=1:numel(new)
    new(k).Units='normalized';
    new(k).Position=[pos(1:2)+rect(k,1:2).*pos(3:4),rect(k,3:4).*pos(3:4)];
    if isa(new(k),'matlab.graphics.axis.Axes'), new(k).PositionConstraint='innerposition'; end
end
end
