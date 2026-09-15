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
        raster('fig_introductory_concepts/Fig_Intro_Concepts.png',[0 .44 1 .54]);
        a = fig_FI_curve('save',false,'visible',false);
        native(a.figs(1),[.045 .05 .48 .365]); close(a.figs);
        pan = panel([.55 .075 .43 .31]); ax=axes(pan); hold(ax,'on');
        th=linspace(0,2*pi,300);
        plot(ax,-.25+cos(th),sin(th),'Color',[.5 .5 .5],'LineWidth',2);
        plot(ax,-.25+.58*cos(th),.58*sin(th),'--','Color',[0 .45 .7],'LineWidth',2);
        plot(ax,-1.25+cos(th),sin(th),'-.','Color',[.85 .33 .1],'LineWidth',2);
        xline(ax,0,':'); yline(ax,0,':'); axis(ax,'equal'); box(ax,'off');
        xlabel(ax,'Re(\lambda), schematic'); ylabel(ax,'Im(\lambda), schematic');
        legend(ax,{'reference','STD: shrink','SFA: left shift'},'Location','southoutside','NumColumns',3,'Box','off');
        title(ax,{'Effective-connectivity intuition','Conceptual approximation'},'FontWeight','normal','FontSize',11);
        notes{end+1}='Disk schematic is conceptual, not an exact transformation of the full adaptive Jacobian. Intro is the archived Sompolinsky illustration; F-I is analytic.';
    case 2
        f=fullfile(cfg.run_dir,'representative_dynamics','representative_dynamics_data.mat');
        if isfile(f)
            D=load(f); draw_dynamics(fig,D); sources{end+1}=f;
        else
            raster('fig_example_timeseries/fig_example_timeseries.png',[0 .035 1 .91]);
            banner('Archived baseline dynamics; midpoint-step simulation pending');
            notes{end+1}='Current archive has baseline input, not a midpoint step. Revised stage has not been run.';
        end
    case 3
        a=fig_eig_heatmap('run_dir',cfg.run_dir,'save',false,'visible',false,'density_scale','loglog');
        native(a.figs(1),[0 .66 1 .32]); sources{end+1}=a.source; close(a.figs);
        f=fullfile(cfg.run_dir,'local_lyapunov','local_lyapunov_data.mat'); D=load(f); sources{end+1}=f;
        pan=panel([0 .045 1 .59]); tl=tiledlayout(pan,3,3,'TileSpacing','compact','Padding','compact');
        st=manuscript_style(); T=D.settings.T;
        for i=1:numel(D.results)
            r=D.results(i); tk=r.topk; t=tk.t_lya(:); col=st.condition_color(r.name);
            ax=nexttile(tl,i); plot(ax,r.t_ex,r.u_ex(1:min(8,size(r.u_ex,1)),:)','LineWidth',.6); xlim(ax,[0 T]); title(ax,r.title,'FontWeight','normal');
            if i==1, ylabel(ax,'input (8 fixed neurons)'); end
            ax=nexttile(tl,3+i); plot(ax,t,tk.finite_LE_spectrum_t(:,1),'Color',col,'LineWidth',1.3); yline(ax,0,':'); xlim(ax,[T/2 T]);
            if i==1, ylabel(ax,'accumulated \lambda_1 (s^{-1})'); end
            ax=nexttile(tl,6+i); h=sum(max(tk.local_LE_spectrum_t,0),2)/log(2);
            plot(ax,t,h,'Color',col,'LineWidth',1.1); xlim(ax,[0 T]); xline(ax,T/2,':'); xlabel(ax,'time (s)');
            if i==1, ylabel(ax,'local positive sum (bit/s)'); end
            title(ax,sprintf('%.2f%% > 0 over [%g,%g] s; K=%d',100*mean(h(t>=T/2 & t<=T)>0),T/2,T,size(tk.local_LE_spectrum_t,2)),'FontSize',10,'FontWeight','normal');
        end
        notes{end+1}='Local positive sum = sum(max(saved local QR rates,0),2)/log(2). Percentage uses identical positive-sum samples over [T/2,T], strict >0, unweighted equal QR segments. It is provisional local h_KS, not an established invariant KS entropy; broader review deferred. Eigenvalue and switching stages are separate examples.';
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
copy_axes(fig,flat,[0 0 1 1]);
close(fig); fig=flat;
set(findall(fig,'Type','axes'),'FontSize',10);
tag=sprintf('Fig_Main%d_Grouped',group);
out=struct('figs',fig,'files',{{}},'source',{sources});
if cfg.save
    if ~isfolder(outdir), mkdir(outdir); end
    % PNG plus editable MATLAB figure: local-rate traces make SVG needlessly large.
    drawnow;
    exportgraphics(fig,fullfile(outdir,[tag '.png']),'Resolution',200);
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
function draw_dynamics(fig,D)
% Saved display arrays only, generated by run_paper_illustrations.
tl=tiledlayout(fig,6,3,'TileSpacing','compact','Padding','compact');
rows={'u','x','r','syn','sfa','std'}; labels={'input u','x','raw rate r','synaptic output','SFA feedback','STD product'};
for c=1:numel(D.results)
    r=D.results(c);
    for j=1:6
        ax=nexttile(tl,(j-1)*3+c); Y=r.(rows{j});
        plot(ax,r.t,Y','LineWidth',.8); xlim(ax,D.settings.display_window); box(ax,'off');
        if c==1, ylabel(ax,labels{j}); end
        if j==1, title(ax,r.title,'FontWeight','normal'); end
        if j==6, xlabel(ax,'time (s)'); else, ax.XTickLabel=[]; end
    end
end
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
