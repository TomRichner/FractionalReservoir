function out=fig_single_neuron_revised(cfg)
% FIG_SINGLE_NEURON_REVISED Replot saved 1TS mechanism examples, no simulation.
arguments
    cfg.run_dir char = ''
    cfg.out_dir char = ''
    cfg.save logical = true
    cfg.visible logical = false
    cfg.verbose char = 'minimal'
end
setup_paths();
f=fullfile(cfg.run_dir,'single_neuron','single_neuron_data.mat');
assert(isfile(f),'fig_single_neuron_revised:MissingStage','Run the paper_illustrations analysis stage first.');
D=load(f); assert(all([D.results.n]==1));
fig=figure('Visible','off','Color','w','Position',[40 40 1000 850]);
tl=tiledlayout(fig,6,3,'TileSpacing','compact','Padding','compact');
fields={'u','x','r','syn','sfa','std'};
labels={'input u','x','raw rate r','synaptic output','SFA feedback','STD product'};
for c=1:3
    r=D.results(c);
    for j=1:6
        ax=nexttile(tl,(j-1)*3+c); plot(ax,r.t,r.(fields{j})','k','LineWidth',1.4);
        xlim(ax,D.settings.display_window); box(ax,'off');
        if c==1, ylabel(ax,labels{j}); end
        if j==1, title(ax,r.title,'FontWeight','normal'); end
        if j==6, xlabel(ax,'time (s)'); else, ax.XTickLabel=[]; end
    end
end
if cfg.visible, fig.Visible='on'; end
out=struct('figs',fig,'files',{{}},'source',f);
if cfg.save
    folder=default_out_dir(cfg.out_dir,mfilename('fullpath'));
    save_figure_stable(folder,'Fig_single_neuron_revised',fig);
    out.files=existing_outputs(folder,'Fig_single_neuron_revised');
end
end
