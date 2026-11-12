function out=fig_single_neuron_revised(cfg)
% FIG_SINGLE_NEURON_REVISED Replot saved 1TS mechanism examples, no simulation.
%
%   out = FIG_SINGLE_NEURON_REVISED('run_dir', d)
%   out = FIG_SINGLE_NEURON_REVISED('run_dir', d, 'conditions', names)
%
% One column per condition of the single_neuron stage (run_paper_illustrations
% on the bundle's single-neuron preset: one unconnected neuron, one mechanism
% per condition), six rows: input, x, raw rate, synaptic output, SFA feedback,
% and the DEPRESSION FACTOR. That row is what the stage saved as `std`: the
% product of the route's b on a product route, the releasable pool b_1 on a
% chained route (never b_1 b_2; the reserve acts only through b_1). Each
% column's depression panel says which, in its corner. Titles come from the
% condition NAME (titles_by_name), not from the title the stage saved.
%
% WHICH COLUMNS ARE DRAWN. cfg.conditions names NETWORK conditions
% (no_adaptation, sfa1_std1, sfa3_stdc2, sfa3_std2), which a single-neuron run
% does not hold, so the selection is mapped onto the single-neuron conditions
% by MECHANISM, not by name or position:
%   * the single-neuron columns isolate one mechanism class each: none
%     (no_adaptation), SFA (sfaX_std0), product depression (sfa0_stdY: one or
%     more independent pools multiplied; one pool is the 1TS depression) and
%     chained depression (sfa0_stdcY);
%   * a column is drawn when its class occurs in at least one SELECTED network
%     condition (sfa1_std1 carries SFA and product depression, sfa3_stdc2
%     carries SFA and chained depression, sfa3_std2 SFA and product
%     depression); the no-mechanism column is always drawn;
%   * columns keep the stage's own order;
%   * a name either side that does not parse as sfaX_stdY / sfaX_stdcY is not
%     dropped: an unparseable single-neuron condition is always drawn, and an
%     unparseable network name contributes no class.
% So in the mu7chain bundle BOTH sets draw all four columns: the chain is the
% paper's multiple-timescale depression, and its single-neuron demonstration
% belongs with it; the product MTS adds no class the single neuron lacks. A
% selection of {no_adaptation} alone draws one column. An empty selection
% draws every column, as before.
%
% Layout: 6 rows x N columns, the figure 1000 px wide per three columns.
arguments
    cfg.run_dir char = ''
    cfg.out_dir char = ''
    cfg.save logical = true
    cfg.visible logical = false
    cfg.verbose char = 'minimal'
    cfg.conditions = {}   % NETWORK condition names; mapped by mechanism (see header)
end
setup_paths();
f=fullfile(cfg.run_dir,'single_neuron','single_neuron_data.mat');
assert(isfile(f),'fig_single_neuron_revised:MissingStage','Run the paper_illustrations analysis stage first.');
D=load(f); assert(all([D.results.n]==1));
keep=single_neuron_columns({D.results.name},cfg.conditions);
R=D.results(keep); N=numel(R);
titles=titles_by_name({R.name},{R.title});
fig=figure('Visible','off','Color','w','Position',[40 40 round(1000*N/3) 850]);
tl=tiledlayout(fig,6,N,'TileSpacing','compact','Padding','compact');
fields={'u','x','r','syn','sfa','std'};
labels={'input u','x','raw rate r','synaptic output','SFA feedback','depression factor'};
for c=1:N
    r=R(c);
    for j=1:6
        ax=nexttile(tl,(j-1)*N+c); plot(ax,r.t,r.(fields{j})','k','LineWidth',1.4);
        xlim(ax,D.settings.display_window); box(ax,'off');
        if c==1, ylabel(ax,labels{j}); end
        if j==1, title(ax,titles{c},'FontWeight','normal'); end
        if j==6
            xlabel(ax,'time (s)');
            note=factor_note(r);
            if ~isempty(note)
                % top centre: the factor falls during the step, so the
                % middle of the panel's top edge stays clear of the trace
                text(ax,0.5,0.95,note,'Units','normalized','HorizontalAlignment','center', ...
                    'VerticalAlignment','top','FontSize',9,'Color',[0.3 0.3 0.3]);
            end
        else
            ax.XTickLabel=[];
        end
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

%% ------------------------------------------------------------------------
function keep=single_neuron_columns(sn_names,requested)
% Logical mask over the single-neuron conditions (see the header's rule).
if isempty(requested) || ((ischar(requested) || isstring(requested)) && all(strlength(requested)==0))
    keep=true(1,numel(sn_names));
    return
end
requested=cellstr(requested);
present=struct('sfa',false,'product',false,'chain',false);
for k=1:numel(requested)
    m=mechanisms(requested{k});
    if isempty(m), continue; end
    present.sfa=present.sfa || m.sfa;
    present.product=present.product || m.product;
    present.chain=present.chain || m.chain;
end
keep=false(1,numel(sn_names));
for k=1:numel(sn_names)
    m=mechanisms(sn_names{k});
    if isempty(m) || ~(m.sfa || m.product || m.chain)
        keep(k)=true;   % unparseable, or no mechanism (no_adaptation): always drawn
    else
        keep(k)=(m.sfa && present.sfa) || (m.product && present.product) || (m.chain && present.chain);
    end
end
end

function m=mechanisms(name)
% Mechanism classes a condition name states; [] when it does not parse.
m=[];
if strcmp(name,'no_adaptation')
    m=struct('sfa',false,'product',false,'chain',false);
    return
end
tok=regexp(name,'^sfa(\d+)_std(c?)(\d+)$','tokens','once');
if isempty(tok), return; end
n_sfa=str2double(tok{1}); is_chain=strcmp(tok{2},'c'); n_std=str2double(tok{3});
m=struct('sfa',n_sfa>0,'product',~is_chain && n_std>0,'chain',is_chain && n_std>0);
end

function s=factor_note(r)
% Which quantity the depression row shows, for this column.
if isfield(r,'std_model')
    model=r.std_model;
elseif ~isempty(regexp(r.name,'_stdc\d','once'))
    model='chain';        % data saved before std_model existed: read the name
elseif ~isempty(regexp(r.name,'_std[1-9]','once'))
    model='product';
else
    model='none';
end
switch model
    case 'chain',   s='b_1 (releasable pool)';
    case 'product'
        if ~isempty(regexp(r.name,'_std1$','once')), s='b (one pool)'; else, s='\Pi_m b_m'; end
    otherwise,      s='';
end
end
