function name = ostruct_write_table(out_dir, fig_tag, caption, F, cols, w_idx, notes)
% OSTRUCT_WRITE_TABLE The _table.md sidecar of an outlier_structures figure: a caption paragraph, notes, and per network and condition the median and range over seeds.
%
%   name = OSTRUCT_WRITE_TABLE(out_dir, fig_tag, caption, F, cols, w_idx)
%   name = OSTRUCT_WRITE_TABLE(..., notes)     notes: cellstr, one paragraph each, after the caption
%
%   F      ostruct_figure_data's output
%   cols   columns of F.T (numeric) to tabulate
%   w_idx  which networks (indices into F.keys), in the order to list them
%   name   the file name written (<fig_tag>_table.md), '' if it could not be opened
%
% Every cell is "median [min, max]" over the finished seeds (one value when
% there is one seed), with the number of seeds in its own column. A last
% section lists the structural check (every weight draw, passed or not) and
% any run the stage skipped or that failed.
%
% See also: ostruct_figure_data, ostruct_mode_measures

if nargin < 7; notes = {}; end
name = [fig_tag '_table.md'];
fid = fopen(fullfile(out_dir, name), 'w');
if fid < 0; name = ''; return; end
closer = onCleanup(@() fclose(fid));
fprintf(fid, '# %s\n\n%s\n\n', fig_tag, caption);
for k = 1:numel(notes)
    fprintf(fid, '%s\n\n', notes{k});
end
fprintf(fid, ['Rates are in 1/s, rotation rates in rad/s, structural eigenvalues (lam0_*_R) in units of the ' ...
    'structural bulk radius R of the drawn matrix; detachments (mode_rho*, right_rho*) are against the ' ...
    'theoretical radius R_eff of each state. Each cell: median [min, max] over the seeds. Definitions: ' ...
    '`src/analysis/outlier_structures/ostruct_mode_measures.m`; the networks and the seed convention: ' ...
    '`src/analysis/outlier_structures_stage_plan.m`. Data: `%s`.\n\n'], F.file);

head = [{'network', 'condition', 'seeds'}, cols];
fprintf(fid, '| %s |\n', strjoin(head, ' | '));
fprintf(fid, '|%s\n', repmat('---|', 1, numel(head)));
for w = w_idx(:)'
    for c = 1:numel(F.cond_names)
        if F.n_seed(w, c) == 0; continue; end
        cells = cell(1, numel(cols));
        for q = 1:numel(cols)
            md = F.med.(cols{q})(w, c); lo = F.lo.(cols{q})(w, c); hi = F.hi.(cols{q})(w, c);
            if F.n_seed(w, c) == 1 || isnan(md) || lo == hi
                cells{q} = sprintf('%.3g', md);
            else
                cells{q} = sprintf('%.3g [%.3g, %.3g]', md, lo, hi);
            end
        end
        fprintf(fid, '| %s | %s | %d | %s |\n', F.keys{w}, F.cond_names{c}, F.n_seed(w, c), strjoin(cells, ' | '));
    end
end

ch = F.S.checks;
% Every network of the groups listed, also one with no finished run (a
% network whose every draw failed the check is in no figure, only here).
nets = F.S.networks;
keys = union(F.keys(w_idx), {nets(ismember({nets.group}, F.group_names(unique(F.group(w_idx))))).key});
ch = ch(ismember(ch.key, keys), :);
fprintf(fid, '\n## The structural check of the weight draws\n\n');
fprintf(fid, '| network | seeds | lam0 / R | expected / R | distance | passed | reason |\n|---|---|---|---|---|---|---|\n');
for k = 1:height(ch)
    ex = 'none';
    if isfinite(ch.expect_re(k)); ex = sprintf('%.3g%+.3gi', ch.expect_re(k), ch.expect_im(k)); end
    fprintf(fid, '| %s | %d-%d | %.3g%+.3gi | %s | %.2g | %d | %s |\n', ch.key{k}, ch.seed1(k), ch.seed2(k), ...
        ch.lam0_re_R(k), ch.lam0_im_R(k), ex, ch.rel_err(k), ch.pass(k), ch.reason{k});
end
r = F.S.runs;
bad = r(~r.ok & ismember(r.key, keys), :);
if height(bad) > 0
    fprintf(fid, '\n## Runs not drawn\n\n| network | condition | seeds | why |\n|---|---|---|---|\n');
    for k = 1:height(bad)
        fprintf(fid, '| %s | %s | %d-%d | %s |\n', bad.key{k}, bad.condition{k}, bad.seed1(k), bad.seed2(k), bad.message{k});
    end
end
end
