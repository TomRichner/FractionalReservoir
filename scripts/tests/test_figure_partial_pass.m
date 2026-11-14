% TEST_FIGURE_PARTIAL_PASS make_all_paper_figures with cfg.only redraws named entries and keeps the record of the rest.
%
% Quick, no simulation: a registry of three toy entries that each write one
% small PNG. Checks that
%   1. a full pass in which one entry fails records it (2 of 3), and
%      read_figure_manifest reads back what make_all_paper_figures returned;
%   2. a partial pass naming the failed entry, once it can succeed, draws
%      that entry alone: the other entries' files are not rewritten, the
%      returned results and the manifest cover the whole set (3 of 3) with
%      the other rows as they were, the manifest has a "Partial passes"
%      section naming the entry, the root's git_provenance.txt is untouched
%      and the pass's own is under redraws/;
%   3. a second partial pass carries the first one's note forward;
%   4. an unknown entry name and a figure root with no manifest are refused,
%      and the refused root is not created.
% Needs a pool allocation in the session, as make_all_paper_figures does
% (pool_set(0) is one). Asserts; prints one PASS line.
%
% See also: make_all_paper_figures, read_figure_manifest, write_figure_manifest

tmp = tempname; mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's'));
run_dir = fullfile(tmp, 'run'); mkdir(run_dir);
fig_root = fullfile(tmp, 'figs');
flag = fullfile(tmp, 'b_can_draw.txt');             % entry b fails until this exists

entry = @(name, in_paper, args) struct('name', name, 'fn', @toy_figure, 'in_paper', in_paper, 'args', {args});
cfg = struct('preset_name', 'toy', 'run_mode', 'fast', 'run_dir', run_dir, 'fig_root', fig_root, ...
    'verbose', 'near-none', 'figures', {{entry('a', true, {'tag', 'A'}), ...
    entry('b', true, {'tag', 'B', 'needs', flag}), entry('c', false, {'tag', 'C'})}});

%% 1. Full pass, entry b fails
[~, r1] = evalc('make_all_paper_figures(cfg)');
assert(isequal({r1.name}, {'a', 'b', 'c'}) && isequal([r1.ok], [true false true]));
M1 = read_figure_manifest(fig_root);
assert(isequal({M1.results.name}, {r1.name}) && isequal([M1.results.ok], [r1.ok]) && ...
    isequal([M1.results.in_paper], [r1.in_paper]) && isequal([M1.results.n_files], [r1.n_files]) && ...
    isequal({M1.results.files}, {r1.files}) && isempty(M1.notes) && ~isempty(M1.generated), ...
    'the manifest must read back as the pass returned it');
assert(contains(M1.results(2).err, 'toy_figure:Needs'));
md1 = fileread(fullfile(fig_root, 'manifest.md'));
assert(contains(md1, '2/3 succeeded') && ~contains(md1, 'Partial passes'));
prov_root = fileread(fullfile(fig_root, 'git_provenance.txt'));
stamp_of = @(name) dir(fullfile(fig_root, name, '*.png')).datenum;
t_a = stamp_of('a'); t_c = stamp_of('c');

%% 2. Partial pass: b alone
fclose(fopen(flag, 'w'));
pause(1.1);                                          % file times have one-second resolution on some volumes
cfg_b = cfg; cfg_b.only = {'b'};
[~, r2] = evalc('make_all_paper_figures(cfg_b)');
assert(isequal({r2.name}, {'a', 'b', 'c'}) && all([r2.ok]), 'the results cover the whole set');
assert(stamp_of('a') == t_a && stamp_of('c') == t_c, 'entries not named must not be redrawn');
assert(isfile(fullfile(fig_root, 'b', 'B.png')));
M2 = read_figure_manifest(fig_root);
assert(isequal(M2.results([1 3]), M1.results([1 3])), 'the other rows are carried as they were');
assert(M2.results(2).ok && M2.results(2).n_files == 1 && isempty(M2.results(2).err));
md2 = fileread(fullfile(fig_root, 'manifest.md'));
assert(contains(md2, '3/3 succeeded') && contains(md2, '## Partial passes') && ~contains(md2, 'FAILED'));
assert(numel(M2.notes) == 2 && contains(M2.notes{1}, M1.generated) && contains(M2.notes{2}, 'redrew `b`'));
assert(strcmp(fileread(fullfile(fig_root, 'git_provenance.txt')), prov_root), ...
    'the root''s git provenance describes the full pass and must stay');
redraws = dir(fullfile(fig_root, 'redraws', '*', 'git_provenance.txt'));
assert(isscalar(redraws), 'the partial pass records its own git state under redraws/');

%% 3. A second partial pass keeps the first one's note
pause(1.1);                                          % its redraws folder is named to the second
cfg_ac = cfg; cfg_ac.only = {'c', 'a'};
[~, r3] = evalc('make_all_paper_figures(cfg_ac)');
assert(isequal({r3.name}, {'a', 'b', 'c'}) && all([r3.ok]));
M3 = read_figure_manifest(fig_root);
assert(numel(M3.notes) == 3 && isequal(M3.notes(1:2), M2.notes) && contains(M3.notes{3}, 'redrew `a`, `c`'));
assert(isequal(M3.results(2), M2.results(2)) && stamp_of('a') > t_a);

%% 4. Refused before anything is drawn
bad = cfg; bad.only = {'no_such_entry'};
assert(throws(@() make_all_paper_figures(bad), 'make_all_paper_figures:UnknownEntry'));
fresh = cfg; fresh.fig_root = fullfile(tmp, 'never_drawn'); fresh.only = {'a'};
assert(throws(@() make_all_paper_figures(fresh), 'read_figure_manifest:NoManifest'));
assert(~isfolder(fresh.fig_root), 'a refused partial pass must not create its figure root');

fprintf('PASS: test_figure_partial_pass -- cfg.only draws the named entries, the manifest keeps the whole set.\n');

%% ------------------------------------------------------------------------
function out = toy_figure(varargin)
% A registry entry's function: one small PNG named by 'tag' into out_dir.
% 'needs': a file that must exist, else it errors (a figure whose stage
% produced nothing).
p = struct(varargin{:});
if isfield(p, 'needs') && ~isfile(p.needs)
    error('toy_figure:Needs', 'Nothing to draw from: %s is missing.', p.needs);
end
if ~isfolder(p.out_dir); mkdir(p.out_dir); end
f = fullfile(p.out_dir, [p.tag '.png']);
imwrite(uint8(255 * ones(4, 4, 3)), f);
out = struct('figs', [], 'files', {{[p.tag '.png']}}, 'source', '');
end

function tf = throws(fn, id)
tf = false;
try
    fn();
catch err
    tf = strcmp(err.identifier, id);
end
end
