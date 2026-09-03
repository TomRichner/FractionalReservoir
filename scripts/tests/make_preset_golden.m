function out_file = make_preset_golden()
% MAKE_PRESET_GOLDEN Freeze srnn_param_preset's outputs as a regression fixture.
%
%   out_file = MAKE_PRESET_GOLDEN()
%
% Captures what srnn_param_preset returns for every valid preset, plus the error
% identifiers for the retired and unknown names, into
% scripts/tests/fixtures/golden_preset_outputs.mat. test_preset_golden compares
% against it.
%
% WHY A FIXTURE RATHER THAN A FROZEN COPY OF THE FUNCTION. The obvious way to
% regression-test a refactor is to copy the old srnn_param_preset.m into tests/
% under a new name and diff the two at run time. That does not work here: the
% refactor DELETES srnn_adaptation_conditions, which the frozen copy would call,
% so the copy would either break or keep a dead function alive to serve its own
% obituary -- along with pairs_input_config, retired_presets and the local
% subfunctions. Freezing the OUTPUTS has no such tail. It cannot rot, and it
% places no constraint on what the refactor is allowed to delete.
%
% THIS MUST BE RUN BEFORE THE REFACTOR STARTS, from untouched code. A fixture
% captured afterwards would compare the new code against itself and prove
% nothing. That is why it is a separate script and a separate commit.
%
% VALIDITY RESTS ON PURITY. srnn_param_preset takes one char argument, reads no
% files and touches no RNG, so its outputs are a function of the name alone --
% verified by calling each preset twice and comparing. If it ever gains state,
% this technique stops being sound.
%
% REGENERATING IS A DELIBERATE ACT. There is no 'force' flag and no overwrite
% prompt: to replace the fixture, delete the .mat and run this again, then say in
% the commit message WHY the expected values moved. A golden file that is
% routinely regenerated to make a test pass is worse than no test, because it
% still reads like evidence.
%
% The .mat needs `git add -f` -- *.mat is gitignored repo-wide (.gitignore:60).
%
% See also: test_preset_golden, srnn_param_preset

setup_paths();

names = { ...
    'default', 'overconnected', ...
    'celltype_pairs_Sc0p2_noise0p025_dualStd_3cond', ...
    'celltype_pairs_Sc0p2_noise0p025_dualStd_4cond', ...
    'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond', ...
    'bursting_pairs', 'sompolinsky_pairs', 'single_neuron_stf', ...
    'single_neuron_dualStd', 'mc_pairs_dualStd'};

% Retired names, from srnn_param_preset's own retired_presets() local. Restated
% here because it is a subfunction and so not reachable from outside the file --
% and because the fixture should record what WAS retired at capture time, not
% track a list that the refactor might edit.
retired = { ...
    'celltype_pairs', ...
    'celltype_pairs_S_c_by_type', 'celltype_pairs_S_c_by_type_n500', ...
    'celltype_pairs_S_c_by_type_n500_fixedF', 'celltype_pairs_all_std_n500', ...
    'celltype_pairs_uniform_std_n500', 'celltype_pairs_uniform_std_n500_mu5p5', ...
    'celltype_pairs_uniform_std_n500_mu5p5_nodrive', ...
    'celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5', ...
    'celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5_noise0p02', ...
    'celltype_pairs_uniform_std_n500_mu5p5_nodrive_sig1p5_noise0p01', ...
    'celltype_pairs_Sc0p2_noise0p025', ...
    'celltype_pairs_Sc0p2_noise0p025_dualStd'};

%% Purity check -- the precondition for the whole technique
fprintf('Checking srnn_param_preset is deterministic...\n');
for k = 1:numel(names)
    [d1, m1, c1] = srnn_param_preset(names{k});
    [d2, m2, c2] = srnn_param_preset(names{k});
    if ~(isequaln(d1, d2) && isequaln(m1, m2) && isequaln(c1, c2))
        error('make_preset_golden:NotDeterministic', ...
            ['Preset ''%s'' returned different values on two consecutive calls. ' ...
             'A golden fixture is only meaningful for a pure function.'], names{k});
    end
end
fprintf('  all %d presets deterministic\n', numel(names));

%% Capture
G = struct();
G.presets = struct('name', {}, 'd', {}, 'model_class', {}, 'conditions', {});
n_cond_total = 0;
for k = 1:numel(names)
    [d, model_class, conditions] = srnn_param_preset(names{k});
    G.presets(end+1) = struct('name', names{k}, 'd', d, ...
        'model_class', model_class, 'conditions', {conditions});
    n_cond_total = n_cond_total + numel(conditions);
    fprintf('  %-48s %-20s %d cond\n', names{k}, model_class, numel(conditions));
end

% Error behaviour is behaviour. A refactor that turned a retired-preset error
% into a silent default would otherwise pass every value comparison above.
G.retired = struct('name', {}, 'identifier', {});
for k = 1:numel(retired)
    G.retired(end+1) = struct('name', retired{k}, ...
        'identifier', error_id(@() srnn_param_preset(retired{k})));
end

G.unknown = struct( ...
    'name', 'definitely_not_a_preset', ...
    'identifier', error_id(@() srnn_param_preset('definitely_not_a_preset')), ...
    'message', error_msg(@() srnn_param_preset('definitely_not_a_preset')));

G.captured = struct('when', char(datetime('now')), 'n_presets', numel(names), ...
    'n_conditions', n_cond_total);

%% Save
fixture_dir = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
if ~isfolder(fixture_dir); mkdir(fixture_dir); end
out_file = fullfile(fixture_dir, 'golden_preset_outputs.mat');
if isfile(out_file)
    error('make_preset_golden:FixtureExists', ...
        ['A fixture already exists:\n  %s\n\n' ...
         'Regenerating is deliberate: delete it and run this again, then say in ' ...
         'the commit message why the expected values moved.'], out_file);
end
save(out_file, 'G');

fprintf('\nWrote %s\n', out_file);
fprintf('  %d presets, %d conditions, %d retired names\n', ...
    numel(names), n_cond_total, numel(retired));
fprintf('  *.mat is gitignored -- commit with:  git add -f %s\n', ...
    strrep(out_file, filesep, '/'));
end

%% ------------------------------------------------------------------------
function id = error_id(fn)
try
    fn();
    id = '';    % did not throw; recorded as such so a lost error is visible
catch ME
    id = ME.identifier;
end
end

function msg = error_msg(fn)
try
    fn();
    msg = '';
catch ME
    msg = ME.message;
end
end
