function ow_plan_compare(tag)
% OW_PLAN_COMPARE The outlier_spectra plan in every run mode, saved under a tag, and compared run by run with the 'before' one.
%
%   ow_plan_compare('before')   save the plans (mu7chain_config of each mode)
%   ow_plan_compare('after')    save them again, then compare with 'before':
%                               per mode, every run with every per-run
%                               setting (network, condition, seed, window_s,
%                               dt_sample, full, file), in order, and the
%                               plan's other settings; print what differs
%
% Written for the change of 2026-11-23 (FD-057): only the full networks
% named in outlier_spectra_long_networks keep the long window. In essential
% and fast the long window equals the short one, so their plans must be
% identical before and after.

setup_paths();
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
dd = fullfile(root, 'data', 'explorations', 'outlier_window');
modes = run_mode_names();
P = struct();
for m = 1:numel(modes)
    P.(modes{m}) = outlier_spectra_stage_plan(mu7chain_config(modes{m}));
end
save(fullfile(dd, ['plans_' tag '.mat']), 'P');
if ~strcmp(tag, 'after'); fprintf('saved the plans (%s)\n', tag); return; end
B = load(fullfile(dd, 'plans_before.mat')); B = B.P;
for m = 1:numel(modes)
    a = B.(modes{m}).runs; b = P.(modes{m}).runs;
    same_runs = isequal(rmfield_if(a), rmfield_if(b));
    fa = setdiff(fieldnames(B.(modes{m})), {'runs'}); fb = setdiff(fieldnames(P.(modes{m})), {'runs'});
    other = setdiff(fb, fa);
    changed = {};
    for f = intersect(fa, fb)'
        if ~isequal(B.(modes{m}).(f{1}), P.(modes{m}).(f{1})); changed{end + 1} = f{1}; end %#ok<AGROW>
    end
    fprintf('%-10s runs %3d -> %3d, per-run settings identical: %d; plan fields changed: %s; added: %s\n', modes{m}, ...
        numel(a), numel(b), same_runs, strjoin_or(changed), strjoin_or(other));
    if ~same_runs && numel(a) == numel(b)
        for k = 1:numel(a)
            if ~isequal(a(k), b(k))
                fprintf('    %s %s [%d %d]: window %g -> %g s, dt_sample %g -> %g s, full %d -> %d\n', b(k).network, ...
                    b(k).condition, b(k).seed, a(k).window_s, b(k).window_s, a(k).dt_sample, b(k).dt_sample, a(k).full, b(k).full);
                if k > 20; fprintf('    ... (first differences only)\n'); break; end
            end
        end
    end
end
end

function r = rmfield_if(r)
% The per-run settings as they are (a hook if a field is ever added).
end

function s = strjoin_or(c)
if isempty(c); s = 'none'; else; s = strjoin(c, ', '); end
end
