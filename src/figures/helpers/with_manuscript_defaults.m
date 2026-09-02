function cleanup = with_manuscript_defaults()
% WITH_MANUSCRIPT_DEFAULTS Set graphics-root defaults, restoring them on exit.
%
%   cleanup = WITH_MANUSCRIPT_DEFAULTS();   % keep `cleanup` alive in scope
%   ... build figures ...
%   % defaults are restored when `cleanup` goes out of scope
%
% Only for plotters that READ the root defaults instead of taking explicit
% arguments. Prefer manuscript_style() and per-object properties; this exists so
% the few places that cannot do that are still safe to batch.
%
% THE FAILURE THIS PREVENTS. plot_memory_capacity_combined and
% Fig_memory_capacity_example each did
%
%   set(0,'DefaultAxesFontSize',14); ... set(0,'DefaultTextInterpreter','none');
%
% and never put them back. Batched, every figure rendered AFTER memory capacity
% inherited them -- and DefaultTextInterpreter = 'none' turns the sensitivity
% sheets' \lambda_1 and \mu_{EE} labels into literal backslash text. Because the
% leak only manifests in a particular ORDER, it survives any amount of testing
% one figure at a time.
%
% Returns an onCleanup object. It must be assigned to a variable that stays in
% scope for as long as the defaults are wanted; dropping it (or calling this
% without capturing the output) restores them immediately, which is useless but
% not harmful.
%
% See also: manuscript_style, onCleanup

s = manuscript_style();

% Note DefaultTextInterpreter is 'tex', NOT the 'none' the old scripts set. The
% manuscript's axis labels are tex (\lambda_1, \mu_{EE}); 'none' was what broke
% them when it leaked. Legends stay 'none' because condition names contain
% underscores that tex would render as subscripts.
cleanup = with_graphics_defaults( ...
    'DefaultAxesFontSize',      s.tick_fs, ...
    'DefaultAxesLineWidth',     s.axis_lw, ...
    'DefaultTextInterpreter',   'tex', ...
    'DefaultLegendInterpreter', 'none');
end
