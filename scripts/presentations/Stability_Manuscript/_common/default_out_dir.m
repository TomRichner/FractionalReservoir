function out_dir = default_out_dir(requested, caller_fullpath)
% DEFAULT_OUT_DIR Where a figure writes: the caller's own folder unless told otherwise.
%
%   out_dir = DEFAULT_OUT_DIR(cfg.out_dir, mfilename('fullpath'))
%
% A figure called STANDALONE writes next to itself, which is convenient when
% iterating on one figure at the MATLAB prompt: the output appears beside the
% code that made it.
%
% A figure called from make_all_paper_figures is given an explicit out_dir under
% cfg.fig_root (default figs/paper/<entry name>), and that is where the paper's
% figures come from. The whole set then lands in one gitignored tree with a
% manifest recording the run_dir and commit that produced it, rather than being
% strewn through the source tree -- which had put ~20 .m files inside 618 MB of
% output, 216 MB of it tracked, and made `.fig` alone 813 MB of git history.
%
% So the default here is the CONVENIENCE path, not the authoritative one. If you
% are wondering where the paper's figures live, the answer is cfg.fig_root.
%
% Note this is a DATA path, not a path bootstrap: mfilename('fullpath') is the
% right tool here, unlike for locating the project root (use
% fileparts(which('setup_paths')) for that -- it is depth-independent).

arguments
    requested       (1,:) char
    caller_fullpath (1,:) char
end

if isempty(requested)
    out_dir = fileparts(caller_fullpath);
else
    out_dir = requested;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end
end
