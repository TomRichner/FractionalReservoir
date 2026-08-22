function out_dir = default_out_dir(requested, caller_fullpath)
% DEFAULT_OUT_DIR Where a figure writes: the caller's own folder unless told otherwise.
%
%   out_dir = DEFAULT_OUT_DIR(cfg.out_dir, mfilename('fullpath'))
%
% Every manuscript figure defaults to writing next to itself, which is what
% makes the folder self-contained: the script, its README, its git provenance
% and its outputs all sit together, and the manuscript references that path.
% Passing a non-empty out_dir overrides it (the master script does this to
% collect a whole set somewhere else).
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
