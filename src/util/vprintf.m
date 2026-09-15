function vprintf(v, needed, fmt, varargin)
%VPRINTF fprintf only when the verbosity setting v is at or above `needed`.
%
%   VPRINTF(v, needed, fmt, ...)
%
% v is the caller's verbosity ('verbose' | 'minimal' | 'near-none', or a
% logical, see verbose_level); needed is the level the message belongs to.
% A per-model or per-job line is 'verbose'; a per-stage or per-batch line is
% 'minimal'; the final outcome of an entry point is 'near-none' (printed at
% every level). Output goes to stdout AND, when an entry point has opened one
% (vlog), to the transcript log. Failures should use vfail, ungated, so they
% are never hidden.
%
% See also: verbose_level, vlog, vfail

if verbose_level(v) >= verbose_level(needed)
    s = sprintf(fmt, varargin{:});
    fprintf('%s', s);
    vlog('append', s);
end
end
