function name = verbose_name(v)
%VERBOSE_NAME Canonical char for a verbosity setting ('verbose'|'minimal'|'near-none').
%
%   name = VERBOSE_NAME(v)
%
% Accepts what verbose_level accepts (the three chars, or a logical) and
% returns the canonical char, so a property setter can store one form.
%
% See also: verbose_level, vprintf

names = {'near-none', 'minimal', 'verbose'};
name = names{verbose_level(v) + 1};
end
