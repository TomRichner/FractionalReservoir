function net = outlier_network(network)
%OUTLIER_NETWORK From a network's name to what outlier_run builds: the model arguments, any added structure, and which mode is the intended one.
%
%   net = outlier_network('EEx1p5')
%
% THE NAMES. A name is one or more tokens joined by '-'. Each token changes
% one thing of the paper's network (mu7chain preset, block means
% [7 -7; 7 -7] in the preset's units, 500 neurons, indegree 100):
%   reference        nothing changed
%   EI3to2           block means [8.4 -5.6; 8.4 -5.6]
%   <XY>x<a>p<b>     the block mean mu_XY times a.b, XY one of EE, EI, IE,
%                    II. The letters are (POSTsynaptic, PREsynaptic), as in
%                    the class's mu_XY_relative: EI is the I-to-E mean, IE
%                    the E-to-I mean. 'EEx1p5' is E-to-E times 1.5;
%                    'EEx1p5-EIx1p3-IEx1p3' strengthens E-to-E and both
%                    arms of the E-I loop.
%   ind<k>           indegree k instead of 100 (the weight unit F follows,
%                    as the class does by default)
%   lrA<a>p<b>W<c>p<d>   low-rank structure added to the drawn matrix
%                    (OutlierNet): an eigenvalue at R (a.b + i c.d), R the
%                    structural bulk radius; rank one and real when the
%                    imaginary part is 0, a rank-two conjugate pair
%                    otherwise. Dense, on random patterns of neurons,
%                    ignores Dale's law.
%
% THE INTENDED MODE (net.target), the eigenvalue of the structural W that
% the measures follow by its eigenvector (outlier_mode_reference):
%   'lowrank'    the eigenvalue of W nearest the one the addition was built
%                to have
%   'blockmean'  the eigenvalue of W nearest the rightmost eigenvalue of
%                the realized block-mean matrix (any <XY>x token, EI3to2)
%   'rightmost'  the rightmost eigenvalue of W (reference, ind<k>: no
%                structure was put in, so the mode of interest is whatever
%                the draw left furthest right; for a seed chosen for its
%                local outlier that is the local outlier)
% With several tokens the first that applies in the order above decides.
%
% Output: struct with fields name, net_args (name-value cell for the model
% constructor), added_structure ([] or the struct OutlierNet takes),
% target, label (for a figure's row or column title, TeX).
%
% See also: outlier_run, OutlierNet, outlier_mode_reference

arguments
    network (1,:) char
end

base = 7;                                   % the preset's |mu_tilde_relative| on every route
mult = ones(2);                             % (post, pre) multipliers of the block means
net_args = {};
added = [];
labels = {};
has_block = false;
route_names = struct('EE', 'E-to-E', 'EI', 'I-to-E', 'IE', 'E-to-I', 'II', 'I-to-I');
idx = struct('E', 1, 'I', 2);

tokens = strsplit(network, '-');
for k = 1:numel(tokens)
    tk = tokens{k};
    route = regexp(tk, '^(EE|EI|IE|II)x(\d+)p(\d+)$', 'tokens', 'once');
    ind   = regexp(tk, '^ind(\d+)$', 'tokens', 'once');
    lr    = regexp(tk, '^lrA(\d+)p(\d+)W(\d+)p(\d+)$', 'tokens', 'once');
    if strcmp(tk, 'reference')
        labels{end + 1} = 'reference'; %#ok<AGROW>
    elseif strcmp(tk, 'EI3to2')
        mult = [1.2 0.8; 1.2 0.8];
        has_block = true;
        labels{end + 1} = 'E:I 3:2'; %#ok<AGROW>
    elseif ~isempty(route)
        f = str2double([route{2} '.' route{3}]);
        mult(idx.(route{1}(1)), idx.(route{1}(2))) = f;
        has_block = true;
        labels{end + 1} = sprintf('%s \\times %s', route_names.(route{1}), num2str(f)); %#ok<AGROW>
    elseif ~isempty(ind)
        net_args = [net_args, {'indegree', str2double(ind{1})}]; %#ok<AGROW>
        labels{end + 1} = sprintf('indegree %s', ind{1}); %#ok<AGROW>
    elseif ~isempty(lr)
        added = struct('a', str2double([lr{1} '.' lr{2}]), 'w', str2double([lr{3} '.' lr{4}]));
        if added.w == 0
            labels{end + 1} = sprintf('rank one, %g R', added.a); %#ok<AGROW>
        else
            labels{end + 1} = sprintf('rank two, (%g \\pm %gi) R', added.a, added.w); %#ok<AGROW>
        end
    else
        error('outlier_network:Name', ['Unknown token ''%s'' in network ''%s'' (reference, EI3to2, ' ...
            '<EE|EI|IE|II>x<a>p<b>, ind<k>, lrA<a>p<b>W<c>p<d>, joined by ''-'').'], tk, network);
    end
end

if has_block
    net_args = [net_args, {'mu_tilde_relative', base .* mult .* [1 -1; 1 -1]}];
end
if ~isempty(added)
    target = 'lowrank';
elseif has_block
    target = 'blockmean';
else
    target = 'rightmost';
end
if numel(labels) > 1
    labels = labels(~strcmp(labels, 'reference'));
end
net = struct('name', network, 'net_args', {net_args}, 'added_structure', added, ...
    'target', target, 'label', strjoin(labels, ', '));
end
