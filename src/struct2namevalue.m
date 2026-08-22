function nv = struct2namevalue(s)
% STRUCT2NAMEVALUE Flatten a scalar struct into a {name, value, ...} cell row.
%
%   args = STRUCT2NAMEVALUE(srnn_param_preset('celltype_pairs'));
%   model = SRNNCellTypePairs(args{:});
%
% The bridge between the struct form a preset (or a condition) is stored in and
% the name-value form the model constructors take. Written out by hand as a
% local subfunction in several scripts before this existed.
%
% Note a preset field whose value is a CELL -- cell_type_names, tau_a -- is
% passed through as that cell, which is what the constructor expects; this is
% why the loop assigns rather than using a comma-separated list expansion.

arguments
    s struct {mustBeScalarOrEmpty}
end

if isempty(s)
    nv = {};
    return
end

f = fieldnames(s);
nv = cell(1, 2*numel(f));
for i = 1:numel(f)
    nv{2*i-1} = f{i};
    nv{2*i}   = s.(f{i});
end
end
