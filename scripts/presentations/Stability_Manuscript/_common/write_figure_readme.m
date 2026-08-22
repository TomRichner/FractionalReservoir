function readme_path = write_figure_readme(out_dir, info)
% WRITE_FIGURE_README Plain-text record of how one manuscript figure was made.
%
%   readme_path = WRITE_FIGURE_README(out_dir, info)
%
% Writes README_<tag>.txt into out_dir. Replaces the hand-rolled
% fprintf(fid, ...) blocks that occupied 18-91 lines in each of ten figure
% scripts (~350 lines in total), all of them re-deriving the same layout, the
% same fopen/onCleanup dance, and the same "FIGURES PRODUCED" list.
%
% Every field of INFO is optional except `tag`; a section is emitted only if its
% field is non-empty, so a conceptual figure with no source run simply has no
% SOURCE section rather than one saying "none".
%
%   info.tag       (required) short name; becomes README_<tag>.txt
%   info.title     one-line headline, underlined with '='
%   info.script    the function that produced the figure
%   info.what      char or cellstr -- WHAT IT SHOWS, free prose
%   info.how       char or cellstr -- HOW IT WAS MADE
%   info.source    struct: field -> value, rendered as an indented list.
%                  Use for the run directory, the preset, the .mat file.
%   info.settings  struct of parameter name -> value, rendered as a table.
%                  Values may be numeric, char, logical or cellstr.
%   info.figures   cellstr of the files written, or a struct array with
%                  .name and .desc
%   info.notes     char or cellstr -- caveats worth reading before believing
%                  the figure. This is where "the histogram range crops 32% of
%                  the distribution" belongs.
%   info.sections  struct array with .heading and .body for anything the fixed
%                  sections above do not cover, appended in order.
%
% Prose fields accept a cellstr; each element becomes a paragraph, and long
% lines are wrapped to 78 columns so the file reads in a terminal.
%
% See also: manuscript_style, save_figure_stable, capture_git_provenance

arguments
    out_dir (1,:) char
    info    struct
end

if ~isfield(info, 'tag') || isempty(info.tag)
    error('write_figure_readme:NoTag', 'info.tag is required.');
end
if ~isfolder(out_dir); mkdir(out_dir); end

readme_path = fullfile(out_dir, sprintf('README_%s.txt', info.tag));
fid = fopen(readme_path, 'w');
if fid < 0
    error('write_figure_readme:CannotOpen', 'Could not open %s', readme_path);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

W = 78;   % wrap column

%% Header
title = getf(info, 'title', sprintf('Stability_Manuscript figure: %s', info.tag));
fprintf(fid, '%s\n', title);
fprintf(fid, '%s\n\n', repmat('=', 1, min(numel(title), W)));

fprintf(fid, 'Generated: %s\n', char(datetime('now')));
if isfield(info, 'script') && ~isempty(info.script)
    fprintf(fid, 'By:        %s\n', info.script);
end
fprintf(fid, '\n');

%% Prose sections
emit_prose(fid, 'WHAT IT SHOWS',  getf(info, 'what', ''), W);
emit_prose(fid, 'HOW IT WAS MADE', getf(info, 'how',  ''), W);

%% Source
if isfield(info, 'source') && ~isempty(info.source) && isstruct(info.source)
    fprintf(fid, 'SOURCE\n');
    f = fieldnames(info.source);
    width = max(cellfun(@numel, f));
    for k = 1:numel(f)
        fprintf(fid, '  %-*s  %s\n', width, f{k}, value_to_str(info.source.(f{k})));
    end
    fprintf(fid, '\n');
end

%% Settings table
if isfield(info, 'settings') && ~isempty(info.settings) && isstruct(info.settings)
    fprintf(fid, 'PARAMETERS AS RUN\n');
    f = fieldnames(info.settings);
    width = max(cellfun(@numel, f));
    for k = 1:numel(f)
        fprintf(fid, '  %-*s  %s\n', width, f{k}, value_to_str(info.settings.(f{k})));
    end
    fprintf(fid, '\n');
end

%% Figures produced
if isfield(info, 'figures') && ~isempty(info.figures)
    fprintf(fid, 'FIGURES PRODUCED (in this folder)\n');
    figs = info.figures;
    if isstruct(figs)
        for k = 1:numel(figs)
            fprintf(fid, '  %s\n', figs(k).name);
            if isfield(figs, 'desc') && ~isempty(figs(k).desc)
                emit_wrapped(fid, figs(k).desc, W, '      ');
            end
        end
    else
        if ischar(figs); figs = {figs}; end
        for k = 1:numel(figs)
            fprintf(fid, '  %s\n', figs{k});
        end
    end
    fprintf(fid, '\n');
end

%% Notes / caveats
emit_prose(fid, 'READING THIS FIGURE', getf(info, 'notes', ''), W);

%% Free-form extra sections
if isfield(info, 'sections') && ~isempty(info.sections)
    for k = 1:numel(info.sections)
        emit_prose(fid, upper(info.sections(k).heading), info.sections(k).body, W);
    end
end
end

%% ------------------------------------------------------------------------
function v = getf(s, name, default)
if isfield(s, name) && ~isempty(s.(name)); v = s.(name); else; v = default; end
end

function emit_prose(fid, heading, body, W)
if isempty(body); return; end
fprintf(fid, '%s\n', heading);
if ischar(body); body = {body}; end
for k = 1:numel(body)
    emit_wrapped(fid, body{k}, W, '  ');
    if k < numel(body); fprintf(fid, '\n'); end
end
fprintf(fid, '\n');
end

function emit_wrapped(fid, text, W, indent)
% Wrap to W columns, preserving explicit newlines as paragraph breaks.
if isempty(text); return; end
paras = strsplit(text, newline);
for p = 1:numel(paras)
    words = strsplit(strtrim(paras{p}));
    words = words(~cellfun(@isempty, words));
    if isempty(words); fprintf(fid, '\n'); continue; end
    line = indent;
    for w = 1:numel(words)
        if numel(line) + numel(words{w}) + 1 > W && ~strcmp(line, indent)
            fprintf(fid, '%s\n', line);
            line = indent;
        end
        if strcmp(line, indent); line = [line words{w}]; %#ok<AGROW>
        else;                    line = [line ' ' words{w}]; %#ok<AGROW>
        end
    end
    fprintf(fid, '%s\n', line);
end
end

function str = value_to_str(v)
% One-line rendering of a settings value. Mirrors what write_run_parameters_md
% does for the sweep directories, so a figure README and a run's parameters.md
% describe the same value the same way.
if ischar(v)
    str = v;
elseif isstring(v) && isscalar(v)
    str = char(v);
elseif islogical(v)
    if isscalar(v); str = mat2str(v); else; str = mat2str(v); end
elseif isnumeric(v)
    if isempty(v);          str = '[]';
    elseif isscalar(v);     str = num2str(v, '%g');
    else;                   str = mat2str(v, 6);
    end
elseif iscellstr(v) %#ok<ISCLSTR>
    str = ['{' strjoin(v, ', ') '}'];
elseif isstruct(v)
    str = ['struct: ' strjoin(fieldnames(v)', ', ')];
elseif isa(v, 'function_handle')
    str = func2str(v);
else
    str = sprintf('<%s>', class(v));
end
end
