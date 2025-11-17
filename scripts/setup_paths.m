function setup_paths()
%SETUP_PATHS Add the repository src directory to the MATLAB path.
%
% Call this helper before running scripts that rely on code located under
% the src directory (sibling of scripts).

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
srcPath = fullfile(projectRoot, 'src');

if ~isfolder(srcPath)
    error('setup_paths:MissingSrc', ...
        'Could not find src directory at %s', srcPath);
end

addpath(genpath(scriptDir));
addpath(genpath(srcPath));
end

