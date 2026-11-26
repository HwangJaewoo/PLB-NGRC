%% Absolute path of the repository root (the folder that contains core/)
function r = Repo_root()
    r = fileparts(fileparts(mfilename('fullpath')));
end
