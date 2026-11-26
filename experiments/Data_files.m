%% Paths of the data file and of the true largest Lyapunov exponent of a model system with N sites (folders are created if missing)
function [data_file, lle_file] = Data_files(sys, N)
    d = fullfile(Repo_root, 'experiments', sys);
    for f = {'data', 'results'}
        if ~isfolder(fullfile(d, f{1}))
            mkdir(fullfile(d, f{1}));
        end
    end
    data_file = fullfile(d, 'data', sprintf('%s_N%03d_data.mat', sys, N));
    lle_file = fullfile(d, 'results', sprintf('%s_N%03d_true_lyapunov.mat', sys, N));
end
