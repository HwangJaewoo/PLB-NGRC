%% Training time of every selected model (median of three fits on one computer) for Table I and the scaling figure; writes timing/training_time.csv
function Measure_training_time()
    n_rep = 3;
    delete(gcp('nocreate'));
    root = Repo_root;
    systems = {'Lorenz96', 'Ikeda', 'Thomas'};
    Ls = cell(numel(systems), 1);
    for a = 1:numel(systems)
        Ls{a} = dir(fullfile(root, 'experiments', systems{a}, 'results', [systems{a}, '_N*.mat']));
    end
    L = vertcat(Ls{:});
    L = L(~contains({L.name}, {'_bo', '_noise', '_true_lyapunov'}));
    n = numel(L);
    files = cell(n + 4, 1);
    times = zeros(n + 4, 1);

    for i = 1:n
        nm = L(i).name;
        sys = nm(1:strfind(nm, '_N') - 1);
        R = load(fullfile(L(i).folder, nm), 'res');
        D = load(fullfile(root, 'experiments', sys, 'data', [nm(1:numel(sys) + 5), '_data.mat']), 'S');
        S_fit = D.S(1:R.res.cfg.n_fit, :);
        if strcmp(R.res.cfg.model, 'RC')
            train = @() PRC_train(S_fit, R.res.best, R.res.cfg.m);
        else
            train = @() Train_NGRC(S_fit, R.res.best, R.res.cfg.m);
        end
        files{i} = nm;
        times(i) = Median_time(train, n_rep);
    end

    D = load(Z500_1D_data(), 'S', 'n_fit');
    for method = {'PNGRC', 'PLBNGRC'}
        nm = ['Z500_1D_', method{1}, '.mat'];
        R = load(fullfile(root, 'experiments', 'Z500_1D', 'results', nm), 'res');
        n = n + 1;
        files{n} = nm;
        times(n) = Median_time(@() Train_NGRC(D.S(1:D.n_fit, :), R.res.best, 1), n_rep);
    end

    cfg = Z500_2D_setup();
    for variant = {'PNGRC', 'PLBNGRC'}
        nm = ['Z500_2D_', variant{1}, '_best.mat'];
        B = load(fullfile(root, 'experiments', 'Z500_2D', 'results', nm), 'p');
        n = n + 1;
        files{n} = nm;
        times(n) = Median_time(@() Z500_2D_train(cfg, B.p), n_rep);
    end

    T = table(files, times, 'VariableNames', {'file', 't_train_s'});
    disp(T);
    writetable(T, fullfile(root, 'timing', 'training_time.csv'));
    fid = fopen(fullfile(root, 'timing', 'computer.txt'), 'w');
    fprintf(fid, 'date: %s\ncpu: %s\nlogical processors: %s\nMATLAB computational threads: %d\nMATLAB: %s\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM'), getenv('PROCESSOR_IDENTIFIER'), getenv('NUMBER_OF_PROCESSORS'), maxNumCompThreads, version);
    fclose(fid);
end

function t = Median_time(train, n_rep)
    tr = zeros(n_rep, 1);
    for r = 1:n_rep
        t0 = tic;
        train();
        tr(r) = toc(t0);
    end
    t = median(tr);
end
