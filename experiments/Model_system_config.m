%% Settings, BO search ranges and result file of one method for the Lorenz-96 system, the coupled Ikeda map or the Thomas system
function [cfg, vars, result_file] = Model_system_config(sys, method, N, noise)
    if nargin < 4
        noise = 0;
    end
    feval([sys, '_data'], N);
    [data_file, lle_file] = Data_files(sys, N);
    D = load(data_file);
    L = load(lle_file, 'lle', 'n_trans', 'n_eval');

    switch sys
        case 'Lorenz96'
            m = 1;
            dt = D.dt;
            T_cl = numel(D.i_test);
            thr = 10;
            vpt_len = round(10 / dt);
            ord = [];
        case 'Ikeda'
            m = 2;
            dt = 1;
            T_cl = 500;
            thr = 20;
            vpt_len = 100;
            ord = 3;
        case 'Thomas'
            m = 1;
            dt = D.dt;
            T_cl = numel(D.i_test);
            thr = 10;
            vpt_len = round(10 / mean(L.lle) / dt);
            ord = 3;
    end

    S = D.S;
    if noise > 0
        S = Add_noise(D, m, noise);
    end

    cfg.N = D.N;
    cfg.m = m;
    cfg.dt = dt;
    if noise > 0
        cfg.sig = noise;
    end
    cfg.n_fit = D.n_fit;
    cfg.i_val = D.i_val;
    cfg.i_test = D.i_test;
    cfg.starts = D.starts;
    cfg.T_cl = T_cl;
    cfg.n_long = 10;
    cfg.thr = thr;
    cfg.metric = @(A, B) Wasserstein1(A, B, 256, m);
    cfg.metric_test = @(A, B) Wasserstein1(A, B, 512, m);
    cfg.max_iter = 100;
    cfg.lle_true = mean(L.lle);
    cfg.lle_trans = L.n_trans;
    cfg.lle_eval = L.n_eval;
    cfg.vpt_thr = 0.3;
    cfg.vpt_len = vpt_len;
    cfg.n_vpt = cfg.n_long;

    v_lam = optimizableVariable('log_lam', [-10, -1]);
    v_nn = optimizableVariable('nn', [1, 3], 'Type', 'integer');
    v_s = optimizableVariable('s', [1, 10], 'Type', 'integer');
    v_eta = optimizableVariable('eta', [0.1, 0.3]);
    switch method
        case 'PNGRC'
            cfg.model = 'PNGRC';
            vars = [v_lam, v_nn];
            cfg.fixed = struct('k', 1, 's', 0);
        case {'PLBNGRC', 'PLBNGRC_inball'}
            cfg.model = 'PLB';
            vars = [v_lam, v_nn, v_s, v_eta];
            cfg.fixed = struct('k', 1);
        case {'LBNGRC', 'LBNGRC_weighted'}
            fit = 'inball';
            if strcmp(method, 'LBNGRC_weighted')
                fit = 'weighted';
            end
            cfg.model = ['LBfull_', fit];
            cfg.N = 1;
            cfg.m = m * D.N;
            vars = [v_lam, v_s, v_eta];
            cfg.fixed = struct('k', 1, 'nn', 0, 'fit', fit);
        case 'PRC'
            cfg.model = 'RC';
            vars = [v_lam, v_nn, optimizableVariable('rho', [0.1, 1.5]), ...
                optimizableVariable('sig_in', [0.1, 2]), optimizableVariable('leak', [0.1, 1])];
            cfg.fixed = struct('Dr', 1200, 'seed', 1, 'n_warm', 500);
    end
    if ~isempty(ord) && ~strcmp(method, 'PRC')
        cfg.fixed.ord = ord;
    end
    if strcmp(method, 'PLBNGRC_inball')
        cfg.fixed.fit = 'inball';
    end
    cfg.S = S;

    name = sprintf('%s_N%03d_%s', sys, N, method);
    if noise > 0
        name = sprintf('%s_noise%02d', name, round(100 * noise));
    end
    result_file = fullfile(Repo_root, 'experiments', sys, 'results', [name, '.mat']);
end

function S = Add_noise(D, m, sig)
    if isfield(D, 'S_raw')
        S_raw = D.S_raw;
    else
        S_raw = D.S * D.sd + D.mu;
    end
    rng(10);
    E = randn(size(S_raw));
    S = S_raw;
    for v = 1:m
        cols = (v - 1) * D.N + (1:D.N);
        sd0 = std(S_raw(1:D.n_fit, cols), 0, 'all');
        S_obs = S_raw(:, cols) + sig * sd0 * E(:, cols);
        mu_o = mean(S_obs(1:D.n_fit, :), 'all');
        sd_o = std(S_obs(1:D.n_fit, :), 0, 'all');
        S(:, cols) = (S_obs - mu_o) / sd_o;
    end
end
