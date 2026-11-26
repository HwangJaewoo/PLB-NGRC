%% Settings, BO search ranges and result file of P-NGRC or PLB-NGRC for the Z500 data along 50N
function [cfg, vars, result_file] = Z500_1D_config(method)
    D = load(Z500_1D_data(), 'S', 'N', 'n_fit', 'i_val', 'i_test');
    cfg.N = D.N;
    cfg.m = 1;
    cfg.n_fit = D.n_fit;
    cfg.i_val = D.i_val;
    cfg.i_test = D.i_test;
    cfg.T_val = numel(D.i_val);
    cfg.gamma = 1;
    cfg.thr = 10;
    cfg.lags = 1:10;
    cfg.max_iter = 100;
    cfg.max_feat = 300;
    cfg.n_real = 10;
    cfg.L = 15;
    cfg.fc_step = 10;
    v_lam = optimizableVariable('log_lam', [-10, -1]);
    v_nn = optimizableVariable('nn', [1, 3], 'Type', 'integer');
    v_k = optimizableVariable('k', [1, 4], 'Type', 'integer');
    v_s = optimizableVariable('s', [1, 4], 'Type', 'integer');
    v_eta = optimizableVariable('eta', [0.1, 0.3]);
    if strcmp(method, 'PNGRC')
        cfg.model = 'PNGRC';
        vars = [v_lam, v_nn, v_k];
        cfg.fixed = struct('s', 0);
    else
        cfg.model = 'PLB';
        vars = [v_lam, v_nn, v_k, v_s, v_eta];
        cfg.fixed = struct();
    end
    cfg.S = D.S;
    d = fullfile(Repo_root, 'experiments', 'Z500_1D', 'results');
    if ~isfolder(d)
        mkdir(d);
    end
    result_file = fullfile(d, ['Z500_1D_', method, '.mat']);
end
