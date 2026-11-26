%% Kuramoto-Sivashinsky benchmark: P-NGRC and PLB-NGRC for L = 22, 44, 100, 200, then the spatial power spectra of the long-term predictions
function KS_run()
    LQ = [22, 64; 44, 128; 100, 256; 200, 512];
    methods = {'PNGRC', 'PLBNGRC'};
    for a = 1:size(LQ, 1)
        for i = 1:numel(methods)
            [cfg, vars, f] = Config(LQ(a, :), methods{i});
            cfg.skip_done = true;
            Run_experiment(cfg, vars, f);
        end
    end
    Spectra(LQ, methods);
end

function [cfg, vars, result_file] = Config(LQ, method)
    [data_file, lle_file] = KS_data(LQ(1), LQ(2));
    D = load(data_file, 'S', 'N', 'dt', 'n_fit', 'i_val', 'i_test', 'starts');
    L = load(lle_file, 'lle', 'n_trans', 'n_eval');
    cfg.N = D.N;
    cfg.m = 1;
    cfg.dt = D.dt;
    cfg.n_fit = D.n_fit;
    cfg.i_val = D.i_val;
    cfg.i_test = D.i_test;
    cfg.starts = D.starts;
    cfg.T_cl = numel(D.i_test);
    cfg.n_long = 10;
    cfg.thr = 10;
    cfg.metric = @(A, B) Wasserstein1(A, B, 256);
    cfg.metric_test = @(A, B) Wasserstein1(A, B, 512);
    cfg.max_iter = 100;
    cfg.lle_true = mean(L.lle);
    cfg.lle_trans = L.n_trans;
    cfg.lle_eval = L.n_eval;
    cfg.vpt_thr = 0.3;
    cfg.vpt_thr_extra = 0.5;
    cfg.vpt_len = round(10 / mean(L.lle) / D.dt);
    cfg.n_vpt = cfg.n_long;
    v_lam = optimizableVariable('log_lam', [-10, -1]);
    v_nn = optimizableVariable('nn', [1, 8], 'Type', 'integer');
    v_s = optimizableVariable('s', [1, 8], 'Type', 'integer');
    v_eta = optimizableVariable('eta', [0.1, 0.3]);
    if strcmp(method, 'PNGRC')
        cfg.model = 'PNGRC';
        vars = [v_lam, v_nn];
        cfg.fixed = struct('k', 1, 's', 0);
    else
        cfg.model = 'PLB';
        vars = [v_lam, v_nn, v_s, v_eta];
        cfg.fixed = struct('k', 1);
    end
    cfg.S = D.S;
    result_file = fullfile(Repo_root, 'experiments', 'KS', 'results', sprintf('KS_L%03d_%s.mat', LQ(1), method));
end

function Spectra(LQ, methods)
    n_long = 10;
    thr = 10;
    nL = size(LQ, 1);
    nM = numel(methods);
    err_mean = nan(nL, nM);
    err_sd = nan(nL, nM);
    err_avgspec = nan(nL, nM);
    n_surv = nan(nL, nM);
    check = nan(nL, nM);
    floor_mean = nan(nL, 1);
    floor_sd = nan(nL, 1);
    floor_avgspec = nan(nL, 1);
    spec = struct('L', cell(1, nL), 'kk', [], 'Pt', [], 'Pfloor', [], 'Pm', []);
    for a = 1:nL
        data_file = KS_data(LQ(a, 1), LQ(a, 2));
        D = load(data_file, 'S', 'n_fit', 'i_test', 'L', 'Q');
        S = D.S;
        T = numel(D.i_test);
        Pt = Spec(S(D.i_test, :));
        nt = norm(Pt);
        starts_floor = round(linspace(1, D.n_fit - T + 1, n_long));
        Pf = zeros(n_long, numel(Pt));
        for j = 1:n_long
            Pf(j, :) = Spec(S(starts_floor(j):starts_floor(j) + T - 1, :));
        end
        ef = vecnorm(Pf - Pt, 2, 2) / nt;
        floor_mean(a) = mean(ef);
        floor_sd(a) = std(ef);
        floor_avgspec(a) = norm(mean(Pf, 1) - Pt) / nt;
        spec(a).L = LQ(a, 1);
        spec(a).kk = (0:D.Q - 1) * 2 * pi / D.L;
        spec(a).Pt = Pt;
        spec(a).Pfloor = mean(Pf, 1);
        spec(a).Pm = nan(nM, numel(Pt));
        for i = 1:nM
            [~, ~, f] = Config(LQ(a, :), methods{i});
            R = load(f, 'res');
            mdl = Train_NGRC(S(1:D.n_fit, :), R.res.best, 1);
            i0 = D.i_test(1) - 1;
            starts = round(linspace(i0, D.i_test(end) - 1, n_long));
            Pm = nan(n_long, numel(Pt));
            for j = 1:n_long
                s0 = starts(j);
                [P, td] = Closed_loop_NGRC(mdl, S(s0 - mdl.k + 1:s0, :), T, thr);
                if j == 1
                    check(a, i) = max(abs(P(:) - R.res.P_long(:)));
                end
                if isnan(td)
                    Pm(j, :) = Spec(P);
                end
            end
            ok = all(isfinite(Pm), 2);
            n_surv(a, i) = sum(ok);
            e = vecnorm(Pm(ok, :) - Pt, 2, 2) / nt;
            err_mean(a, i) = mean(e);
            err_sd(a, i) = std(e);
            Pavg = mean(Pm(ok, :), 1);
            err_avgspec(a, i) = norm(Pavg - Pt) / nt;
            spec(a).Pm(i, :) = Pavg;
            fprintf('KS L = %d %s: spectrum error %.3f +- %.3f, error of the mean spectrum %.3f\n', ...
                LQ(a, 1), methods{i}, err_mean(a, i), err_sd(a, i), err_avgspec(a, i));
        end
    end
    save(fullfile(Repo_root, 'experiments', 'KS', 'results', 'KS_spectrum.mat'), 'LQ', 'methods', 'err_mean', 'err_sd', ...
        'err_avgspec', 'n_surv', 'check', 'floor_mean', 'floor_sd', 'floor_avgspec', 'spec');
end

function P = Spec(X)
    Xc = X - mean(X, 'all');
    P = mean(abs(fft(Xc, [], 2)).^2, 1) / size(X, 2);
end
