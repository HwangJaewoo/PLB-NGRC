%% Z500 experiment: BO on the long-term statistics of the validation set, then long-term statistics and short-term correlation on the test set
function res = Z500_1D_experiment(cfg, vars, result_file)
    S = cfg.S;
    cfg.Pk_val = Spectrum(S(cfg.i_val, :));
    cfg.acf_val = Acf(S(cfg.i_val, :), cfg.lags);
    base = result_file(1:end - 4);
    diary([base, '_log.txt']);
    c_diary = onCleanup(@() diary('off'));
    fprintf('==== %s | %s ====\n', datestr(now, 'yyyy-mm-dd HH:MM'), result_file);

    if isfile(result_file)
        load(result_file, 'res');
        old = res.cfg;
        same = strcmp(old.model, cfg.model) && old.max_iter == cfg.max_iter ...
            && old.n_fit == cfg.n_fit && old.T_val == cfg.T_val && isequal(old.fixed, cfg.fixed);
        if ~same
            error('Saved result was made with different settings. Delete %s and rerun.', result_file);
        end
        best = res.best;
        fprintf('loaded %s (val %.4f)\n', result_file, res.val_obj);
    else
        t0 = tic;
        feas = @(T) Feasible(T, cfg.max_feat);
        cfg_saved = rmfield(cfg, 'S');
        bo = Run_BO(@(p) Objective(p, S, cfg), vars, cfg.max_iter, [base, '_bo.mat'], cfg_saved, ...
            {'XConstraintFcn', feas});
        n_use = min(cfg.max_iter, numel(bo.ObjectiveTrace));
        [f_best, i_best] = min(bo.ObjectiveTrace(1:n_use));
        best = Set_fixed(table2struct(bo.XTrace(i_best, :)), cfg.fixed);
        res = struct();
        res.best = best;
        res.val_obj = f_best;
        res.i_best = i_best;
        res.bo_x = bo.XTrace;
        res.bo_f = bo.ObjectiveTrace;
        res.t_bo = toc(t0);
        fprintf('BO done (val %.4f, %.1f min this session)\n', res.val_obj, res.t_bo / 60);
    end
    disp(best);

    mdl = Train_NGRC(S(1:cfg.n_fit, :), best, cfg.m);
    R = Residuals(mdl, S(1:cfg.n_fit, :));
    k = mdl.k;
    i0 = cfg.i_test(1) - 1;
    res.res_frac = var(R(:)) / var(S(1:cfg.n_fit, :), 0, 'all');

    X_test = S(cfg.i_test, :);
    n_half = floor(numel(cfg.i_test) / 2);
    XA = X_test(1:n_half, :);
    XB = X_test(n_half + 1:2 * n_half, :);
    res.floor = Long_metrics(XA, XB, cfg.lags);

    LM = nan(cfg.n_real, 4);
    ndiv = 0;
    res.P_show = [];
    for r = 1:cfg.n_real
        rng(100 + r, 'twister');
        [P, td] = Closed_loop(mdl, S(i0 - k + 1:i0, :), 2 * n_half, cfg.thr, R, cfg.gamma);
        if ~isnan(td)
            ndiv = ndiv + 1;
            fprintf('  realization %2d/%d: diverged\n', r, cfg.n_real);
            continue
        end
        if isempty(res.P_show)
            res.P_show = P;
        end
        LM(r, :) = 0.5 * (Long_metrics(XA, P(1:n_half, :), cfg.lags) ...
            + Long_metrics(XA, P(n_half + 1:end, :), cfg.lags));
        fprintf('  realization %2d/%d: D %.4f\n', r, cfg.n_real, LM(r, 1));
    end
    res.LM = LM;
    res.ndiv = ndiv;

    L = cfg.L;
    starts = cfg.i_test(1):cfg.fc_step:(cfg.i_test(end) - L);
    A = nan(numel(starts), L);
    E = nan(numel(starts), L);
    Ap = nan(numel(starts), L);
    Ep = nan(numel(starts), L);
    for j = 1:numel(starts)
        s0 = starts(j);
        [Pf, tdf] = Closed_loop(mdl, S(s0 - k + 1:s0, :), L, cfg.thr, [], 0);
        for l = 1:L
            if isnan(tdf) || l < tdf
                A(j, l) = Corr(Pf(l, :), S(s0 + l, :));
                E(j, l) = mean((Pf(l, :) - S(s0 + l, :)).^2);
            end
            Ap(j, l) = Corr(S(s0, :), S(s0 + l, :));
            Ep(j, l) = mean((S(s0, :) - S(s0 + l, :)).^2);
        end
    end
    res.acc = mean(A, 1, 'omitnan');
    res.acc_p = mean(Ap, 1);
    res.ss = 1 - mean(E, 1, 'omitnan') ./ mean(Ep, 1);
    res.hor = Horizon(res.acc, 0.6);
    res.hor_p = Horizon(res.acc_p, 0.6);

    fprintf('long-term: D %.4f (floor %.4f), spec %.3f (%.3f), acf %.3f (%.3f), std ratio %.2f | div %d/%d\n', ...
        mean(LM(:, 1), 'omitnan'), res.floor(1), mean(LM(:, 2), 'omitnan'), res.floor(2), ...
        mean(LM(:, 3), 'omitnan'), res.floor(3), mean(LM(:, 4), 'omitnan'), ndiv, cfg.n_real);
    fprintf('short-term: ACC 1/3/5 d %.2f/%.2f/%.2f, horizon %.2f d (persistence %.2f d)\n', ...
        res.acc(1), res.acc(3), res.acc(5), res.hor, res.hor_p);

    res.cfg = rmfield(cfg, 'S');
    save(result_file, 'res');
end

function f = Objective(p, S, cfg)
    if istable(p)
        p = table2struct(p);
    end
    p = Set_fixed(p, cfg.fixed);

    mdl = Train_NGRC(S(1:cfg.n_fit, :), p, cfg.m);
    R = Residuals(mdl, S(1:cfg.n_fit, :));

    rng(1, 'twister');
    s0 = cfg.n_fit;
    [P, td] = Closed_loop(mdl, S(s0 - mdl.k + 1:s0, :), cfg.T_val, cfg.thr, R, cfg.gamma);

    if isnan(td)
        Xv = S(cfg.i_val, :);
        f = Wasserstein1(Xv, P, 256) ...
            + norm(Spectrum(P) - cfg.Pk_val) / norm(cfg.Pk_val) ...
            + norm(Acf(P, cfg.lags) - cfg.acf_val) / norm(cfg.acf_val);
    else
        f = 3 + (1 - td / cfg.T_val);
    end
end

function ok = Feasible(T, max_feat)
    d = (2 * T.nn + 1) .* T.k;
    ok = (d + 1) .* (d + 2) / 2 <= max_feat;
end

function [P, td] = Closed_loop(mdl, H, T, thr, R, gamma)
    P = nan(T, size(H, 2));
    td = NaN;
    nR = size(R, 1);

    for t = 1:T
        y = Predict_NGRC(mdl, Delay_input(H, mdl.nn, mdl.m, mdl.k));
        x = reshape(y, 1, []);
        if gamma > 0
            x = x + gamma * R(randi(nR), :);
        end

        if any(~isfinite(x)) || max(abs(x)) > thr
            td = t;
            return
        end

        P(t, :) = x;
        H = [H(2:end, :); x];
    end
end

function R = Residuals(mdl, S)
    N = size(S, 2);
    U = Delay_input(S(1:end-1, :), mdl.nn, mdl.m, mdl.k);
    Y = reshape(S(mdl.k+1:end, :), [], 1);

    Yhat = zeros(size(Y));
    cs = 20000;
    for i0 = 1:cs:numel(Y)
        ii = i0:min(i0 + cs - 1, numel(Y));
        Yhat(ii) = Predict_NGRC(mdl, U(ii, :));
    end

    R = reshape(Y - Yhat, [], N);
end

function v = Long_metrics(Ref, Sim, lags)
    Pr = Spectrum(Ref);
    Ps = Spectrum(Sim);
    ar = Acf(Ref, lags);
    as = Acf(Sim, lags);
    v = [Wasserstein1(Ref, Sim, 512), norm(Ps - Pr) / norm(Pr), norm(as - ar) / norm(ar), std(Sim(:)) / std(Ref(:))];
end

function Pk = Spectrum(X)
    Xc = X - mean(X, 'all');
    Pk = mean(abs(fft(Xc, [], 2)).^2, 1) / size(X, 2);
end

function r = Acf(X, lags)
    r = zeros(1, numel(lags));
    for i = 1:numel(lags)
        l = lags(i);
        A = X(1:end - l, :);
        B = X(1 + l:end, :);
        A = A - mean(A);
        B = B - mean(B);
        r(i) = mean(sum(A .* B) ./ sqrt(sum(A.^2) .* sum(B.^2)));
    end
end

function c = Corr(a, b)
    a = a - mean(a);
    b = b - mean(b);
    c = sum(a .* b) / sqrt(sum(a.^2) * sum(b.^2));
end

function h = Horizon(acc, thr)
    a = [1, acc];
    idx = find(a < thr, 1);
    if isempty(idx)
        h = inf;
        return
    end
    h = (idx - 2) + (a(idx - 1) - thr) / (a(idx - 1) - a(idx));
end
