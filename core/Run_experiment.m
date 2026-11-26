%% Hyperparameter selection by BO on the validation set, then evaluation on the test set (VPT, W1, largest Lyapunov exponent, training time)
function res = Run_experiment(cfg, vars, result_file)
    S = cfg.S;
    E = Engine(cfg.model);
    sigma = std(S(1:cfg.n_fit, :), 0, 1);
    base = result_file(1:end - 4);
    diary([base, '_log.txt']);
    c_diary = onCleanup(@() diary('off'));
    fprintf('==== %s | %s ====\n', datestr(now, 'yyyy-mm-dd HH:MM'), result_file);

    if isfile(result_file)
        load(result_file, 'res');
        Check_settings(res.cfg, cfg, result_file);
        best = res.best;
        fprintf('loaded %s (val %.4f)\n', result_file, res.val_obj);
        if isfield(cfg, 'skip_done') && cfg.skip_done && isfield(res, 'lle_all') && isfield(res, 't_train')
            fprintf('evaluation already complete, skipped\n');
            return
        end
    else
        cfg_saved = rmfield(cfg, 'S');
        t0 = tic;
        bo = Run_BO(@(p) Objective(p, S, cfg, E), vars, cfg.max_iter, [base, '_bo.mat'], cfg_saved);
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
        res.var_name = arrayfun(@(v) v.Name, vars, 'UniformOutput', false);
        res.var_range = arrayfun(@(v) v.Range, vars, 'UniformOutput', false);
        res.var_type = arrayfun(@(v) v.Type, vars, 'UniformOutput', false);
        fprintf('BO done (val %.4f, %.1f min this session)\n', res.val_obj, res.t_bo / 60);
    end
    disp(best);

    n_rep = 3;
    t_rep = zeros(n_rep, 1);
    for r = 1:n_rep
        t1 = tic;
        mdl = E.train(S(1:cfg.n_fit, :), best, cfg.m);
        t_rep(r) = toc(t1);
    end
    res.t_train = median(t_rep);
    fprintf('training time: %.3f s (median of %d)\n', res.t_train, n_rep);
    nh = E.n_hist(mdl);
    i0 = cfg.i_test(1) - 1;

    starts_long = round(linspace(i0, cfg.i_test(end) - 1, cfg.n_long));
    X_test = S(cfg.i_test, :);
    div = false(cfg.n_long, 1);
    D_all = nan(cfg.n_long, 1);
    for j = 1:cfg.n_long
        s0 = starts_long(j);
        [P, td] = E.loop(mdl, S(s0 - nh + 1:s0, :), cfg.T_cl, cfg.thr);
        if j == 1
            res.P_long = P;
        end
        div(j) = ~isnan(td);
        if ~div(j)
            D_all(j) = cfg.metric_test(X_test, P);
            fprintf('  long-term %2d/%d: survived, D = %.4f\n', j, cfg.n_long, D_all(j));
        else
            fprintf('  long-term %2d/%d: diverged at step %d\n', j, cfg.n_long, td);
        end
    end
    res.diverged = div(1);
    res.div_all = div;
    res.surv = mean(~div);
    res.D = mean(D_all, 'omitnan');
    res.D_all = D_all;
    fprintf('long-term: survived %d / %d runs of %d steps | D = %.4f\n', sum(~div), cfg.n_long, cfg.T_cl, res.D);

    starts_floor = round(linspace(1, cfg.n_fit - cfg.T_cl + 1, cfg.n_long));
    D_fl = zeros(cfg.n_long, 1);
    for j = 1:cfg.n_long
        blk = S(starts_floor(j):starts_floor(j) + cfg.T_cl - 1, :);
        D_fl(j) = cfg.metric_test(X_test, blk);
    end
    res.D_floor = mean(D_fl);

    starts = round(linspace(i0, cfg.i_test(end) - cfg.vpt_len, cfg.n_vpt));
    vpt = zeros(cfg.n_vpt, 1);
    thr_x = [];
    if isfield(cfg, 'vpt_thr_extra')
        thr_x = cfg.vpt_thr_extra(:)';
    end
    vpt_x = zeros(cfg.n_vpt, numel(thr_x));
    for j = 1:cfg.n_vpt
        s0 = starts(j);
        Pv = E.loop(mdl, S(s0 - nh + 1:s0, :), cfg.vpt_len, cfg.thr);
        Pv(isnan(Pv)) = 1e3;
        n = VPT(Pv, S(s0 + 1:s0 + cfg.vpt_len, :), sigma, cfg.vpt_thr);
        vpt(j) = n * cfg.dt * cfg.lle_true;
        for q = 1:numel(thr_x)
            vpt_x(j, q) = VPT(Pv, S(s0 + 1:s0 + cfg.vpt_len, :), sigma, thr_x(q)) * cfg.dt * cfg.lle_true;
        end
    end
    res.vpt = vpt;
    fprintf('short-term: VPT = %.2f +- %.2f Lyapunov times\n', mean(vpt), std(vpt));
    if ~isempty(thr_x)
        res.vpt_extra = vpt_x;
        res.vpt_thr_extra = thr_x;
        for q = 1:numel(thr_x)
            fprintf('short-term: VPT (NRMSE %.2f) = %.2f +- %.2f Lyapunov times\n', thr_x(q), mean(vpt_x(:, q)), std(vpt_x(:, q)));
        end
    end

    lle_all = nan(cfg.n_long, 1);
    for j = 1:cfg.n_long
        if div(j)
            continue
        end
        s0 = starts_long(j);
        rng(j);
        lam = E.lyapunov(mdl, S(s0 - nh + 1:s0, :), cfg.lle_trans, cfg.lle_eval);
        lle_all(j) = lam / cfg.dt;
        fprintf('  model LLE %2d/%d: %.4f\n', j, cfg.n_long, lle_all(j));
    end
    res.lle_all = lle_all;
    res.lle = mean(lle_all, 'omitnan');
    fprintf('model LLE = %.3f +- %.3f over %d starts (true %.3f)\n', ...
        res.lle, std(lle_all, 'omitnan'), sum(~isnan(lle_all)), cfg.lle_true);

    cfg_saved = cfg;
    for f = {'S', 'S_ref'}
        if isfield(cfg_saved, f{1})
            cfg_saved = rmfield(cfg_saved, f{1});
        end
    end
    res.cfg = cfg_saved;
    save(result_file, 'res');
end

function E = Engine(model)
    if strcmp(model, 'RC')
        E.train = @PRC_train;
        E.n_hist = @(mdl) mdl.n_warm;
        E.loop = @Closed_loop_PRC;
        E.lyapunov = @(mdl, H, n_trans, n_eval) Largest_lyapunov(@(R) PRC_step(mdl, R, PRC_readout(mdl, R)), ...
            PRC_synchronize(mdl, H), n_trans, n_eval);
    else
        E.train = @Train_NGRC;
        E.n_hist = @(mdl) mdl.k;
        E.loop = @Closed_loop_NGRC;
        E.lyapunov = @(mdl, H, n_trans, n_eval) Largest_lyapunov(@(Hs) Step_history(mdl, Hs), H, n_trans, n_eval);
    end
end

function f = Objective(p, S, cfg, E)
    if istable(p)
        p = table2struct(p);
    end
    p = Set_fixed(p, cfg.fixed);
    mdl = E.train(S(1:cfg.n_fit, :), p, cfg.m);
    nh = E.n_hist(mdl);
    f_all = zeros(numel(cfg.starts), 1);
    for j = 1:numel(cfg.starts)
        s0 = cfg.starts(j);
        [P, td] = E.loop(mdl, S(s0 - nh + 1:s0, :), cfg.T_cl, cfg.thr);
        if isnan(td)
            f_all(j) = cfg.metric(S(cfg.i_val, :), P);
        else
            f_all(j) = 1 + (1 - td / cfg.T_cl);
        end
    end
    f = mean(f_all);
end

function n = VPT(P, X, sigma, thr)
    err = sqrt(mean(((P - X) ./ sigma).^2, 2));
    n = find(err > thr, 1);
    if isempty(n)
        n = size(P, 1);
    else
        n = n - 1;
    end
end

function H = Step_history(mdl, H)
    y = Predict_NGRC(mdl, Delay_input(H, mdl.nn, mdl.m, mdl.k));
    H = [H(2:end, :); reshape(y, 1, [])];
end

function Check_settings(old, cfg, result_file)
    same = strcmp(old.model, cfg.model) && old.max_iter == cfg.max_iter ...
        && old.n_fit == cfg.n_fit && old.T_cl == cfg.T_cl && old.thr == cfg.thr ...
        && isequal(old.starts, cfg.starts) && isequal(old.fixed, cfg.fixed);
    if isfield(cfg, 'sig')
        same = same && isfield(old, 'sig') && old.sig == cfg.sig;
    end
    if ~same
        error('Saved result was made with different settings. Delete %s and rerun.', result_file);
    end
end
