%% BO of one variant (PLBNGRC, PNGRC, season_only, latitude_only, neither, in_features) on the long-term statistics of the validation set; resumes from its checkpoint
function Z500_2D_BO(variant)
    n_eval = 100;
    cfg = Z500_2D_setup();
    cfg.n_seed = 3;
    cfg.T_stab = 3650;
    cfg.res_stride = 1;
    cfg.clip_q = 0.99;
    cfg.obj_div = 3;
    cfg.variant = variant;
    d = fullfile(Repo_root, 'experiments', 'Z500_2D', 'results');
    if ~isfolder(d)
        mkdir(d);
    end
    bo_file = fullfile(d, ['Z500_2D_', variant, '_bo.mat']);
    cfg.best_file = fullfile(d, ['Z500_2D_', variant, '_best.mat']);
    diary(fullfile(d, ['Z500_2D_', variant, '_log.txt']));
    c_diary = onCleanup(@() diary('off'));
    fprintf('==== %s | Z500 2D %s ====\n', datestr(now, 'yyyy-mm-dd HH:MM'), variant);

    if strcmp(variant, 'PNGRC')
        vars = [optimizableVariable('la', [-10, -1]), optimizableVariable('k', [1, 4], 'Type', 'integer')];
    else
        vars = [optimizableVariable('la', [-10, -1]), optimizableVariable('eta', [0.1, 0.3]), ...
            optimizableVariable('s', [1, 4], 'Type', 'integer'), optimizableVariable('k', [1, 4], 'Type', 'integer')];
    end
    switch variant
        case 'latitude_only'
            vars = [vars, optimizableVariable('blat', [0, 20])];
        case 'season_only'
            vars = [vars, optimizableVariable('bsea', [0, 20])];
        case 'PLBNGRC'
            vars = [vars, optimizableVariable('blat', [0, 20]), optimizableVariable('bsea', [0, 20])];
    end

    t0 = tic;
    if isfile(bo_file)
        Lb = load(bo_file, 'BayesoptResults');
        if numel(Lb.BayesoptResults.ObjectiveTrace) == 0
            delete(bo_file);
        end
    end
    if isfile(bo_file)
        n_done = numel(Lb.BayesoptResults.ObjectiveTrace);
        fprintf('resuming %s: %d / %d evaluations done\n', bo_file, n_done, n_eval);
        if n_done >= n_eval
            BayesoptResults = Lb.BayesoptResults;
        else
            Monitor([], 'reset', cfg, t0, n_eval, n_done, Lb.BayesoptResults.MinObjective);
            BayesoptResults = resume(Lb.BayesoptResults, 'MaxObjectiveEvaluations', n_eval - n_done, ...
                'OutputFcn', {@saveToFile, @(r, st) Monitor(r, st, cfg, t0, n_eval, n_done)}, 'SaveFileName', bo_file);
        end
    else
        rng(0);
        Monitor([], 'reset', cfg, t0, n_eval, 0, Inf);
        BayesoptResults = bayesopt(@(x) Objective(x, cfg), vars, 'MaxObjectiveEvaluations', n_eval, ...
            'IsObjectiveDeterministic', true, 'AcquisitionFunctionName', 'expected-improvement-plus', ...
            'PlotFcn', [], 'Verbose', 0, 'OutputFcn', {@saveToFile, @(r, st) Monitor(r, st, cfg, t0, n_eval, 0)}, ...
            'SaveFileName', bo_file);
    end

    n_use = min(n_eval, numel(BayesoptResults.ObjectiveTrace));
    [obj, i_best] = min(BayesoptResults.ObjectiveTrace(1:n_use));
    p = Make_p(BayesoptResults.XTrace(i_best, :), variant);
    n_done = i_best;
    save(cfg.best_file, 'p', 'obj', 'n_done');
    fprintf('[%s] best: s = %d, k = %d, alpha = %.2e, eta = %.2f, beta_lat = %.1f, beta_sea = %.1f, obj = %.4f\n', ...
        variant, p.s, p.k, p.alpha, p.eta, p.bl, p.bs, obj);
end

function p = Make_p(x, variant)
    bl = 0;
    bs = 0;
    if ismember('blat', x.Properties.VariableNames)
        bl = x.blat;
    end
    if ismember('bsea', x.Properties.VariableNames)
        bs = x.bsea;
    end
    if strcmp(variant, 'PNGRC')
        s = 0;
        eta = 1;
    else
        s = x.s;
        eta = x.eta;
    end
    p = struct('s', s, 'k', x.k, 'alpha', 10^x.la, 'eta', eta, 'bl', bl, 'bs', bs, 'fa', strcmp(variant, 'in_features'));
end

function f = Objective(x, cfg)
    p = Make_p(x, cfg.variant);
    mdl = Z500_2D_train(cfg, p);
    R = Z500_2D_residuals(mdl, cfg, cfg.res_stride, cfg.clip_q);
    Xv = cfg.S(cfg.i_val, :);
    nv = numel(cfg.i_val);
    T = max(nv, cfg.T_stab);
    f = 0;
    for sd = 1:cfg.n_seed
        [Xl, td] = Z500_2D_closed_loop(mdl, cfg, R, cfg.i_val(1) - 1, T, sd, 1);
        if ~isnan(td)
            f = cfg.obj_div;
            return
        end
        f = f + Z500_2D_metric(Xl(1:nv, :), Xv, cfg);
    end
    f = f / cfg.n_seed;
end

function stop = Monitor(results, state, cfg, t0, n_total, n_prev, best0)
    persistent best
    stop = false;
    if ischar(state) && strcmp(state, 'reset')
        best = best0;
        return
    end
    if isempty(best)
        best = Inf;
    end
    if ~strcmp(state, 'iteration')
        return
    end
    n_now = numel(results.ObjectiveTrace);
    if n_now < n_prev
        n_now = n_prev + n_now;
    end
    el = toc(t0) / 60;
    eta = el / max(n_now - n_prev, 1) * (n_total - n_now);
    fprintf('[BO %3d/%3d] f = %8.4f | best = %8.4f | elapsed %6.1f min | ETA %6.1f min | %s\n', ...
        n_now, n_total, results.ObjectiveTrace(end), results.MinObjective, el, eta, datestr(now, 'HH:MM'));
    if results.MinObjective < best
        best = results.MinObjective;
        p = Make_p(results.XAtMinObjective, cfg.variant);
        obj = best;
        n_done = n_now;
        save(cfg.best_file, 'p', 'obj', 'n_done');
    end
end
