%% Bayesian optimization (rng(0)) with a checkpoint after every evaluation; an interrupted run resumes
function bo = Run_BO(fun, vars, n_total, ckpt, cfg_saved, extra)
    if nargin < 6
        extra = {};
    end
    cfg_file = [ckpt(1:end - 4), '_cfg.mat'];
    if isempty(gcp('nocreate'))
        parpool;
    end
    t0 = tic;
    if isfile(ckpt)
        L = load(ckpt, 'BayesoptResults');
        if numel(L.BayesoptResults.ObjectiveTrace) == 0
            delete(ckpt);
            fprintf('removed empty checkpoint %s\n', ckpt);
        end
    end
    if isfile(ckpt)
        if ~isfile(cfg_file)
            error('%s exists without %s. Delete the checkpoint to start over.', ckpt, cfg_file);
        end
        C = load(cfg_file, 'cfg_saved');
        if ~isequaln(Strip(C.cfg_saved), Strip(cfg_saved))
            error('Checkpoint %s was made with different settings. Delete it (and %s) to start over.', ckpt, cfg_file);
        end
        L = load(ckpt, 'BayesoptResults');
        bo0 = L.BayesoptResults;
        n_done = numel(bo0.ObjectiveTrace);
        fprintf('resuming BO from %s: %d / %d evaluations done (best so far %.4f)\n', ckpt, n_done, n_total, bo0.MinObjective);
        if n_done >= n_total
            bo = bo0;
            return
        end
        Monitor([], 'reset', t0, n_total, n_done);
        bo = resume(bo0, 'MaxObjectiveEvaluations', n_total - n_done, ...
            'OutputFcn', {@saveToFile, @(r, s) Monitor(r, s, t0, n_total, n_done)}, 'SaveFileName', ckpt);
    else
        save(cfg_file, 'cfg_saved');
        rng(0);
        Monitor([], 'reset', t0, n_total, 0);
        bo = bayesopt(fun, vars, ...
            'MaxObjectiveEvaluations', n_total, ...
            'IsObjectiveDeterministic', true, ...
            'AcquisitionFunctionName', 'expected-improvement-plus', ...
            'UseParallel', true, ...
            'Verbose', 0, 'PlotFcn', [], ...
            'OutputFcn', {@saveToFile, @(r, s) Monitor(r, s, t0, n_total, 0)}, ...
            'SaveFileName', ckpt, extra{:});
    end
end

function c = Strip(c)
    f = fieldnames(c);
    for i = 1:numel(f)
        if isa(c.(f{i}), 'function_handle')
            c = rmfield(c, f{i});
        end
    end
end

function stop = Monitor(results, state, t0, n_total, n_prev)
    persistent best
    stop = false;
    if ischar(state) && strcmp(state, 'reset')
        best = Inf;
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
        x = results.XAtMinObjective;
        nm = x.Properties.VariableNames;
        txt = cell(1, numel(nm));
        for i = 1:numel(nm)
            v = x.(nm{i});
            if iscategorical(v)
                v = double(v);
            end
            txt{i} = sprintf('%s = %.4g', nm{i}, v);
        end
        fprintf('    >> new best %.4f: %s\n', best, strjoin(txt, ', '));
    end
end
