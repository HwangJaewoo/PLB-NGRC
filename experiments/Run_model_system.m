%% All experiments on one model system: main comparison, ablation, observational noise and scaling with N
function Run_model_system(sys)
    switch sys
        case 'Lorenz96'
            N_list = [20, 40, 100, 200];
            noise_list = [];
        otherwise
            N_list = [10, 20, 50, 100];
            noise_list = [0.01, 0.05, 0.10];
    end
    N0 = N_list(1);

    for method = {'PNGRC', 'PLBNGRC', 'PRC', 'PLBNGRC_inball', 'LBNGRC', 'LBNGRC_weighted'}
        [cfg, vars, f] = Model_system_config(sys, method{1}, N0);
        cfg.skip_done = true;
        Run_experiment(cfg, vars, f);
    end

    for noise = noise_list
        for method = {'PNGRC', 'PLBNGRC'}
            [cfg, vars, f] = Model_system_config(sys, method{1}, N0, noise);
            cfg.skip_done = true;
            Run_experiment(cfg, vars, f);
        end
    end

    for method = {'PNGRC', 'PLBNGRC', 'PRC'}
        [~, ~, f0] = Model_system_config(sys, method{1}, N0);
        R0 = load(f0, 'res');
        b0 = R0.res.best;
        for N = N_list(2:end)
            [cfg, vars, f] = Model_system_config(sys, method{1}, N);
            cfg.skip_done = true;
            names = arrayfun(@(v) v.Name, vars, 'UniformOutput', false);
            if strcmp(method{1}, 'PRC')
                cfg.fixed = struct('Dr', b0.Dr, 'seed', b0.seed, 'n_warm', b0.n_warm, 'nn', b0.nn);
                vars = vars(~strcmp(names, 'nn'));
            else
                cfg.fixed = struct('k', b0.k);
                if isfield(b0, 'ord')
                    cfg.fixed.ord = b0.ord;
                end
                cfg.fixed.s = b0.s;
                cfg.fixed.nn = b0.nn;
                vars = vars(ismember(names, {'log_lam', 'eta'}));
            end
            Run_experiment(cfg, vars, f);
        end
    end
end
