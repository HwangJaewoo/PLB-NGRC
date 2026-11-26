%% Z500 over the latitude band: PLB-NGRC and P-NGRC, and the ablation of the conditioning with the other hyperparameters fixed and selected again
function Z500_2D_run()
    d = fullfile(Repo_root, 'experiments', 'Z500_2D', 'results');
    for variant = {'PLBNGRC', 'PNGRC'}
        Z500_2D_BO(variant{1});
        Z500_2D_test(variant{1});
    end

    B = load(fullfile(d, 'Z500_2D_PLBNGRC_best.mat'), 'p', 'n_done');
    for variant = {'season_only', 'latitude_only', 'neither', 'in_features'}
        p = B.p;
        p.bl = p.bl * strcmp(variant{1}, 'latitude_only');
        p.bs = p.bs * strcmp(variant{1}, 'season_only');
        p.fa = strcmp(variant{1}, 'in_features');
        obj = NaN;
        n_done = B.n_done;
        name = ['fixed_', variant{1}];
        save(fullfile(d, ['Z500_2D_', name, '_best.mat']), 'p', 'obj', 'n_done');
        Z500_2D_test(name);
    end

    for variant = {'season_only', 'latitude_only', 'neither', 'in_features'}
        Z500_2D_BO(variant{1});
        Z500_2D_test(variant{1});
    end
end
