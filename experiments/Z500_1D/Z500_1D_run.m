%% Z500 along 50N: P-NGRC and PLB-NGRC selected by BO on long-term statistics of the validation set, then evaluated on the test set
function Z500_1D_run()
    for method = {'PNGRC', 'PLBNGRC'}
        [cfg, vars, f] = Z500_1D_config(method{1});
        Z500_1D_experiment(cfg, vars, f);
    end
end
