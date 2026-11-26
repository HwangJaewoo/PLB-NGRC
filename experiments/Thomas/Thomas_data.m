%% Extended Thomas system time series with N sites (training, validation and test sets) and the largest Lyapunov exponent of the true system
function Thomas_data(N)
    [data_file, lle_file] = Data_files('Thomas', N);
    if isfile(data_file) && isfile(lle_file)
        return
    end
    rng(0);
    b = 0.10;
    h = 0.01;
    dt = 0.5;
    T_warm = 2000;
    n_fit = 10000;
    n_val = 2500;
    n_test = 5000;
    n_starts = 10;

    n_sub = round(dt / h);
    n_keep = n_fit + n_val + n_test;
    x = 4 * randn(1, N);
    x = Thomas_flow(x, b, h, round(T_warm / h));
    S_raw = zeros(n_keep, N);
    for t = 1:n_keep
        x = Thomas_flow(x, b, h, n_sub);
        S_raw(t, :) = x;
    end
    mu = mean(S_raw(1:n_fit, :), 'all');
    sd = std(S_raw(1:n_fit, :), 0, 'all');
    S = (S_raw - mu) / sd;
    i_val = (n_fit + 1):(n_fit + n_val);
    i_test = (n_fit + n_val + 1):n_keep;
    starts = round(linspace(i_val(1) - 1, i_val(end) - 1, n_starts));
    save(data_file, 'S', 'mu', 'sd', 'N', 'b', 'h', 'dt', 'n_fit', 'i_val', 'i_test', 'starts');

    rng(1);
    n_trans = round(50 / dt);
    n_eval = round(1000 / dt);
    lle = zeros(numel(starts), 1);
    parfor j = 1:numel(starts)
        lle(j) = Largest_lyapunov(@(z) Thomas_flow(z, b, h, n_sub), S(starts(j), :) * sd + mu, n_trans, n_eval) / dt;
    end
    fprintf('Thomas N = %d: true largest Lyapunov exponent %.4f +- %.4f\n', N, mean(lle), std(lle));
    save(lle_file, 'lle', 'n_trans', 'n_eval');
end
