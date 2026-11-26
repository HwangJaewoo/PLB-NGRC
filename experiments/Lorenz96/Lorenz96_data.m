%% Lorenz-96 time series with N sites (training, validation and test sets) and the largest Lyapunov exponent of the true system
function Lorenz96_data(N)
    [data_file, lle_file] = Data_files('Lorenz96', N);
    if isfile(data_file) && isfile(lle_file)
        return
    end
    rng(0);
    F = 8;
    h = 0.001;
    dt = 0.01;
    T_warm = 100;
    T_fit = 100;
    T_val = 25;
    T_test = 50;
    n_starts = 10;

    n_sub = round(dt / h);
    n_warm = round(T_warm / dt);
    n_fit = round(T_fit / dt);
    n_val = round(T_val / dt);
    n_test = round(T_test / dt);
    n_keep = n_fit + n_val + n_test;

    x = F + 0.01 * randn(1, N);
    S_raw = zeros(n_keep, N);
    for t = 1:(n_warm + n_keep)
        x = Lorenz96_flow(x, h, n_sub, F);
        if t > n_warm
            S_raw(t - n_warm, :) = x;
        end
    end
    mu = mean(S_raw(1:n_fit, :), 'all');
    sd = std(S_raw(1:n_fit, :), 0, 'all');
    S = (S_raw - mu) / sd;
    i_val = (n_fit + 1):(n_fit + n_val);
    i_test = (n_fit + n_val + 1):n_keep;
    starts = round(linspace(i_val(1) - 1, i_val(end) - 1, n_starts));
    save(data_file, 'S', 'mu', 'sd', 'N', 'F', 'h', 'dt', 'n_fit', 'i_val', 'i_test', 'starts');

    rng(1);
    n_trans = round(10 / dt);
    n_eval = round(200 / dt);
    step_fun = @(x) Lorenz96_flow(x, h, n_sub, F);
    lle = zeros(numel(starts), 1);
    for j = 1:numel(starts)
        lle(j) = Largest_lyapunov(step_fun, S(starts(j), :) * sd + mu, n_trans, n_eval) / dt;
    end
    fprintf('Lorenz96 N = %d: true largest Lyapunov exponent %.3f +- %.3f\n', N, mean(lle), std(lle));
    save(lle_file, 'lle', 'n_trans', 'n_eval');
end
