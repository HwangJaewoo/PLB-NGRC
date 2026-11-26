%% Coupled Ikeda map time series with N sites (training, validation and test sets) and the largest Lyapunov exponent of the true map
function Ikeda_data(N)
    [data_file, lle_file] = Data_files('Ikeda', N);
    if isfile(data_file) && isfile(lle_file)
        return
    end
    rng(0);
    u = 0.9;
    eps_c = 0.3;
    T_warm = 2000;
    n_fit = 3000;
    n_val = 750;
    n_test = 2000;
    n_starts = 10;
    n_keep = n_fit + n_val + n_test;

    z0 = 0.1 * (randn(N, 1) + 1i * randn(N, 1));
    x = [real(z0).', imag(z0).'];
    S_raw = zeros(n_keep, 2 * N);
    for t = 1:(T_warm + n_keep)
        x = Ikeda_map(x, N, u, eps_c);
        if t > T_warm
            S_raw(t - T_warm, :) = x;
        end
    end
    cx = 1:N;
    cy = N + (1:N);
    mu = [mean(S_raw(1:n_fit, cx), 'all'), mean(S_raw(1:n_fit, cy), 'all')];
    sd = [std(S_raw(1:n_fit, cx), 0, 'all'), std(S_raw(1:n_fit, cy), 0, 'all')];
    S = S_raw;
    S(:, cx) = (S_raw(:, cx) - mu(1)) / sd(1);
    S(:, cy) = (S_raw(:, cy) - mu(2)) / sd(2);
    i_val = (n_fit + 1):(n_fit + n_val);
    i_test = (n_fit + n_val + 1):n_keep;
    starts = round(linspace(i_val(1) - 1, i_val(end) - 1, n_starts));
    save(data_file, 'S', 'S_raw', 'mu', 'sd', 'N', 'u', 'eps_c', 'n_fit', 'i_val', 'i_test', 'starts');

    rng(1);
    n_trans = 50;
    n_eval = 450;
    step_fun = @(x) Ikeda_map(x, N, u, eps_c);
    lle = zeros(numel(starts), 1);
    for j = 1:numel(starts)
        lle(j) = Largest_lyapunov(step_fun, S_raw(starts(j), :), n_trans, n_eval);
    end
    fprintf('Ikeda N = %d: true largest Lyapunov exponent %.3f +- %.3f\n', N, mean(lle), std(lle));
    save(lle_file, 'lle', 'n_trans', 'n_eval');
end
