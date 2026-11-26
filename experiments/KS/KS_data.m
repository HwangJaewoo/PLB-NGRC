%% Kuramoto-Sivashinsky time series for domain length L with Q grid points and the largest Lyapunov exponent of the true system
function [data_file, lle_file] = KS_data(L, Q)
    d = fullfile(Repo_root, 'experiments', 'KS');
    for f = {'data', 'results'}
        if ~isfolder(fullfile(d, f{1}))
            mkdir(fullfile(d, f{1}));
        end
    end
    data_file = fullfile(d, 'data', sprintf('KS_L%03d_data.mat', round(L)));
    lle_file = fullfile(d, 'results', sprintf('KS_L%03d_true_lyapunov.mat', round(L)));
    if isfile(data_file) && isfile(lle_file)
        return
    end
    rng(0);
    dt = 0.25;
    n_sub = 1;
    T_warm = 10000;
    n_fit = 10000;
    n_val = 2500;
    n_test = 5000;
    n_starts = 10;

    P = Setup(L, Q, dt / n_sub);
    n_keep = n_fit + n_val + n_test;
    u = 0.1 * randn(1, Q);
    u = u - mean(u);
    u = KS_flow(u, P, round(T_warm / dt) * n_sub);
    S_raw = zeros(n_keep, Q);
    for t = 1:n_keep
        u = KS_flow(u, P, n_sub);
        S_raw(t, :) = u;
    end
    if any(~isfinite(S_raw(:)))
        error('KS integration blew up');
    end
    mu = mean(S_raw(1:n_fit, :), 'all');
    sd = std(S_raw(1:n_fit, :), 0, 'all');
    S = (S_raw - mu) / sd;
    N = Q;
    i_val = (n_fit + 1):(n_fit + n_val);
    i_test = (n_fit + n_val + 1):n_keep;
    starts = round(linspace(i_val(1) - 1, i_val(end) - 1, n_starts));
    save(data_file, 'S', 'mu', 'sd', 'N', 'L', 'Q', 'dt', 'n_sub', 'n_fit', 'i_val', 'i_test', 'starts');

    rng(1);
    n_trans = round(100 / dt);
    n_eval = round(2000 / dt);
    lle = zeros(numel(starts), 1);
    parfor j = 1:numel(starts)
        lle(j) = Largest_lyapunov(@(z) KS_flow(z, P, n_sub), S(starts(j), :) * sd + mu, n_trans, n_eval) / dt;
    end
    fprintf('KS L = %d: true largest Lyapunov exponent %.4f +- %.4f\n', L, mean(lle), std(lle));
    save(lle_file, 'lle', 'n_trans', 'n_eval');
end

function P = Setup(L, Q, h)
    k = [0:Q / 2 - 1, 0, -Q / 2 + 1:-1]' * (2 * pi / L);
    Lk = k.^2 - k.^4;
    E = exp(h * Lk);
    E2 = exp(h * Lk / 2);
    M = 32;
    r = exp(1i * pi * ((1:M) - 0.5) / M);
    LR = h * Lk(:, ones(M, 1)) + r(ones(Q, 1), :);
    P.Qc = h * real(mean((exp(LR / 2) - 1) ./ LR, 2));
    P.f1 = h * real(mean((-4 - LR + exp(LR) .* (4 - 3 * LR + LR.^2)) ./ LR.^3, 2));
    P.f2 = h * real(mean((2 + LR + exp(LR) .* (-2 + LR)) ./ LR.^3, 2));
    P.f3 = h * real(mean((-4 - 3 * LR - LR.^2 + exp(LR) .* (4 - LR)) ./ LR.^3, 2));
    P.E = E;
    P.E2 = E2;
    P.g = -0.5i * k;
    P.L = L;
    P.Q = Q;
    P.h = h;
end
