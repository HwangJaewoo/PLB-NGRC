%% Training of parallel RC (Pathak et al. 2018): one reservoir per site, readout shared by all sites, ridge regression
function mdl = PRC_train(S, p, m)
    mdl = Build(p, m);
    [T, Nm] = size(S);
    N = Nm / m;
    Dr = mdl.Dr;
    G = zeros(Dr, Dr);
    B = zeros(Dr, m);
    nr = 0;
    R = zeros(N, Dr);
    cs = max(1, floor(5000 / N));
    Qc = zeros(cs * N, Dr);
    Yc = zeros(cs * N, m);
    c = 0;
    for t = 1:T - 1
        R = PRC_step(mdl, R, S(t, :));
        if t > mdl.n_warm
            Q = R;
            Q(:, 2:2:end) = R(:, 2:2:end).^2;
            Qc(c * N + (1:N), :) = Q;
            Yc(c * N + (1:N), :) = reshape(S(t + 1, :), N, m);
            c = c + 1;
            if c == cs || t == T - 1
                G = G + Qc(1:c * N, :)' * Qc(1:c * N, :);
                B = B + Qc(1:c * N, :)' * Yc(1:c * N, :);
                nr = nr + c * N;
                c = 0;
            end
        end
    end
    mdl.W = (G / nr + mdl.lam * eye(Dr)) \ (B / nr);
end

function mdl = Build(p, m)
    if istable(p)
        p = table2struct(p);
    end
    rng(p.seed, 'twister');
    Dr = p.Dr;
    d0 = (2 * p.nn + 1) * m;
    A = sprand(Dr, Dr, 3 / Dr);
    [ia, ja, va] = find(A);
    A = sparse(ia, ja, 2 * va - 1, Dr, Dr);
    rho0 = max(abs(eig(full(A))));
    A = A * (p.rho / rho0);
    mdl.A_t = A.';
    mdl.idx = mod((0:Dr - 1)', d0) + 1;
    mdl.w = p.sig_in * (2 * rand(1, Dr) - 1);
    mdl.leak = p.leak;
    mdl.nn = p.nn;
    mdl.m = m;
    mdl.Dr = Dr;
    mdl.n_warm = p.n_warm;
    mdl.lam = 10^p.log_lam;
end
