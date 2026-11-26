%% Training of PLB-NGRC for the Z500 data over the latitude band: partition conditioned on latitude and season, weighted fitting
function mdl = Z500_2D_train(cfg, p)
    k = p.k;
    NB = cfg.NB;
    P = cfg.P;
    d1 = size(NB, 2) * k;
    fa = isfield(p, 'fa') && p.fa;
    na = 3 * fa;
    [I2, J2] = find(triu(ones(d1 + na)));
    D = 1 + d1 + na + numel(I2);
    bas = struct('fa', fa, 'I2', I2, 'J2', J2);
    J = 2^p.s;
    tt_all = (k:cfg.n_fit - 1)';
    if p.s == 0
        C = zeros(1, d1 + 3);
        sig = Inf;
    else
        rng(0);
        ts = tt_all(randi(numel(tt_all), ceil(cfg.n_sub / P), 1));
        Z = Coords(Z500_2D_local_states(cfg.Saug, ts, NB, k), ts, cfg, p);
        idx = Ball_tree(Z, p.s);
        C = zeros(J, d1 + 3);
        r = zeros(J, 1);
        for j = 1:J
            C(j, :) = mean(Z(idx == j, :), 1);
            r(j) = max(vecnorm(Z(idx == j, :) - C(j, :), 2, 2));
        end
        sig = p.eta * r;
    end
    G = zeros(D, D, J);
    B = zeros(1, D, J);
    sw = zeros(J, 1);
    for c0 = 1:cfg.chunk:numel(tt_all)
        tt = tt_all(c0:min(c0 + cfg.chunk - 1, numel(tt_all)));
        L = Z500_2D_local_states(cfg.Saug, tt, NB, k);
        y = reshape(cfg.S(tt + 1, :), [], 1);
        w = Z500_2D_blend_weights(Coords(L, tt, cfg, p), C, sig);
        for j = 1:J
            wj = w(:, j);
            wj(wj <= 1e-3) = 0;
            O = Z500_2D_features(L - C(j, 1:d1), bas, cfg, tt);
            G(:, :, j) = G(:, :, j) + O' * (O .* wj);
            B(:, :, j) = B(:, :, j) + (y .* wj)' * O;
            sw(j) = sw(j) + sum(wj);
        end
    end
    W = zeros(J, D);
    for j = 1:J
        W(j, :) = (B(:, :, j) / sw(j)) / (G(:, :, j) / sw(j) + p.alpha * eye(D));
    end
    mdl = struct('p', p, 'k', k, 'd1', d1, 'fa', fa, 'I2', I2, 'J2', J2, 'C', C, 'sig', sig, 'W', W);
end

function Z = Coords(L, tt, cfg, p)
    nt = numel(tt);
    Z = [L, p.bl * kron(cfg.latP, ones(nt, 1)), p.bs * repmat(cfg.ph(tt, :), cfg.P, 1)];
end

function idx = Ball_tree(Z, s)
    n = size(Z, 1);
    idx = ones(n, 1);
    for lev = 1:s
        new = zeros(n, 1);
        for g = 1:2^(lev - 1)
            m = find(idx == g);
            [~, a] = max(var(Z(m, :), 0, 1));
            [~, o] = sort(Z(m, a));
            h = floor(numel(m) / 2);
            new(m(o(1:h))) = 2 * g - 1;
            new(m(o(h + 1:end))) = 2 * g;
        end
        idx = new;
    end
end
