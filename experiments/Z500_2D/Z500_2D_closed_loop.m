%% Autonomous prediction of T days from day s0 with randomly drawn one-step errors (gamma = 1) or without (gamma = 0)
function [X, td] = Z500_2D_closed_loop(mdl, cfg, R, s0, T, seed, gamma)
    rng(seed);
    H = cfg.S(s0 - mdl.k + 1:s0, :);
    X = nan(T, cfg.P);
    td = NaN;
    for t = 1:T
        y = Step(mdl, cfg, H, s0 + t - 1);
        if gamma > 0
            y = y + gamma * R(randi(size(R, 1)), :);
        end
        if any(~isfinite(y)) || max(abs(y)) > cfg.div_tol
            td = t;
            return
        end
        X(t, :) = y;
        H = [H(2:end, :); y];
    end
end

function y = Step(mdl, cfg, H, t)
    k = size(H, 1);
    NB = cfg.NB;
    Ha = [H, zeros(k, 1)];
    nnb = size(NB, 2);
    L = zeros(cfg.P, nnb * k);
    c = 0;
    for q = 0:k - 1
        row = Ha(end - q, :);
        for n = 1:nnb
            c = c + 1;
            L(:, c) = row(NB(:, n))';
        end
    end
    Z = [L, mdl.p.bl * cfg.latP, mdl.p.bs * repmat(cfg.ph(t, :), cfg.P, 1)];
    w = Z500_2D_blend_weights(Z, mdl.C, mdl.sig);
    y = zeros(cfg.P, 1);
    for j = 1:size(mdl.C, 1)
        y = y + w(:, j) .* (Z500_2D_features(L - mdl.C(j, 1:mdl.d1), mdl, cfg, t) * mdl.W(j, :)');
    end
    y = y';
end
