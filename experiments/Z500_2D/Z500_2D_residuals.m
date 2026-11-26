%% One-step prediction errors on the training set (largest 1% excluded), used as added errors in long-term predictions
function R = Z500_2D_residuals(mdl, cfg, stride, clip_q)
    tt = (mdl.k:stride:cfg.n_fit - 1)';
    R = zeros(numel(tt), cfg.P);
    for c0 = 1:cfg.chunk:numel(tt)
        ic = c0:min(c0 + cfg.chunk - 1, numel(tt));
        R(ic, :) = cfg.S(tt(ic) + 1, :) - Predict(mdl, cfg, tt(ic));
    end
    rm = max(abs(R), [], 2);
    R = R(rm <= quantile(rm, clip_q), :);
end

function Y = Predict(mdl, cfg, tt)
    nt = numel(tt);
    L = Z500_2D_local_states(cfg.Saug, tt, cfg.NB, mdl.k);
    Z = [L, mdl.p.bl * kron(cfg.latP, ones(nt, 1)), mdl.p.bs * repmat(cfg.ph(tt, :), cfg.P, 1)];
    w = Z500_2D_blend_weights(Z, mdl.C, mdl.sig);
    y = zeros(size(L, 1), 1);
    for j = 1:size(mdl.C, 1)
        y = y + w(:, j) .* (Z500_2D_features(L - mdl.C(j, 1:mdl.d1), mdl, cfg, tt) * mdl.W(j, :)');
    end
    Y = reshape(y, nt, []);
end
