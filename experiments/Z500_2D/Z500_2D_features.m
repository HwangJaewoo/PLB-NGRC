%% Polynomial features of the centered local states (latitude and seasonal phase appended only for the in-features variant)
function O = Z500_2D_features(Lc, mdl, cfg, tt)
    if isfield(mdl, 'fa') && mdl.fa
        nt = numel(tt);
        Lc = [Lc, kron(cfg.latP, ones(nt, 1)), repmat(cfg.ph(tt, :), cfg.P, 1)];
    end
    O = [ones(size(Lc, 1), 1), Lc, Lc(:, mdl.I2) .* Lc(:, mdl.J2)];
end
