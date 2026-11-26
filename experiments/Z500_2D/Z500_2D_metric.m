%% Objective of the Z500 data over the latitude band: W1 + zonal spectrum error + autocorrelation error
function [f, parts] = Z500_2D_metric(X, Y, cfg)
    w1 = W1(X, Y);
    Px = Spectrum(X, cfg);
    Py = Spectrum(Y, cfg);
    Cx = Acf(X, cfg.max_lag);
    Cy = Acf(Y, cfg.max_lag);
    es = norm(Px - Py) / norm(Py);
    ea = norm(Cx - Cy) / norm(Cy);
    f = w1 + es + ea;
    parts = [w1, es, ea];
end

function Pw = Spectrum(X, cfg)
    T = size(X, 1);
    Z = reshape(X', cfg.Nlon, cfg.Nlat * T);
    F = fft(Z, [], 1);
    Pw = mean(abs(F(2:cfg.Nlon / 2 + 1, :)).^2, 2);
end

function c = Acf(X, max_lag)
    X = X - mean(X, 1);
    v = sum(X.^2, 1);
    c = zeros(max_lag, 1);
    for d = 1:max_lag
        c(d) = mean(sum(X(1 + d:end, :) .* X(1:end - d, :), 1) ./ v);
    end
end

function d = W1(a, b)
    p = linspace(0.0005, 0.9995, 2000);
    qa = quantile(a(:), p);
    qb = quantile(b(:), p);
    d = mean(abs(qa - qb)) / (quantile(b(:), 0.99) - quantile(b(:), 0.01));
end
