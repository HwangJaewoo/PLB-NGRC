%% Test-set evaluation of one variant: correlation of 15-day forecasts and long-term statistics of ten-year predictions with added errors
function Z500_2D_test(variant)
    n_seed = 10;
    clip_q = 0.99;
    fc_len = 15;
    fc_step = 10;
    thr = 0.6;
    d = fullfile(Repo_root, 'experiments', 'Z500_2D', 'results');
    out_file = fullfile(d, ['Z500_2D_', variant, '_test.mat']);
    if isfile(out_file)
        fprintf('%s exists, skipped\n', out_file);
        return
    end
    cfg = Z500_2D_setup();
    load(fullfile(d, ['Z500_2D_', variant, '_best.mat']), 'p', 'obj');
    if ~isfield(p, 'fa')
        p.fa = false;
    end
    fprintf('[%s] s = %d, k = %d, alpha = %.2e, eta = %.2f, beta_lat = %.1f, beta_sea = %.1f, fa = %d (val obj %.4f)\n', ...
        variant, p.s, p.k, p.alpha, p.eta, p.bl, p.bs, p.fa, obj);
    mdl = Z500_2D_train(cfg, p);
    R = Z500_2D_residuals(mdl, cfg, 1, clip_q);
    fprintf('one-step error variance / data variance = %.3f\n', var(R(:)) / var(cfg.S(1:cfg.n_fit, :), 0, 'all'));

    i_test = cfg.i_test;
    S = cfg.S;
    starts = i_test(1):fc_step:(i_test(end) - fc_len);
    acc = zeros(numel(starts), fc_len);
    acc_p = zeros(numel(starts), fc_len);
    for j = 1:numel(starts)
        s0 = starts(j);
        Xf = Z500_2D_closed_loop(mdl, cfg, R, s0, fc_len, 0, 0);
        for l = 1:fc_len
            acc(j, l) = Corr(Xf(l, :), S(s0 + l, :));
            acc_p(j, l) = Corr(S(s0, :), S(s0 + l, :));
        end
    end
    acc_m = mean(acc, 1);
    accp_m = mean(acc_p, 1);
    fprintf('correlation after 1/2/3/5 days: %.2f / %.2f / %.2f / %.2f (persistence %.2f / %.2f / %.2f / %.2f)\n', ...
        acc_m([1, 2, 3, 5]), accp_m([1, 2, 3, 5]));
    fprintf('horizon (correlation %.1f): %.2f days (persistence %.2f days)\n', thr, Horizon(acc_m, thr), Horizon(accp_m, thr));

    n_half = floor(numel(i_test) / 2);
    X1 = S(i_test(1:n_half), :);
    X2 = S(i_test(n_half + 1:2 * n_half), :);
    [~, ref] = Z500_2D_metric(X2, X1, cfg);
    parts = nan(n_seed, 3);
    sdr = nan(n_seed, 1);
    td_all = nan(n_seed, 1);
    X_first = [];
    for sd = 1:n_seed
        [Xl, td] = Z500_2D_closed_loop(mdl, cfg, R, i_test(1) - 1, numel(i_test), sd, 1);
        td_all(sd) = td;
        if isnan(td)
            [~, pa] = Z500_2D_metric(Xl(1:n_half, :), X1, cfg);
            [~, pb] = Z500_2D_metric(Xl(n_half + 1:2 * n_half, :), X1, cfg);
            parts(sd, :) = 0.5 * (pa + pb);
            sdr(sd) = std(Xl(:)) / std(X1(:));
            if isempty(X_first)
                X_first = Xl;
            end
        end
    end
    fprintf('diverged: %d / %d\n', sum(~isnan(td_all)), n_seed);
    fprintf('W1 %.4f +- %.4f (reference %.4f), spectrum error %.4f +- %.4f (%.4f), autocorrelation error %.4f +- %.4f (%.4f), std ratio %.3f\n', ...
        mean(parts(:, 1), 'omitnan'), std(parts(:, 1), 'omitnan'), ref(1), mean(parts(:, 2), 'omitnan'), std(parts(:, 2), 'omitnan'), ref(2), ...
        mean(parts(:, 3), 'omitnan'), std(parts(:, 3), 'omitnan'), ref(3), mean(sdr, 'omitnan'));
    save(out_file, 'p', 'acc_m', 'accp_m', 'parts', 'ref', 'sdr', 'td_all', 'X_first', 'starts');
end

function c = Corr(xp, xt)
    a = xp - mean(xp);
    b = xt - mean(xt);
    c = sum(a .* b) / sqrt(sum(a.^2) * sum(b.^2));
end

function h = Horizon(acc, thr)
    a = [1, acc];
    i = find(a < thr, 1);
    if isempty(i)
        h = numel(acc);
    else
        h = (i - 2) + (a(i - 1) - thr) / (a(i - 1) - a(i));
    end
end
