%% Kuramoto-Sivashinsky: space-time plots, mean power spectra, and VPT, W1 and largest Lyapunov exponent versus L (output figure9.png)
clear; clc; close all;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'core'), genpath(fullfile(root, 'experiments')));
img_dir = fullfile(root, 'figures', 'img');
if ~isfolder(img_dir)
    mkdir(img_dir);
end
K = load(fullfile(root, 'experiments', 'KS', 'results', 'KS_spectrum.mat'), 'LQ', 'spec', 'err_avgspec', 'n_surv');

L_st = 200;
T_LT = 5;
fig_w = 22;
fig_h = 19;
fs = 10;
x_tick = 0:2:8;
x_lim = [-0.2, 8];
col = [0, 0.2, 0.8; 0.85, 0, 0];
lab = {'True', 'P-NGRC', 'PLB-NGRC'};

nL = size(K.LQ, 1);
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, fig_w, fig_h]);
x_L = 0.085;
x_cb = 0.940;  w_cb = 0.015;
g_cb = 0.012;
g_t = 0.045;  g_b = 0.06;
y_t = 0.72;  h_t = 0.23;
y_b = 0.42;  h_b = 0.19;
y_m = 0.12;  h_m = 0.18;
g_m = 0.08;
w_t = (x_cb - g_cb - x_L - 2 * g_t) / 3;
w_b = (x_cb + w_cb - x_L - 3 * g_b) / 4;
pos_t = arrayfun(@(k) [x_L + (k - 1) * (w_t + g_t), y_t, w_t, h_t], 1:3, 'UniformOutput', false);
pos_b = arrayfun(@(k) [x_L + (k - 1) * (w_b + g_b), y_b, w_b, h_b], 1:4, 'UniformOutput', false);
w_m = (x_cb + w_cb - x_L - 2 * g_m) / 3;
pos_m = arrayfun(@(k) [x_L + (k - 1) * (w_m + g_m), y_m, w_m, h_m], 1:3, 'UniformOutput', false);

tag = sprintf('KS_L%03d', L_st);
D = load(fullfile(root, 'experiments', 'KS', 'data', [tag, '_data.mat']), 'S', 'mu', 'sd', 'L', 'Q', 'dt', 'i_test');
q = load(fullfile(root, 'experiments', 'KS', 'results', [tag, '_true_lyapunov.mat']), 'lle');
lam = mean(q.lle);
n_st = round(T_LT / (lam * D.dt));
t_lt = (1:n_st) * D.dt * lam;
x = (0:D.Q - 1) * D.L / D.Q;
U = {D.S(D.i_test(1:n_st), :), [], []};
meth = {'PNGRC', 'PLBNGRC'};
for i = 1:2
    R = load(fullfile(root, 'experiments', 'KS', 'results', [tag, '_', meth{i}, '.mat']), 'res');
    U{i + 1} = R.res.P_long(1:n_st, :);
end
for k = 1:3
    U{k} = U{k} .* D.sd + D.mu;
end
c_max = ceil(max(cellfun(@(Z) max(abs(Z(:))), U)));
st_lab = {'True', 'P-NGRC', 'PLB-NGRC'};
ax_st = gobjects(1, 3);
for k = 1:3
    ax_st(k) = axes(fig, 'Position', pos_t{k});
    imagesc(ax_st(k), t_lt, x, U{k}');
    set(ax_st(k), 'YDir', 'normal', 'FontSize', fs, 'Box', 'on');
    set(ax_st(k), 'CLim', [-c_max, c_max]);
    xlabel(ax_st(k), '$\lambda_{\max} t$', 'Interpreter', 'latex', 'FontSize', fs + 2);
    if k == 1
        ylabel(ax_st(k), '$x$', 'Interpreter', 'latex', 'FontSize', fs + 2);
    end
    xlim(ax_st(k), [0, T_LT]);
    ylim(ax_st(k), [0, D.L]);
    xticks(ax_st(k), 0:1:T_LT);
    yticks(ax_st(k), 0:D.L / 2:D.L);
    title(ax_st(k), st_lab{k}, 'FontWeight', 'normal', 'FontSize', fs);
    text(ax_st(k), 0, 1.02, sprintf('(%s)', char('a' + k - 1)), 'Units', 'normalized', 'FontWeight', 'bold', ...
        'FontSize', fs + 2, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end
colormap(fig, parula);
cb = colorbar(ax_st(3));
cb.Ticks = [-c_max, 0, c_max];
cb.Position = [x_cb, y_t, w_cb, h_t];
ax_st(3).Position = pos_t{3};
cb.Label.String = '$u$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = fs + 2;
cb.FontSize = fs;
ax_all = gobjects(1, nL);
h_leg = gobjects(1, 3);
y_rng = [inf, -inf];
for a = 1:nL
    sp = K.spec(a);
    Q = numel(sp.Pt);
    iq = 2:Q / 2;
    kk = sp.kk(iq);
    ax = axes(fig, 'Position', pos_b{a});
    ax_all(a) = ax;
    h_leg(1) = semilogy(ax, kk, sp.Pt(iq), 'k-', 'LineWidth', 1.2);
    hold(ax, 'on');
    y_rng = [min(y_rng(1), min(sp.Pt(iq))), max(y_rng(2), max(sp.Pt(iq)))];
    for i = 1:2
        h_leg(i + 1) = semilogy(ax, kk, sp.Pm(i, iq), '-', 'Color', col(i, :), 'LineWidth', 1);
        y_rng = [min(y_rng(1), min(sp.Pm(i, iq))), max(y_rng(2), max(sp.Pm(i, iq)))];
        fprintf('L = %3d %-8s: error of the averaged spectrum %.3f, %d predictions\n', ...
            sp.L, lab{i + 1}, K.err_avgspec(a, i), K.n_surv(a, i));
    end
    xlim(ax, x_lim);
    xticks(ax, x_tick);
    set(ax, 'FontSize', fs, 'Box', 'on', 'YMinorTick', 'off', 'XMinorTick', 'off');
    grid(ax, 'on');
    set(ax, 'YMinorGrid', 'off', 'XMinorGrid', 'off');
    title(ax, sprintf('$L = %d$', sp.L), 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', fs + 2);
    text(ax, 0, 1.02, sprintf('(%s)', char('a' + a + 2)), 'Units', 'normalized', 'FontWeight', 'bold', ...
        'FontSize', fs + 2, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    xlabel(ax, 'Wavenumber $q$', 'Interpreter', 'latex', 'FontSize', fs + 2);
    if a == 1
        ylabel(ax, 'Power', 'FontSize', fs);
    end
end
e_lo = floor(log10(y_rng(1)));
e_hi = ceil(log10(y_rng(2)));
st = ceil((e_hi - e_lo) / 3);
e_lo = e_hi - st * ceil((e_hi - e_lo) / st);
for a = 1:nL
    ylim(ax_all(a), 10.^[e_lo, e_hi]);
    yticks(ax_all(a), 10.^(e_lo:st:e_hi));
end
models = {'PNGRC', 'PLBNGRC'};
M = struct('vpt', zeros(nL, 2), 'vpt_s', zeros(nL, 2), 'w1', zeros(nL, 2), 'w1_s', zeros(nL, 2), ...
    'lle', zeros(nL, 2), 'lle_s', zeros(nL, 2), 'ref', zeros(nL, 1), 'true', zeros(nL, 1));
for a = 1:nL
    tg = sprintf('KS_L%03d', K.LQ(a, 1));
    for i = 1:2
        R = load(fullfile(root, 'experiments', 'KS', 'results', [tg, '_', models{i}, '.mat']), 'res');
        r = R.res;
        v = r.vpt_extra(:, find(r.vpt_thr_extra == 0.5, 1));
        M.vpt(a, i) = mean(v);                  M.vpt_s(a, i) = std(v);
        M.w1(a, i) = mean(r.D_all, 'omitnan');  M.w1_s(a, i) = std(r.D_all, 'omitnan');
        M.lle(a, i) = mean(r.lle_all, 'omitnan'); M.lle_s(a, i) = std(r.lle_all, 'omitnan');
        M.ref(a) = r.D_floor;
    end
    q = load(fullfile(root, 'experiments', 'KS', 'results', [tg, '_true_lyapunov.mat']), 'lle');
    M.true(a) = mean(q.lle);
    fprintf('L = %3d: VPT(0.5) %.2f / %.2f | W1 %.4f / %.4f (ref %.4f) | lambda %.3f / %.3f (true %.3f)\n', K.LQ(a, 1), ...
        M.vpt(a, :), M.w1(a, :), M.ref(a), M.lle(a, :), M.true(a));
end
Lname = arrayfun(@(x) sprintf('%d', x), K.LQ(:, 1)', 'UniformOutput', false);
ax_m = gobjects(1, 3);
ylab_m = {'VPT', '$\mathcal{W}_1$', '$\lambda_{\max}$'};
mm = {M.vpt, M.w1, M.lle};
ss = {M.vpt_s, M.w1_s, M.lle_s};
rr = {[], M.ref, M.true};
for k = 1:3
    ax_m(k) = axes(fig, 'Position', pos_m{k});
    hold(ax_m(k), 'on');
    b = bar(ax_m(k), 1:nL, mm{k}, 0.8, 'EdgeColor', 'none');
    for i = 1:2
        b(i).FaceColor = col(i, :);
        errorbar(ax_m(k), b(i).XEndPoints, mm{k}(:, i), ss{k}(:, i), 'k', 'LineStyle', 'none', 'LineWidth', 1.0, 'CapSize', 4);
    end
    top = max(mm{k}(:) + ss{k}(:));
    if ~isempty(rr{k})
        xe = [0.5, (1:nL - 1) + 0.5, nL + 0.5];
        plot(ax_m(k), reshape([xe(1:end - 1); xe(2:end)], 1, []), reshape([rr{k}'; rr{k}'], 1, []), '--', ...
            'Color', 'k', 'LineWidth', 1.2);
        top = max(top, max(rr{k}));
    end
    xlim(ax_m(k), [0.5, nL + 0.5]);
    set(ax_m(k), 'XTick', 1:nL, 'XTickLabel', Lname, 'FontSize', fs, 'Box', 'on');
    set(fig, 'CurrentAxes', ax_m(k));
    Tight_lim(0, top);
    grid(ax_m(k), 'on');
    xlabel(ax_m(k), '$L$', 'Interpreter', 'latex', 'FontSize', fs + 2);
    if k == 1
        ylabel(ax_m(k), ylab_m{k}, 'FontSize', fs);
    else
        ylabel(ax_m(k), ylab_m{k}, 'Interpreter', 'latex', 'FontSize', fs + 2);
    end
    text(ax_m(k), 0, 1.02, sprintf('(%s)', char('g' + k)), 'Units', 'normalized', 'FontWeight', 'bold', ...
        'FontSize', fs + 2, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end

lgd = legend(ax_all(1), h_leg, lab, 'Orientation', 'horizontal', 'FontSize', fs, 'Box', 'off');
lgd.Units = 'normalized';
lgd.Position(1:2) = [(x_L + x_cb + w_cb) / 2 - lgd.Position(3) / 2, 0.01];
for a = 1:nL
    ax_all(a).Position = pos_b{a};
end
for k = 1:3
    ax_m(k).Position = pos_m{k};
end
exportgraphics(fig, fullfile(img_dir, 'figure9.png'), 'Resolution', 600);

function Tight_lim(lo, hi)
    pad = 0.02 * (hi - lo);
    lo = lo - (lo < 0) * pad;
    hi = hi + (hi > 0) * pad;
    e = 10^floor(log10((hi - lo) / 4));
    best = [inf, 0, 0, 0];
    for c = [1, 2, 5, 10, 20, 50] * e
        y = [floor(lo / c), ceil(hi / c)] * c;
        n = round((y(2) - y(1)) / c);
        if n >= 2 && n <= 4 && (y(2) - y(1) < best(1) - 1e-12 || abs(y(2) - y(1) - best(1)) < 1e-12)
            best = [y(2) - y(1), y, c];
        end
    end
    ylim(best(2:3));
    yticks(best(2):best(4):best(3));
end
