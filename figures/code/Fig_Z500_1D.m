%% Z500 along 50N: map, forecast correlation, one forecast, and long-term statistics of P-NGRC and PLB-NGRC (output figure10.png)
clear; clc; close all;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'core'), genpath(fullfile(root, 'experiments')));
img_dir = fullfile(root, 'figures', 'img');
if ~isfolder(img_dir)
    mkdir(img_dir);
end
z_dir = fullfile(root, 'experiments', 'Z500_1D');

D = load(fullfile(z_dir, 'data', 'Z500_1D_data.mat'), 'S', 'lon36');
S = D.S;
lon36 = double(D.lon36);
A = load(fullfile(z_dir, 'results', 'Z500_1D_PNGRC.mat'), 'res');
B = load(fullfile(z_dir, 'results', 'Z500_1D_PLBNGRC.mat'), 'res');
r = {A.res, B.res};
cfg = B.res.cfg;
N = size(S, 2);

fig_w = 22;
fig_h = 15;
fs = 11;
fs_lab = 12;
lw = 1.2;
ms = 4;
leads = [1, 2, 3];
lat0 = 50;
lat_min = 20;
thr = 0.6;
col = [0, 0.2, 0.8; 0.85, 0, 0];
c_thr = [0, 0.6, 0];

s0 = cfg.i_test(1);
i0 = cfg.i_test(1) - 1;
n_half = floor(numel(cfg.i_test) / 2);
X1 = S(cfg.i_test(1:n_half), :);
X2 = S(cfg.i_test(n_half + 1:2 * n_half), :);
Pf = cell(1, 2);
Xr = cell(1, 2);
for i = 1:2
    mdl = Train_NGRC(S(1:cfg.n_fit, :), r{i}.best, cfg.m);
    Pf{i} = Z500_closed_loop(mdl, S(s0 - mdl.k + 1:s0, :), max(leads), cfg.thr, [], 0);
    R = Z500_residuals(mdl, S(1:cfg.n_fit, :));
    Xr{i} = cell(1, cfg.n_real);
    for q = 1:cfg.n_real
        rng(100 + q, 'twister');
        [Pq, td] = Z500_closed_loop(mdl, S(i0 - mdl.k + 1:i0, :), 2 * n_half, cfg.thr, R, cfg.gamma);
        if isnan(td)
            Xr{i}{q} = Pq;
        end
    end
    Xr{i} = Xr{i}(~cellfun(@isempty, Xr{i}));
end
names = {'P-NGRC', 'PLB-NGRC'};
for i = 1:2
    fprintf('%-8s: ACC 1 d %.2f (persistence %.2f), horizon %.2f d (persistence %.2f d), long-term runs %d / %d\n', ...
        names{i}, r{i}.acc(1), r{i}.acc_p(1), r{i}.hor, r{i}.hor_p, numel(Xr{i}), cfg.n_real);
    for l = leads
        fprintf('   first day of the test set, lead %d d: correlation %.2f\n', l, Corr(Pf{i}(l, :), S(s0 + l, :)));
    end
end

f = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, fig_w, fig_h]);

ax_a = axes(f);
Draw_map(ax_a, lat0, lon36, lat_min, fs);
Label('(a)');

ax_b = axes(f);
hold on;
L = numel(r{2}.acc);
yline(thr, '--', 'Color', c_thr, 'LineWidth', lw, 'Alpha', 1, 'HandleVisibility', 'off');
h_p = plot(1:L, r{2}.acc_p, '--s', 'Color', 'k', 'LineWidth', lw, 'MarkerSize', ms);
h_2 = plot(1:L, r{2}.acc, '-o', 'Color', col(2, :), 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', col(2, :));
h_1 = plot(1:L, r{1}.acc, '--o', 'Color', col(1, :), 'LineWidth', lw, 'MarkerSize', ms);
xlim([0, 15]);
xticks(0:5:15);
Corr_lim(min([r{1}.acc, r{2}.acc, r{2}.acc_p]));
xlabel('Lead time (days)');
ylabel('Correlation');
leg = legend([h_2, h_1, h_p], {'PLB-NGRC', 'P-NGRC', 'Persistence'}, 'Box', 'off');
Label('(b)');
grid on;
box on;

lon_p = [lon36(:); 360];
v = [reshape(S(s0 + leads, :), 1, []), reshape(Pf{1}(leads, :), 1, []), reshape(Pf{2}(leads, :), 1, [])];
ym = ceil(max(abs(v)));
labels_p = {'(c)', '(d)', '(e)'};
ax_p = gobjects(1, numel(leads));
for q = 1:numel(leads)
    l = leads(q);
    ax_p(q) = axes(f);
    hold on;
    yline(0, '-', 'Color', [0.7, 0.7, 0.7], 'Alpha', 1);
    xt = S(s0 + l, :);
    plot(lon_p, [xt, xt(1)], '-', 'Color', 'k', 'LineWidth', lw);
    sty = {'--', '-'};
    for i = [2, 1]
        xp = Pf{i}(l, :);
        plot(lon_p, [xp, xp(1)], sty{i}, 'Color', col(i, :), 'LineWidth', lw);
    end
    xlim([0, 360]);
    ylim([-ym, ym]);
    yticks([-ym, 0, ym]);
    Lon_ticks();
    xlabel('Longitude');
    if q == 1
        ylabel('Z500 anomaly');
    end
    Label(labels_p{q});
    grid on;
    box on;
end

edges = -4.5:0.25:4.5;
xc = edges(1:end - 1) + diff(edges) / 2;
ax_s = gobjects(1, 3);
ax_s(1) = axes(f);
hold on;
pk = 0;
for i = 1:2
    Y = zeros(numel(Xr{i}), numel(xc));
    for q = 1:numel(Xr{i})
        Y(q, :) = histcounts(Xr{i}{q}(:), edges, 'Normalization', 'pdf');
    end
    Band(xc, Y, col(i, :), false);
    pk = max(pk, max(mean(Y, 1) + std(Y, 0, 1)));
end
p1 = histcounts(X1(:), edges, 'Normalization', 'pdf');
p2 = histcounts(X2(:), edges, 'Normalization', 'pdf');
plot(xc, p1, '-', 'Color', 'k', 'LineWidth', lw);
plot(xc, p2, '--', 'Color', 'k', 'LineWidth', lw);
Mean_lines(xc, Xr, col, lw, @(X) histcounts(X(:), edges, 'Normalization', 'pdf'));
xlim([-4, 4]);
xticks(-4:2:4);
Tight_lim(0, max([pk, p1, p2]));
xlabel('$X_i(t)$', 'Interpreter', 'latex');
ylabel('Probability density');
Label('(f)');
grid on;
box on;

kk = 1:N / 2;
ax_s(2) = axes(f);
hold on;
for i = 1:2
    Band(kk, Stack(Xr{i}, @(X) Spectrum(X, N)), col(i, :), true);
end
s1 = Spectrum(X1, N);
s2 = Spectrum(X2, N);
plot(kk, s1, '-', 'Color', 'k', 'LineWidth', lw);
plot(kk, s2, '--', 'Color', 'k', 'LineWidth', lw);
Mean_lines(kk, Xr, col, lw, @(X) Spectrum(X, N));
set(gca, 'YScale', 'log', 'YMinorGrid', 'off', 'YMinorTick', 'off');
xlim([0, 20]);
xticks(0:5:20);
allv = [s1, s2, mean(Stack(Xr{1}, @(X) Spectrum(X, N)), 1), mean(Stack(Xr{2}, @(X) Spectrum(X, N)), 1)];
ylim(10.^[floor(log10(min(allv))), ceil(log10(max(allv)))]);
yticks(10.^(floor(log10(min(allv))):ceil(log10(max(allv)))));
xlabel('Zonal wavenumber');
ylabel('Power');
Label('(g)');
grid on;
box on;

max_lag = numel(cfg.lags);
lag = 0:max_lag;
ax_s(3) = axes(f);
hold on;
yline(0, '-', 'Color', [0.7, 0.7, 0.7], 'Alpha', 1);
for i = 1:2
    Band(lag, Stack(Xr{i}, @(X) Acf(X, max_lag)), col(i, :), false);
end
plot(lag, Acf(X1, max_lag), '-', 'Color', 'k', 'LineWidth', lw);
plot(lag, Acf(X2, max_lag), '--', 'Color', 'k', 'LineWidth', lw);
Mean_lines(lag, Xr, col, lw, @(X) Acf(X, max_lag));
xlim([0, max_lag]);
xticks(0:5:max_lag);
ylim([0, 1]);
yticks(0:0.5:1);
xlabel('Lag (days)');
ylabel('Autocorrelation');
Label('(h)');
grid on;
box on;

set(findall(f, '-property', 'FontSize'), 'FontSize', fs);
for ax = [ax_a, ax_b, ax_p, ax_s]
    set(ax.Title, 'FontSize', fs_lab, 'FontWeight', 'bold');
end
L0 = 0.07;
R0 = 0.985;
y_1 = 0.72;
h_1 = 0.24;
w_a = h_1 * fig_h / fig_w;
set(ax_a, 'Position', [L0 + 0.01, y_1, w_a, h_1]);
x_b = 0.36;
set(ax_b, 'Position', [x_b, y_1, R0 - x_b, h_1]);
Legend_top(leg, ax_b);
y_2 = 0.41;
h_2 = 0.17;
gap = 0.07;
w_2 = (R0 - L0 - 2 * gap) / 3;
for q = 1:3
    set(ax_p(q), 'Position', [L0 + (q - 1) * (w_2 + gap), y_2, w_2, h_2], 'XTickLabelRotation', 0);
end
y_3 = 0.075;
h_3 = 0.19;
for q = 1:3
    set(ax_s(q), 'Position', [L0 + (q - 1) * (w_2 + gap), y_3, w_2, h_3]);
end
exportgraphics(f, fullfile(img_dir, 'figure10.png'), 'Resolution', 600);

function Legend_top(leg, ax)
    set(leg, 'Orientation', 'horizontal', 'Units', 'normalized', 'ItemTokenSize', [20, 18]);
    drawnow;
    lp = leg.Position;
    ap = ax.Position;
    leg.Position = [ap(1) + ap(3) - lp(3) - 0.005, ap(2) + ap(4) - lp(4) - 0.01, lp(3), lp(4)];
end

function Corr_lim(lo)
    if lo >= 0
        ylim([0, 1]);
        yticks(0:0.5:1);
    else
        ylim([-0.5, 1]);
        yticks(-0.5:0.5:1);
    end
end

function Label(s)
    title(s, 'Units', 'normalized', 'Position', [0, 1.02], 'HorizontalAlignment', 'left');
end

function Y = Stack(Xr, fun)
    if isempty(Xr)
        Y = [];
        return
    end
    Y = cell2mat(cellfun(@(X) fun(X), Xr(:), 'UniformOutput', false));
end

function Band(x, Y, c, is_log)
    m = mean(Y, 1);
    s = std(Y, 0, 1);
    lo = m - s;
    if is_log
        lo = max(lo, m / 10);
    end
    fill([x, fliplr(x)], [lo, fliplr(m + s)], c, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function Mean_lines(x, Xr, col, lw, fun)
    sty = {'--', '-'};
    for i = [2, 1]
        plot(x, mean(Stack(Xr{i}, fun), 1), sty{i}, 'Color', col(i, :), 'LineWidth', lw);
    end
end

function P = Spectrum(X, Nlon)
    F = fft(X, [], 2);
    P = mean(abs(F(:, 2:Nlon / 2 + 1)).^2, 1) / Nlon;
end

function c = Acf(X, max_lag)
    X = X - mean(X, 1);
    v = sum(X.^2, 1);
    c = ones(1, max_lag + 1);
    for d = 1:max_lag
        c(d + 1) = mean(sum(X(1 + d:end, :) .* X(1:end - d, :), 1) ./ v);
    end
end

function c = Corr(a, b)
    a = a - mean(a);
    b = b - mean(b);
    c = sum(a .* b) / sqrt(sum(a.^2) * sum(b.^2));
end

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

function Lon_ticks()
    xticks(0:90:360);
    xticklabels({sprintf('0%s', char(176)), sprintf('90%sE', char(176)), sprintf('180%s', char(176)), ...
        sprintf('90%sW', char(176)), sprintf('0%s', char(176))});
end

function [P, td] = Z500_closed_loop(mdl, H, T, thr, R, gamma)
    P = nan(T, size(H, 2));
    td = NaN;
    nR = size(R, 1);
    for t = 1:T
        y = Predict_NGRC(mdl, Delay_input(H, mdl.nn, mdl.m, mdl.k));
        x = reshape(y, 1, []);
        if gamma > 0
            x = x + gamma * R(randi(nR), :);
        end
        if any(~isfinite(x)) || max(abs(x)) > thr
            td = t;
            return
        end
        P(t, :) = x;
        H = [H(2:end, :); x];
    end
end

function R = Z500_residuals(mdl, S)
    N = size(S, 2);
    U = Delay_input(S(1:end-1, :), mdl.nn, mdl.m, mdl.k);
    Y = reshape(S(mdl.k+1:end, :), [], 1);
    Yhat = zeros(size(Y));
    cs = 20000;
    for i0 = 1:cs:numel(Y)
        ii = i0:min(i0 + cs - 1, numel(Y));
        Yhat(ii) = Predict_NGRC(mdl, U(ii, :));
    end
    R = reshape(Y - Yhat, [], N);
end

function Draw_map(ax, lat0, lon_pts, lat_min, fs)
    axes(ax);
    hold on;
    axis equal off;
    c_coast = [0, 0, 0];
    c_grid = [0.85, 0.85, 0.85];
    c_lat0 = [0, 0, 0];
    lon_c = 0:1:360;
    for la = [30, 70]
        [x, y] = Proj(la * ones(size(lon_c)), lon_c);
        plot(x, y, '-', 'Color', c_grid, 'LineWidth', 0.6);
    end
    for lo = 0:30:330
        [x, y] = Proj(linspace(lat_min, 90, 50), lo * ones(1, 50));
        plot(x, y, '-', 'Color', c_grid, 'LineWidth', 0.6);
    end
    [x, y] = Proj(lat_min * ones(size(lon_c)), lon_c);
    plot(x, y, '-', 'Color', c_coast, 'LineWidth', 0.8);
    [la, lo] = Coast_lines();
    if ~isempty(la)
        la(la < lat_min) = NaN;
        [x, y] = Proj(la, lo);
        plot(x, y, '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);
    end
    [x, y] = Proj(lat0 * ones(size(lon_c)), lon_c);
    plot(x, y, '-', 'Color', c_lat0, 'LineWidth', 1.2);
    [x, y] = Proj(lat0 * ones(size(lon_pts)), lon_pts);
    plot(x, y, 'o', 'Color', c_lat0, 'MarkerFaceColor', c_lat0, 'MarkerSize', 2.5);
    [x, y] = Proj(lat0 - 14, 195);
    text(x, y, sprintf('%d%sN', lat0, char(176)), 'Color', c_lat0, 'FontSize', fs, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    labels = {sprintf('0%s', char(176)), sprintf('90%sE', char(176)), ...
        sprintf('180%s', char(176)), sprintf('90%sW', char(176))};
    lons = [0, 90, 180, 270];
    ha = {'center', 'left', 'center', 'right'};
    va = {'top', 'middle', 'bottom', 'middle'};
    for q = 1:4
        [x, y] = Proj(lat_min - 2, lons(q));
        text(x, y, labels{q}, 'FontSize', fs, 'Color', 'k', 'HorizontalAlignment', ha{q}, 'VerticalAlignment', va{q});
    end
    r_out = Rad(lat_min);
    xlim([-r_out, r_out]);
    ylim([-r_out, r_out]);
    set(ax, 'Clipping', 'off');
end

function [x, y] = Proj(lat, lon)
    x = Rad(lat) .* sind(lon);
    y = -Rad(lat) .* cosd(lon);
end

function r = Rad(lat)
    r = tand((90 - lat) / 2);
end

function [la, lo] = Coast_lines()
    la = [];
    lo = [];
    try
        C = load('coastlines');
        la = C.coastlat;
        lo = C.coastlon;
        return
    catch
    end
    try
        T = load('topo');
        M = contourc(0.5:1:359.5, -89.5:1:89.5, T.topo, [0, 0]);
        hdr = false(1, size(M, 2));
        k = 1;
        while k < size(M, 2)
            hdr(k) = true;
            k = k + M(2, k) + 1;
        end
        lo = M(1, :);
        la = M(2, :);
        lo(hdr) = NaN;
        la(hdr) = NaN;
    catch
    end
end
