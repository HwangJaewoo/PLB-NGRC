%% Z500 over the latitude band: map, fields, forecast correlation and long-term statistics of P-NGRC and PLB-NGRC (output figure11.png)
clear; clc; close all;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'core'), genpath(fullfile(root, 'experiments')));
img_dir = fullfile(root, 'figures', 'img');
if ~isfolder(img_dir)
    mkdir(img_dir);
end
z_dir = fullfile(root, 'experiments', 'Z500_2D');

cfg = Z500_2D_setup();
D = load(Z500_2D_data(), 'lon', 'lat');
lon = double(D.lon(:)');
lat = double(D.lat(:)');
A = load(fullfile(z_dir, 'results', 'Z500_2D_PNGRC_test.mat'));
B = load(fullfile(z_dir, 'results', 'Z500_2D_PLBNGRC_test.mat'));
r = {A, B};
S = cfg.S;
Nlon = cfg.Nlon;
Nlat = cfg.Nlat;

fig_w = 22;
fig_h = 22.6;
fs = 11;
fs_lab = 12;
lw = 1.2;
ms = 4;
leads = [1, 2, 3];
lat_min = 20;
thr = 0.6;
cl = [-3, 3];
col = [0, 0.2, 0.8; 0.85, 0, 0];
c_thr = [0, 0.6, 0];

s0 = cfg.i_test(1);
n_half = floor(numel(cfg.i_test) / 2);
X1 = S(cfg.i_test(1:n_half), :);
X2 = S(cfg.i_test(n_half + 1:2 * n_half), :);
Pf = cell(1, 2);
Xr = cell(1, 2);
for i = 1:2
    p = r{i}.p;
    if ~isfield(p, 'fa')
        p.fa = false;
    end
    mdl = Z500_2D_train(cfg, p);
    [Pf{i}, ~] = Z500_2D_closed_loop(mdl, cfg, [], s0, max(leads), 0, 0);
    R = Z500_2D_residuals(mdl, cfg, 1, 0.99);
    Xr{i} = cell(1, 10);
    for sd = 1:10
        [Xl, td] = Z500_2D_closed_loop(mdl, cfg, R, cfg.i_test(1) - 1, numel(cfg.i_test), sd, 1);
        if isnan(td)
            Xr{i}{sd} = Xl;
        end
    end
    Xr{i} = Xr{i}(~cellfun(@isempty, Xr{i}));
end
names = {'P-NGRC', 'PLB-NGRC'};
for i = 1:2
    a = r{i}.acc_m;
    ap = r{i}.accp_m;
    j = find(a < ap, 1);
    fprintf('%-8s: ACC 1 d %.2f (persistence %.2f), horizon %.2f d (persistence %.2f d), long-term runs %d / 10\n', ...
        names{i}, a(1), ap(1), Horizon(a, thr), Horizon(ap, thr), numel(Xr{i}));
    fprintf('   ACC 1-15 d: %s\n   persistence: %s\n', mat2str(round(a, 2)), mat2str(round(ap, 2)));
    if isempty(j)
        fprintf('   stays above persistence up to 15 d\n');
    else
        fprintf('   first lead below persistence: %d d\n', j);
    end
    for l = leads
        fprintf('   first day of the test set, lead %d d: correlation %.2f\n', l, Corr(Pf{i}(l, :), S(s0 + l, :)));
    end
end

f = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, fig_w, fig_h]);

ax_a = axes(f);
Draw_band(ax_a, lat, lon, lat_min, fs);
Label('(a)');

ax_b = axes(f);
hold on;
L = numel(r{2}.acc_m);
yline(thr, '--', 'Color', c_thr, 'LineWidth', lw, 'Alpha', 1, 'HandleVisibility', 'off');
h_p = plot(1:L, r{2}.accp_m, '--s', 'Color', 'k', 'LineWidth', lw, 'MarkerSize', ms);
h_2 = plot(1:L, r{2}.acc_m, '-o', 'Color', col(2, :), 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', col(2, :));
h_1 = plot(1:L, r{1}.acc_m, '--o', 'Color', col(1, :), 'LineWidth', lw, 'MarkerSize', ms);
xlim([0, 15]);
xticks(0:5:15);
Corr_lim(min([r{1}.acc_m, r{2}.acc_m, r{2}.accp_m]));
xlabel('Lead time (days)');
ylabel('Correlation');
leg = legend([h_2, h_1, h_p], {'PLB-NGRC', 'P-NGRC', 'Persistence'}, 'Box', 'off');
Label('(b)');
grid on;
box on;

lab_t = {'(c)', '(d)', '(e)'};
lab_n = {'(f)', '(g)', '(h)'};
lab_q = {'(i)', '(j)', '(k)'};
ax_t = gobjects(1, numel(leads));
ax_n = gobjects(1, numel(leads));
ax_q = gobjects(1, numel(leads));
for q = 1:numel(leads)
    l = leads(q);
    ax_t(q) = axes(f);
    Field(S(s0 + l, :), lon, lat, Nlon, Nlat, cl);
    Field_style(lab_t{q}, sprintf('True, %d day%s', l, Plural(l)), q == 1, false);
    ax_n(q) = axes(f);
    Field(Pf{1}(l, :), lon, lat, Nlon, Nlat, cl);
    Field_style(lab_n{q}, sprintf('P-NGRC, %d day%s', l, Plural(l)), q == 1, false);
    ax_q(q) = axes(f);
    Field(Pf{2}(l, :), lon, lat, Nlon, Nlat, cl);
    Field_style(lab_q{q}, sprintf('PLB-NGRC, %d day%s', l, Plural(l)), q == 1, true);
end
cb = colorbar(ax_q(end));
cb.Ticks = [cl(1), 0, cl(2)];

edges = -4.5:0.25:4.5;
xc = edges(1:end - 1) + diff(edges) / 2;
ax_s = gobjects(1, 3);
ax_s(1) = axes(f);
hold on;
pk = 0;
for i = 1:2
    Y = Stack(Xr{i}, @(X) histcounts(X(:), edges, 'Normalization', 'pdf'));
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
Label('(l)');
grid on;
box on;

kk = 1:Nlon / 2;
ax_s(2) = axes(f);
hold on;
for i = 1:2
    Band(kk, Stack(Xr{i}, @(X) Spectrum(X, Nlon)), col(i, :), true);
end
s1 = Spectrum(X1, Nlon);
s2 = Spectrum(X2, Nlon);
plot(kk, s1, '-', 'Color', 'k', 'LineWidth', lw);
plot(kk, s2, '--', 'Color', 'k', 'LineWidth', lw);
Mean_lines(kk, Xr, col, lw, @(X) Spectrum(X, Nlon));
set(gca, 'YScale', 'log', 'YMinorGrid', 'off', 'YMinorTick', 'off');
xlim([0, 20]);
xticks(0:5:20);
allv = [s1, s2, mean(Stack(Xr{1}, @(X) Spectrum(X, Nlon)), 1), mean(Stack(Xr{2}, @(X) Spectrum(X, Nlon)), 1)];
ylim(10.^[floor(log10(min(allv))), ceil(log10(max(allv)))]);
yticks(10.^(floor(log10(min(allv))):ceil(log10(max(allv)))));
xlabel('Zonal wavenumber');
ylabel('Power');
Label('(m)');
grid on;
box on;

max_lag = cfg.max_lag;
lag = 0:max_lag;
ax_s(3) = axes(f);
hold on;
yline(0, '-', 'Color', [0.7, 0.7, 0.7], 'Alpha', 1);
lo_k = 1;
for i = 1:2
    Y = Stack(Xr{i}, @(X) Acf(X, max_lag));
    Band(lag, Y, col(i, :), false);
    lo_k = min(lo_k, min(mean(Y, 1) - std(Y, 0, 1)));
end
plot(lag, Acf(X1, max_lag), '-', 'Color', 'k', 'LineWidth', lw);
plot(lag, Acf(X2, max_lag), '--', 'Color', 'k', 'LineWidth', lw);
Mean_lines(lag, Xr, col, lw, @(X) Acf(X, max_lag));
xlim([0, max_lag]);
xticks(0:5:max_lag);
Corr_lim(lo_k);
xlabel('Lag (days)');
ylabel('Autocorrelation');
Label('(n)');
grid on;
box on;

set(findall(f, '-property', 'FontSize'), 'FontSize', fs);
for ax = [ax_a, ax_b, ax_t, ax_n, ax_q, ax_s]
    set(ax.Title, 'FontSize', fs_lab, 'FontWeight', 'bold');
end
L0 = 0.07;
R0 = 0.955;
cb_w = 0.012;
cb_gap = 0.01;
R1 = R0 - cb_gap - cb_w;
y_1 = 0.832;
h_1 = 0.143;
w_a = h_1 * fig_h / fig_w;
set(ax_a, 'Position', [L0 + 0.01, y_1, w_a, h_1]);
x_b = 0.37;
set(ax_b, 'Position', [x_b, y_1, R0 - x_b, h_1]);
Legend_top(leg, ax_b);
y_2 = 0.6305;
y_n = 0.4708;
y_3 = 0.311;
h_f = 0.1093;
gap_f = 0.03;
w_f = (R1 - L0 - 2 * gap_f) / 3;
for q = 1:numel(leads)
    set(ax_t(q), 'Position', [L0 + (q - 1) * (w_f + gap_f), y_2, w_f, h_f]);
    set(ax_n(q), 'Position', [L0 + (q - 1) * (w_f + gap_f), y_n, w_f, h_f]);
    set(ax_q(q), 'Position', [L0 + (q - 1) * (w_f + gap_f), y_3, w_f, h_f]);
end
set([ax_t, ax_n, ax_q], 'XTickLabelRotation', 0);
set(cb, 'Position', [R1 + cb_gap, y_3, cb_w, y_2 + h_f - y_3]);
y_4 = 0.0673;
h_4 = 0.1513;
gap_4 = 0.07;
w_4 = (R0 - L0 - 2 * gap_4) / 3;
for q = 1:3
    set(ax_s(q), 'Position', [L0 + (q - 1) * (w_4 + gap_4), y_4, w_4, h_4]);
end
exportgraphics(f, fullfile(img_dir, 'figure11.png'), 'Resolution', 600);

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

function s = Plural(l)
    s = '';
    if l > 1
        s = 's';
    end
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
    Z = reshape(X', Nlon, []);
    F = fft(Z, [], 1);
    P = mean(abs(F(2:Nlon / 2 + 1, :)).^2, 2)' / Nlon;
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

function h = Horizon(acc, thr)
    a = [1, acc];
    i = find(a < thr, 1);
    if isempty(i)
        h = numel(acc);
    else
        h = (i - 2) + (a(i - 1) - thr) / (a(i - 1) - a(i));
    end
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

function Field(v, lon, lat, Nlon, Nlat, cl)
    Z = reshape(double(v), Nlon, Nlat)';
    [LO, LA] = meshgrid([lon, 360], lat);
    pcolor(LO, LA, [Z, Z(:, 1)]);
    shading interp;
    set(gca, 'CLim', cl, 'Layer', 'top', 'Box', 'on');
    colormap(gca, jet);
end

function Field_style(label, name, show_ylabel, show_xlabel)
    xlim([0, 360]);
    Lon_ticks();
    yticks([30, 50, 70]);
    yticklabels({sprintf('30%sN', char(176)), sprintf('50%sN', char(176)), sprintf('70%sN', char(176))});
    if ~show_ylabel
        yticklabels({});
    end
    if show_xlabel
        xlabel('Longitude');
    else
        xticklabels({});
    end
    title(label, 'Units', 'normalized', 'Position', [0, 1.02], 'HorizontalAlignment', 'left');
    text(1, 1.02, name, 'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
end

function Lon_ticks()
    xticks(0:90:360);
    xticklabels({sprintf('0%s', char(176)), sprintf('90%sE', char(176)), sprintf('180%s', char(176)), ...
        sprintf('90%sW', char(176)), sprintf('0%s', char(176))});
end

function Draw_band(ax, lat, lon, lat_min, fs)
    axes(ax);
    hold on;
    axis equal off;
    c_coast = [0.5, 0.5, 0.5];
    c_grid = [0.85, 0.85, 0.85];
    c_band = [0, 0, 0];
    lon_c = 0:1:360;
    for lo = 0:30:330
        [x, y] = Proj(linspace(lat_min, 90, 50), lo * ones(1, 50));
        plot(x, y, '-', 'Color', c_grid, 'LineWidth', 0.6);
    end
    [x, y] = Proj(lat_min * ones(size(lon_c)), lon_c);
    plot(x, y, '-', 'Color', [0, 0, 0], 'LineWidth', 0.8);
    [la, lo] = Coast_lines();
    if ~isempty(la)
        la(la < lat_min) = NaN;
        [x, y] = Proj(la, lo);
        plot(x, y, '-', 'Color', c_coast, 'LineWidth', 0.5);
    end
    for la0 = [min(lat), max(lat)]
        [x, y] = Proj(la0 * ones(size(lon_c)), lon_c);
        plot(x, y, '-', 'Color', c_band, 'LineWidth', 1.2);
    end
    [LO, LA] = meshgrid(lon, lat);
    [x, y] = Proj(LA(:), LO(:));
    plot(x, y, '.', 'Color', c_band, 'MarkerSize', 4);
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
