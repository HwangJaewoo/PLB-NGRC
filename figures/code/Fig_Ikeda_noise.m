%% Coupled Ikeda map with observational noise: observed and PLB-NGRC trajectories, and VPT, W1 and largest Lyapunov exponent for 1%, 5% and 10% noise (output figure6.png)
clear; clc; close all;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'core'), genpath(fullfile(root, 'experiments')));
img_dir = fullfile(root, 'figures', 'img');
if ~isfolder(img_dir)
    mkdir(img_dir);
end
dir_ik = fullfile(root, 'experiments', 'Ikeda');

sig_list = [0.01, 0.05, 0.10];
ig = 1;
fig_w = 22;
fig_h = 11;
fs = 10;
ms_true = 6;
ms_pred = 6;
c_plb = [0.85, 0, 0];
c_bar = [0.9, 0.3, 0.3];
bw = 0.6;
lw_err = 1.2;
lw_ref = 1.2;
nl = numel(sig_list);
lev = arrayfun(@(s) sprintf('%d%%', round(100 * s)), sig_list, 'UniformOutput', false);

load(fullfile(dir_ik, 'data', 'Ikeda_N010_data.mat'), 'S_raw', 'N', 'n_fit', 'i_test');
load(fullfile(dir_ik, 'results', 'Ikeda_N010_true_lyapunov.mat'), 'lle');
lle_true = mean(lle);
cx = 1:N;
cy = N + (1:N);

[vpt_m, vpt_s, w1_m, w1_s, w1_ref, lle_m, lle_s] = deal(zeros(1, nl));
X_att = cell(1, nl);
P_att = cell(1, nl);
for q = 1:nl
    sig = sig_list(q);
    R = load(fullfile(dir_ik, 'results', sprintf('Ikeda_N010_PLBNGRC_noise%02d.mat', round(100 * sig))), 'res');
    r = R.res;
    vpt_m(q) = mean(r.vpt);                vpt_s(q) = std(r.vpt);
    w1_m(q) = mean(r.D_all, 'omitnan');    w1_s(q) = std(r.D_all, 'omitnan');
    lle_m(q) = mean(r.lle_all, 'omitnan'); lle_s(q) = std(r.lle_all, 'omitnan');
    w1_ref(q) = r.D_floor;

    rng(10);
    E = randn(size(S_raw));
    sd0 = [std(S_raw(1:n_fit, cx), 0, 'all'), std(S_raw(1:n_fit, cy), 0, 'all')];
    S_obs = S_raw;
    S_obs(:, cx) = S_raw(:, cx) + sig * sd0(1) * E(:, cx);
    S_obs(:, cy) = S_raw(:, cy) + sig * sd0(2) * E(:, cy);
    mu = [mean(S_obs(1:n_fit, cx), 'all'), mean(S_obs(1:n_fit, cy), 'all')];
    sd = [std(S_obs(1:n_fit, cx), 0, 'all'), std(S_obs(1:n_fit, cy), 0, 'all')];
    S = S_obs;
    S(:, cx) = (S_obs(:, cx) - mu(1)) / sd(1);
    S(:, cy) = (S_obs(:, cy) - mu(2)) / sd(2);

    T_cl = size(r.P_long, 1);
    fprintf('noise %4s: VPT %.2f, W1 %.4f, lambda %.3f, diverged %d/10\n', ...
        lev{q}, vpt_m(q), w1_m(q), lle_m(q), sum(r.div_all));

    X_att{q} = S(i_test(1:T_cl), :);
    P_att{q} = r.P_long;
end

v = cell2mat(cellfun(@(Z) Z(:, [ig, N + ig]), [X_att, P_att], 'UniformOutput', false)');
[xl, xt] = Nice_lim(min(v(:, 1)), max(v(:, 1)));
[yl, yt] = Nice_lim(min(v(:, 2)), max(v(:, 2)));

f = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, fig_w, fig_h]);

lab_top = {'(a)', '(b)', '(c)'};
ax_top = gobjects(1, nl);
for q = 1:nl
    ax_top(q) = axes(f);
    hold on;
    plot(X_att{q}(:, ig), X_att{q}(:, N + ig), '.', 'Color', 'k', 'MarkerSize', ms_true);
    plot(P_att{q}(:, ig), P_att{q}(:, N + ig), '.', 'Color', c_plb, 'MarkerSize', ms_pred);
    xlim(xl);  ylim(yl);  xticks(xt);  yticks(yt);
    xlabel(sprintf('$x_{%d}$', ig), 'Interpreter', 'latex');
    ylabel(sprintf('$y_{%d}$', ig), 'Interpreter', 'latex');
    Label(lab_top{q});
    grid on;
    box on;
end

ax_d = axes(f);
Bar_err(vpt_m, vpt_s, bw, lw_err, c_bar);
Nice_ylim(max(vpt_m + vpt_s));
ylabel('VPT');
Bar_style('(d)', lev);

ax_e = axes(f);
Bar_err(w1_m, w1_s, bw, lw_err, c_bar);
x_edge = [0.4, (1:nl - 1) + 0.5, nl + 0.6];
x_ref = reshape([x_edge(1:end - 1); x_edge(2:end)], 1, []);
y_ref = reshape([w1_ref; w1_ref], 1, []);
plot(x_ref, y_ref, '--', 'Color', 'k', 'LineWidth', lw_ref);
Nice_ylim(max([w1_m + w1_s, w1_ref]));
ylabel('$\mathcal{W}_1$', 'Interpreter', 'latex');
Bar_style('(e)', lev);

ax_f = axes(f);
Bar_err(lle_m, lle_s, bw, lw_err, c_bar);
yline(lle_true, '--', 'Color', 'k', 'LineWidth', lw_ref, 'Alpha', 1);
Nice_ylim2(min(lle_m - lle_s), max([lle_m + lle_s, lle_true]));
ylabel('$\lambda_{\max}$', 'Interpreter', 'latex');
Bar_style('(f)', lev);

set(findall(f, '-property', 'FontSize'), 'FontSize', fs);
for ax = [ax_top, ax_d, ax_e, ax_f]
    set(ax.Title, 'FontSize', fs + 2, 'FontWeight', 'bold');
    if strcmp(ax.XLabel.Interpreter, 'latex'), set(ax.XLabel, 'FontSize', fs + 2); end
    if strcmp(ax.YLabel.Interpreter, 'latex'), set(ax.YLabel, 'FontSize', fs + 2); end
end
L = 0.06;  R = 0.99;  gap = 0.08;
w = (R - L - 2 * gap) / 3;
y_1 = 0.58;  y_2 = 0.10;  h = 0.36;
x_col = L + (0:2) * (w + gap);
for q = 1:nl
    set(ax_top(q), 'Position', [x_col(q), y_1, w, h]);
end
set(ax_d, 'Position', [x_col(1), y_2, w, h]);
set(ax_e, 'Position', [x_col(2), y_2, w, h]);
set(ax_f, 'Position', [x_col(3), y_2, w, h]);
exportgraphics(f, fullfile(img_dir, 'figure6.png'), 'Resolution', 600);

function Bar_err(m, s, bw, lw_err, c)
    hold on;
    b = bar(1:numel(m), m, bw, 'FaceColor', c, 'EdgeColor', 'none');
    b.BaseLine.Visible = 'off';
    errorbar(1:numel(m), m, s, 'k', 'LineStyle', 'none', 'LineWidth', lw_err, 'CapSize', 6);
end

function Bar_style(label, names)
    xlim([0.4, numel(names) + 0.6]);
    set(gca, 'XTick', 1:numel(names), 'XTickLabel', names);
    xlabel('Noise level');
    Label(label);
    grid on;
    box on;
end

function Label(s)
    title(s, 'Units', 'normalized', 'Position', [0, 1.02], 'HorizontalAlignment', 'left');
end

function step = Nice_step(span, n_max)
    e = 10^floor(log10(span / n_max));
    c = [1, 2, 5, 10] * e;
    step = c(find(span ./ c <= n_max, 1));
end

function [lim, ticks] = Nice_lim(lo, hi)
    n_max = 6;
    pad = 0.03 * (hi - lo);
    lo = lo - pad;
    hi = hi + pad;
    step = Nice_step(hi - lo, n_max);
    lim = [floor(lo / step), ceil(hi / step)] * step;
    if (lim(2) - lim(1)) / step > n_max
        step = Nice_step((lim(2) - lim(1)) * 1.01, n_max);
        lim = [floor(lo / step), ceil(hi / step)] * step;
    end
    ticks = lim(1):step:lim(2);
end

function Nice_ylim(top)
    Tight_lim(0, top);
end

function Nice_ylim2(lo, hi)
    Tight_lim(min(0, lo), max(0, hi));
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
