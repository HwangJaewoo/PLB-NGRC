%% Coupled Ikeda map: true and predicted trajectories, and VPT, W1 and largest Lyapunov exponent of the four methods (output figure4.png)
clear; clc; close all;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'core'), genpath(fullfile(root, 'experiments')));
img_dir = fullfile(root, 'figures', 'img');
if ~isfolder(img_dir)
    mkdir(img_dir);
end
dir_ik = fullfile(root, 'experiments', 'Ikeda');

load(fullfile(dir_ik, 'data', 'Ikeda_N010_data.mat'), 'S', 'N', 'i_test');
load(fullfile(dir_ik, 'results', 'Ikeda_N010_true_lyapunov.mat'), 'lle');
A = load(fullfile(dir_ik, 'results', 'Ikeda_N010_PNGRC.mat'), 'res');
B = load(fullfile(dir_ik, 'results', 'Ikeda_N010_PLBNGRC.mat'), 'res');
r = {A.res, B.res};

ig = 1;
fig_w = 22;
fig_h = 11.6;
fs = 10;
ms = 6;
col = [0, 0.2, 0.8; 0.85, 0, 0];
bw = 0.6;
lw_err = 1.2;
lw_ref = 1.2;

T = size(r{1}.P_long, 1);
X = S(i_test(1:T), :);
traj = {X, r{1}.P_long, r{2}.P_long};
v = cell2mat(cellfun(@(Z) Z(:, [ig, N + ig]), traj, 'UniformOutput', false)');
[xl, xt] = Nice_lim(min(v(:, 1)), max(v(:, 1)));
[yl, yt] = Nice_lim(min(v(:, 2)), max(v(:, 2)));
c_traj = [0, 0, 0; col];

names_b = {'P-NGRC', 'LB-NGRC', 'P-RC', 'PLB-NGRC'};
files_b = {'Ikeda_N010_PNGRC.mat', 'Ikeda_N010_LBNGRC.mat', 'Ikeda_N010_PRC.mat', 'Ikeda_N010_PLBNGRC.mat'};
col_b = [0, 0.2, 0.8; 0.93, 0.55, 0; 0, 0.55, 0.25; 0.85, 0, 0];
nb = numel(files_b);
[vpt_m, vpt_s, w1_m, w1_s, lle_m, lle_s] = deal(zeros(1, nb));
for i = 1:nb
    tmp = load(fullfile(dir_ik, 'results', files_b{i}), 'res');
    rb = tmp.res;
    vpt_m(i) = mean(rb.vpt);
    vpt_s(i) = std(rb.vpt);
    w1_m(i) = mean(rb.D_all, 'omitnan');
    w1_s(i) = std(rb.D_all, 'omitnan');
    lle_m(i) = mean(rb.lle_all, 'omitnan');
    lle_s(i) = std(rb.lle_all, 'omitnan');
    fprintf('%-12s: VPT %.2f +- %.2f | W1 %.4f +- %.4f | lambda %.3f +- %.3f | diverged %d/10\n', names_b{i}, ...
        vpt_m(i), vpt_s(i), w1_m(i), w1_s(i), lle_m(i), lle_s(i), sum(rb.div_all));
    if i == 1
        w1_ref = rb.D_floor;
    end
end
lle_true = mean(lle);
fprintf('W1 reference %.4f, true lambda %.3f\n', w1_ref, lle_true);

f = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, fig_w, fig_h]);

lab_top = {'(a)', '(b)', '(c)'};
ax_top = gobjects(1, 3);
for q = 1:3
    ax_top(q) = axes(f);
    plot(traj{q}(:, ig), traj{q}(:, N + ig), '.', 'Color', c_traj(q, :), 'MarkerSize', ms);
    xlim(xl);
    ylim(yl);
    xticks(xt);
    yticks(yt);
    xlabel(sprintf('$x_{%d}$', ig), 'Interpreter', 'latex');
    ylabel(sprintf('$y_{%d}$', ig), 'Interpreter', 'latex');
    Label(lab_top{q});
    grid on;
    box on;
end

ax_d = axes(f);
hb = Bar_err(vpt_m, vpt_s, bw, lw_err, col_b, false);
Nice_ylim(max(vpt_m + vpt_s));
ylabel('VPT');
Panel_style('(d)', names_b);

ax_e = axes(f);
Bar_err(w1_m, w1_s, bw, lw_err, col_b, true);
yline(w1_ref, '--', 'Color', 'k', 'LineWidth', lw_ref, 'Alpha', 1);
Log_ylim(min([w1_m, w1_ref]), max(w1_m + w1_s));
ylabel('$\mathcal{W}_1$', 'Interpreter', 'latex');
Panel_style('(e)', names_b);

ax_f = axes(f);
Bar_err(lle_m, lle_s, bw, lw_err, col_b, false);
yline(lle_true, '--', 'Color', 'k', 'LineWidth', lw_ref, 'Alpha', 1);
yline(0, '-', 'Color', 'k', 'LineWidth', 0.5, 'Alpha', 1);
Nice_ylim2(min(lle_m - lle_s), max([lle_m + lle_s, lle_true]));
ylabel('$\lambda_{\max}$', 'Interpreter', 'latex');
Panel_style('(f)', names_b);

set(findall(f, '-property', 'FontSize'), 'FontSize', fs);
axs = [ax_top, ax_d, ax_e, ax_f];
for ax = axs
    set(ax.Title, 'FontSize', fs + 2, 'FontWeight', 'bold');
    if strcmp(ax.XLabel.Interpreter, 'latex'), set(ax.XLabel, 'FontSize', fs + 2); end
    if strcmp(ax.YLabel.Interpreter, 'latex'), set(ax.YLabel, 'FontSize', fs + 2); end
end
L = 0.06;  R = 0.99;  gap = 0.08;
w = (R - L - 2 * gap) / 3;
y_1 = 0.60;  y_2 = 0.13;  h = 0.34;
x_col = L + (0:2) * (w + gap);
for q = 1:3
    set(ax_top(q), 'Position', [x_col(q), y_1, w, h]);
end
set(ax_d, 'Position', [x_col(1), y_2, w, h]);
set(ax_e, 'Position', [x_col(2), y_2, w, h]);
set(ax_f, 'Position', [x_col(3), y_2, w, h]);
lg = legend(ax_d, hb, names_b, 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', fs);
lg.Units = 'normalized';
lg.Position(1:2) = [0.5 - lg.Position(3) / 2, 0.005];
exportgraphics(f, fullfile(img_dir, 'figure4.png'), 'Resolution', 600);

function hb = Bar_err(m, s, bw, lw_err, col, is_log)
    hold on;
    n = numel(m);
    hb = gobjects(1, n);
    base = 0;
    if is_log
        set(gca, 'YScale', 'log');
        base = 10^(floor(log10(min(m(m > 0)))) - 1);
    end
    for i = 1:n
        hb(i) = bar(i, m(i), bw, 'FaceColor', col(i, :), 'EdgeColor', 'none', 'BaseValue', base);
        hb(i).BaseLine.Visible = 'off';
    end
    s_lo = s;
    if is_log
        s_lo = min(s, 0.9 * m);
    end
    errorbar(1:n, m, s_lo, s, 'k', 'LineStyle', 'none', 'LineWidth', lw_err, 'CapSize', 5, 'HandleVisibility', 'off');
end

function Log_ylim(lo, hi)
    e = floor(log10(lo)):ceil(log10(hi));
    ylim(10.^[e(1), e(end)]);
    yticks(10.^e);
    yticklabels(arrayfun(@(k) sprintf('10^{%d}', k), e, 'UniformOutput', false));
    set(gca, 'YMinorTick', 'off', 'YMinorGrid', 'off');
end

function Panel_style(label, names)
    xlim([0.4, numel(names) + 0.6]);
    set(gca, 'XTick', 1:numel(names), 'XTickLabel', {}, 'XGrid', 'off');
    Label(label);
    grid on;
    box on;
end

function Label(s)
    title(s, 'Units', 'normalized', 'Position', [0, 1.02], 'HorizontalAlignment', 'left');
end

function step = Nice_step(span, n_max)
    if nargin < 2
        n_max = 4;
    end
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
