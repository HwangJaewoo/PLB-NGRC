%% Lorenz-96: VPT, W1 and largest Lyapunov exponent of P-NGRC, LB-NGRC, P-RC and PLB-NGRC (output figure3.png)
clear; clc; close all;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'core'), genpath(fullfile(root, 'experiments')));
img_dir = fullfile(root, 'figures', 'img');
if ~isfolder(img_dir)
    mkdir(img_dir);
end
res_dir = fullfile(root, 'experiments', 'Lorenz96', 'results');

load(fullfile(res_dir, 'Lorenz96_N020_true_lyapunov.mat'), 'lle');
files = {'Lorenz96_N020_PNGRC.mat', 'Lorenz96_N020_LBNGRC.mat', 'Lorenz96_N020_PRC.mat', 'Lorenz96_N020_PLBNGRC.mat'};
nm = numel(files);
r = cell(1, nm);
for i = 1:nm
    tmp = load(fullfile(res_dir, files{i}), 'res');
    r{i} = tmp.res;
end

[vpt_m, vpt_s, w1_m, w1_s, lle_m, lle_s] = deal(zeros(1, nm));
for i = 1:nm
    vpt_m(i) = mean(r{i}.vpt);
    vpt_s(i) = std(r{i}.vpt);
    w1_m(i) = mean(r{i}.D_all, 'omitnan');
    w1_s(i) = std(r{i}.D_all, 'omitnan');
    lle_m(i) = mean(r{i}.lle_all, 'omitnan');
    lle_s(i) = std(r{i}.lle_all, 'omitnan');
end
w1_ref = r{1}.D_floor;
lle_true = mean(lle);
names = {'P-NGRC', 'LB-NGRC', 'P-RC', 'PLB-NGRC'};
for i = 1:nm
    fprintf('%-12s VPT %.2f +- %.2f | W1 %.4f +- %.4f | lambda %.3f +- %.3f | survived %d / %d\n', names{i}, ...
        vpt_m(i), vpt_s(i), w1_m(i), w1_s(i), lle_m(i), lle_s(i), sum(~r{i}.div_all), numel(r{i}.div_all));
end
fprintf('reference W1 %.4f, true lambda %.3f\n', w1_ref, lle_true);

fig_w = 22;
fig_h = 5.6;
fs = 9;
col = [0, 0.2, 0.8; 0.93, 0.55, 0; 0, 0.55, 0.25; 0.85, 0, 0];
bw = 0.7;
lw_err = 1.2;
lw_ref = 1.2;

L = 0.06;
R = 0.99;
gap = 0.08;
y0 = 0.22;
h = 0.66;
w = (R - L - 2 * gap) / 3;

f = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, fig_w, fig_h]);

ax_a = axes(f);
hb = Bar_err(vpt_m, vpt_s, bw, lw_err, col, false);
Nice_ylim(max(vpt_m + vpt_s));
ylabel('VPT');
Panel_style('(a)', names);

ax_b = axes(f);
Bar_err(w1_m, w1_s, bw, lw_err, col, true);
yline(w1_ref, '--', 'Color', 'k', 'LineWidth', lw_ref, 'Alpha', 1);
Log_ylim(min([w1_m, w1_ref]), max(w1_m + w1_s));
ylabel('$\mathcal{W}_1$', 'Interpreter', 'latex');
Panel_style('(b)', names);

ax_c = axes(f);
Bar_err(lle_m, lle_s, bw, lw_err, col, false);
yline(lle_true, '--', 'Color', 'k', 'LineWidth', lw_ref, 'Alpha', 1);
Nice_ylim(max([lle_m + lle_s, lle_true]));
ylabel('$\lambda_{\max}$', 'Interpreter', 'latex');
Panel_style('(c)', names);

set(findall(f, '-property', 'FontSize'), 'FontSize', fs);
for ax = [ax_a, ax_b, ax_c]
    set(ax.Title, 'FontSize', fs + 2, 'FontWeight', 'bold');
end
set([ax_b.YLabel, ax_c.YLabel], 'FontSize', fs + 2);
set(ax_a, 'Position', [L, y0, w, h]);
set(ax_b, 'Position', [L + w + gap, y0, w, h]);
set(ax_c, 'Position', [L + 2 * (w + gap), y0, w, h]);
lg = legend(ax_a, hb, names, 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', fs);
lg.Units = 'normalized';
lg.Position(1:2) = [0.5 - lg.Position(3) / 2, 0.01];
exportgraphics(f, fullfile(img_dir, 'figure3.png'), 'Resolution', 600);

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
    n = numel(names);
    xlim([0.4, n + 0.6]);
    set(gca, 'XTick', 1:n, 'XTickLabel', {}, 'XGrid', 'off');
    title(label, 'Units', 'normalized', 'Position', [0, 1.02], 'HorizontalAlignment', 'left');
    grid on;
    box on;
end

function Nice_ylim(top)
    Tight_lim(0, top);
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
