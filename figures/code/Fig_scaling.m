%% Scaling with the number of sites N: VPT, W1, largest Lyapunov exponent and training time of P-NGRC, P-RC and PLB-NGRC (output figure8.png)
clear; clc; close all;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'core'), genpath(fullfile(root, 'experiments')));
img_dir = fullfile(root, 'figures', 'img');
if ~isfolder(img_dir)
    mkdir(img_dir);
end

sys_list = {'Lorenz96', 'Ikeda', 'Thomas'};
N_lists = {[20, 40, 100, 200], [10, 20, 50, 100], [10, 20, 50, 100]};
fs = 10;
sys_name = {'Lorenz-96', 'Ikeda', 'Thomas'};
ylab = {'VPT', '$\mathcal{W}_1$', '$\lambda_{\max}$', 'Training time (s)'};
col = struct('PNGRC', [0, 0.2, 0.8], 'PRC', [0, 0.55, 0.25], 'PLBNGRC', [0.85, 0, 0]);
sty = {'LineWidth', 1.2, 'MarkerSize', 4};
sty_eb = [sty, {'CapSize', 4}];

T = readtable(fullfile(root, 'timing', 'training_time.csv'), 'TextType', 'char');
methods = {'PNGRC', 'PRC', 'PLBNGRC'};
L = cell(1, 3);
fld = {'surv', 'D', 'D_floor', 'vpt', 'vpt_sd', 'lle', 'lle_sd', 'lle_true', 't_train'};
for r = 1:3
    N_list = N_lists{r};
    Z = nan(numel(N_list), numel(methods));
    S = struct('N_list', N_list, 'models', {methods});
    for q = 1:numel(fld)
        S.summ.(fld{q}) = Z;
    end
    d = fullfile(root, 'experiments', sys_list{r}, 'results');
    for j = 1:numel(N_list)
        Q = load(fullfile(d, sprintf('%s_N%03d_true_lyapunov.mat', sys_list{r}, N_list(j))), 'lle');
        for i = 1:numel(methods)
            nm = sprintf('%s_N%03d_%s.mat', sys_list{r}, N_list(j), methods{i});
            R = load(fullfile(d, nm), 'res');
            S.summ.surv(j, i) = R.res.surv;
            S.summ.D(j, i) = R.res.D;
            S.summ.D_floor(j, i) = R.res.D_floor;
            S.summ.vpt(j, i) = mean(R.res.vpt);
            S.summ.vpt_sd(j, i) = std(R.res.vpt);
            S.summ.lle(j, i) = R.res.lle;
            S.summ.lle_sd(j, i) = std(R.res.lle_all, 'omitnan');
            S.summ.lle_true(j, i) = mean(Q.lle);
            S.summ.t_train(j, i) = T.t_train_s(strcmp(T.file, nm));
        end
    end
    L{r} = S;
end

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 22, 16]);
tiledlayout(4, 3, 'TileSpacing', 'loose', 'Padding', 'compact');
ax_leg = [];
h_leg = [];
k = 0;
for c = 1:4
    for r = 1:3
        ax = nexttile;
        k = k + 1;
        if isempty(L{r})
            axis(ax, 'off');
            continue
        end
        S = L{r};
        N = S.N_list;
        hold(ax, 'on');
        h = gobjects(1, numel(S.models));
        div_lab = cell(numel(N) * numel(S.models), 3);
        n_lab = 0;
        for i = 1:numel(S.models)
            cc = col.(S.models{i});
            switch c
                case 1
                    h(i) = errorbar(ax, N, S.summ.vpt(:, i), S.summ.vpt_sd(:, i), 'o-', 'Color', cc, 'MarkerFaceColor', cc, sty_eb{:});
                case 2
                    h(i) = plot(ax, N, S.summ.D(:, i), 'o-', 'Color', cc, 'MarkerFaceColor', cc, sty{:});
                    if i == 1
                        plot(ax, N, S.summ.D_floor(:, i), 'k--', 'LineWidth', 1);
                    end
                case 3
                    h(i) = errorbar(ax, N, S.summ.lle(:, i), S.summ.lle_sd(:, i), 'o-', 'Color', cc, 'MarkerFaceColor', cc, sty_eb{:});
                    if i == 1
                        plot(ax, N, S.summ.lle_true(:, i), 'k--', 'LineWidth', 1);
                    end
                case 4
                    h(i) = plot(ax, N, S.summ.t_train(:, i), 'o-', 'Color', cc, 'MarkerFaceColor', cc, sty{:});
            end
            if c == 2 || c == 3
                if c == 2
                    v = S.summ.D(:, i);
                else
                    v = S.summ.lle(:, i);
                end
                for j = find(S.summ.surv(:, i) == 0)'
                    jf = find(~isnan(v(1:j - 1)), 1, 'last');
                    if ~isempty(jf)
                        n_lab = n_lab + 1;
                        div_lab(n_lab, :) = {N(j), v(jf), cc};
                    end
                end
            end
        end
        set(ax, 'XScale', 'log', 'XTick', N, 'FontSize', fs, 'Box', 'on');
        grid(ax, 'on');
        xlim(ax, [N(1) / 1.2, N(end) * 1.2]);
        if c == 2 || c == 4
            set(ax, 'YScale', 'log');
        end
        Tick_align(ax, c == 2 || c == 4);
        yl = ylim(ax);
        for q = 1:n_lab
            y = div_lab{q, 2};
            y = y + 0.15 * (yl(2) - yl(1)) * sign(mean(yl) - y);
            text(ax, div_lab{q, 1}, y, 'div.', 'Color', div_lab{q, 3}, 'FontSize', fs, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        if c == 1
            title(ax, sys_name{r}, 'FontWeight', 'normal');
            text(ax, 0, 1.02, sprintf('(%s)', char('a' + r - 1)), 'Units', 'normalized', 'FontWeight', 'bold', ...
                'FontSize', fs + 2, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
        end
        if c == 4
            xlabel(ax, '$N$', 'Interpreter', 'latex', 'FontSize', fs + 2);
        end
        if r == 1
            if any(ylab{c} == '$')
                ylabel(ax, ylab{c}, 'Interpreter', 'latex', 'FontSize', fs + 2);
            else
                ylabel(ax, ylab{c}, 'FontSize', fs);
            end
        end
        if isempty(ax_leg)
            ax_leg = ax;
            h_leg = h;
        end
    end
end
if ~isempty(ax_leg)
    lgd = legend(ax_leg, h_leg, {'P-NGRC', 'P-RC', 'PLB-NGRC'}, 'Orientation', 'horizontal', 'FontSize', fs, 'Box', 'off');
    lgd.Layout.Tile = 'north';
end
exportgraphics(fig, fullfile(img_dir, 'figure8.png'), 'Resolution', 600);

function Tick_align(ax, is_log)
    yl = ylim(ax);
    if is_log
        lo = 10^floor(log10(yl(1)) + 1e-9);
        hi = 10^ceil(log10(yl(2)) - 1e-9);
        t = 10.^(round(log10(lo)):round(log10(hi)));
    else
        steps = kron(10.^(-6:3), [1, 2, 5]);
        for st = steps
            lo = st * floor(yl(1) / st + 1e-9);
            hi = st * ceil(yl(2) / st - 1e-9);
            if round((hi - lo) / st) <= 3
                break
            end
        end
        t = lo:st:hi;
    end
    ylim(ax, [lo, hi]);
    yticks(ax, t);
    if is_log
        yticklabels(ax, arrayfun(@(e) sprintf('10^{%d}', e), round(log10(t)), 'UniformOutput', false));
    end
    set(ax, 'YMinorTick', 'off', 'XMinorTick', 'off', 'YMinorGrid', 'off', 'XMinorGrid', 'off');
end
