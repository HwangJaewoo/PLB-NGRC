%% Lorenz-96: true and PLB-NGRC time series of one site and trajectories projected onto (x1, x2) (output figure2.png)
clear; clc; close all;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'core'), genpath(fullfile(root, 'experiments')));
img_dir = fullfile(root, 'figures', 'img');
if ~isfolder(img_dir)
    mkdir(img_dir);
end

load(fullfile(root, 'experiments', 'Lorenz96', 'data', 'Lorenz96_N020_data.mat'), 'S', 'dt', 'i_test');
B = load(fullfile(root, 'experiments', 'Lorenz96', 'results', 'Lorenz96_N020_PLBNGRC.mat'), 'res');
X = S(i_test, :);
P = B.res.P_long;

T_a = 10;
i0 = 10;
j1 = 1;
j2 = 2;
fig_w = 17.8;
fig_h = 5.6;
fs = 9;
lw_ts = 0.8;
lw_pp = 0.3;
c_plb = [0.85, 0, 0];

n_a = round(T_a / dt);
t_a = (0:n_a) * dt;
n_pp = min(size(X, 1), size(P, 1));
v = [X(1:n_pp, [j1, j2]); P(1:n_pp, [j1, j2])];
lim = 3;
tk = [-lim, 0, lim];

f = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, fig_w, fig_h]);

ax_a = axes(f);
plot(t_a, X(1:n_a + 1, i0), '-', 'Color', 'k', 'LineWidth', lw_ts);
hold on;
plot(t_a, P(1:n_a + 1, i0), '-', 'Color', c_plb, 'LineWidth', lw_ts);
xlim([0, T_a]);
xticks(0:2:T_a);
ya = [X(1:n_a + 1, i0); P(1:n_a + 1, i0)];
ylim([-lim, lim]);
yticks(tk);
xlabel('$t$', 'Interpreter', 'latex');
ylabel(sprintf('$x_{%d}$', i0), 'Interpreter', 'latex');
Label('(a)');
box on;
grid on;

ax_b = axes(f);
plot(X(1:n_pp, j1), X(1:n_pp, j2), '-', 'Color', 'k', 'LineWidth', lw_pp);
Pp_style('(b)', j1, j2, lim, tk);

ax_c = axes(f);
plot(P(1:n_pp, j1), P(1:n_pp, j2), '-', 'Color', c_plb, 'LineWidth', lw_pp);
Pp_style('(c)', j1, j2, lim, tk);

set(findall(f, '-property', 'FontSize'), 'FontSize', fs);
for ax = [ax_a, ax_b, ax_c]
    set(ax.Title, 'FontSize', fs + 2, 'FontWeight', 'bold');
end

bot = 0.19;
h = 0.70;
w_pp = h * fig_h / fig_w;
x_a = 0.06; w_a = 0.33;
gap = 0.075;
x_b = x_a + w_a + gap;
x_c = x_b + w_pp + gap;
set(ax_a, 'Position', [x_a, bot, w_a, h]);
set(ax_b, 'Position', [x_b, bot, w_pp, h]);
set(ax_c, 'Position', [x_c, bot, w_pp, h]);

exportgraphics(f, fullfile(img_dir, 'figure2.png'), 'Resolution', 600);

function Pp_style(label, j1, j2, lim, tk)
    xlim([-lim, lim]);
    ylim([-lim, lim]);
    xticks(tk);
    yticks(tk);
    xlabel(sprintf('$x_{%d}$', j1), 'Interpreter', 'latex');
    ylabel(sprintf('$x_{%d}$', j2), 'Interpreter', 'latex');
    Label(label);
    box on;
    grid on;
end

function Label(s)
    title(s, 'Units', 'normalized', 'Position', [0, 1.02], 'HorizontalAlignment', 'left');
end
