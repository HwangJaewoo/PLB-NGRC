%% Training of NGRC: ball tree partition of the local states and one ridge-regression polynomial per ball (weighted or in-ball fitting)
function mdl = Train_NGRC(S, p, m)
    if istable(p)
        p = table2struct(p);
    end
    if ~isfield(p, 'k')
        p.k = 1;
    end
    if ~isfield(p, 's')
        p.s = 0;
    end
    if ~isfield(p, 'eta')
        p.eta = 0.3;
    end
    if ~isfield(p, 'ord')
        p.ord = 2;
    end
    if ~isfield(p, 'fit')
        p.fit = 'weighted';
    end
    if ~any(strcmp(p.fit, {'weighted', 'inball'}))
        error('p.fit must be weighted or inball');
    end
    lam = 10^p.log_lam;

    U = Delay_input(S(1:end - 1, :), p.nn, m, p.k);
    N = size(S, 2) / m;
    Y = zeros(size(U, 1), m);
    for v = 1:m
        Y(:, v) = reshape(S(p.k + 1:end, (v - 1) * N + (1:N)), [], 1);
    end

    d = size(U, 2);
    if p.s == 0
        C = zeros(1, d);
        sig = 1;
        lab = ones(size(U, 1), 1);
    else
        [C, Rb, lab] = Ball_tree(U, p.s);
        sig = max(p.eta * Rb, 1e-6);
    end
    B = size(C, 1);

    nf = 1 + d + d * (d + 1) / 2;
    if p.ord >= 3
        nf = nf + d * (d + 1) * (d + 2) / 6;
    end
    W = zeros(nf, m, B);
    nrow = size(U, 1);
    max_el = 2e8;
    if nrow * max(B, nf) <= max_el
        Pi = Blend_weights(U, C, sig);
        for b = 1:B
            if strcmp(p.fit, 'inball')
                sel = lab == b;
                wb = ones(nnz(sel), 1);
            else
                w = Pi(:, b);
                sel = w > 1e-3;
                wb = w(sel);
            end
            O = Poly_features(U(sel, :) - C(b, :), p.ord);
            Ow = O .* wb;
            sw = sum(wb);
            W(:, :, b) = (O' * Ow / sw + lam * eye(nf)) \ (Ow' * Y(sel, :) / sw);
        end
    else
        G = zeros(nf, nf, B);
        Bm = zeros(nf, m, B);
        sw = zeros(B, 1);
        cs = max(1, floor(2e7 / max(B, nf)));
        for c0 = 1:cs:nrow
            ii = (c0:min(c0 + cs - 1, nrow))';
            if strcmp(p.fit, 'inball')
                Pc = double(lab(ii) == (1:B));
            else
                Pc = Blend_weights(U(ii, :), C, sig);
                Pc(Pc <= 1e-3) = 0;
            end
            for b = 1:B
                w = Pc(:, b);
                sel = w > 0;
                if ~any(sel)
                    continue
                end
                O = Poly_features(U(ii(sel), :) - C(b, :), p.ord);
                Ow = O .* w(sel);
                G(:, :, b) = G(:, :, b) + O' * Ow;
                Bm(:, :, b) = Bm(:, :, b) + Ow' * Y(ii(sel), :);
                sw(b) = sw(b) + sum(w(sel));
            end
        end
        for b = 1:B
            W(:, :, b) = (G(:, :, b) / sw(b) + lam * eye(nf)) \ (Bm(:, :, b) / sw(b));
        end
    end

    mdl.nn = p.nn;
    mdl.k = p.k;
    mdl.m = m;
    mdl.s = p.s;
    mdl.ord = p.ord;
    mdl.fit = p.fit;
    mdl.C = C;
    mdl.sig = sig;
    mdl.W = W;
end

function [C, R, lab] = Ball_tree(U, s)
    if s == 0
        C = mean(U, 1);
        R = max(sqrt(sum((U - C).^2, 2)));
        lab = ones(size(U, 1), 1);
        return
    end
    [~, a] = max(var(U, 0, 1));
    [~, order] = sort(U(:, a));
    h = floor(size(U, 1) / 2);
    [C1, R1, l1] = Ball_tree(U(order(1:h), :), s - 1);
    [C2, R2, l2] = Ball_tree(U(order(h + 1:end), :), s - 1);
    C = [C1; C2];
    R = [R1; R2];
    lab = zeros(size(U, 1), 1);
    lab(order(1:h)) = l1;
    lab(order(h + 1:end)) = l2 + size(C1, 1);
end
