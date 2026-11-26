%% Largest Lyapunov exponent per step of the map step_fun (Benettin method with central differences)
function lam = Largest_lyapunov(step_fun, x, n_trans, n_eval)
    eps_fd = 1e-7;
    for t = 1:n_trans
        x = step_fun(x);
    end
    v = randn(size(x));
    v = v / norm(v(:));
    lyap_sum = 0;
    for t = 1:n_eval
        xp = step_fun(x + eps_fd * v);
        xm = step_fun(x - eps_fd * v);
        w = (xp - xm) / (2 * eps_fd);
        nw = norm(w(:));
        lyap_sum = lyap_sum + log(nw);
        v = w / nw;
        x = step_fun(x);
        if any(~isfinite(x(:)))
            lam = NaN;
            return
        end
    end
    lam = lyap_sum / n_eval;
end
