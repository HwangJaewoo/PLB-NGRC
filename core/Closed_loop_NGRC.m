%% Autonomous prediction of T steps from the history H; td is the step of divergence (NaN if none)
function [P, td] = Closed_loop_NGRC(mdl, H, T, thr)
    P = nan(T, size(H, 2));
    td = NaN;
    for t = 1:T
        y = Predict_NGRC(mdl, Delay_input(H, mdl.nn, mdl.m, mdl.k));
        x = reshape(y, 1, []);
        if any(~isfinite(x)) || max(abs(x)) > thr
            td = t;
            return
        end
        P(t, :) = x;
        H = [H(2:end, :); x];
    end
end
