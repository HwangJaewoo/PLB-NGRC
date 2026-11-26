%% Autonomous prediction of parallel RC for T steps after synchronization with Sw; td is the step of divergence (NaN if none)
function [P, td] = Closed_loop_PRC(mdl, Sw, T, thr)
    R = PRC_synchronize(mdl, Sw);
    P = nan(T, size(Sw, 2));
    td = NaN;
    for t = 1:T
        x = PRC_readout(mdl, R);
        if any(~isfinite(x)) || max(abs(x)) > thr
            td = t;
            return
        end
        P(t, :) = x;
        R = PRC_step(mdl, R, x);
    end
end
