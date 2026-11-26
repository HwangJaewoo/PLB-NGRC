%% Drives the reservoirs of parallel RC with the true states Sw before an autonomous prediction
function R = PRC_synchronize(mdl, Sw)
    N = size(Sw, 2) / mdl.m;
    R = zeros(N, mdl.Dr);
    for t = 1:size(Sw, 1)
        R = PRC_step(mdl, R, Sw(t, :));
    end
end
