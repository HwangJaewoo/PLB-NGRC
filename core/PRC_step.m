%% One reservoir update of all sites of parallel RC driven by the state x
function R = PRC_step(mdl, R, x)
    U = Delay_input(x, mdl.nn, mdl.m, 1);
    R = (1 - mdl.leak) * R + mdl.leak * tanh(R * mdl.A_t + U(:, mdl.idx) .* mdl.w);
end
