%% Readout of parallel RC shared by all sites (every second reservoir node squared)
function x = PRC_readout(mdl, R)
    Q = R;
    Q(:, 2:2:end) = R(:, 2:2:end).^2;
    Y = Q * mdl.W;
    x = reshape(Y, 1, []);
end
