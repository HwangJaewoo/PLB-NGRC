%% Lorenz-96 system integrated by n_sub fourth-order Runge-Kutta steps of size h
function x = Lorenz96_flow(x, h, n_sub, F)
    for j = 1:n_sub
        k1 = Rhs(x, F);
        k2 = Rhs(x + 0.5 * h * k1, F);
        k3 = Rhs(x + 0.5 * h * k2, F);
        k4 = Rhs(x + h * k3, F);
        x = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    end
end

function dx = Rhs(x, F)
    dx = (circshift(x, -1) - circshift(x, 2)) .* circshift(x, 1) - x + F;
end
