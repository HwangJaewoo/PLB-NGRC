%% Extended Thomas system with N sites integrated by n_sub fourth-order Runge-Kutta steps of size h
function x = Thomas_flow(x, b, h, n_sub)
    for j = 1:n_sub
        k1 = Rhs(x, b);
        k2 = Rhs(x + 0.5 * h * k1, b);
        k3 = Rhs(x + 0.5 * h * k2, b);
        k4 = Rhs(x + h * k3, b);
        x = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    end
end

function dx = Rhs(x, b)
    dx = sin(circshift(x, -1)) - b * x;
end
