%% One step of the coupled Ikeda map; x = [real parts, imaginary parts] of the N sites
function x = Ikeda_map(x, N, u, eps_c)
    z = x(1:N) + 1i * x(N + 1:2 * N);
    g = 1 + u * z .* exp(1i * (0.4 - 6 ./ (1 + abs(z).^2)));
    z = (1 - eps_c) * g + 0.5 * eps_c * (circshift(g, 1) + circshift(g, -1));
    x = [real(z), imag(z)];
end
