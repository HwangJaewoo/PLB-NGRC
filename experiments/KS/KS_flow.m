%% n_sub ETDRK4 steps of the Kuramoto-Sivashinsky equation (coefficients P from KS_data); u is a row vector in physical space
function u = KS_flow(u, P, n_sub)
    v = fft(u(:));
    for j = 1:n_sub
        Nv = P.g .* fft(real(ifft(v)).^2);
        a = P.E2 .* v + P.Qc .* Nv;
        Na = P.g .* fft(real(ifft(a)).^2);
        b = P.E2 .* v + P.Qc .* Na;
        Nb = P.g .* fft(real(ifft(b)).^2);
        c = P.E2 .* a + P.Qc .* (2 * Nb - Nv);
        Nc = P.g .* fft(real(ifft(c)).^2);
        v = P.E .* v + Nv .* P.f1 + 2 * (Na + Nb) .* P.f2 + Nc .* P.f3;
        v = fft(real(ifft(v)));
    end
    u = real(ifft(v)).';
end
