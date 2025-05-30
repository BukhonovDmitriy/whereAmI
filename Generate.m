function moved_points = Generate(R, angle_offset, N, Kx, Ky, dxy, zero_offset_var, noise_var)
    x0 = normrnd(0, zero_offset_var);
    y0 = normrnd(0, zero_offset_var);
    M = [Kx 0; 0 Ky] * [1 dxy; dxy 1];

    points = R * [
        cos((0:N-1) * 2 * pi / N + angle_offset);
        sin((0:N-1) * 2 * pi / N + angle_offset)];

    moved_points = M * points + ...
        [ones(1, N) * x0; ones(1, N) * y0] + ...
        normrnd(0, noise_var, 2, N);
end