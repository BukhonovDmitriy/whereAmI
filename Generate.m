function moved_points = Generate(R, angle_offset, N, Kx, Ky, dxy, x0, y0, noise_sigma)
    M = [Kx 0; 0 Ky] * [1 dxy; dxy 1];

    points = R * [
        cos((0:N-1) * 2 * pi / N + angle_offset);
        sin((0:N-1) * 2 * pi / N + angle_offset)];

    moved_points = M * points + ...
        [ones(1, N) * x0; ones(1, N) * y0] + ...
        normrnd(0, noise_sigma, 2, N);
end