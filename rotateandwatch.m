function [var, max_err]=rotateandwatch(options)
    arguments
        options.R   double = 1
        options.ang_amplitude double = pi/12
        options.N   double = 8
        options.Kx  double = 1.1
        options.Ky  double = 0.9
        options.dxy double = 0.05
        options.zero_offset_var double = 0.1 * pi/180  % 0.1 градусов в секунду * радиан в градусах
        options.noise_var double = 0.03
        options.iters int32 = 1000
        options.ang_freq int32 = 100
    end
    R   = options.R;
    ang_amplitude = options.ang_amplitude;
    N   = options.N;
    Kx  = options.Kx;
    Ky  = options.Ky;
    dxy = options.dxy;
    zero_offset_var = options.zero_offset_var;
    noise_var = options.noise_var;
    iters = options.iters;
    ang_freq = options.ang_freq;

    diffs = zeros(iters, ang_freq * 2 + 1);
    fprintf("####################\n");
    for i = 1:iters
        for j = -ang_freq:ang_freq
            ang = ang_amplitude * (double(j) / double(ang_freq));
            mp = Generate(R, ang, N, Kx, Ky, dxy, zero_offset_var, noise_var);
            angle_eval = ApproximateAngle(mp);
            diffs(i, j + 1001) = abs(angle_eval - ang);
        end
        
        if mod (i, iters/20) == 0
            fprintf(".");
        end
    end
    fprintf("\n");
    
    var = sum(diffs.^2, "all") / (double(iters) * double(ang_freq * 2 + 1) - 1);
    max_err = max(diffs, [], "all");
end

function recovered_points = Approximate(moved_points)
    sz = size(moved_points);
    N = sz(2);

    K = [moved_points(2,1:end).^2;...
        moved_points(1, 1:end) .* moved_points(2, 1:end);...
        moved_points(1,1:end);...
        moved_points(2,1:end);...
        ones(1,N)]';

    V = (K' * K)^-1 * K' * (-moved_points(1, 1:end)'.^2);

    [x0, y0] = approx_x0y0(V);
    Q = approx_Q(V, x0, y0);
    M = approx_M(Q);

    recovered_points = M^-1 * (moved_points - [ones(1, N) * x0; ones(1, N) * y0]);
end

function recovered_angle = ApproximateAngle(moved_points)
    recovered_points = Approximate(moved_points);
    recovered_angle = EvaluateAngle(recovered_points);
end

function [x0, y0] = approx_x0y0(V)
    A = V(1);
    B = V(2);
    C = V(3);
    D = V(4);

    x0 = (B*D - 2*A*C) / (4*A-B^2);
    y0 = (B*C - 2*D) / (4*A-B^2);
end

function Q = approx_Q(V,x0, y0) %будем "возвращать" точки на единичную окружность
    A = V(1);
    B = V(2);
    E = V(5);
    Q11 = - 1 / (E - x0^2 - A*y0^2 - B*x0*y0);

    Q = [Q11 B*Q11/2; B*Q11/2 A*Q11];
end

function M = approx_M(Q)
    O = Q^-1;
    Kx = ((O(1,1)*O(2,2) + (O(1,1)*O(2,2)*(O(1,1)*O(2,2) - O(1,2)^2))^0.5)...
        / (2*O(2,2)))^0.5;
    Ky = ((O(1,1)*O(2,2) + (O(1,1)*O(2,2)*(O(1,1)*O(2,2) - O(1,2)^2))^0.5)...
        / (2*O(1,1)))^0.5;
    dxy = O(1,2) / (2 * Kx * Ky);
    M = [Kx 0; 0 Ky] * [1 dxy; dxy 1];
end