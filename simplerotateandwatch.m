function [var, max_err]=rotateandwatch(options)
    arguments
        options.R   double = 1
        options.ang_amplitude double = pi/12
        options.N   double = 8
        options.Kx  double = 1.1
        options.Ky  double = 0.9
        options.dxy double = 0.05
        options.x0  double = 0.1
        options.y0  double = 0.15
        options.noise_sigma double = 0.03
        options.iters int32 = 1000
        options.ang_freq int32 = 100
    end
    R   = options.R;
    ang_amplitude = options.ang_amplitude;
    N   = options.N;
    Kx  = options.Kx;
    Ky  = options.Ky;
    dxy = options.dxy;
    x0  = options.x0;
    y0  = options.y0;
    noise_sigma = options.noise_sigma;
    iters = options.iters;
    ang_freq = options.ang_freq;

    diffs = zeros(iters, ang_freq * 2 + 1);
    fprintf("####################\n");
    for i = 1:iters
        for j = -ang_freq:ang_freq
            ang = ang_amplitude * (double(j) / double(ang_freq));
            mp = Generate(R, ang, N, Kx, Ky, dxy, x0, y0, noise_sigma);
            angle_eval = ApproximateAngle(R, mp);
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

% здесь никаких эллипсов, просто сдвигаем. R для похожести на функцию из
% rotateandwatch.m
function recovered_points = Approximate(R, moved_points)
    sz = size(moved_points);
    N = sz(2);
    
    offset = sum(moved_points, 2) / N;
    
    recovered_points = moved_points - offset .* ones(sz);
end

function recovered_angle = ApproximateAngle(R, moved_points)
    recovered_points = Approximate(R, moved_points);
    recovered_angle = EvaluateAngle(recovered_points);
end
