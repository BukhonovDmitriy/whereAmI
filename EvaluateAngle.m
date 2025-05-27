function angle = EvaluateAngle(points)
    angle = atan2(points(2,1), points(1,1));
%return
    sz = size(points);
    N = sz(2);
    angles = zeros(N);

    for i = 1:N
        angles(i) = atan2(points(2,i), points(1,i)) - (i-1) * 2*pi/N;
        if angles(i) < -pi
            angles(i) = angles(i) + 2 * pi;
        end
    end

    angle = sum(angles, "all") / N;
end