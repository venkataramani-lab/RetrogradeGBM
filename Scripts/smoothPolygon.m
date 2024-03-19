function smoothedPoints = smoothPolygon(points)
    smoothedPoints = points;
    alpha = 0.8; % can be adjusted

    for i = 1:size(points, 1)
        smoothedPoints(i, :) = (1 - alpha) * points(i, :) + alpha * (points(mod(i, size(points, 1)) + 1, :) + points(mod(i - 2, size(points, 1)) + 1, :)) / 2;
    end
end