function distance = pointToPolylineDistance(point, polyline)
    % Calculate the distance between a point and a polyline

    % Check if the polyline is closed (if so, remove the last point)
    if size(polyline, 1) > 1 && isequal(polyline(1, :), polyline(end, :))
        polyline = polyline(1:end-1, :);
    end

    % Check if there are enough points in the polyline
    if size(polyline, 1) < 2
        error('Polyline must have at least two points.');
    end

    % Calculate the distance to the nearest point in the polyline
    distances = sqrt((point(1) - polyline(:, 1)).^2 + (point(2) - polyline(:, 2)).^2);
    distance = min(distances);
end