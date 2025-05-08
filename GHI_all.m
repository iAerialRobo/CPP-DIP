
smallRadius = 50; 
largeRadius = 175; 
minDistance = 80; 
nearRadius = largeRadius - smallRadius;

template = imread('xxx.JPG');
figure;
imshow(template);
hold on; 

filename = 'palm_tree_coordinates.txt';
fileID = fopen(filename, 'r');

headerLine = fgetl(fileID); 
dims = sscanf(headerLine, 'Image Dimensions: Width=%d, Height=%d');
areaWidth = dims(1);
areaHeight = dims(2);

fgetl(fileID); 
data = textscan(fileID, '%f %f', 'Delimiter', ',', 'CollectOutput', true);
points = data{1};
fclose(fileID);

circleCenters = []; 
circleCount = 1; 

covered = false(1, numPoints); 
startPoint = points(1, :); 

while ~all(covered)
    
    uncoveredPoints = points(~covered, :);
    distances = zeros(size(uncoveredPoints, 1), 1);
    
    for i = 1:size(uncoveredPoints, 1)
        distances(i) = norm(uncoveredPoints(i, :) - startPoint);
    end
    
    [~, idx] = min(distances);
    nearestPoint = uncoveredPoints(idx, :);
   
    bestCenter = nearestPoint;
    maxCovered = 0;
    bestCoveredIndices = [];
nearRadius = largeRadius - smallRadius;
    
    for dx = -nearRadius:1:nearRadius
        for dy = -nearRadius:1:nearRadius
            testCenter = nearestPoint + [dx, dy];
            distances = vecnorm(points - testCenter, 2, 2);
            coveredIndices = find(distances + smallRadius <= largeRadius);
            numCovered = sum(~covered(coveredIndices)); 
         
            if numCovered > maxCovered && norm(testCenter - nearestPoint) <= largeRadius
                bestCenter = testCenter;
                maxCovered = numCovered;
                bestCoveredIndices = coveredIndices;
            end
        end
    end
    
    if ~isempty(bestCoveredIndices)
        bestCenter = mean(points(bestCoveredIndices, :), 1); 
    end

        distancesToCenter = vecnorm(points(bestCoveredIndices, :) - bestCenter, 2, 2);
    maxDistanceToSmallCircle = max(distancesToCenter);
    feasibleRadius = largeRadius - maxDistanceToSmallCircle - smallRadius;
    
    circleCenters = [circleCenters; bestCenter, feasibleRadius]; 
    
    theta = linspace(0, 2*pi, 100);
    x = bestCenter(1) + largeRadius * cos(theta);
    y = bestCenter(2) + largeRadius * sin(theta);
    plot(x, y, 'r-', 'LineWidth', 0.5); 
       text(bestCenter(1), bestCenter(2), num2str(circleCount), 'Color', 'g', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10);

    
    covered(bestCoveredIndices) = true;

    
    startPoint = bestCenter;

    
    circleCount = circleCount + 1;
end
hold on; 

function tsp = solveTSP(circleCenters)  
    
    n = size(circleCenters, 1);
    start_idx = randi(n);  
    path = start_idx;  
    remaining_indices = setdiff(1:n, path);  

    nearest_points_idx = nearestTwoPoints(circleCenters, remaining_indices, start_idx);
    path = [path, nearest_points_idx];  
    remaining_indices = setdiff(remaining_indices, nearest_points_idx);  

    while ~isempty(remaining_indices) 
        new_point_idx = remaining_indices(randi(length(remaining_indices)));
        remaining_indices = setdiff(remaining_indices, new_point_idx);   
        path = insertPoint(circleCenters, path, new_point_idx);
    end
    
    tsp = path;  
end

function nearest_points_idx = nearestTwoPoints(circleCenters, remaining_indices, first_point_idx)
    
    distances = zeros(length(remaining_indices), 1);
    for i = 1:length(remaining_indices)
        current_point_idx = remaining_indices(i);
        distances(i) = norm(circleCenters(first_point_idx, :) - circleCenters(current_point_idx, :));
    end
    
    [~, idx] = sort(distances);  
    nearest_points_idx = remaining_indices(idx(1:2));  
end

function path = insertPoint(circleCenters, path, new_point_idx)
    
    distance_changes = [];
    angle_changes = [];
    intersection_counts = [];
    
    for i = 1:length(path)
        
        next_idx = mod(i, length(path)) + 1;         
        
        distance_change = calculateDistanceChange(circleCenters, path(i), path(next_idx), new_point_idx);
        
        angle_change = calculateAngleChange(circleCenters, path(i), path(next_idx), new_point_idx);
        intersection_count = countIntersections(circleCenters, path, path(i), path(next_idx), new_point_idx);
               
        distance_changes = [distance_changes, distance_change];
        angle_changes = [angle_changes, angle_change];
        intersection_counts = [intersection_counts, intersection_count];
       
    end
      
    normalized_distances = (distance_changes - min(distance_changes)) / (max(distance_changes) - min(distance_changes));
    normalized_angles = (angle_changes - min(angle_changes)) / (max(angle_changes) - min(angle_changes));      
    combined_scores = 0.3*normalized_distances + 0.7*normalized_angles + intersection_counts;  
    [~, min_idx] = min(combined_scores);
    next_idx = mod(min_idx, length(path)) + 1;  
    best_path = [path(1:min_idx), new_point_idx, path(next_idx:end)]; 
    path = best_path;
end

function distance_change = calculateDistanceChange(circleCenters, p1_idx, p2_idx, new_point_idx) 
    p1 = circleCenters(p1_idx, :);
    p2 = circleCenters(p2_idx, :);
    new_point = circleCenters(new_point_idx, :);
    original_distance = norm(p1 - p2);  
    new_distance1 = norm(p1 - new_point);
    new_distance2 = norm(p2 - new_point);   
    distance_change = (new_distance1 + new_distance2) - original_distance;
end

function angle_change = calculateAngleChange(circleCenters, p1_idx, p2_idx, new_point_idx)
    p1 = circleCenters(p1_idx, :);
    p2 = circleCenters(p2_idx, :);
    new_point = circleCenters(new_point_idx, :);   
    angle_change = calculateAngle(p1, new_point, p2);
end

function count = countIntersections(circleCenters, path, p1_idx, p2_idx, new_point_idx)
    count = 0;
    new_segment = [circleCenters(p1_idx, :); circleCenters(new_point_idx, :)]; 
    new_segment2 = [circleCenters(new_point_idx, :); circleCenters(p2_idx, :)]; 
    
    
    for i = 1:length(path)
        next_idx = mod(i, length(path)) + 1; 
        seg = [circleCenters(path(i), :); circleCenters(path(next_idx), :)];
              
        if any([p1_idx, p2_idx] == path(i)) || any([p1_idx, p2_idx] == path(next_idx))
            continue;
        end
                
        if isIntersecting(new_segment, seg) || isIntersecting(new_segment2, seg)
            count = count + 1;
        end
    end
end

function intersect = isIntersecting(seg1, seg2)
    A = seg1(1, :); B = seg1(2, :);
    C = seg2(1, :); D = seg2(2, :);   
    d1 = sign((B(1) - A(1)) * (C(2) - A(2)) - (B(2) - A(2)) * (C(1) - A(1)));
    d2 = sign((B(1) - A(1)) * (D(2) - A(2)) - (B(2) - A(2)) * (D(1) - A(1)));
    d3 = sign((D(1) - C(1)) * (A(2) - C(2)) - (D(2) - C(2)) * (A(1) - C(1)));
    d4 = sign((D(1) - C(1)) * (B(2) - C(2)) - (D(2) - C(2)) * (B(1) - C(1))); 
    intersect = (d1 * d2 < 0) && (d3 * d4 < 0);

end

function dist = totalDistance(circleCenters, path)    
    dist = 0;
    for i = 1:length(path)
        next_idx = mod(i, length(path)) + 1;  
        dist = dist + norm(circleCenters(path(i), :) - circleCenters(path(next_idx), :));  
    end
end

function angle = calculateAngle(p1, p2, p3)
    v1 = p2 - p1;
    v2 = p3 - p2;
    cosTheta = dot(v1, v2) / (norm(v1) * norm(v2) + eps); 
    angle = acosd(cosTheta); 
    
end

bestRoute = solveTSP(circleCenters);
segmentLengths = [];  

for i = 1:length(bestRoute)-1
    startCenter = circleCenters(bestRoute(i), 1:2);
    endCenter = circleCenters(bestRoute(i+1), 1:2);
    segmentLength = norm(startCenter - endCenter);  
    
    segmentLengths = [segmentLengths, segmentLength];  
end
totalPathLength = sum(segmentLengths);
totalAngle = 0;
for i = 2:(length(bestRoute) - 1)
    
    p1 = circleCenters(bestRoute(i - 1), :);
    p2 = circleCenters(bestRoute(i), :);
    p3 = circleCenters(bestRoute(i + 1), :);  
    angle = calculateAngle(p1, p2, p3);
    totalAngle = totalAngle + angle;
end
disp([ ' Best Route Length: ', num2str(totalPathLength), ...
          ' Cumulative Angle: ', num2str(totalAngle)]);
for i = 1:(size(bestRoute, 2) - 1)
    plot([circleCenters(bestRoute(i), 1), circleCenters(bestRoute(i+1), 1)], ...
         [circleCenters(bestRoute(i), 2), circleCenters(bestRoute(i+1), 2)], 'r-', 'LineWidth', 2);
end
plot([circleCenters(bestRoute(end), 1), circleCenters(bestRoute(1), 1)], ...
     [circleCenters(bestRoute(end), 2), circleCenters(bestRoute(1), 2)], 'r-', 'LineWidth', 2);
title('Optimal Path with Ant Colony Optimization');
hold off;


