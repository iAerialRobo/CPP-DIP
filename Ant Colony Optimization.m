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

numPoints = size(points, 1); 
for i = 1:numPoints
theta = linspace(0, 2*pi, 100);
x = points(i, 1) + smallRadius * cos(theta);
y = points(i, 2) + smallRadius * sin(theta);
plot(x, y, 'b-'); 
plot(points(i, 1), points(i, 2), 'bo', 'MarkerFaceColor', 'b');
end


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

    
    x_feasible = bestCenter(1) + feasibleRadius * cos(theta);
    y_feasible = bestCenter(2) + feasibleRadius * sin(theta);
    fill(x_feasible, y_feasible, 'g', 'FaceAlpha', 0.3); 
    plot(bestCenter(1), bestCenter(2), 'go', 'MarkerFaceColor', 'g'); 

    
    text(bestCenter(1), bestCenter(2), num2str(circleCount), 'Color', 'g', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10);

    
    covered(bestCoveredIndices) = true;

    
    startPoint = bestCenter;

    
    circleCount = circleCount + 1;
end
hold on; 


numAnts = 100;      
numIterations = 1000; 
alpha = 1;         
beta = 1;          
gamma = 1;         
delta = 4;       
evaporationRate = 0.3; 
Q = 1000;            
numCenters = size(circleCenters, 1); 

distMatrix = squareform(pdist(circleCenters));

minDist = min(distMatrix(:)); 
maxDist = max(distMatrix(:)); 

normalizedDistMatrix = (distMatrix - minDist) / (maxDist - minDist);
numPoints = size(circleCenters, 1);  
rewardHistory = zeros(1, numIterations); 
distanceHistory = zeros(1, numIterations);
angHistory = zeros(1, numIterations); 
internHistory = zeros(1, numIterations); 


pheromone = ones(numPoints, numPoints);


function angle = calculateAngle(p1, p2, p3)
    v1 = p2 - p1;
    v2 = p3 - p2;
    cosTheta = dot(v1, v2) / (norm(v1) * norm(v2) + eps); 
    angle = acosd(cosTheta); 
    
end


distances = zeros(numCenters, numCenters);
for i = 1:numCenters
    for j = 1:numCenters
        if i ~= j
            distances(i, j) = norm(circleCenters(i, 1:2) - circleCenters(j, 1:2));
        else
            distances(i, j) = Inf; 
        end
    end
end

function isIntersecting = checkLineSegmentIntersection(p1, q1, p2, q2, circleCenters)
    
    d1 = direction(p1, q1, p2, circleCenters);
    d2 = direction(p1, q1, q2, circleCenters);
    d3 = direction(p2, q2, p1, circleCenters);
    d4 = direction(p2, q2, q1, circleCenters);
    
    
    if (d1 * d2 < 0) && (d3 * d4 < 0)
        isIntersecting = true;
    else
        isIntersecting = false;
    end
end

function d = direction(p, q, r, circleCenters)
    

    d = (circleCenters(q, 1) - circleCenters(p, 1)) * (circleCenters(r, 2) - circleCenters(p, 2)) - (circleCenters(r, 1) - circleCenters(p, 1)) * (circleCenters(q, 2) - circleCenters(p, 2));
end


bestRoute = [];

numPoints = size(circleCenters,1);
epsilon_max = 1.0; 
epsilon_min = 0.001;
epsilon_max1 = 0.9; 
epsilon_min1 = 0.1;
targetEpsilon = 0.2;
lambda = -log(targetEpsilon) / (0.80 * numIterations);
lambda1 = -log(targetEpsilon) / (0.80 * numIterations);
minLength = Inf;
minAngle = Inf;
minIntern = Inf;
bestRouteLength = inf;
bestAngle = inf;
bestIntern = inf;
for iteration = 1:numIterations
    fprintf('Iteration: %.2f\n', iteration);
    routes = zeros(numAnts, numPoints+1); 
    routeLengths = zeros(numAnts, 1);  
    routeRewards = zeros(numAnts, 1);  
    angleRewards = zeros(numAnts, 1);
    internRewards = zeros(numAnts, 1);
    routesQuality = zeros(numAnts, 1);
    cumulativeAngles = zeros(numAnts, 1); 
    cumulativeIntern = zeros(numAnts, 1); 
    total_pheromone = 0;
    epsilon = max(epsilon_min, epsilon_max * exp(-lambda * iteration)); 
    evaporationRate = max(epsilon_min1, epsilon_max1 * exp(-lambda1 * iteration));
    for ant = 1:numAnts
        counts = 0;
        
        currentPoint = randi(numPoints);
        visited = false(1, numPoints);
        visited(currentPoint) = true;
        route = [];
        route = [route, currentPoint];  
        
        
        for step = 2:numPoints
            probabilities = zeros(1, numPoints);
            totalProbability = 0;
            
            
            for nextPoint = 1:numPoints
                if ~visited(nextPoint)  
                    pheromoneVal = pheromone(currentPoint, nextPoint)^alpha;
                    distanceVal = (1/normalizedDistMatrix(currentPoint, nextPoint))^beta;

                    
                    if step > 2 && step < numPoints
                        prevPoint = route(end-1); 
                        angle = calculateAngle(circleCenters(prevPoint, :), ...
                                               circleCenters(currentPoint, :), ...
                                               circleCenters(nextPoint, :));
                        anglePenalty = (1/(1+angle))^gamma;

                    else
                        anglePenalty = 1; 
                    end
                    Intercount = 0;
                    if size(route, 2) > 2
                        latestStart = route(end); 
                        latestEnd = nextPoint;   
                        for k = 1:step-2  
                            if k == 1 && step == size(route, 1)
                                continue;
                            end           
                            
                            existingStart = route(k);
                            existingEnd = route(k+1);

                            
                            if checkLineSegmentIntersection(latestStart, latestEnd, existingStart, existingEnd, circleCenters)
                                Intercount = Intercount + 1 ;
                            end
                        end
                        InterReward = (1/(Intercount+1))^delta;
                        probabilities(nextPoint) = pheromoneVal*distanceVal*anglePenalty*InterReward;
                    else
                        probabilities(nextPoint) = pheromoneVal*distanceVal*anglePenalty;
                    end
                    
                    
 
                    totalProbability = totalProbability + probabilities(nextPoint);
                end
            end
            
            
            if totalProbability > 0
                
                testprobabilities = probabilities;
                probabilities = probabilities / totalProbability;
            else
                
                probabilities(:) = 1; 
                probabilities(visited) = 0; 
                probabilities = probabilities / numPoints; 
            end
            

            if rand < epsilon            
                
                nonZeroIndices = find(probabilities > 0);
                
                nextPoint = nonZeroIndices(randi(length(nonZeroIndices)));
            else
                
                [~, nextPoint] = max(probabilities);
            end

            
            route = [route, nextPoint];  
            visited(nextPoint) = true;
            if step > 1
                startCenter = circleCenters(route(end-1),1:2);
                endCenter = circleCenters(nextPoint,1:2);
                segmentLength = norm(startCenter - endCenter);
                
                routeReward = (1/normalizedDistMatrix(route(end), route(end-1)))^beta;
                routeLengths(ant) = routeLengths(ant) + segmentLength;
                routeRewards(ant) = routeRewards(ant) + routeReward;
                
            end
            if step > 2 && step < numPoints
                prevPoint = route(end-2); 
                segmentAngle = calculateAngle(circleCenters(prevPoint, :), ...
                                       circleCenters(currentPoint, :), ...
                                       circleCenters(nextPoint, :));
                cumulativeAngles(ant) = cumulativeAngles(ant) + abs(segmentAngle); 
                angleRewards(ant) = angleRewards(ant) + 3*(1/((1+segmentAngle)/180))^gamma;
                
            end
            
            if step> 3
                latestStart = route(end-1); 
                latestEnd = route(end);   
                segmentIntern = 0;
                for k = 1:step-2  
                    if k == 1 && step == size(route, 1)
                        continue;
                    end           
                    
                    existingStart = route(k);
                    existingEnd = route(k+1);

                    
                    if checkLineSegmentIntersection(latestStart, latestEnd, existingStart, existingEnd, circleCenters)
                        cumulativeIntern(ant) = cumulativeIntern(ant) + 1 ;
                        segmentIntern = segmentIntern + 1;
                    end
                end
            end

            currentPoint = nextPoint;
        end

        
        route = [route, route(1)];
        routes(ant, :) = route;
        
        startCenter = circleCenters(route(end-1),1:2);
        endCenter = circleCenters(route(end),1:2);
        segmentLength = norm(startCenter - endCenter);
        routeReward = (1/normalizedDistMatrix(route(end), route(end-1)))^beta;
        routeRewards(ant) = routeRewards(ant) + routeReward;
        routeLengths(ant) = routeLengths(ant) + segmentLength;
        
        p1 = circleCenters(route(end-2),1:2);
        p2 = circleCenters(route(end-1),1:2);
        p3 =circleCenters(route(end),1:2);
        angle = calculateAngle(p1, p2, p3);
        cumulativeAngles(ant) = cumulativeAngles(ant) + angle;
        angleRewards(ant) = angleRewards(ant) + (1/((1+segmentAngle)/180))^gamma;
        
        latestStart = route(end-1); 
        latestEnd = route(end);   
        segmentIntern = 0;
        for k = 2:size(route, 2)-2  
             if k == 1 && i == size(route, 2)
                continue;
            end           
            
            existingStart = route(k);
            existingEnd = route(k+1);

            
            if checkLineSegmentIntersection(latestStart, latestEnd, existingStart, existingEnd, circleCenters)
                cumulativeIntern(ant) = cumulativeIntern(ant) + 1 ;
                segmentIntern = segmentIntern + 1;
            end
        end     

        
        internRewards(ant) = exp(-cumulativeIntern(ant) / max(cumulativeIntern + 1e-6));
        routesQuality(ant) = routeLengths(ant)*0.1*0.3 + cumulativeAngles(ant)*0.7 + cumulativeIntern(ant)*1500;
  
        
         if cumulativeIntern(ant) < minIntern
            minIntern = cumulativeIntern(ant);
         end

        if cumulativeIntern(ant) == minIntern
            if routeLengths(ant)*0.1*0.3 + cumulativeAngles(ant)*0.7 < bestRouteLength*0.1*0.3 + bestAngle*0.7
               bestRouteLength = routeLengths(ant);
                bestRoute = route;
                bestAngle = cumulativeAngles(ant);   
                bestIntern = cumulativeIntern(ant);
                AntID = i;
            end
         elseif routesQuality(ant) < bestRouteLength*0.1*0.3 + bestAngle*0.7 + bestIntern*1500

            bestRouteLength = routeLengths(ant);
            bestRoute = route;
            bestAngle = cumulativeAngles(ant);   
            bestIntern = cumulativeIntern(ant);
            AntID = i;
         end
          
    end

    
    [sortedQuality, sortedIndices] = sort(routesQuality, 'descend');  
    numTopAnts = floor(numAnts * 0.05);  

    
    pheromone = (1 - evaporationRate) * pheromone;  

    
    for ant = 1:numTopAnts
        angleFactor = exp(-cumulativeAngles(sortedIndices(ant)) / max(cumulativeAngles + 1e-6));
        InterReward = exp(-cumulativeIntern(sortedIndices(ant)) / max(cumulativeIntern + 1e-6));
        normLength = (routeLengths(sortedIndices(ant)) - min(routeLengths)) / (max(routeLengths) - min(routeLengths) + 1e-6);
        for i = 1:(numPoints - 1)
            
            pheromone(routes(sortedIndices(ant), i), routes(sortedIndices(ant), i+1)) = ...
                pheromone(routes(sortedIndices(ant), i), routes(sortedIndices(ant), i+1)) + ...
                + (1 - normLength)*0.7 + 0.3 * angleFactor + 2* InterReward;
                
        end
    end
    
    
    for i = 1:length(bestRoute)-1
        start = bestRoute(i);       
        end_ = bestRoute(i+1);      
        
        
        total_pheromone = total_pheromone + pheromone(start, end_);
    end

    
    rewardHistory(iteration) = total_pheromone;
    distanceHistory(iteration) = bestRouteLength;
    angHistory(iteration) = bestAngle;
    internHistory(iteration) = minIntern;
end



figure; 
subplot(4, 1, 1); 
plot(1:iteration, rewardHistory, 'b-', 'LineWidth', 2); 
xlabel('Episode');
ylabel('Total Reward');
title('Reward per Episode');
grid on;


subplot(4, 1, 2); 
plot(1:iteration, distanceHistory, 'r-', 'LineWidth', 2); 
xlabel('Episode');
ylabel('Total Distance');
title('Distance per Episode');
grid on;


subplot(4, 1, 3); 
plot(1:iteration, angHistory, 'r-', 'LineWidth', 2); 
xlabel('Episode');
ylabel('Total ang');
title('Ang per Episode');
grid on;


subplot(4, 1, 4); 
plot(1:iteration, internHistory, 'r-', 'LineWidth', 2); 
xlabel('Episode');
ylabel('Total Intersection');
title('Intersection per Episode');
grid on;

numPoints = size(circleCenters, 1);
distanceMatrix = zeros(numPoints);

for i = 1:numPoints
    for j = i+1:numPoints
        distanceMatrix(i, j) = norm(circleCenters(i, :) - circleCenters(j, :));
        distanceMatrix(j, i) = distanceMatrix(i, j); 
    end
end


fprintf('Best path length visiting all regions: %.f units\n', bestRouteLength);
fprintf('Best angle visiting all regions: %.f degree\n', bestAngle);
fprintf('Best Intern visiting all regions: %.f \n', bestIntern);


totalAngle = 0;

for i = 2:(length(bestRoute) - 1)
    
    p1 = circleCenters(bestRoute(i - 1), :);
    p2 = circleCenters(bestRoute(i), :);
    p3 = circleCenters(bestRoute(i + 1), :);

    
    angle = calculateAngle(p1, p2, p3);
    totalAngle = totalAngle + angle;
end

    disp(['Iteration ', num2str(iteration), ' Best Route Length: ', num2str(bestRouteLength), ...
          ' Cumulative Angle: ', num2str(totalAngle)]);


figure;
imshow(template);
hold on; 
numPoints = size(points, 1); 
for i = 1:numPoints
    theta = linspace(0, 2*pi, 100);
    x = points(i, 1) + smallRadius * cos(theta);
    y = points(i, 2) + smallRadius * sin(theta);
    plot(x, y, 'b-'); 
    plot(points(i, 1), points(i, 2), 'bo', 'MarkerFaceColor', 'b'); 
end

for i = 1:size(circleCenters, 1)
    bestCenter = circleCenters(i,1:2);
    theta = linspace(0, 2*pi, 100);
    x = bestCenter(1) + largeRadius * cos(theta);
    y = bestCenter(2) + largeRadius * sin(theta);
    plot(x, y, 'r-','LineWidth', 0.5); 
    plot(bestCenter(1), bestCenter(2), 'ro', 'MarkerFaceColor', 'r'); 
    text(bestCenter(1), bestCenter(2), num2str(i), 'Color', 'r', ...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10); 
end




for i = 1:size(circleCenters, 1)
    plot([circleCenters(bestRoute(i), 1), circleCenters(bestRoute(i+1), 1)], ...
         [circleCenters(bestRoute(i), 2), circleCenters(bestRoute(i+1), 2)], 'r-', 'LineWidth', 2);
end

title('Optimal Path with Ant Colony Optimization');
hold off;
