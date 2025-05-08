

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

global circleCenters;
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

    [~, idx1] = min(distances);
    nearestPoint = uncoveredPoints(idx1, :);

    
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

function points = getPointsInRegion(center, radius)
    
    gridResolution = 30;    
    [X, Y] = meshgrid(center(1)-radius:gridResolution:center(1)+radius, ...
                      center(2)-radius:gridResolution:center(2)+radius);  
    distances = sqrt((X - center(1)).^2 + (Y - center(2)).^2); 
    mask = distances <= radius;   
    x = X(mask);
    y = Y(mask);   
    points = [x, y];
end

numCenters = size(circleCenters, 1); 
gamma = 0.2; 
epsilon = 0.5; 
episodes = 40000; 
beta = 2; 
Alpha = 1; 
omega = 100;
epsilon_max = 1.0; 
epsilon_min = 0.001;
lambda = 1.0000e-04; 

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

R = zeros(numCenters, numCenters); 

for state = 1:numCenters   
    localDistances = distances(state, :);   
    finiteLocalDistances = localDistances(isfinite(localDistances));   
    localMinDistance = min(finiteLocalDistances);
    localMaxDistance = max(finiteLocalDistances);   
    normalizedLocalDistances = (localDistances - localMinDistance) / ...
                               (localMaxDistance - localMinDistance);  
    normalizedLocalDistances = 1 - normalizedLocalDistances;   
    normalizedLocalDistances(~isfinite(localDistances)) = Inf;   
    R(state, :) = normalizedLocalDistances;
end

function angle = calculateAngle(p1, p2, p3)
    v1 = p2 - p1;
    v2 = p3 - p2;
    cosTheta = dot(v1, v2) / (norm(v1) * norm(v2) + eps); 
    angle = acosd(cosTheta); 
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


Q = zeros(numCenters, numCenters);
returns = zeros(numCenters, numCenters);
visitCounts = zeros(numCenters, numCenters);


state = randi(numCenters); 
visited = false(1, numCenters);
visited(state) = true;
path = [];
bestroute = [];
shortestdistance = Inf;
minAng = inf;
minIntercount = inf;
inistate = state;

while ~all(visited)   
    qValues = R(state, :);    
    validActions = find(~visited & (1:numCenters) ~= state);
    visited(state) = true;   
    if isempty(validActions)
        error('No valid actions available.');
    end   
    [~, maxIdx] = max(qValues(validActions));
    action = validActions(maxIdx);   
    path = [path; state, action];
    nextState = action;    
    lastState = state;
    state = nextState;
    visited(state) = true;
end

firstState = path(1, 1); 
lastState = path(end, 2); 
path = [path; lastState, firstState]; 
Intercount = 0;
G = 0; 
for t = size(path, 1):-1:1
    state = path(t, 1);
    action = path(t, 2); 
    reward = R(state, action);
    if t>1 
       p1 = circleCenters(path(t-1, 1), 1:2); 
       p2 = circleCenters(path(t, 1), 1:2);
       p3 = circleCenters(path(t, 2), 1:2);
       angle = calculateAngle(p1, p2, p3);
    end
    anglePenalty = 1-angle/180;
    
    if t>2
        latestStart = path(t, 1); 
        latestEnd = path(t, 2);   
        for k = 1:t-2  
             if k == 1 && t == size(path, 1)
                continue;
            end                    
            existingStart = path(k, 1);
            existingEnd = path(k, 2);    
            if checkLineSegmentIntersection(latestStart, latestEnd, existingStart, existingEnd, circleCenters)
                Intercount = Intercount + 1 ;
            end
        end        
    end
    
    G = reward + anglePenalty*beta +gamma * G;
    fprintf('reward,anglepenalty is: %.2f -> %.2f -> %.2f-> %.2f\n', reward, anglePenalty, G, Intercount);
    visitCounts(state, action) = visitCounts(state, action) + 1;
    returns(state, action) = returns(state, action) + (G - returns(state, action)) / visitCounts(state, action);
    Q(state, action) = returns(state, action);
end

rewardHistory = zeros(1, episodes); 
distanceHistory = zeros(1, episodes);
angHistory = zeros(1, episodes); 
internHistory = zeros(1, episodes); 

for ep = 1:episodes
    state = randi(numCenters); 
    visited = false(1, numCenters);
    visited(state) = true;
    reward = 0;
    distance = 0;
    path = [];
    epsilon = max(epsilon_min, epsilon_max * exp(-lambda * ep)); 
    Intercount = 0;

    while ~all(visited)
        qValues = Q(state, :);
        validActions = find(~visited & (1:numCenters) ~= state);
        if isempty(validActions)
            error('No valid actions available.');
        end
        if rand < epsilon            
            action = validActions(randi(length(validActions)));
        else
            
            [~, maxIdx] = max(qValues(validActions));
            action = validActions(maxIdx);
        end        
        path = [path; state, action];
        nextState = action;
        lastState = state;
        state = nextState;
        visited(state) = true;
    end
    
    firstState = path(1, 1); 
    lastState = path(end, 2); 
    path = [path; lastState, firstState]; 
    
    G = 0; 
    for t = size(path, 1):-1:1
        state = path(t, 1);
        action = path(t, 2);
        
        reward = R(state, action);
        if isinf(reward)
            fprintf('reward is -Inf\n, %.2f -> %.2f', state, action);
            break;
        end

        if t>1 
           p1 = circleCenters(path(t-1, 1), 1:2); 
           p2 = circleCenters(path(t, 1), 1:2);
           p3 = circleCenters(path(t, 2), 1:2);
           angle = calculateAngle(p1, p2, p3);
        end
        anglePenalty = (1-angle/180)^2;

        if t>2
            latestStart = path(t, 1); 
            latestEnd = path(t, 2);   
            for k = 1:t-2  
                 if k == 1 && t == size(path, 1)
                    continue;
                end               
                existingStart = path(k, 1);
                existingEnd = path(k, 2);
    
                if checkLineSegmentIntersection(latestStart, latestEnd, existingStart, existingEnd, circleCenters)
                    Intercount = Intercount + 1 ;
                end
            end        
        end
        
        InterReward = exp(-10000 * Intercount /  sum(1:size(path, 1) - 2));
        G = reward*Alpha + anglePenalty*beta + InterReward*omega + gamma * G;        
        visitCounts(state, action) = visitCounts(state, action) + 1;
        returns(state, action) = returns(state, action) + (G - returns(state, action)) / visitCounts(state, action);
        Q(state, action) = returns(state, action);       
    end
    
    G = 0;
    Dis = 0;
    Ang = 0;

    for t = size(path, 1):-1:1
    state = path(t, 1);
    action = path(t, 2);

    if t>2 
       p1 = circleCenters(path(t-1, 1), 1:2); 
       p2 = circleCenters(path(t, 1), 1:2);
       p3 = circleCenters(path(t, 2), 1:2);
       angle = calculateAngle(p1, p2, p3);
    end
    
    Ang = Ang + angle;

    reward = Q(state, action);
    G = G + reward; 
    
    distance = distances(state, action);
    Dis = Dis + distance;
    end
    
    route = [path(:,1); path(size(path, 1),2)];

    if 0.1 * Dis * 0.3 + Ang * 0.7 + Intercount * 500 < 0.1 * shortestdistance * 0.3 + minAng * 0.7 + minIntercount * 500
        shortestdistance = Dis;
        bestroute = route;
        minAng = Ang;
        minIntercount = Intercount;
    end
    
    internHistory(ep) = minIntercount;
    rewardHistory(ep) = G;
    distanceHistory(ep) = Dis;
    angHistory(ep) = Ang;
end



figure; 
subplot(4, 1, 1); 
plot(1:episodes, rewardHistory, 'b-', 'LineWidth', 2); 
xlabel('Episode');
ylabel('Total Reward');
title('Reward per Episode');
grid on;

subplot(4, 1, 2); 
plot(1:episodes, distanceHistory, 'r-', 'LineWidth', 2); 
xlabel('Episode');
ylabel('Total Distance');
title('Distance per Episode');
grid on;

subplot(4, 1, 3); 
plot(1:episodes, angHistory, 'r-', 'LineWidth', 2); 
xlabel('Episode');
ylabel('Total ang');
title('Ang per Episode');
grid on;

subplot(4, 1, 4); 
plot(1:episodes, internHistory, 'r-', 'LineWidth', 2); 
xlabel('Episode');
ylabel('Total Intersection');
title('Intersection per Episode');
grid on;

optimalPath = zeros(length(bestroute), 2);
optimalPath(1,:) = circleCenters(bestroute(1), 1:2); 
optimalPath(length(bestroute),:) = circleCenters(bestroute(1), 1:2); 
lastCenter = circleCenters(bestroute(1), 1:2);

for i = 2:length(bestroute)-1
    startCenter = circleCenters(bestroute(i), 1:2); 
    endCenter = circleCenters(bestroute(i+1), 1:2); 
    radius = circleCenters(bestroute(i), 3);
    OriginalDistance = norm(lastCenter - startCenter) + norm(startCenter - endCenter);
    OriginalAngle = calculateAngle(lastCenter, startCenter, endCenter);
    OriginalCost = 0.1*OriginalDistance * 0.4 + OriginalAngle * 0.6;
    bestCost = OriginalCost;
    bestPoint = startCenter;
    
    candidatePoints = getPointsInRegion(startCenter, radius);
    
    for j = 1:size(candidatePoints, 1)
        candidateCenter = candidatePoints(j, :);
        
        dist1 = norm(lastCenter - candidateCenter);
        dist2 = norm(candidateCenter - endCenter);
        totalDistance = dist1 + dist2;

        angle = calculateAngle(lastCenter, candidateCenter, endCenter);

        cost = 0.1*totalDistance * 0.4 + angle * 0.6;
        
        if cost < bestCost
            bestCost = cost;
            bestPoint = candidateCenter;
        end   
    end

    optimalPath(i,:) = bestPoint;
    lastCenter = bestPoint;
end

totalAngle = 0;
Angle = 0;

segmentLengths = [];  
totalPathLength = 0;  

for i = 1:length(bestroute)-1

    startCenter = optimalPath(i);
    endCenter = optimalPath(i+1);
    segmentLength = norm(startCenter - endCenter);  
    
    segmentLengths = [segmentLengths, segmentLength];  
end
totalPathLength = sum(segmentLengths);
fprintf('Shortest path length visiting all regions: %.2f units\n', totalPathLength);
totalAngle = 0;
angles = [];
    for i = 2:length(bestroute)-1             
        p1 = optimalPath(i-1,:);
        p2 = optimalPath(i,:);
        p3 = optimalPath(i+1,:);         
        angle = calculateAngle(p1, p2, p3);
        angles = [angles; bestroute(i), angle];
        totalAngle = totalAngle + angle;
    end
fprintf('Shortest angle visiting all regions: %.2f degree\n', totalAngle);

figure;
h = imshow(template);
hold on;
set(h, 'AlphaData', 1);

hold on; 

for i= 1:length(circleCenters)  
    theta = linspace(0, 2*pi, 100);   
    fill_color = [117/255, 159/255, 204/255]; 
    theta = linspace(0, 2*pi, 100);
    x = circleCenters(i, 1) + largeRadius * cos(theta);
    y = circleCenters(i, 2) + largeRadius * sin(theta);
    fill(x, y, fill_color, 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.64); 
    x_feasible = circleCenters(i, 1) + circleCenters(i, 3) * cos(theta);
    y_feasible = circleCenters(i, 2) + circleCenters(i, 3) * sin(theta);
    fill_color = [47/255, 245/255, 41/255]; 
    fill(x_feasible, y_feasible, fill_color,'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.43); 
    
end
for i = 1:length(bestroute)-1
    skipOuterLoop = false;
    startCenter = optimalPath(i,:);
    endCenter = optimalPath(i+1,:);  
    plot([startCenter(1), endCenter(1)], [startCenter(2), endCenter(2)], 'r-', 'LineWidth', 4.5);
    plot(startCenter(1), startCenter(2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
    plot(endCenter(1), endCenter(2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4); 
end

title('Optimal Path with Ant Colony Optimization');
hold off;


