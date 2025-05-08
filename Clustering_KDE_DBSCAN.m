inputFile = 'expanded_polygon_equal_distance.txt'; 
boundary = readmatrix(inputFile, 'Delimiter', ' '); 
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

k = 10; 
density = zeros(numPoints, 1);

for i = 1:numPoints
    distances = sqrt((points(i,1) - points(:,1)).^2 + (points(i,2) - points(:,2)).^2);
    [~, idx] = sort(distances);
    density(i) = sum(distances(idx(2:k+1)) < 270); 
end

customColors = [
    218/255, 236/255, 209/255;   
    165/255, 213/255, 158/255;  
    90/255, 182/255, 106/255; 
    40/255, 141/255, 73/255  
];

customColormap = interp1(linspace(0, 1, size(customColors, 1)), customColors, linspace(0, 1, 256));
gridX = linspace(min(points(:,1)), max(points(:,1)), 200);
gridY = linspace(min(points(:,2)), max(points(:,2)), 200);
[xx, yy] = meshgrid(gridX, gridY);
bandwidth = 80; 
densityValues = ksdensity(points, [xx(:), yy(:)], 'Bandwidth', bandwidth);
densityGrid = reshape(densityValues, length(gridY), length(gridX));

figure;
imagesc(gridX, gridY, densityGrid);
colorbar;
set(gca, 'YDir', 'reverse');
title('Kernel Density Estimation');
xlabel('X');
ylabel('Y');

densityPerPoint = zeros(size(points, 1), 1); 
for i = 1:size(points, 1)
    
    [~, xIdx] = min(abs(gridX - points(i,1)));
    [~, yIdx] = min(abs(gridY - points(i,2)));
    
    densityPerPoint(i) = densityGrid(yIdx, xIdx);
end

threshold = 2.8; 
highDensityPoints = points(densityPerPoint*10^7 > threshold, :);
lowDensityPoints = points(densityPerPoint*10^7 <= threshold, :);
epsilon = 300;
minPts = 11;
idx = dbscan(highDensityPoints, epsilon, minPts);

figure;
hold on;

numClusters = max(idx); 

colors = lines(numClusters); 
scatter(highDensityPoints(:,1), highDensityPoints(:,2),  100, idx, 'filled');
 
highDensityPoints = vertcat(highDensityPoints, lowDensityPoints);
idx = [idx; -1 * ones(size(lowDensityPoints, 1), 1)];

colors = lines(numClusters); 
for i = 1:numClusters
    clusterPoints = highDensityPoints(idx == i, :);
    if size(clusterPoints, 1) >= 3
        k = convhull(clusterPoints(:,1), clusterPoints(:,2));
        polygonPoints = clusterPoints(k, :);  
        fill(clusterPoints(k,1), clusterPoints(k,2), colors(i, :),'LineWidth', 1, 'FaceAlpha', 0.3, 'EdgeColor', colors(i, :));
    end
end

plot(boundary(:,1), boundary(:,2), 'k-', 'LineWidth', 1);

hScatter = scatter(points(:,1), points(:,2), 50, density, 'filled');

colormap(customColormap);

colorbar;
title('Voronoi Diagram within Polygon');
xlabel('X');
ylabel('Y');
set(gca, 'YDir', 'reverse');
axis equal;

hold off;
