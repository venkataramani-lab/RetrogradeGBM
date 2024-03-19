% INPUT:
% 1. Results-csv of input cells of the last time point (X AND Y OF THE CENTROIDS OF THE CELLS MUST BE IN COLUMNS 3 AND 4, INVERT Y MUST BE SELECTED IN SET MEASUREMENTS IN FIJI)
% 2. determine clustering parameters (epsilon, minPoints) and enter x and y of the image
% 3. determine output folder for clusters of the last time point 
% 4. enter Results-csv folder with input cells of all time points (X AND Y OF THE CENTROIDS OF THE CELLS MUST BE IN COLUMNS 3 AND 4, INVERT Y MUST BE SELECTED IN SET MEASUREMENTS IN FIJI)
% 5. enter Results-csv-folder with starter cells of all time points (X AND Y OF THE CENTROIDS OF THE CELLS MUST BE IN COLUMNS 3 AND 4, INVERT Y MUST BE SELECTED IN SET MEASUREMENTS IN FIJI)
% 6. output starter information csv
% 7. output input information csv

% Smoothing factor alpha can be adjusted in smoothPolygon.m
% NEEDS THE FUNCTIONS pointToPolylineDistance.m AND smoothPolygon.m

csv_path = uigetfile('Select csv file with coordinates of input cells.');
raw = readtable(csv_path, 'NumHeaderLines', 1);
T = table2array(raw(:, 3:4));
x = T(:, 1);
y = T(:, 2);

cluster_quality = "n";
while cluster_quality == "n"
    prompt = 'epsilon: ';
    epsilon = input(prompt); %800
    prompt = 'minpoints: ';
    minpoints = input(prompt); %10
    
    % DBSCAN
    idx = dbscan(T, epsilon, minpoints);
    
    % Visualization of the clusters
    h = figure;
    gscatter(x, y, idx);
    title(sprintf('DBSCAN-Clustering-epsilon%d-minpoints%d', epsilon, minpoints));
    pause(5);

    prompt = 'Happy with the clustering? (y/n) ';
    cluster_quality = input(prompt, 's');
end

prompt = 'Enter x of the image: ';
numCols = input(prompt);

prompt = 'Enter x of the image: ';
numRows = input(prompt);

% Asking for output folder
outputFolder = uigetdir('Choose an output folder.');

if outputFolder == 0
    error('Selection of the storage location canceled. The program is terminated.');
end

% save DBScan png
outputFilename = fullfile(outputFolder, sprintf('DBSCAN_Clustering_epsilon%d_minpoints%d.png', epsilon, minpoints));
saveas(h, outputFilename, 'png');

uniqueClusters = unique(idx);

roiImage = zeros(numRows, numCols);

for i = 1:length(uniqueClusters)
    clusterIdx = uniqueClusters(i);
    
    % Select cells in the current cluster
    clusterCells = T(idx == clusterIdx, :);
    
    % Find the outermost points of the cluster
    minX = min(clusterCells(:, 1));
    maxX = max(clusterCells(:, 1));
    minY = min(clusterCells(:, 2));
    maxY = max(clusterCells(:, 2));
    
    % Export ROI as tif
    currentROI_positive = poly2mask([minX, maxX, maxX, minX], [minY, minY, maxY, maxY], numRows, numCols);
    currentROI = ~poly2mask([minX, maxX, maxX, minX], [minY, minY, maxY, maxY], numRows, numCols);

    imwrite(currentROI, fullfile(outputFolder, sprintf('Cluster_%d_ROI.tif', clusterIdx)));

    % If ClusterIdx is -1, reverse the values
    if clusterIdx ~= -1
        roiImage = roiImage | currentROI_positive;
    end
end

roiImageClusterMinus1 = roiImage;

% save tif
imwrite(roiImageClusterMinus1, fullfile(outputFolder, 'Cluster_-1_ROI.tif'));

uniqueClusters = unique(idx);

roiImage = zeros(numRows, numCols);

for i = 1:length(uniqueClusters)
    clusterIdx = uniqueClusters(i);
    
    % Select cells in the current cluster
    clusterCells = T(idx == clusterIdx, :);
    
    % If there are less than 3 points in the cluster, skip the ConvexHull creation
    if size(clusterCells, 1) < 3
        continue;
    end

    % create ConvexHull 
    k = convhull(clusterCells(:, 1), clusterCells(:, 2));

    % Export ROI as tif
    currentROI_positive = poly2mask(clusterCells(k, 1), clusterCells(k, 2), numRows, numCols);
    currentROI = ~poly2mask(clusterCells(k, 1), clusterCells(k, 2), numRows, numCols);

    imwrite(currentROI, fullfile(outputFolder, sprintf('Cluster_%d_ConvexHull.tif', clusterIdx)));

    % If ClusterIdx is -1, reverse the values
    if clusterIdx ~= -1
        roiImage = roiImage | currentROI_positive;
    end
end

roiImageClusterMinus1 = roiImage;

% save tif
imwrite(roiImageClusterMinus1, fullfile(outputFolder, 'Cluster_-1_ConvexHull.tif'));

% Ask user for the folder with the csv files
csvFolderPath = uigetdir('Select the folder with the CSV files of the other points in time.');

if csvFolderPath == 0
    error('Selection of the folder canceled. The program is terminated.');
end

% Choose all csv files in the folder
csvFiles = dir(fullfile(csvFolderPath, '*.csv'));

new_T = cell(length(csvFiles), 1);

% Iterate through each CSV file in the folder
for j = 1:length(csvFiles)
    currentCsvPath = fullfile(csvFiles(j).folder, csvFiles(j).name);
    
    new_raw = readtable(currentCsvPath, 'NumHeaderLines', 1);
    new_T{j} = table2array(new_raw(:, 3:4));
    new_x = new_T{j}(:, 1);
    new_y = new_T{j}(:, 2);
    
    timepointImage = zeros(numRows, numCols);
    timepointImage_smoothed = zeros(numRows, numCols);

    for i = 1:length(uniqueClusters)
        clusterIdx = uniqueClusters(i);
        
        if clusterIdx == -1
            continue;
        end
        
        clusterCells = T(idx == clusterIdx, :);
        
        if size(clusterCells, 1) < 3
            fprintf('Cluster %d: Skipping convex hull computation - Insufficient points.\n', clusterIdx);
            continue;
        end
        
        % create ConvexHull
        k = convhull(clusterCells(:, 1), clusterCells(:, 2));
        
        hullPoints = clusterCells(k, :);
        
        % Select cells from the new CSV file that are within the Convex Hull
        new_cells_inside_hull = new_T{j}(inpolygon(new_T{j}(:, 1), new_T{j}(:, 2), hullPoints(:, 1), hullPoints(:, 2)), :);

        % create new ConvexHulls
        try
            new_k = convhull(new_cells_inside_hull(:, 1), new_cells_inside_hull(:, 2));
            convexHullPoints = [new_cells_inside_hull(new_k, 1), new_cells_inside_hull(new_k, 2)];
            % Smoothing of Convex Hull
            smoothedPoints = smoothPolygon(convexHullPoints);
            new_k_smoothed = convhull(smoothedPoints(:, 1), smoothedPoints(:, 2));

        catch exception
            fprintf('Cluster %d, Timepoint %s: Skipping convex hull computation - Insufficient points.\n', clusterIdx, csvFiles(j).name);
            zeroImage = zeros(numRows, numCols);
            new_outputFilename = fullfile(outputFolder, sprintf('Cluster_%d_New_ConvexHull_%s.tif', clusterIdx, csvFiles(j).name));
            imwrite(zeroImage, new_outputFilename);
            timepointImage = timepointImage | zeroImage;

            new_outputFilename_smoothed = fullfile(outputFolder, sprintf('Smoothed_Cluster_%d_New_ConvexHull_%s.tif', clusterIdx, csvFiles(j).name));
            imwrite(zeroImage, new_outputFilename_smoothed);
            timepointImage_smoothed = timepointImage_smoothed | zeroImage;
            
            continue;
        end

        % save tif
        new_outputFilename = fullfile(outputFolder, sprintf('Cluster_%d_New_ConvexHull_%s.tif', clusterIdx, csvFiles(j).name));
        imwrite(~poly2mask(new_cells_inside_hull(new_k, 1), new_cells_inside_hull(new_k, 2), numRows, numCols), new_outputFilename);

        % save tif
        new_outputFilename_smoothed = fullfile(outputFolder, sprintf('Smoothed_Cluster_%d_New_ConvexHull_%s.tif', clusterIdx, csvFiles(j).name));
        imwrite(~poly2mask(smoothedPoints(new_k_smoothed, 1), smoothedPoints(new_k_smoothed, 2), numRows, numCols), new_outputFilename_smoothed);

        timepointImage = timepointImage | poly2mask(new_cells_inside_hull(new_k, 1), new_cells_inside_hull(new_k, 2), numRows, numCols);
        timepointImage_smoothed = timepointImage_smoothed | poly2mask(smoothedPoints(new_k_smoothed, 1), smoothedPoints(new_k_smoothed, 2), numRows, numCols);
    end
    
    outputFilenameTimepoint = fullfile(outputFolder, sprintf('All_ConvexHulls_Timepoint_%s.tif', csvFiles(j).name));
    imwrite(~timepointImage, outputFilenameTimepoint);

    outputFilenameTimepoint_smoothed = fullfile(outputFolder, sprintf('Smoothed_All_ConvexHulls_Timepoint_%s.tif', csvFiles(j).name));
    imwrite(~timepointImage_smoothed, outputFilenameTimepoint_smoothed);
end

csvFolderPath_starter = uigetdir('Select the folder with the starter CSV files.');

csvFiles_starter = dir(fullfile(csvFolderPath_starter, '*.csv'));

[outputCsvFile, outputCsvPath] = uiputfile('*.csv', 'Select a storage location for the output CSV file.');

if outputCsvFile == 0
    error('Selection of the storage location for the output CSV file canceled. The program is terminated.');
end

outputCsvFile = fopen(fullfile(outputCsvPath, outputCsvFile), 'w');
fprintf(outputCsvFile, 'X,Y,Cluster,OriginalFile,NextCluster,DistanceToNextCluster\n');

for j = 1:length(csvFiles_starter)
    currentCsvPath_starter = fullfile(csvFiles_starter(j).folder, csvFiles_starter(j).name);
    
    new_raw_starter = readtable(currentCsvPath_starter, 'NumHeaderLines', 1);
    new_T_starter{j} = table2array(new_raw_starter(:, 3:4));
    new_x_starter = new_T_starter{j}(:, 1);
    new_y_starter = new_T_starter{j}(:, 2);

    fprintf('Total points in new_T_starter: %d\n', size(new_T_starter{j}, 1));

    existingCoordinates = [];

    for i = 1:length(uniqueClusters)
        clusterIdx = uniqueClusters(i);
        
        if clusterIdx == -1
            continue;
        end
        
        clusterCells = T(idx == clusterIdx, :);
        
        fprintf('Cluster %d: Total points in clusterCells: %d\n', clusterIdx, size(clusterCells, 1));

        k = convhull(clusterCells(:, 1), clusterCells(:, 2));

        hullPoints = clusterCells(k, :);

        inConvexHull = inpolygon(new_x_starter, new_y_starter, hullPoints(:, 1), hullPoints(:, 2));
        cellsInCluster = new_T_starter{j}(inConvexHull, :);

        fprintf('Cluster %d: Total points in cellsInCluster: %d\n', clusterIdx, size(cellsInCluster, 1));

        existingCoordinates = [existingCoordinates; cellsInCluster];

        for k = 1:size(cellsInCluster, 1)
            
            fprintf(outputCsvFile, '%f,%f,%d,%s,0,0\n', cellsInCluster(k, 1), cellsInCluster(k, 2), clusterIdx, csvFiles_starter(j).name);
        end
    end

    missingCoordinates = setdiff(new_T_starter{j}, existingCoordinates, 'rows');
    existingCoordinates_timepoint{j} = existingCoordinates;

    for k = 1:size(missingCoordinates, 1)
        currentPoint = missingCoordinates(k, :);

        minDistance = Inf;
        closestHull = 100;

        for i = 1:length(uniqueClusters)
            clusterIdx = uniqueClusters(i);

            if clusterIdx == -1
                continue;
            end

            clusterCells = T(idx == clusterIdx, :);

            if size(clusterCells, 1) < 3
                continue;
            end

            l = convhull(clusterCells(:, 1), clusterCells(:, 2));

            hullPoints = clusterCells(l, :);
            new_cells_inside_hull_neu = new_T{j}(inpolygon(new_T{j}(:, 1), new_T{j}(:, 2), hullPoints(:, 1), hullPoints(:, 2)), :);

            if size(new_cells_inside_hull_neu, 1) < 3
                continue;
            end

            distance = pointToPolylineDistance(currentPoint, new_cells_inside_hull_neu);

            if distance < minDistance
                minDistance = distance;
                closestHull = clusterIdx;
            end
        end

        fprintf(outputCsvFile, '%f,%f,-1,%s,%d,%f\n', currentPoint(1), currentPoint(2), csvFiles_starter(j).name, closestHull, minDistance);
    end
end

fclose(outputCsvFile);
disp('Output CSV file was successfully created.');

[outputCsvFile_inputcells, outputCsvPath_inputcells] = uiputfile('*.csv', 'Select a storage location for the output CSV file.');

if outputCsvFile_inputcells == 0
    error('Selection of the storage location for the output CSV file canceled. The program is terminated.');
end

outputCsvFile_inputcells = fopen(fullfile(outputCsvPath_inputcells, outputCsvFile_inputcells), 'w');
fprintf(outputCsvFile_inputcells, 'X,Y,Cluster,OriginalFile,DistanceToStarter,DistanceToStarterInCluster\n');

for j = 1:length(csvFiles)
    currentCsvPath = fullfile(csvFiles(j).folder, csvFiles(j).name);
    
    new_raw = readtable(currentCsvPath, 'NumHeaderLines', 1);
    new_T = table2array(new_raw(:, 3:4));
    new_x = new_T(:, 1);
    new_y = new_T(:, 2);

    fprintf('Total points in new_T: %d\n', size(new_T, 1));

    existingCoordinates = [];

    for i = 1:length(uniqueClusters)
        clusterIdx = uniqueClusters(i);
        
        if clusterIdx == -1
            continue;
        end
        
        clusterCells = T(idx == clusterIdx, :);
        
        fprintf('Cluster %d: Total points in clusterCells: %d\n', clusterIdx, size(clusterCells, 1));

        k = convhull(clusterCells(:, 1), clusterCells(:, 2));

        hullPoints = clusterCells(k, :);

        inConvexHull = inpolygon(new_x, new_y, hullPoints(:, 1), hullPoints(:, 2));
        cellsInCluster = new_T(inConvexHull, :);

        fprintf('Cluster %d: Total points in cellsInCluster: %d\n', clusterIdx, size(cellsInCluster, 1));

        existingCoordinates = [existingCoordinates; cellsInCluster];

        for k = 1:size(cellsInCluster, 1)
            currentPoint = cellsInCluster(k, :);
            minDistance = Inf;
            minDistance_Cluster = Inf;
    
            for m = 1:length(new_T_starter{j})
                currentPoint_starter = new_T_starter{j}(m, :);
    
                distance = pdist2(currentPoint, currentPoint_starter);
    
                if distance < minDistance
                    minDistance = distance;
                end
            end

            for m = 1:length(existingCoordinates_timepoint{j})
                currentPoint_starter = existingCoordinates_timepoint{j}(m, :);
    
                distance = pdist2(currentPoint, currentPoint_starter);
    
                if distance < minDistance_Cluster
                    minDistance_Cluster = distance;
                end
            end
            fprintf(outputCsvFile_inputcells, '%f,%f,%d,%s,%f,%f\n', cellsInCluster(k, 1), cellsInCluster(k, 2), clusterIdx, csvFiles(j).name, minDistance, minDistance_Cluster);
        end
    end

    missingCoordinates = setdiff(new_T, existingCoordinates, 'rows');

    for k = 1:size(missingCoordinates, 1)
        currentPoint = missingCoordinates(k, :);

        minDistance = Inf;
        minDistance_Cluster = Inf;

        for i = 1:length(new_T_starter{j})
            currentPoint_starter = new_T_starter{j}(i, :);

            distance = pdist2(currentPoint, currentPoint_starter);

            if distance < minDistance
                minDistance = distance;
            end
        end

        for m = 1:length(existingCoordinates_timepoint{j})
            currentPoint_starter = existingCoordinates_timepoint{j}(m, :);
    
            distance = pdist2(currentPoint, currentPoint_starter);
    
            if distance < minDistance_Cluster
                minDistance_Cluster = distance;
            end
        end
        fprintf(outputCsvFile_inputcells, '%f,%f,-1,%s,%f,%f\n', currentPoint(1), currentPoint(2), csvFiles(j).name, minDistance, minDistance_Cluster);
    end
end

fclose(outputCsvFile_inputcells);
disp('Output CSV file was successfully created.');


% THEN THE CONVEX HULLS CAN BE EVALUATED IN FIJI WITH ANALYZE PARTICLE (FOR EXAMPLE AREA AND CIRCULARITY)