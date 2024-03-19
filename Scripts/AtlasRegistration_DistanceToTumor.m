%% REQUIRES THE FOLLOWING JSON FILES FROM NUTIL IN THE FOLDER:
%% 1. cell_coordinates_TumorHemisphere.json
%% 2. cell_coordinates_NonTumorHemisphere.json
%% 3. tumor_coordinates_TumorHemisphere.json
%% 4. tumor_coordinates_NonTumorHemisphere.json

folderPath = uigetdir('Select the folder in which the files are located.');

cellFilePath1 = fullfile(folderPath, 'cell_coordinates_TumorHemisphere.json');
cellFilePath2 = fullfile(folderPath, 'cell_coordinates_NonTumorHemisphere.json');

tumorFilePath1 = fullfile(folderPath, 'tumor_coordinates_TumorHemisphere.json');
tumorFilePath2 = fullfile(folderPath, 'tumor_coordinates_NonTumorHemisphere.json');

cellData1 = jsondecode(fileread(cellFilePath1));
cellData2 = jsondecode(fileread(cellFilePath2));
tumorData1 = jsondecode(fileread(tumorFilePath1));
tumorData2 = jsondecode(fileread(tumorFilePath2));
cellData = [cellData1; cellData2];
tumorData = [tumorData1; tumorData2];

numIndices = numel(cellData);
Sum_numEntries = 0;

% cellData
for i = 1:numIndices

    if isempty(cellData(i).triplets)
        continue;
    end

    triplets = cellData(i).triplets;
    
    numEntries(i) = numel(triplets);
    
    Sum_numEntries = Sum_numEntries + numEntries(i);

end

resultArray_cells = zeros(Sum_numEntries/3, 3); % Each line in JSON file contains 3 consecutive numbers, therefore *3

resultIndex = 1;

for i = 1:numIndices

    if isempty(cellData(i).triplets)
        continue;
    end

    triplets = cellData(i).triplets;

    for j = 1:3:numEntries(i)
        tripletValues = triplets(j:j+2);
         
        resultArray_cells(resultIndex, 1:3) = tripletValues;
        
        resultIndex = resultIndex + 1;
    end
end

numIndices = numel(tumorData);
Sum_numEntries = 0;

% tumorData
for i = 1:numIndices

    if isempty(tumorData(i).triplets)
        continue;
    end

    triplets = tumorData(i).triplets;
    
    numEntries(i) = numel(triplets);
    
    Sum_numEntries = Sum_numEntries + numEntries(i);

end

resultArray_tumor = zeros(Sum_numEntries/3, 3); % Each line in JSON file contains 3 consecutive numbers, therefore *3

resultIndex = 1;

for i = 1:numIndices

    if isempty(tumorData(i).triplets)
        continue;
    end

    triplets = tumorData(i).triplets;

    for j = 1:3:numEntries(i)
        tripletValues = triplets(j:j+2);
         
        resultArray_tumor(resultIndex, 1:3) = tripletValues;
        
        resultIndex = resultIndex + 1;
    end
end

numRows_cells = size(resultArray_cells, 1);
numRows_tumor = size(resultArray_tumor, 1);

% create csv
[outputCsvFile, outputCsvPath] = uiputfile('*.csv', 'Select a location for the output CSV file.');
outputCsvFile = fopen(fullfile(outputCsvPath, outputCsvFile), 'w');
fprintf(outputCsvFile, 'Number,Cell_x,Cell_y,Cell_z,NextTumorPixel_x,NextTumorPixel_y,NextTumorPixel_z,Distance,Distance_um\n');

for i = 1:numRows_cells

    minDistance = Inf;
    nearestTumorPoint = [];

    x_cell = resultArray_cells(i,1);
    y_cell = resultArray_cells(i,2);
    z_cell = resultArray_cells(i,3);

    for j = 1:numRows_tumor
        x_tumor = resultArray_tumor(j,1);
        y_tumor = resultArray_tumor(j,2);
        z_tumor = resultArray_tumor(j,3);

        distance = sqrt((x_cell - x_tumor)^2 + (y_cell - y_tumor)^2 + (z_cell - z_tumor)^2);

        if distance < minDistance
            minDistance = distance;
            minDistance_um = minDistance * 25; % 1 QuickNII-Coordinate Unit = 25 um in ABA Mouse (https://www.nitrc.org/plugins/mwiki/index.php?title=quicknii:Coordinate_systems)
            nearest_x_tumor = x_tumor;
            nearest_y_tumor = y_tumor;
            nearest_z_tumor = z_tumor;
        end
    end

    fprintf(outputCsvFile, '%d,%f,%f,%f,%f,%f,%f,%f,%f\n', i, x_cell, y_cell, z_cell, nearest_x_tumor, nearest_y_tumor, nearest_z_tumor, minDistance, minDistance_um);
end