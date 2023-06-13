% This script generates the default HARM data using UK-SSP5 as a template.
% It takes the original cell indexes produced by Katie Jenkins (UEA) and
% transforms them into 2D (or 3D or 4D) based arrays.

% Load Katie's original cell arrays
load('heat_impact_data.mat')
load('Mortality_Threshold_RR_Cold.mat')

% For temporary data (awaiting Katie to upload final version 22/2/23):
building_threshold = building_threshold_retro_shading;

% Load grid data for transforming back to 2D
grid_idx = 1:82;
grid_idy = 1:112;
load('PreProcessedData/LSM12.mat')
LSM = LSM12;
load('PreProcessedData/lat_UK_RCM.mat')
lats = lat_UK_RCM;
load('PreProcessedData/long_UK_RCM.mat')
lons = long_UK_RCM;

%% Find cell ID for converting to 2D
disp('Converting input data to 2D')
% Create new empty arrays for data
population_SSP5_2D = nan(length(grid_idx),length(grid_idy),5,5); % lon x lat x age groups x decades
Mortality_RR_2D = nan(length(grid_idx),length(grid_idy),5); % lon x lat x age groups
dailyDeathRate_2D = nan(length(grid_idx),length(grid_idy),5); % lon x lat x age groups
Mortality_MMT_2D = nan(length(grid_idx),length(grid_idy)); % lon x lat
HARM_regs = nan(length(grid_idx),length(grid_idy)); % lon x lat
building_numbers_SSP5_2D = nan(length(grid_idx),length(grid_idy),4,5); % lon x lat x building types x decades
building_threshold_2D = nan(length(grid_idx),length(grid_idy),4,2,5); % lon x lat x building types x adaptation/no adaptation x decades
    
% Go through each grid cell
for i = 1:length(grid_idx)
    for j = 1:length(grid_idy)
        
        % Only process if it is a land point
        if LSM(i,j) == 1
            
            % Pull out corresponding lat-long
            latID = lats(i,j);
            lonID = lons(i,j);
            X = dailyDeathRate(:,1);
            Y = dailyDeathRate(:,2);
            
            % Adjust lat and long to avoid rounding errors
            latID = fix(latID*1e5)/1e5;
            lonID = fix(lonID*1e5)/1e5;
            X = fix(X*1e5)/1e5;
            Y = fix(Y*1e5)/1e5;
            
            % Find row of Katie's data that matches 2D data
            ii = find(X==lonID);
            jj = find(Y==latID);
            
            % Double check that the row has correctly been found
            if ii ~= jj
                disp('An error has occurred identifying the correct grid ID')
                return
            end
            
            % Take Katie's data and save to correct part of output array
            for times = 1:5
                population_SSP5_2D(i,j,:,times) = population_SSP5{times,1}(ii,3:7);
                building_numbers_SSP5_2D(i,j,:,times) = building_numbers_SSP5{times,1}(ii,3:6);
                building_threshold_2D(i,j,:,1,times) = building_threshold{times,1}(ii,4:7);
                building_threshold_2D(i,j,:,2,times) = building_threshold{times,1}(ii,8:11);
            end
            Mortality_RR_2D(i,j,:) = Mortality_Threshold_RR(ii,5:9);
            Mortality_MMT_2D(i,j) = Mortality_Threshold_RR(ii,4);
            dailyDeathRate_2D(i,j,:) = dailyDeathRate(ii,3:7);
            HARM_regs(i,j) = Mortality_Threshold_RR(ii,3);

        end
    end
end
disp(' ')


% Remove RoI
population_SSP5_2D(repmat(Mortality_MMT_2D == 0,1,1,5,5)) = nan;
building_numbers_SSP5_2D(repmat(Mortality_MMT_2D == 0,1,1,4,5)) = nan;
building_threshold_2D(repmat(Mortality_MMT_2D == 0,1,1,4,2,5)) = nan;
dailyDeathRate_2D(repmat(Mortality_MMT_2D == 0,1,1,5)) = nan;
Mortality_RR_2D(repmat(Mortality_MMT_2D == 0,1,1,5)) = nan;
Mortality_RR_2D_cold(repmat(Mortality_MMT_2D == 0,1,1,5)) = nan;
Mortality_MMT_2D_cold(Mortality_MMT_2D == 0) = nan;
Mortality_MMT_2D(Mortality_MMT_2D == 0) = nan;

% Save the data
save('PreProcessedData/population_SSP5_2D.mat','population_SSP5_2D')
save('PreProcessedData/building_numbers_SSP5_2D.mat','building_numbers_SSP5_2D')
save('PreProcessedData/building_threshold_2D.mat','building_threshold_2D')
save('PreProcessedData/dailyDeathRate_2D.mat','dailyDeathRate_2D')
save('PreProcessedData/Mortality_RR_2D.mat','Mortality_RR_2D')
save('PreProcessedData/Mortality_MMT_2D.mat','Mortality_MMT_2D')
save('PreProcessedData/Mortality_RR_2D_cold.mat','Mortality_RR_2D_cold')
save('PreProcessedData/Mortality_MMT_2D_cold.mat','Mortality_MMT_2D_cold')
