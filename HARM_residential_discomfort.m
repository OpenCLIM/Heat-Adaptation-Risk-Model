%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   HARM model (Heat - Adaptation - Risk - Model)
%
%   Heat related Residential discomfort calculations
%
%   This version uses post-processed UKCP18 TMean data provided from HEAT model 
%   (Alan Kennedy-Asser, University of Bristol); building and population
%   data from UDM (Craig Robson, NCL); and residential building overheating 
%   data from UCL (Lauren Ferguson and Anna Mavrogianni, UCL)
%   
%   27/03/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT DATA
%   heat_impact_data.mat
%       
%   1.  Gridded population data {nClimScen}[nGridCells, 7]
%       Column 1&2 lat and long
%       column 3 {1} = Baseline 2011 
%       column 4-7 {1} = gridded poulation per age group
%       Columns 4-7 {2-5} = For 2020, 2030, 2050, 2080 UK from UDM (based
%       on UK-SSP projections.
%       NOTE: Only use total population for residential discomfort
%       calculations.
%
%   2. 12km gridded TMean daily data {nGridCells}(12, 14800). .csv files of 
%       daily TMean (one per grid cell) from HEAT model for baseline and warming levels 1.5°C
%       2.0°C, 3.0°C and 4.0°C
%
%       Grid coordinates also extracted from file names to provide grid ID for 
%       mapping based on lat_UK_RCM.mat and long_UK_RCM.mat.
%
%   3/4. lat_UK_RCM.mat and long_UK_RCM.mat: [82,112]. From Heat outputs, these are read 
%      in and provide latitude and longitude based on coordinate given in file 
%      names from 7 above. Saved in gridID.
%
%   5. deploymentPercent.(5,1). Percentage of total buildings where adaptation
%     (energy and shading retrofit) is installed. From UCL model. 
%     rows = 5 timeperiods (BL, 2020, 2030, 2050, 2080).
%
%   6. building_thresholds. {5,1}(2205, 11). TMean theresholds above which
%      residential discomfort would occur overnight in bedrooms. Based on UCL
%      model data. For 5 timeperiods (no data for 2080/2100 at present). Columns 1-2 lon and
%      lat, column 3 = GOR, column 4-7 = threholds for 4 building types with
%      NO adaptation; columns 8-11 = thresholds for 4 building types WITH
%      adaptation. Adaptation reflects energy retrofit and shading devices.
%
%   7. building_numbers. Total buildings per 12km gridcell for 4 building
%      types: Flats, detached, semi-detached, Terrace from UDM. Assumed
%      1.2ppl per dwelling in estimating future buildings based on density
%      of population.
%
%   8. building_name. String text for building names (in order) - for plots.
%      Flats, Detached, Semi-Detached, Terraced
%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CURRENT MODEL OUTPUTS:
%   1. Intermediate results: ResultsExceedDays.
%   Average annual number of days that TMean thresholds are exceeded by 
%   building type. Note - these values are based only on TMean data and 
%   threshold data, and do not account for the population or number of 
%   dwellings.
%
%   2. Final result: ResultsAnnualFrequencyEvents
%   Impact-based definition. Average annual frequency of overheating events 
%   (can be classified as medium, severe or extreme) where residents facing 
%   discomfort from dwellings overheating. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp ('Metric: UK heat related residential discomfort')
%% 1. Search for available files; Read Mean Temperature data from HEAT .csv files; set up arrays; read coordinates
%  Read each gridcell array to a cell Array inc. the cell ID.
%  Each TMean array represents daily data e.g. [1-360] per year [1-30][nCols=10,800] * 12 GCMS [nRows = 12]
%  Only need to do once - if files exist then skip.

% RUN EITHER BIAS CORRECTED OR RAW DEPENDENT ON USER DEFINED PARAMETER %
% To account for different number of days and leap years in files

if BiasCorrected == 0 %False, reads raw data which has 30 day months, no leap years
   disp ('Using non-bias corrected data')

    if (ClimateScenario == 1 || ClimateScenario == 6) %% past or all climate scenarios
        if exist('ARCADIA-Tmean-raw-past', 'file')
            disp ('loading input data Past TMean')

            dd = dir('ARCADIA-Tmean-raw-past\*.csv');
            fileNames = {dd.name}; 
            data_past_tmean = cell(numel(fileNames),2);
            data_past_tmean(:,1) = regexprep(fileNames, '_UKCP18-12km_T.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)   
                fname = fullfile('ARCADIA-Tmean-raw-past', fileNames(ii)); %allocate each .csv data file to cell array
                data_past_tmean{ii,2} = dlmread(fname{1},',',0,0);
            end

            save('input_data_past_tmean.mat', 'data_past_tmean', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together

        else
            disp ('No Past files exist')
        end
    end

    if (ClimateScenario == 2 || ClimateScenario == 6) %% 1.5degree or all climate scenarios
        if exist('ARCADIA_Tmean_raw_1_5C', 'file')
            disp ('loading input data 1.5degree TMean')

            dd = dir('ARCADIA_Tmean_raw_1_5C\*.csv');
            fileNames = {dd.name}; 
            data_1_5C_tmean = cell(numel(fileNames),2);
            data_1_5C_tmean(:,1) = regexprep(fileNames, '_UKCP18-12km_T.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)
                fname = fullfile('ARCADIA_Tmean_raw_1_5C', fileNames(ii)); %allocate each .csv data file to cell array
                 data_1_5C_tmean{ii,2} = dlmread(fname{1},',',0,0);
            end

            save('input_data_1_5C_tmean.mat', 'data_1_5C_tmean', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together

        else
            disp ('No 1.5degreeC files exist')
        end
    end
        %  2degreeC
    if (ClimateScenario == 3 || ClimateScenario == 6) %% 2 degree or all climate scenarios
        if exist('ARCADIA_Tmean_raw_2_0C', 'file')
            disp ('loading input data 2.0degree TMean')

            dd = dir('ARCADIA_Tmean_raw_2_0C\*.csv');
            fileNames = {dd.name}; 
            data_2_0C_tmean = cell(numel(fileNames),2);
            data_2_0C_tmean(:,1) = regexprep(fileNames, '_UKCP18-12km_T.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)   
                fname = fullfile('ARCADIA_Tmean_raw_2_0C', fileNames(ii)); %allocate each .csv data file to cell array
                data_2_0C_tmean{ii,2} = dlmread(fname{1},',',0,0);
            end

            save('input_data_2_0C_tmean.mat', 'data_2_0C_tmean', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together

        else
            disp ('No 2.0DegreeC files exist')
        end
    end

        %  3degreeC 
     if (ClimateScenario == 4 || ClimateScenario == 6) %% 3degreeC or all climate scenarios
        if exist('ARCADIA_Tmean_raw_3_0C', 'file')
            disp ('loading input data 3.0degree TMean')

            dd = dir('ARCADIA_Tmeanraw_3_0C\*.csv');
            fileNames = {dd.name}; 
            data_3_0C_tmean = cell(numel(fileNames),2);
            data_3_0C_tmean(:,1) = regexprep(fileNames, '_UKCP18-12km_T.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)   
                fname = fullfile('ARCADIA_Tmean_raw_3_0C', fileNames(ii)); %allocate each .csv data file to cell array
                data_3_0C_tmean{ii,2} = dlmread(fname{1},',',0,0);
            end

            save('input_data_3_0C_tmean.mat', 'data_3_0C_tmean', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together

        else
            disp ('No 3.0DegreeC files exist')
        end
     end
         %  4degreeC 
     if (ClimateScenario == 5 || ClimateScenario == 6) %% 4degreeC or all climate scenarios
        if exist('ARCADIA_Tmean_raw_4_0C', 'file')
            disp ('loading input data 4.0degree TMean')

            dd = dir('ARCADIA_Tmean_raw_4_0C\*.csv');
            fileNames = {dd.name}; 
            data_4_0C_tmean = cell(numel(fileNames),2);
            data_4_0C_tmean(:,1) = regexprep(fileNames, '_UKCP18-12km_T.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)   
                fname = fullfile('ARCADIA_Tmean_raw_4_0C', fileNames(ii)); %allocate each .csv data file to cell array
                data_4_0C_tmean{ii,2} = dlmread(fname{1},',',0,0);
            end

            save('input_data_4_0C_tmean.mat', 'data_4_0C_tmean', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together
        else
            disp ('No 4.0DegreeC files exist')
        end
     end
 
elseif BiasCorrected == 1 %reads bias corrected data and removes leap years.
       disp ('Using bias corrected data')
       
    %  Past / baseline
    if (ClimateScenario == 1 || ClimateScenario == 6) %% baseline or all climate scenarios  
        if exist ('input_data_past_tmean.mat', 'file') %%I use this to speed up
        %if running offline from DAFNI - checks if file already exists.
           load ('input_data_past_tmean.mat')
           disp ('loading input data Past TMean')
           nClimScen = nClimScen+1; %counts # climate scenarios being run together
    
    elseif exist('ARCADIA-TMean-BC-past', 'file')
            disp ('Read in daily TMean data (Past) for mortality estimates')

            dd = dir('ARCADIA-TMean-BC-past\*.csv');
            fileNames = {dd.name}; 
            data_past_tmean = cell(numel(fileNames),2);
            data_past_tmean(:,1) = regexprep(fileNames, '_UKCP18-biascorr_past.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)   
                fname = fullfile('ARCADIA-TMean-BC-past', fileNames(ii)); %allocate each .csv data file to cell array
                data_past_tmean{ii,2} = dlmread(fname{1},',',0,0);
                data_past_tmean{ii,2}(:,1461:1461:end) = []; %delete leap year, every 1461 columns. Bias corrected data uses Gregorian Calendar.
                data_past_tmean{ii,2}(:,10950) = 0; %Add back in dummy value for 31st December so years are even.
            end

            save('input_data_past_tmean.mat', 'data_past_tmean', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together

        else
            disp ('No Past files exist')
        end
    end
    
    if (ClimateScenario == 2 || ClimateScenario == 6) %% 1.5 degree or all climate scenarios
        %  1.5degreeC 
        if exist ('input_data_1_5C_tmean.mat', 'file')  %%remove for running on DAFNI
           load ('input_data_1_5C_tmean.mat')           %
           disp ('loading input data 1.5C')             %
           nClimScen = nClimScen+1; %counts # climate scenarios being run together
    
        elseif exist('ARCADIA-Tmean-BC-15', 'file')
            disp ('Read in daily Tmean data (1.5DegreeC) for mortality estimates')

            dd = dir('ARCADIA-Tmean-BC-15\*.csv'); 
            fileNames = {dd.name}; 
            data_1_5C_tmean = cell(numel(fileNames),2);
            data_1_5C_tmean(:,1) = regexprep(fileNames, '_UKCP18-12km_T.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)
                fname = fullfile('ARCADIA-Tmean-BC-15', fileNames(ii)); %allocate each .csv data file to cell array
                data_1_5C_tmean{ii,2} = dlmread(fname{1},',',0,0);
                data_1_5C_tmean{ii,2}(:,1461:1461:end) = []; %delete leap year, every 1461 columns. Bias corrected data uses Gregorian Calendar.
                data_1_5C_tmean{ii,2}(:,10950) = 0; %Add back in dummy value for 31st December so years are even.
                
            end
            save('input_data_1_5C_tmean.mat', 'data_1_5C_tmean', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together

        else
            disp ('No 1.5degreeC files exist')
        end
    end
    
    %  2degreeC
    if (ClimateScenario == 3 || ClimateScenario == 6) %% 2 degree or all climate scenarios
    
        if exist ('input_data_2_0C_tmean.mat', 'file')  %%remove for running on DAFNI
           load ('input_data_2_0C_tmean.mat')           %
           disp ('loading input data 2.0C')             %
           nClimScen = nClimScen+1; %counts # climate scenarios being run together
           
        elseif exist('ARCADIA-Tmean-BC-20', 'file')
            disp ('Read in daily Tmean data (2.0DegreeC) for mortality estimates')

            dd = dir('ARCADIA-Tmean-BC-20\*.csv');
            fileNames = {dd.name}; 
            data_2_0C_tmean = cell(numel(fileNames),2);
            data_2_0C_tmean(:,1) = regexprep(fileNames, '_UKCP18-12km_T.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)   
                fname = fullfile('ARCADIA-Tmean-BC-20', fileNames(ii)); %allocate each .csv data file to cell array
                data_2_0C_tmean{ii,2} = dlmread(fname{1},',',0,0);
                data_2_0C_tmean{ii,2}(:,1461:1461:end) = []; %delete leap year, every 1461 columns. Bias corrected data uses Gregorian Calendar.
                data_2_0C_tmean{ii,2}(:,10950) = 0; %Add back in dummy value for 31st December so years are even.
            end
            save('input_data_2_0C_tmean.mat', 'data_2_0C_tmean', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together

        else
            disp ('No 2degreeC files exist')
        end
    end
    
    %  3degreeC
    if (ClimateScenario == 4 || ClimateScenario == 6) %% 3 degree or all climate scenarios
       if exist ('input_data_3_0C_tmean.mat', 'file')  %%remove for running on DAFNI
           load ('input_data_3_0C_tmean.mat')           %
           disp ('loading input data 3.0C')             %
           nClimScen = nClimScen+1; %counts # climate scenarios being run together 
           
        elseif exist('ARCADIA-Tmean-BC-30', 'file')
            disp ('Read in daily Tmean data (3.0DegreeC) for mortality estimates')

            dd = dir('ARCADIA-Tmean-BC-30\*.csv');
            fileNames = {dd.name}; 
            data_3_0C_tmean = cell(numel(fileNames),2);
            data_3_0C_tmean(:,1) = regexprep(fileNames, '_UKCP18-12km_T.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)   
                fname = fullfile('ARCADIA-Tmean-BC-30', fileNames(ii)); %allocate each .csv data file to cell array
                data_3_0C_tmean{ii,2} = dlmread(fname{1},',',0,0);
                data_3_0C_tmean{ii,2}(:,1461:1461:end) = []; %delete leap year, every 1461 columns. Bias corrected data uses Gregorian Calendar.
                data_3_0C_tmean{ii,2}(:,10950) = 0; %Add back in dummy value for 31st December so years are even.
            end
            
            save('input_data_3_0C_tmean.mat', 'data_3_0C_tmean', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together

        else
            disp ('No 3.0DegreeC files exist')
        end
    end
    
    % 4 degree
    if (ClimateScenario == 5 || ClimateScenario == 6) %% 4 degree or all climate scenarios
         if exist ('input_data_4_0C_tmean.mat', 'file')  %%remove for running on DAFNI
           load ('input_data_4_0C_tmean.mat')           %
           disp ('loading input data 4.0C')             %
           nClimScen = nClimScen+1; %counts # climate scenarios being run together 
        
         elseif exist('ARCADIA-Tmean-BC-40', 'file')
            disp ('Read in daily Tmean data (4.0DegreeC) for mortality estimates')

            dd = dir('ARCADIA-Tmean-BC-40\*.csv');
            fileNames = {dd.name}; 
            data_4_0C_tmean = cell(numel(fileNames),2);
            data_4_0C_tmean(:,1) = regexprep(fileNames, '_UKCP18-12km_T.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)   
                fname = fullfile('ARCADIA-Tmean-BC-40', fileNames(ii)); %allocate each .csv data file to cell array
                data_4_0C_tmean{ii,2} = dlmread(fname{1},',',0,0);
                data_4_0C_tmean{ii,2}(:,1461:1461:end) = []; %delete leap year, every 1461 columns. Bias corrected data uses Gregorian Calendar.
                data_4_0C_tmean{ii,2}(:,10950) = 0; %Add back in dummy value for 31st December so years are even.
            end
            
            save('input_data_4_0C_tmean.mat', 'data_4_0C_tmean', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together
        else
            disp ('No 4.0DegreeC files exist')
        end
    end
end

%% Set lat and long for use in mapping outputs and linking population and building data to TMean data
%(format DS --> currently convert to DMS in ArcGIS for plotting)

disp('calculating gridID data')
load lat_UK_RCM
load long_UK_RCM

%set array size based on the climate scenario selected.
if(ClimateScenario == 1 || ClimateScenario == 6)
    [nGridCells,nVar] = size(data_past_tmean);
    HEAT_Data = data_past_tmean;
elseif (ClimateScenario == 2)
    [nGridCells,nVar] = size(data_1_5C_tmean); 
    HEAT_Data = data_1_5C_tmean;
elseif (ClimateScenario == 3)
    [nGridCells,nVar] = size(data_2_0C_tmean); 
    HEAT_Data = data_2_0C_tmean;
elseif (ClimateScenario == 4)
    [nGridCells,nVar] = size(data_3_0C_tmean); 
    HEAT_Data = data_3_0C_tmean;
elseif (ClimateScenario == 5)
    [nGridCells,nVar] = size(data_4_0C_tmean); 
    HEAT_Data = data_4_0C_tmean;
end
    
    cellIndex = split(HEAT_Data(:,1),"_"); % remove hyphen. cellIndex used to pull out corresponding land-based lon_lat from UK & ROI grid
    cellIndex = str2double(cellIndex);
    gridID = zeros(nGridCells,2); % sets corresponding lon & lat for gridCells for plotting

for i = 1:nGridCells
    x = cellIndex(i,1);
    y = cellIndex(i,2);
    gridID(i,2) = lat_UK_RCM(x,y); %corresponding lon_lat for each cell read in. Corresponds to rows of [data_climScen] and [population].
    gridID(i,1) = long_UK_RCM(x,y);
end
gridID = fix(gridID*1e5)/1e5;

lonIDIndex = fix(building_thresholds{1}(:,1)*1e5)/1e5; %truncate decimals to allow match to gridID - different order to IDs in gridID
latIDIndex = fix(building_thresholds{1}(:,2)*1e5)/1e5; %truncate decimals to allow match to gridID - different order to IDs in gridID

save('gridID', 'gridID', 'lonIDIndex', 'latIDIndex')
clear cellIndex %not used past this point - use gridID

%% Define size of variables / Cell Arrays for below calculations
    [nRows,nCols] = size(HEAT_Data{1,2}); %set up TMean data array size
    years = nCols/daysPerYear; %to check years correspond to 30
 
%% CALCULATE RESIDENTIAL OVERHEATING
%loop through gridded TMean data to calculate overheating risk in bedrooms
%with and without socio-eocnomic change, with and without adaptation, and for all nClimScens.

h=1; %h used to count through nGridCells.
avDaysExceed = zeros(nRows,nGridCells); % days that exceed TMean thresholds by build type (intermediate output)
AnnualFrequencyEvent = zeros(nRows,nGridCells); % probability of residential overheating risk (final output)
ResultsExceedDays = zeros(nGridCells,14); %not cell for BL as no adaptation/SE change. Intermediate outputs  - averaged across years and RCMS. Av, Median, 10th and 90th percentile.
ResultsAnnualFrequency = zeros(nGridCells,6); %not cell for BL as no adaptation/SE change - across all build types. Final outputs  - averaged across years and RCMS. Av, Median, 10th and 90th percentile.
nBuilds = 4;  %Flats, detached, semi-detached, Terrace - pre-determined by UDM categories.
adaptScenarios = 2; %results run with and without adaptation (based on option selected by user) qunderpinned by UCL data.
consecutive_days = 0;

% Based on user set up - dependent on the year socio-economic data is selected for, the deployment of adaptation action (uptake) will differ(from UCL)
if populationYear15degree == 2020 || populationYear2degree == 2020 || populationYear3degree == 2020 || populationYear4degree == 2020
    deploymentYear = 2;
elseif populationYear15degree == 2030 || populationYear2degree == 2030 || populationYear3degree == 2030 || populationYear4degree == 2030
    deploymentYear = 3;
elseif populationYear15degree == 2050 || populationYear2degree == 2050 || populationYear3degree == 2050 || populationYear4degree == 2050
    deploymentYear = 4;
elseif populationYear15degree == 2080 || populationYear2degree == 2080 || populationYear3degree == 2080 || populationYear4degree == 2080
    deploymentYear = 5;
else
    deploymentYear = 1; %past/basleine
end

% select correct building data - can remove when on DAFNI as will be defined in user parameter set-up?
if UK_SSP == 1
    building_numbers = building_numbers_SSP1;
elseif UK_SSP == 2
    building_numbers = building_numbers_SSP2;
elseif UK_SSP == 3
    building_numbers = building_numbers_SSP3;
elseif UK_SSP == 4
    building_numbers = building_numbers_SSP4;
elseif UK_SSP == 5
    building_numbers = building_numbers_SSP5;
end

%%Calculate residential discomfort: BASELINE - climate change only (no adaptation or SE change in BL run)
if (ClimateScenario == 1 || ClimateScenario == 6) %% baseline or all climate scenarios
    if exist ('input_data_past_tmean.mat', 'file') 
        disp ('Calculating residential overheating risk - Near Past')

        %1. THRESHOLD DETECTOR - pulls out DAYS where TMean thresholds are exceeded based on building_threshold_.mat underpinned by data
        % from UCL model for different time periods {past, 2030s, 2050s, 2080s} without and with adaptation for each build type
        % building_thresholds = [lon, lat, GOR, flat, detached, semi-detached, terrace (no adaptation), flat, detached, semi-detached, terrace (with adaptation)]
        ppl_per_dwellings = zeros(nGridCells, 1); %average ppl per household - UDM uses 2.5ppl per dwelling. Calculated here for each gridcell based on total dwellings and population
        cnt = 3; %to set correct columns for recording data in table for each building type.

        %pre-allocate arrays
        annualEvents = cell(nGridCells,1);   %store intermediate value of days per year TMax threshold exceeded.    
        for n = 1:nGridCells
            annualEvents{n} = zeros(nRows, nCols);
        end

        ppl_per_event_thres = cell(nGridCells, 1); %ppl affected after x-conscutive days totaled across all build types affected, per gridcell
        for n = 1:nGridCells
            ppl_per_event_thres{n} = zeros(nRows, nCols);
        end 

        ppl_per_event = cell(nGridCells, 1); %ppl affected after 5-conscutive days totalled across all build types affected, for all gridcells
        for n = 1:nGridCells
            ppl_per_event{n} = zeros(nRows, nCols);
        end
        
        events = cell(nGridCells, 1); %ppl affected after 5-conscutive days totalled across all build types affected, for all gridcells
        for n = 1:nGridCells
            events{n} = zeros(nRows, nCols);
        end
        
        %Go through each build type/threshold in turn...
        for b = 1:nBuilds 

            %pre-allocated arrays --> reset for each building type.
            days_exceeded = cell(nGridCells, 1); %days exceeded, for all gridcells by build types
            for n = 1:nGridCells
                days_exceeded{n} = zeros(nRows, nCols);
            end

            annDaysExceed = cell(nGridCells,1);   %store intermediate value of days per year TMax threshold exceeded.    
            for n = 1:nGridCells
                annDaysExceed{n} = zeros(nRows, nCols);
            end

            i=1;
            while (h <= nGridCells)
                if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) %looks through cellID to match to data_past_tmean to building/pop/threshold data as some grid cells may be missing, or may be in differet order in future files etc.

                    if b==1 % just calculate once. Calculates average ppl per dwellings across all dwelling types.
                        if population{1}(h,3)> 0 && sum(building_numbers{1}(h,3:6)) > 0 %avoid dividing by zeros.
                        ppl_per_dwellings(h,1) = population{1}(h,3)/sum(building_numbers{1}(h,3:6));  %equally distributed across build types. UDM uses uniform value of 2.5ppl per dwelling.    
                        end
                    end
                    %Threshold detector
                    for r = 1:nRows
                        for c = 1:nCols
                            if data_past_tmean{i,2}(r,c) < building_thresholds{1}(h,b+3)
                            % do nothing - speeds up if statement by specifying most likey option first..
                            elseif data_past_tmean{i,2}(r,c) >= building_thresholds{1}(h,b+3)
                                days_exceeded{h}(r,c) = 1;
                                ppl_per_event{h}(r,c) = ppl_per_event{h}(r,c)+(ppl_per_dwellings(h,1)*building_numbers{1}(h,b+2)); %cumulative ppl at risk per gridcell, increments across building types.
                            end
                        end
                    end

                    h=h+1;
                    i=1; %reset counter
                    r=1; %reset counter
                    c=1; %reset counter

                else
                    i=i+1;
                end
            end

            %2. CALCULATE days per year thresholds exceeded for each building type (per grid cell)[Temp threshold output that can
            %be plotted/mapped per buidling type (intermediate output).
            for h = 1:nGridCells
                annDaysExceed{h} = squeeze(nansum(reshape(days_exceeded{h}.',daysPerYear,size(days_exceeded{h},2)./daysPerYear,[]))).'; %re-shape array
                for r = 1:nRows
                    avDaysExceed(r,h)= mean(annDaysExceed{h}(r,:),[1,2]); %mean across years for each GCM
                end
                    percentile_low = prctile(avDaysExceed,(10),1);%10th percentile across years/GCMs
                    percentile_high = prctile(avDaysExceed,(90),1);%90th percentile across years/GCMs
                    av = mean(avDaysExceed,1);%average
            end

            %3. Save threshold results for plotting: gridded - by building type (b) across columns
            if b == 1
                ResultsExceedDays(1:h,cnt:cnt+2) = [av.', percentile_low.', percentile_high.']; %3 variables[nGridCells*3]
                ResultsExceedDays(1:h,1) = lonIDIndex(:,1); %longitude
                ResultsExceedDays(1:h,2) = latIDIndex(:,1); %latitude
                cnt = cnt+3; %increment columns to record results for each building type.

            elseif b > 1 && b <=nBuilds
                ResultsExceedDays(1:h,cnt:cnt+2) = [av.', percentile_low.', percentile_high.']; %3 variables[nGridCells*3]
                cnt = cnt+3;
            end
            h=1; % re-set counter for looping through next building type    
        end

        % 4. When all building types have run and ppl_per_event calculated, identify where residents affected by discomfort meet overheating
        % event criteria based on user-defined variable (RDProbabilityLevel).
        if b == nBuilds % check all build types have run;
            for h = 1:nGridCells
                for r = 1:nRows
                    for c = 1:nCols
                        if ppl_per_event{h}(r,c) > 0 && ppl_per_event{h}(r,c)/population{1}(h,3) >= RDProbabilityLevel    % checks if % of total population per grid cell affected
                        ppl_per_event_thres{h}(r,c) = 1; % Is % population affected threshold met? 0 = False, 1 = True.
                        end
                    end
                end
            end
        end    
        
        %5. Identify occurences where the number of consecutive days where ppl_per_event_thres = 1 >= RDNumberConsecutiveDays. This defines the overheating 'event'
        for h = 1:nGridCells
            for r = 1:nRows
                for c = 1:nCols
                     if ppl_per_event_thres{h}(r,c) == 0
                        % do nothing - speeds up if statement by specifying most likely option first..
                        consecutive_days = 0;
                     elseif ppl_per_event_thres{h}(r,c) == 1
                        consecutive_days = consecutive_days+1;
                        if consecutive_days == RDNumberConsecutiveDays
                            events{h}(r,c) = 1;
                            consecutive_days = 0;
                        end
                     end
                end
                consecutive_days = 0;
            end
        end
        
        %6. CALCULATE annual probability of 'events' occuring (per grid cell).
        for h = 1:nGridCells  %[12*30]
            annualEvents{h} = squeeze(nansum(reshape(events{h}.',daysPerYear,size(events{h},2)./daysPerYear,[]))).'; %re-shape array

            for r = 1:nRows
                AnnualFrequencyEvent(r,h)= mean(annualEvents{h}(r,:),[1,2]); %mean across 30 years per RCM.
            end
                percentile_low = prctile(AnnualFrequencyEvent,(10),1);%10th percentile across years/GCMs
                percentile_high = prctile(AnnualFrequencyEvent,(90),1);%90th percentile across years/GCMs
                median = prctile(AnnualFrequencyEvent,(50),1);%50th percentile across years/GCMs
                av = mean(AnnualFrequencyEvent,1);%average
        end

        %7. Save results for plotting: gridded
        ResultsAnnualFrequency(:,1) = lonIDIndex(:,1); %longitude
        ResultsAnnualFrequency(:,2) = latIDIndex(:,1); %latitude
        ResultsAnnualFrequency(:,3:6) = [av.', median.', percentile_low.', percentile_high.'];

         % Files to SAVE
         fname = sprintf('RD_ResultsPast_%d_%0.2f_%d.mat', populationYearbaseline, RDProbabilityLevel, RDNumberConsecutiveDays);
         ResultsExceedDays_past = ResultsExceedDays(:,:);
         ResultsAnnualFrequencyEvents = ResultsAnnualFrequency(:,:);
         save(fname, 'ResultsAnnualFrequencyEvents', 'ResultsExceedDays_past')

         h=1; % re-set here if ClimateScenario == 6 and code continues...
    else
        disp ('no input data - past')

    end
end

%% Results for 1.5 degree for climate change only AND climate change and SE change (both with and without adaptation)
if (ClimateScenario == 2 || ClimateScenario == 6) %% 1.5degree or all climate scenarios
    if exist ('input_data_1_5C_tmean.mat', 'file')
    disp ('Calculating residential overheating risk - 1.5DegreeC')
    
    cnt = 3; %to set correct columns for recording data in table for each build type.
    ppl_per_dwellings = zeros(nGridCells, 1); %average ppl per household - UDM uses 2.5ppl per dwelling. Calculated here for each gridcell based on total dwellings and population
    p=0; %counter to make sure the correct population cell array is read for years when warming set to occur by user. First run climate change only so p=0.
    
        for l = 1:socioEcScen % Compute for CC-only and then SE+CC (population and dwelling numbers)

            %pre-allocate arrays - cell arrays as with and without adaptation
            ResultsExceedDays = cell(adaptScen,1); %Intermediate outputs  - averaged across years and RCMS. Av, Median, 10th and 90th percentile.
            for n = 1:adaptScen
                ResultsExceedDays{n} = zeros(nGridCells, 14);
            end
            
            ResultsAnnualFrequency = cell(adaptScen,1); %Final outputs  - averaged across years and RCMS. Av, Median, 10th and 90th percentile.
            for n = 1:adaptScen
                ResultsAnnualFrequency{n} = zeros(nGridCells, 6);
            end
            
            ppl_per_event = cell(nGridCells, 1); %ppl affected after 5-conscutive days totalled across all build types affected, for all gridcells
            for n = 1:nGridCells
                ppl_per_event{n} = zeros(nRows, nCols);
            end

            ppl_per_event_adapt = cell(nGridCells, 1); %ppl affected after 5-conscutive days totalled across all build types affected, for all gridcells
            for n = 1:nGridCells
                ppl_per_event_adapt{n} = zeros(nRows, nCols);
            end

            for a = 1:adaptScen % run without and then with adapation (different threshold values)
                if a == 1 %no adaptation
                    adaptCnt = 3; %select correct thresholds from columns
                elseif a == 2 %adaptation
                    adaptCnt = 7; %select correct thresholds from columns
                end
                
                % pre-allocate cell arrays for intermediate values        
                events = cell(nGridCells, 1); %number of events, caulculated where x% population affected for x-number consecutive days.
                for n = 1:nGridCells
                    events{n} = zeros(nRows, nCols);
                end
                
                annualEvents = cell(nGridCells,1);   %store intermediate value of events per year.    
                for n = 1:nGridCells
                    annualEvents{n} = zeros(nRows, nCols);
                end           

                ppl_per_event_thres = cell(nGridCells, 1); %ppl affected after 5-conscutive days totalled across all build types affected, for all gridcells
                for n = 1:nGridCells
                    ppl_per_event_thres{n} = zeros(nRows, nCols);
                end

            %1. THRESHOLD DETECTOR - pulls out DAYS where thresholds are exceeded based on building_thresholds.mat from UCL model team

                for b = 1:nBuilds
                    %pre-allocated arrays --> reset for each building type.
                    days_exceeded = cell(nGridCells, 1); %days exceeded, for all gridcells by build types
                    for n = 1:nGridCells
                    days_exceeded{n} = zeros(nRows, nCols);
                    end
                    
                    annDaysExceed = cell(nGridCells,1);   %store intermediate value of days per year TMax threshold exceeded.    
                    for n = 1:nGridCells
                        annDaysExceed{n} = zeros(nRows, nCols);
                    end
                    
                    i=1;
                    while (h <= nGridCells)
                        if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) %looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                           
                            if b==1 % just calculate once. Calculates average ppl per dwellings across all dwelling types.
                                if population{1+p}(h,3)> 0 && sum(building_numbers{1+p}(h,3:6)) > 0 %avoid dividing by zeros.
                                    ppl_per_dwellings(h,1) = population{1+p}(h,3)/sum(building_numbers{1+p}(h,3:6));  %equally distributed across build types. UDM uses uniform value of 2.5ppl per dwelling.    
                                end
                            end
                            %Threshold detector
                            for r = 1:nRows
                                for c = 1:nCols
                                    if data_1_5C_tmean{i,2}(r,c) < building_thresholds{1+p}(h,b+adaptCnt) %no adaptation and adaptation thresholds
                                    % do nothing - speeds up if statement by specifying most likey option first..
                                    elseif data_1_5C_tmean{i,2}(r,c) >= building_thresholds{1+p}(h,b+adaptCnt)
                                        days_exceeded{h}(r,c) = 1;
                                        if a == 1 % no adaptation of housing stock
                                            ppl_per_event{h}(r,c) = ppl_per_event{h}(r,c)+(ppl_per_dwellings(h,1)*building_numbers{1+p}(h,b+2)); %cumulative ppl at risk per gridcell, increments across building types.
                                        elseif a == 2 % adaptation of housing stock
                                            ppl_per_event_adapt{h}(r,c) = ppl_per_event_adapt{h}(r,c)+(ppl_per_dwellings(h,1)*building_numbers{1+p}(h,b+2)); %assumes 100% deployment across population and build type
                                        end
                                    end
                                end
                            end

                            h=h+1;
                            i=1; %reset counter
                            r=1; %reset counter
                            c=1; %reset counter

                        else
                            i=i+1;
                        end
                    end

                    %2. CALCULATE days per year thresholds exceeded for each building type (per grid cell)[Temp threshold output that can
                    %be plotted/mapped per buidling type (intermediate output).
                    for h = 1:nGridCells
                        annDaysExceed{h} = squeeze(nansum(reshape(days_exceeded{h}.',daysPerYear,size(days_exceeded{h},2)./daysPerYear,[]))).'; %re-shape array
                        for r = 1:nRows
                            avDaysExceed(r,h)= mean(annDaysExceed{h}(r,:),[1,2]); %mean across 30 years for each RCM
                        end
                            percentile_low = prctile(avDaysExceed,(10),1);%10th percentile across years/GCMs
                            percentile_high = prctile(avDaysExceed,(90),1);%90th percentile across years/GCMs
                            av = mean(avDaysExceed,1);%average
                    end

                    %3. Save threshold results for plotting: gridded - by building type (b) across columns
                    if b == 1
                        ResultsExceedDays{a}(1:h,cnt:cnt+2) = [av.', percentile_low.', percentile_high.']; %3 variables[nGridCells*3]
                        ResultsExceedDays{a}(1:h,1) = lonIDIndex(:,1); %longitude
                        ResultsExceedDays{a}(1:h,2) = latIDIndex(:,1); %latitude
                        cnt = cnt+3; %increment columns to record results for each building type.

                    elseif b > 1 && b <=nBuilds
                        ResultsExceedDays{a}(1:h,cnt:cnt+2) = [av.', percentile_low.', percentile_high.']; %3 variables[nGridCells*3]
                        cnt = cnt+3;
                    end
                    h=1; % re-set counter for looping through next building type
                end

                % 4. When all building types have run and ppl_per_event calculated, identify where residents affected by discomfort meet overheating
                % event criteria based on user-defined variable (RDProbabilityLevel).
                if b == nBuilds && a == 1 % once all build types have run; No adaptation.
                    for h = 1:nGridCells
                        for r = 1:nRows
                            for c = 1:nCols
                                if ppl_per_event{h}(r,c) > 0 && ppl_per_event{h}(r,c)/population{1+p}(h,3) >= RDProbabilityLevel    % RDProbabilityLevel set by user.
                                    ppl_per_event_thres{h}(r,c) = 1; %is % population affected threshold met? 0 = False, 1 = True.
                                end
                            end
                        end
                    end
                elseif b == nBuilds && a == 2 % once all build types have run; With adaptation. Account for deployment levels here (from UCL - 100% of houses not adapated in grid cell, a proportion will still be affected under thresholds for no adaptation).
                    for h = 1:nGridCells
                        for r = 1:nRows
                            for c = 1:nCols
                                if ppl_per_event{h}(r,c) > 0 && ppl_per_event_adapt{h}(r,c) > 0 && ppl_per_event{h}(r,c) == ppl_per_event_adapt{h}(r,c) && ppl_per_event_adapt{h}(r,c)/population{1+p}(h,3) >= RDProbabilityLevel    % RDProbabilityLevel set by user.
                                    ppl_per_event_thres{h}(r,c) = 1; % adaptation has not mitigated ANY risk.
                                elseif ppl_per_event{h}(r,c) > 0 && ppl_per_event_adapt{h}(r,c) == 0 && (ppl_per_event{h}(r,c)*(1-deploymentPercent(deploymentYear,1)/100))/population{1+p}(h,3) >= RDProbabilityLevel % where adaptation has mitigated event across all building types
                                    ppl_per_event_thres{h}(r,c) = 1; % adaptation mitigated ALL risk but 100% of houses are not adapted - calculates residual risk.
                                elseif ppl_per_event{h}(r,c) > 0 && ppl_per_event_adapt{h}(r,c) > 0 && (ppl_per_event{h}(r,c)*(1-((ppl_per_event_adapt{h}(r,c)/ppl_per_event{h}(r,c))/(100*deploymentPercent(deploymentYear,1)))))/population{1+p}(h,3) >= RDProbabilityLevel
                                    ppl_per_event_thres{h}(r,c) = 1; % where adaptation has mitigated event across ONLY SOME building types BUT NOT ALL so proportional
                                end       
                            end
                        end
                    end 
                end

                %5. Identify occurences where the number of consecutive days where ppl_per_event_thres = 1 >= RDNumberConsecutiveDays. This defines the overheating 'event'
                for h = 1:nGridCells
                    for r = 1:nRows
                        for c = 1:nCols
                             if ppl_per_event_thres{h}(r,c) == 0
                                % do nothing - speeds up if statement by specifying most likely option first..
                                consecutive_days = 0;
                             elseif ppl_per_event_thres{h}(r,c) == 1
                                consecutive_days = consecutive_days+1;
                                if consecutive_days == RDNumberConsecutiveDays
                                    events{h}(r,c) = 1;
                                    consecutive_days = 0;
                                end
                             end
                        end
                        consecutive_days = 0;
                    end
                end
                
                %6. CALCULATE annual probability of 'events' occuring (per grid cell).
                for h = 1:nGridCells  %[12*30]
                    annualEvents{h} = squeeze(nansum(reshape(events{h}.',daysPerYear,size(events{h},2)./daysPerYear,[]))).'; %re-shape array

                    for r = 1:nRows
                        AnnualFrequencyEvent(r,h)= mean(annualEvents{h}(r,:),[1,2]); %mean across years for each GCM
                    end
                        percentile_low = prctile(AnnualFrequencyEvent,(10),1);%10th percentile across years/GCMs
                        percentile_high = prctile(AnnualFrequencyEvent,(90),1);%90th percentile across years/GCMs
                        median = prctile(AnnualFrequencyEvent,(50),1);%50th percentile across years/GCMs
                        av = mean(AnnualFrequencyEvent,1);%average
                end

                %7. Save results for plotting: gridded
                ResultsAnnualFrequency{a}(:,1) = lonIDIndex(:,1); %longitude
                ResultsAnnualFrequency{a}(:,2) = latIDIndex(:,1); %latitude
                ResultsAnnualFrequency{a}(:,3:6) = [av.', median.', percentile_low.', percentile_high.'];
                
                cnt = 3; % re-set to calculate with adaptation.
                h=1; % re-set to calculate with adaptation.
            end    

     % Files to SAVE
            if l==1 %climate change only (with and without adaptation)
            fname = sprintf('RD_Results_1_5C_%s_%0.2f_%d_%0.2f.mat', '_CCOnly', RDProbabilityLevel, RDNumberConsecutiveDays, deploymentPercent(deploymentYear,1));
            ResultsAnnualFrequencyEvents_1_5_CC = ResultsAnnualFrequency;
            ResultsExceedDays_1_5_CC = ResultsExceedDays; % reflects exceedance of TMean and role of adaptation on changing thrsholds. SE change has no effect here.
            save(fname, 'ResultsAnnualFrequencyEvents_1_5_CC', 'ResultsExceedDays_1_5_CC')
            
            elseif l==2 %Climate change and population change (with and without adaptation)
            fname = sprintf('RD_Results_1_5C_%s_%s_%d_SSP%d_%0.2f_%d_%0.2f.mat', '_SE', 'Adapt', populationYear15degree, UK_SSP, RDProbabilityLevel, RDNumberConsecutiveDays, deploymentPercent(deploymentYear,1));
            ResultsAnnualFrequencyEvents_1_5_SE = ResultsAnnualFrequency; %ResultsExceedDays doesn't change with SE scenario.
            save(fname, 'ResultsAnnualFrequencyEvents_1_5_SE')
            end
            
            if populationYear15degree == 2020
                p=1;
            elseif populationYear15degree == 2030
                p=2;
            elseif populationYear15degree == 2050
                p=3;
            elseif populationYear15degree == 2080 || populationYear15degree == 2100 %2080, or 2100 can be used for end-century. TBC
                p=4;
            end
        end
        h=1; % re-set here if ClimateScenario == 6 and code continues...
    else
        disp ('no input data - 1.5C')
    end
end

%% Below calculates results for 2 degree for climate change only AND climate change and SE change (both with and without adaptation)
if (ClimateScenario == 3 || ClimateScenario == 6) %% 2degree or all climate scenarios
    if exist ('input_data_2_0C_tmean.mat', 'file')
    disp ('Calculating residential overheating risk - 2DegreeC')
    
    cnt = 3; %to set correct columns for recording data in table for each build type.
    ppl_per_dwellings = zeros(nGridCells, 1); %average ppl per household - UDM uses 2.5ppl per dwelling. Calculated here for each gridcell based on total dwellings and population
    p=0; %counter to make sure the correct population cell array is read for years when warming set to occur by user. First run climate change only so p=0.
      
        for l = 1:socioEcScen % Compute for CC-only and then SE+CC (population and dwelling numbers)

            %pre-allocate arrays - cell arrays as with and without adaptation
            ResultsExceedDays = cell(adaptScen,1); %Intermediate outputs  - averaged across years and RCMS. Av, Median, 10th and 90th percentile.
            for n = 1:adaptScen
                ResultsExceedDays{n} = zeros(nGridCells, 14);
            end
            
            ResultsAnnualFrequency = cell(adaptScen,1); %Final outputs  - averaged across years and RCMS. Av, Median, 10th and 90th percentile.
            for n = 1:adaptScen
                ResultsAnnualFrequency{n} = zeros(nGridCells, 6);
            end
            
            ppl_per_event = cell(nGridCells, 1); %ppl affected after 5-conscutive days totalled across all build types affected, for all gridcells
            for n = 1:nGridCells
                ppl_per_event{n} = zeros(nRows, nCols);
            end

            ppl_per_event_adapt = cell(nGridCells, 1); %ppl affected after 5-conscutive days totalled across all build types affected, for all gridcells
            for n = 1:nGridCells
                ppl_per_event_adapt{n} = zeros(nRows, nCols);
            end

            for a = 1:adaptScen % run without and then with adapation (different threshold values)
                if a == 1 %no adaptation
                    adaptCnt = 3; %select correct thresholds from columns
                elseif a == 2 %adaptation
                    adaptCnt = 7; %select correct thresholds from columns
                end
                
                % pre-allocate cell arrays for intermediate values        
                events = cell(nGridCells, 1); %number of events, caulculated where x% population affected for x-number consecutive days.
                for n = 1:nGridCells
                    events{n} = zeros(nRows, nCols);
                end
                
                annualEvents = cell(nGridCells,1);   %store intermediate value of days per year TMax threshold exceeded.    
                for n = 1:nGridCells
                    annualEvents{n} = zeros(nRows, nCols);
                end           

                ppl_per_event_thres = cell(nGridCells, 1); %ppl affected after 5-conscutive days totalled across all build types affected, for all gridcells
                for n = 1:nGridCells
                    ppl_per_event_thres{n} = zeros(nRows, nCols);
                end

            %1. THRESHOLD DETECTOR - pulls out DAYS where thresholds are exceeded based on building_thresholds.mat from UCL model team
                for b = 1:nBuilds
                    %pre-allocated arrays --> reset for each building type.
                    days_exceeded = cell(nGridCells, 1); %days exceeded, for all gridcells by build types
                    for n = 1:nGridCells
                    days_exceeded{n} = zeros(nRows, nCols);
                    end
                    
                    annDaysExceed = cell(nGridCells,1);   %store intermediate value of days per year TMax threshold exceeded.    
                    for n = 1:nGridCells
                        annDaysExceed{n} = zeros(nRows, nCols);
                    end
                    
                    i=1;
                    while (h <= nGridCells)
                        if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) %looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.

                        if b==1 % just calculate once. Calculates average ppl per dwellings across all dwelling types.
                           if population{1+p}(h,3)> 0 && sum(building_numbers{1+p}(h,3:6)) > 0 %avoid dividing by zeros.
                                ppl_per_dwellings(h,1) = population{1+p}(h,3)/sum(building_numbers{1+p}(h,3:6));  %equally distributed across build types. UDM uses uniform value of 2.5ppl per dwelling.    
                           end
                        end
                        %Threshold detector
                        for r = 1:nRows
                            for c = 1:nCols
                                if data_2_0C_tmean{i,2}(r,c) < building_thresholds{1+p}(h,b+adaptCnt) %no adaptation and adaptation thresholds
                                % do nothing - speeds up if statement by specifying most likey option first..
                                elseif data_2_0C_tmean{i,2}(r,c) >= building_thresholds{1+p}(h,b+adaptCnt)
                                days_exceeded{h}(r,c) = 1;
                                    if a == 1 % no adaptation of housing stock
                                        ppl_per_event{h}(r,c) = ppl_per_event{h}(r,c)+(ppl_per_dwellings(h,1)*building_numbers{1+p}(h,b+2)); %cumulative ppl at risk per gridcell, increments across building types.
                                    elseif a == 2 % 100% adaptation of housing stock
                                        ppl_per_event_adapt{h}(r,c) = ppl_per_event_adapt{h}(r,c)+(ppl_per_dwellings(h,1)*building_numbers{1+p}(h,b+2)); %assumes 100% deployment across population and build type
                                    end
                                end
                            end
                        end

                        h=h+1;
                        i=1; %reset counter
                        r=1; %reset counter
                        c=1; %reset counter

                        else
                        i=i+1;
                        end
                    end

                    %2. CALCULATE days per year thresholds exceeded for each building type (per grid cell)[Temp threshold output that can
                    %be plotted/mapped per buidling type (intermediate output).
                    for h = 1:nGridCells
                        annDaysExceed{h} = squeeze(nansum(reshape(days_exceeded{h}.',daysPerYear,size(days_exceeded{h},2)./daysPerYear,[]))).'; %re-shape array
                        for r = 1:nRows
                            avDaysExceed(r,h)= mean(annDaysExceed{h}(r,:),[1,2]); %mean across years for each GCM
                        end
                            percentile_low = prctile(avDaysExceed,(10),1);%10th percentile across years/GCMs
                            percentile_high = prctile(avDaysExceed,(90),1);%90th percentile across years/GCMs
                            av = mean(avDaysExceed,1);%average
                    end

                    %3. Save threshold results for plotting: gridded - by building type (b) across columns
                    if b == 1
                        ResultsExceedDays{a}(1:h,cnt:cnt+2) = [av.', percentile_low.', percentile_high.']; %3 variables[nGridCells*3]
                        ResultsExceedDays{a}(1:h,1) = lonIDIndex(:,1); %longitude
                        ResultsExceedDays{a}(1:h,2) = latIDIndex(:,1); %latitude
                        cnt = cnt+3; %increment columns to record results for each building type.

                    elseif b > 1 && b <=nBuilds
                        ResultsExceedDays{a}(1:h,cnt:cnt+2) = [av.', percentile_low.', percentile_high.']; %3 variables[nGridCells*3]
                        cnt = cnt+3;
                    end
                    h=1; % re-set counter for looping through next building type
                end

                % 4. When all building types have run and ppl_per_event calculated, identify where residents affected by discomfort meet overheating
                % event criteria based on user-defined variable (RDProbabilityLevel).
                if b == nBuilds && a == 1 % once all build types have run; No adaptation.
                    for h = 1:nGridCells
                        for r = 1:nRows
                            for c = 1:nCols
                                if ppl_per_event{h}(r,c) > 0 && ppl_per_event{h}(r,c)/population{1+p}(h,3) >= RDProbabilityLevel    % RDProbabilityLevel set by user.
                                ppl_per_event_thres{h}(r,c) = 1; %is % population affected threshold met? 0 = False, 1 = True.
                                end
                            end
                        end
                    end
                elseif b == nBuilds && a == 2 % once all build types have run; With adaptation. Account for deployment levels here (from UCL - 100% of houses not adapated in grid cell, a proportion will still be affected under thresholds for no adaptation). 
                    for h = 1:nGridCells
                        for r = 1:nRows
                            for c = 1:nCols
                                if ppl_per_event{h}(r,c) > 0 && ppl_per_event_adapt{h}(r,c) > 0 && ppl_per_event{h}(r,c) == ppl_per_event_adapt{h}(r,c) && ppl_per_event_adapt{h}(r,c)/population{1+p}(h,3) >= RDProbabilityLevel    % RDProbabilityLevel set by user.
                                    ppl_per_event_thres{h}(r,c) = 1; % adaptation has not mitigated ANY risk.
                                elseif ppl_per_event{h}(r,c) > 0 && ppl_per_event_adapt{h}(r,c) == 0 && (ppl_per_event{h}(r,c)*(1-deploymentPercent(deploymentYear,1)/100))/population{1+p}(h,3) >= RDProbabilityLevel % where adaptation has mitigated event across all building types
                                    ppl_per_event_thres{h}(r,c) = 1; % adaptation mitigated ALL risk but 100% of houses are not adapted - calculates residual risk.
                                elseif ppl_per_event{h}(r,c) > 0 && ppl_per_event_adapt{h}(r,c) > 0 && (ppl_per_event{h}(r,c)*(1-((ppl_per_event_adapt{h}(r,c)/ppl_per_event{h}(r,c))/(100*deploymentPercent(deploymentYear,1)))))/population{1+p}(h,3) >= RDProbabilityLevel
                                    ppl_per_event_thres{h}(r,c) = 1; % where adaptation has mitigated event across ONLY SOME building types BUT NOT ALL so proportional
                                end           
                            end
                        end
                    end 
                end

                %5. Identify occurences where the number of consecutive days where ppl_per_event_thres = 1 >= RDNumberConsecutiveDays. This defines the overheating 'event'
                for h = 1:nGridCells
                    for r = 1:nRows
                        for c = 1:nCols
                             if ppl_per_event_thres{h}(r,c) == 0
                                % do nothing - speeds up if statement by specifying most likely option first..
                                consecutive_days = 0;
                             elseif ppl_per_event_thres{h}(r,c) == 1
                                consecutive_days = consecutive_days+1;
                                if consecutive_days == RDNumberConsecutiveDays
                                    events{h}(r,c) = 1;
                                    consecutive_days = 0;
                                end
                             end
                        end
                        consecutive_days = 0;
                    end
                end
                
                %6. CALCULATE annual probability of 'events' occuring (per grid cell).
                for h = 1:nGridCells  %[12*30]
                    annualEvents{h} = squeeze(nansum(reshape(events{h}.',daysPerYear,size(events{h},2)./daysPerYear,[]))).'; %re-shape array

                    for r = 1:nRows
                        AnnualFrequencyEvent(r,h)= mean(annualEvents{h}(r,:),[1,2]); %mean across years for each GCM. Number of recorded events/total possible number of events
                    end
                        percentile_low = prctile(AnnualFrequencyEvent,(10),1);%10th percentile across years/GCMs
                        percentile_high = prctile(AnnualFrequencyEvent,(90),1);%90th percentile across years/GCMs
                        median = prctile(AnnualFrequencyEvent,(50),1);%50th percentile across years/GCMs
                        av = mean(AnnualFrequencyEvent,1);%average
                end

                %7. Save results for plotting: gridded
                ResultsAnnualFrequency{a}(:,1) = lonIDIndex(:,1); %longitude
                ResultsAnnualFrequency{a}(:,2) = latIDIndex(:,1); %latitude
                ResultsAnnualFrequency{a}(:,3:6) = [av.', median.', percentile_low.', percentile_high.'];
            
                cnt = 3; % re-set to calculate with adaptation.
                h=1; % re-set to calculate with adaptation.
            end    

     % Files to SAVE
            if l==1 %climate change only (with and without adaptation)
            fname = sprintf('RD_Results_2_0C_%s_%0.2f_%d_%0.2f.mat', '_CCOnly', RDProbabilityLevel, RDNumberConsecutiveDays, deploymentPercent(deploymentYear,1));
            ResultsAnnualFrequencyEvents_2_0_CC = ResultsAnnualFrequency;
            ResultsExceedDays_2_0_CC = ResultsExceedDays; % reflects exceedance of TMean and role of adaptation on changing thrsholds. SE change has no effect here.
            save(fname, 'ResultsAnnualFrequencyEvents_2_0_CC', 'ResultsExceedDays_2_0_CC')
            
            elseif l==2 %Climate change and population change (with and without adaptation)
            fname = sprintf('RD_Results_2_0C_%s_%s_%d_SSP%d_%0.2f__%d_%0.2f.mat', '_SE', 'Adapt', populationYear2degree, UK_SSP, RDProbabilityLevel, RDNumberConsecutiveDays, deploymentPercent(deploymentYear,1));
            ResultsAnnualFrequencyEvents_2_0_SE = ResultsAnnualFrequency;
            save(fname, 'ResultsAnnualFrequencyEvents_2_0_SE')
            end
            
            if populationYear2degree == 2020
                p=1;
            elseif populationYear2degree == 2030
                p=2;
            elseif populationYear2degree == 2050
                p=3;
            elseif populationYear2degree == 2080 || populationYear2degree == 2100 %2080, or 2100 can be used for end-century. TBC
                p=4;
            end
        end
        h=1; % re-set here if ClimateScenario == 6 and code continues...
    else
        disp ('no input data - 2.0C')
    end
end

%% Below calculates results for 3 degree for climate change only AND climate change and SE change (both with and without adaptation)
if (ClimateScenario == 4 || ClimateScenario == 6) %% 3degree or all climate scenarios
    if exist ('input_data_3_0C_tmean.mat', 'file')
    disp ('Calculating residential overheating risk - 3DegreeC')
    
    cnt = 3; %to set correct columns for recording data in table for each build type.
    ppl_per_dwellings = zeros(nGridCells, 1); %average ppl per household - UDM uses 2.5ppl per dwelling. Calculated here for each gridcell based on total dwellings and population
    p=0; %counter to make sure the correct population cell array is read for years when warming set to occur by user. First run climate change only so p=0.
        
        for l = 1:socioEcScen % Compute for CC-only and then SE+CC (population and dwelling numbers)

            %pre-allocate arrays - cell arrays as with and without adaptation
            ResultsExceedDays = cell(adaptScen,1); %Intermediate outputs  - averaged across years and RCMS. Av, Median, 10th and 90th percentile.
            for n = 1:adaptScen
                ResultsExceedDays{n} = zeros(nGridCells, 14);
            end
            
            ResultsAnnualFrequency = cell(adaptScen,1); %Final outputs  - averaged across years and RCMS. Av, Median, 10th and 90th percentile.
            for n = 1:adaptScen
                ResultsAnnualFrequency{n} = zeros(nGridCells, 6);
            end
            
            ppl_per_event = cell(nGridCells, 1); %ppl affected after 5-conscutive days totalled across all build types affected, for all gridcells
            for n = 1:nGridCells
                ppl_per_event{n} = zeros(nRows, nCols);
            end

            ppl_per_event_adapt = cell(nGridCells, 1); %ppl affected after 5-conscutive days totalled across all build types affected, for all gridcells
            for n = 1:nGridCells
                ppl_per_event_adapt{n} = zeros(nRows, nCols);
            end

            for a = 1:adaptScen % run without and then with adapation (different threshold values)
                if a == 1 %no adaptation
                    adaptCnt = 3; %select correct thresholds from columns
                elseif a == 2 %adaptation
                    adaptCnt = 7; %select correct thresholds from columns
                end
                
                % pre-allocate cell arrays for intermediate values    
                events = cell(nGridCells, 1); %number of events, caulculated where x% population affected for x-number consecutive days.
                for n = 1:nGridCells
                    events{n} = zeros(nRows, nCols);
                end
    
                annualEvents = cell(nGridCells,1);   %store intermediate value of days per year TMax threshold exceeded.    
                for n = 1:nGridCells
                    annualEvents{n} = zeros(nRows, nCols);
                end           

                ppl_per_event_thres = cell(nGridCells, 1); %ppl affected after 5-conscutive days totalled across all build types affected, for all gridcells
                for n = 1:nGridCells
                    ppl_per_event_thres{n} = zeros(nRows, nCols);
                end

            %1. THRESHOLD DETECTOR - pulls out DAYS where thresholds are exceeded based on building_thresholds.mat from UCL model team
                for b = 1:nBuilds
                    %pre-allocated arrays --> reset for each building type.
                    days_exceeded = cell(nGridCells, 1); %days exceeded, for all gridcells by build types
                    for n = 1:nGridCells
                    days_exceeded{n} = zeros(nRows, nCols);
                    end
                    
                    annDaysExceed = cell(nGridCells,1);   %store intermediate value of days per year TMax threshold exceeded.    
                    for n = 1:nGridCells
                        annDaysExceed{n} = zeros(nRows, nCols);
                    end
                    
                    i=1;
                    while (h <= nGridCells)
                        if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) %looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                           
                            if b==1 % just calculate once. Calculates average ppl per dwellings across all dwelling types.
                                if population{1+p}(h,3)> 0 && sum(building_numbers{1+p}(h,3:6)) > 0 %avoid dividing by zeros.
                                    ppl_per_dwellings(h,1) = population{1+p}(h,3)/sum(building_numbers{1+p}(h,3:6));  %equally distributed across build types. UDM uses uniform value of 2.5ppl per dwelling.    
                                end
                            end
                            %Threshold detector
                            for r = 1:nRows
                                for c = 1:nCols
                                    if data_3_0C_tmean{i,2}(r,c) < building_thresholds{1+p}(h,b+adaptCnt) %no adaptation and adaptation thresholds
                                    % do nothing - speeds up if statement by specifying most likey option first..
                                    elseif data_3_0C_tmean{i,2}(r,c) >= building_thresholds{1+p}(h,b+adaptCnt)
                                    days_exceeded{h}(r,c) = 1;
                                        if a == 1 % no adaptation of housing stock
                                            ppl_per_event{h}(r,c) = ppl_per_event{h}(r,c)+(ppl_per_dwellings(h,1)*building_numbers{1+p}(h,b+2)); %cumulative ppl at risk per gridcell, increments across building types.
                                        elseif a == 2 % 100% adaptation of housing stock
                                            ppl_per_event_adapt{h}(r,c) = ppl_per_event_adapt{h}(r,c)+(ppl_per_dwellings(h,1)*building_numbers{1+p}(h,b+2)); %assumes 100% deployment across population and build type
                                        end
                                    end
                                end
                            end

                            h=h+1;
                            i=1; %reset counter
                            r=1; %reset counter
                            c=1; %reset counter
                        else
                            i=i+1;
                        end
                    end

                    %2. CALCULATE days per year thresholds exceeded for each building type (per grid cell)[Temp threshold output that can
                    %be plotted/mapped per buidling type (intermediate output).
                    for h = 1:nGridCells
                        annDaysExceed{h} = squeeze(nansum(reshape(days_exceeded{h}.',daysPerYear,size(days_exceeded{h},2)./daysPerYear,[]))).'; %re-shape array
                        for r = 1:nRows
                            avDaysExceed(r,h)= mean(annDaysExceed{h}(r,:),[1,2]); %mean across years for each GCM
                        end
                            percentile_low = prctile(avDaysExceed,(10),1);%10th percentile across years/GCMs
                            percentile_high = prctile(avDaysExceed,(90),1);%90th percentile across years/GCMs
                            av = mean(avDaysExceed,1);%average
                    end

                    %3. Save threshold results for plotting: gridded - by building type (b) across columns
                    if b == 1
                        ResultsExceedDays{a}(1:h,cnt:cnt+2) = [av.', percentile_low.', percentile_high.']; %3 variables[nGridCells*3]
                        ResultsExceedDays{a}(1:h,1) = lonIDIndex(:,1); %longitude
                        ResultsExceedDays{a}(1:h,2) = latIDIndex(:,1); %latitude
                        cnt = cnt+3; %increment columns to record results for each building type.

                    elseif b > 1 && b <=nBuilds
                        ResultsExceedDays{a}(1:h,cnt:cnt+2) = [av.', percentile_low.', percentile_high.']; %3 variables[nGridCells*3]
                        cnt = cnt+3;
                    end
                    h=1; % re-set counter for looping through next building type
                end

                % 4. When all building types have run and ppl_per_event calculated, identify where residents affected by discomfort meet overheating
                % event criteria based on user-defined variable (RDProbabilityLevel).
                if b == nBuilds && a == 1 % once all build types have run; No adaptation.
                    for h = 1:nGridCells
                        for r = 1:nRows
                            for c = 1:nCols
                                if ppl_per_event{h}(r,c) > 0 && ppl_per_event{h}(r,c)/population{1+p}(h,3) >= RDProbabilityLevel    % RDProbabilityLevel set by user.
                                    ppl_per_event_thres{h}(r,c) = 1;
                                end
                            end
                        end
                    end
                 elseif b == nBuilds && a == 2 % once all build types have run; With adaptation. Account for deployment levels here (from UCL - 100% of houses not adapated in grid cell, a proportion will still be affected under thresholds for no adaptation).
                    for h = 1:nGridCells
                        for r = 1:nRows
                            for c = 1:nCols
                                if ppl_per_event{h}(r,c) > 0 && ppl_per_event_adapt{h}(r,c) > 0 && ppl_per_event{h}(r,c) == ppl_per_event_adapt{h}(r,c) && ppl_per_event_adapt{h}(r,c)/population{1+p}(h,3) >= RDProbabilityLevel    % RDProbabilityLevel set by user.
                                    ppl_per_event_thres{h}(r,c) = 1; % adaptation has not mitigated ANY risk.
                                elseif ppl_per_event{h}(r,c) > 0 && ppl_per_event_adapt{h}(r,c) == 0 && (ppl_per_event{h}(r,c)*(1-deploymentPercent(deploymentYear,1)/100))/population{1+p}(h,3) >= RDProbabilityLevel % where adaptation has mitigated event across all building types
                                    ppl_per_event_thres{h}(r,c) = 1; % adaptation mitigated ALL risk but 100% of houses are not adapted - calculates residual risk.
                                elseif ppl_per_event{h}(r,c) > 0 && ppl_per_event_adapt{h}(r,c) > 0 && (ppl_per_event{h}(r,c)*(1-((ppl_per_event_adapt{h}(r,c)/ppl_per_event{h}(r,c))/(100*deploymentPercent(deploymentYear,1)))))/population{1+p}(h,3) >= RDProbabilityLevel
                                    ppl_per_event_thres{h}(r,c) = 1; % where adaptation has mitigated event across ONLY SOME building types BUT NOT ALL so proportional
                                end       
                            end
                        end
                    end 
                end
                 
                %5. Identify occurences where the number of consecutive days where ppl_per_event_thres = 1 >= RDNumberConsecutiveDays. This defines the overheating 'event'
                for h = 1:nGridCells
                    for r = 1:nRows
                        for c = 1:nCols
                             if ppl_per_event_thres{h}(r,c) == 0
                                % do nothing - speeds up if statement by specifying most likely option first..
                                consecutive_days = 0;
                             elseif ppl_per_event_thres{h}(r,c) == 1
                                consecutive_days = consecutive_days+1;
                                if consecutive_days == RDNumberConsecutiveDays
                                    events{h}(r,c) = 1;
                                    consecutive_days = 0;
                                end
                             end
                        end
                        consecutive_days = 0;
                    end
                end

                %6. CALCULATE annual probability of 'events' occuring (per grid cell).
                for h = 1:nGridCells  %[12*30]
                    annualEvents{h} = squeeze(nansum(reshape(events{h}.',daysPerYear,size(events{h},2)./daysPerYear,[]))).'; %re-shape array

                    for r = 1:nRows
                        AnnualFrequencyEvent(r,h)= mean(annualEvents{h}(r,:),[1,2]); %mean across years for each GCM
                    end
                        percentile_low = prctile(AnnualFrequencyEvent,(10),1);%10th percentile across years/GCMs
                        percentile_high = prctile(AnnualFrequencyEvent,(90),1);%90th percentile across years/GCMs
                        median = prctile(AnnualFrequencyEvent,(50),1);%50th percentile across years/GCMs
                        av = mean(AnnualFrequencyEvent,1);%average
                end

                %7. Save results for plotting: gridded
                ResultsAnnualFrequency{a}(:,1) = lonIDIndex(:,1); %longitude
                ResultsAnnualFrequency{a}(:,2) = latIDIndex(:,1); %latitude
                ResultsAnnualFrequency{a}(:,3:6) = [av.', median.', percentile_low.', percentile_high.'];
                
                cnt = 3; % re-set to calculate with adaptation.
                h=1; % re-set to calculate with adaptation.
            end    

     % Files to SAVE
            if l==1 %climate change only (with and without adaptation)
            fname = sprintf('RD_Results_3_0C_%s_%0.2f_%d_%0.2f.mat', '_CCOnly', RDProbabilityLevel, RDNumberConsecutiveDays, deploymentPercent(deploymentYear,1));
            ResultsAnnualFrequencyEvents_3_0_CC = ResultsAnnualFrequency;
            ResultsExceedDays_3_0_CC = ResultsExceedDays; % reflects exceedance of TMean and role of adaptation on changing thrsholds. SE change has no effect here.
            save(fname, 'ResultsAnnualFrequencyEvents_3_0_CC', 'ResultsExceedDays_3_0_CC')
            
            elseif l==2 %Climate change and population change (with and without adaptation)
            fname = sprintf('RD_Results_3_0C_%s_%s_%d_SSP%d_%0.2f__%d_%0.2f.mat', '_SE', 'Adapt', populationYear3degree, UK_SSP, RDProbabilityLevel, RDNumberConsecutiveDays, deploymentPercent(deploymentYear,1));
            ResultsAnnualFrequencyEvents_3_0_SE = ResultsAnnualFrequency;
            save(fname, 'ResultsAnnualFrequencyEvents_3_0_SE')
            end
            
            if populationYear3degree == 2020
                p=1;
            elseif populationYear3degree == 2030
                p=2;
            elseif  populationYear3degree == 2050
                p=3;
            elseif populationYear3degree == 2080 || populationYear3degree == 2100 %2080, or 2100 can be used for end-century. TBC
                p=4;
            end
        end
        h=1; % re-set here if ClimateScenario == 6 and code continues...
    else
        disp ('no input data - 3.0C')
    end
end

%% Below calculates results for 4 degree for climate change only AND climate change and SE change (both with and without adaptation)
if (ClimateScenario == 5 || ClimateScenario == 6) %% 4degree or all climate scenarios
    if exist ('input_data_4_0C_tmean.mat', 'file')
    disp ('Calculating residential overheating risk - 4DegreeC')
    
    cnt = 3; %to set correct columns for recording data in table for each build type.
    ppl_per_dwellings = zeros(nGridCells, 1); %average ppl per household - UDM uses 2.5ppl per dwelling. Calculated here for each gridcell based on total dwellings and population
    p=0; %counter to make sure the correct population cell array is read for years when warming set to occur by user. First run climate change only so p=0.
    
        for l = 1:socioEcScen % Compute for CC-only and then SE+CC (population and dwelling numbers)

            %pre-allocate arrays - cell arrays as with and without adaptation
            ResultsExceedDays = cell(adaptScen,1); %Intermediate outputs  - averaged across years and RCMS. Av, Median, 10th and 90th percentile.
            for n = 1:adaptScen
                ResultsExceedDays{n} = zeros(nGridCells, 14);
            end
            
            ResultsAnnualFrequency = cell(adaptScen,1); %Final outputs  - averaged across years and RCMS. Av, Median, 10th and 90th percentile.
            for n = 1:adaptScen
                ResultsAnnualFrequency{n} = zeros(nGridCells, 6);
            end
            
            ppl_per_event = cell(nGridCells, 1); %ppl affected after 5-conscutive days totalled across all build types affected, for all gridcells
            for n = 1:nGridCells
                ppl_per_event{n} = zeros(nRows, nCols);
            end

            ppl_per_event_adapt = cell(nGridCells, 1); %ppl affected after 5-conscutive days totalled across all build types affected, for all gridcells
            for n = 1:nGridCells
                ppl_per_event_adapt{n} = zeros(nRows, nCols);
            end

            for a = 1:adaptScen % run without and then with adapation (different threshold values)
                if a == 1 %no adaptation
                    adaptCnt = 3; %select correct thresholds from columns
                elseif a == 2 %adaptation
                    adaptCnt = 7; %select correct thresholds from columns
                end
                
                % pre-allocate cell arrays for intermediate values        
                events = cell(nGridCells, 1); %number of events, caulculated where x% population affected for x-number consecutive days.
                for n = 1:nGridCells
                    events{n} = zeros(nRows, nCols);
                end
                
                annualEvents = cell(nGridCells,1);   %store intermediate value of days per year TMax threshold exceeded.    
                for n = 1:nGridCells
                    annualEvents{n} = zeros(nRows, nCols);
                end           

                ppl_per_event_thres = cell(nGridCells, 1); %ppl affected after 5-conscutive days totalled across all build types affected, for all gridcells
                for n = 1:nGridCells
                    ppl_per_event_thres{n} = zeros(nRows, nCols);
                end

            %1. THRESHOLD DETECTOR - pulls out DAYS where thresholds are exceeded based on building_thresholds.mat from UCL model team

                for b = 1:nBuilds
                    %pre-allocated arrays --> reset for each building type.
                    days_exceeded = cell(nGridCells, 1); %days exceeded, for all gridcells by build types
                    for n = 1:nGridCells
                    days_exceeded{n} = zeros(nRows, nCols);
                    end
                    
                    annDaysExceed = cell(nGridCells,1);   %store intermediate value of days per year TMax threshold exceeded.    
                    for n = 1:nGridCells
                        annDaysExceed{n} = zeros(nRows, nCols);
                    end
                    
                    i=1;
                    while (h <= nGridCells)
                        if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) %looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                           
                            if b==1 % just calculate once. Calculates average ppl per dwellings across all dwelling types.
                                if population{1+p}(h,3)> 0 && sum(building_numbers{1+p}(h,3:6)) > 0 %avoid dividing by zeros.
                                    ppl_per_dwellings(h,1) = population{1+p}(h,3)/sum(building_numbers{1+p}(h,3:6));  %equally distributed across build types. UDM uses uniform value of 2.5ppl per dwelling.    
                                end
                            end
                            %Threshold detector
                            for r = 1:nRows
                                for c = 1:nCols
                                    if data_4_0C_tmean{i,2}(r,c) < building_thresholds{1+p}(h,b+adaptCnt) %no adaptation and adaptation thresholds
                                    % do nothing - speeds up if statement by specifying most likey option first..
                                    elseif data_4_0C_tmean{i,2}(r,c) >= building_thresholds{1+p}(h,b+adaptCnt)
                                    days_exceeded{h}(r,c) = 1;
                                        if a == 1 % no adaptation of housing stock
                                            ppl_per_event{h}(r,c) = ppl_per_event{h}(r,c)+(ppl_per_dwellings(h,1)*building_numbers{1+p}(h,b+2)); %cumulative ppl at risk per gridcell, increments across building types.
                                        elseif a == 2 % 100% adaptation of housing stock
                                            ppl_per_event_adapt{h}(r,c) = ppl_per_event_adapt{h}(r,c)+(ppl_per_dwellings(h,1)*building_numbers{1+p}(h,b+2)); %assumes 100% deployment across population and build type
                                        end
                                    end
                                end
                            end

                            h=h+1;
                            i=1; %reset counter
                            r=1; %reset counter
                            c=1; %reset counter

                        else
                            i=i+1;
                        end
                    end

                    %2. CALCULATE days per year thresholds exceeded for each building type (per grid cell)[Temp threshold output that can
                    %be plotted/mapped per buidling type (intermediate output).
                    for h = 1:nGridCells
                        annDaysExceed{h} = squeeze(nansum(reshape(days_exceeded{h}.',daysPerYear,size(days_exceeded{h},2)./daysPerYear,[]))).'; %re-shape array
                        for r = 1:nRows
                            avDaysExceed(r,h)= mean(annDaysExceed{h}(r,:),[1,2]); %mean across years for each GCM
                        end
                            percentile_low = prctile(avDaysExceed,(10),1);%10th percentile across years/GCMs
                            percentile_high = prctile(avDaysExceed,(90),1);%90th percentile across years/GCMs
                            av = mean(avDaysExceed,1);%average
                    end

                    %3. Save threshold results for plotting: gridded - by building type (b) across columns
                    if b == 1
                        ResultsExceedDays{a}(1:h,cnt:cnt+2) = [av.', percentile_low.', percentile_high.']; %3 variables[nGridCells*3]
                        ResultsExceedDays{a}(1:h,1) = lonIDIndex(:,1); %longitude
                        ResultsExceedDays{a}(1:h,2) = latIDIndex(:,1); %latitude
                        cnt = cnt+3; %increment columns to record results for each building type.

                    elseif b > 1 && b <=nBuilds
                        ResultsExceedDays{a}(1:h,cnt:cnt+2) = [av.', percentile_low.', percentile_high.']; %3 variables[nGridCells*3]
                        cnt = cnt+3;
                    end
                    h=1; % re-set counter for looping through next building type
                end

                % 4. When all building types have run and ppl_per_event calculated, identify where residents affected by discomfort meet overheating
                % event criteria based on user-defined variable (RDProbabilityLevel).
                if b == nBuilds && a == 1 % once all build types have run; No adaptation.
                    for h = 1:nGridCells
                        for r = 1:nRows
                            for c = 1:nCols
                                if ppl_per_event{h}(r,c) > 0 && ppl_per_event{h}(r,c)/population{1+p}(h,3) >= RDProbabilityLevel    % RDProbabilityLevel set by user.
                                    ppl_per_event_thres{h}(r,c) = 1;
                                end
                            end
                        end
                    end
                 elseif b == nBuilds && a == 2 % once all build types have run; With adaptation. Account for deployment levels here (from UCL - 100% of houses not adapated in grid cell, a proportion will still be affected under thresholds for no adaptation).
                    for h = 1:nGridCells
                        for r = 1:nRows
                            for c = 1:nCols
                                if ppl_per_event{h}(r,c) > 0 && ppl_per_event_adapt{h}(r,c) > 0 && ppl_per_event{h}(r,c) == ppl_per_event_adapt{h}(r,c) && ppl_per_event_adapt{h}(r,c)/population{1+p}(h,3) >= RDProbabilityLevel    % RDProbabilityLevel set by user.
                                    ppl_per_event_thres{h}(r,c) = 1; % adaptation has not mitigated ANY risk.
                                elseif ppl_per_event{h}(r,c) > 0 && ppl_per_event_adapt{h}(r,c) == 0 && (ppl_per_event{h}(r,c)*(1-deploymentPercent(deploymentYear,1)/100))/population{1+p}(h,3) >= RDProbabilityLevel % where adaptation has mitigated event across all building types
                                    ppl_per_event_thres{h}(r,c) = 1; % adaptation mitigated ALL risk but 100% of houses are not adapted - calculates residual risk.
                                elseif ppl_per_event{h}(r,c) > 0 && ppl_per_event_adapt{h}(r,c) > 0 && (ppl_per_event{h}(r,c)*(1-((ppl_per_event_adapt{h}(r,c)/ppl_per_event{h}(r,c))/(100*deploymentPercent(deploymentYear,1)))))/population{1+p}(h,3) >= RDProbabilityLevel
                                    ppl_per_event_thres{h}(r,c) = 1; % where adaptation has mitigated event across ONLY SOME building types BUT NOT ALL so proportional
                                end       
                            end
                        end
                    end 
                end
                 
                
                %5. Identify occurences where the number of consecutive days where ppl_per_event_thres = 1 >= RDNumberConsecutiveDays. This defines the overheating 'event'
                for h = 1:nGridCells
                    for r = 1:nRows
                        for c = 1:nCols
                             if ppl_per_event_thres{h}(r,c) == 0
                                % do nothing - speeds up if statement by specifying most likely option first..
                                consecutive_days = 0;
                             elseif ppl_per_event_thres{h}(r,c) == 1
                                consecutive_days = consecutive_days+1;
                                if consecutive_days == RDNumberConsecutiveDays
                                    events{h}(r,c) = 1;
                                    consecutive_days = 0;
                                end
                             end
                        end
                        consecutive_days = 0;
                    end
                end

                %6. CALCULATE annual probability of 'events' occuring (per grid cell).
                for h = 1:nGridCells  %[12*30]
                    annualEvents{h} = squeeze(nansum(reshape(events{h}.',daysPerYear,size(events{h},2)./daysPerYear,[]))).'; %re-shape array

                    for r = 1:nRows
                        AnnualFrequencyEvent(r,h)= mean(annualEvents{h}(r,:),[1,2]); %mean across years for each GCM
                    end
                        percentile_low = prctile(AnnualFrequencyEvent,(10),1);%10th percentile across years/GCMs
                        percentile_high = prctile(AnnualFrequencyEvent,(90),1);%90th percentile across years/GCMs
                        median = prctile(AnnualFrequencyEvent,(50),1);%50th percentile across years/GCMs
                        av = mean(AnnualFrequencyEvent,1);%average
                end

                %7. Save results for plotting: gridded
                ResultsAnnualFrequency{a}(:,1) = lonIDIndex(:,1); %longitude
                ResultsAnnualFrequency{a}(:,2) = latIDIndex(:,1); %latitude
                ResultsAnnualFrequency{a}(:,3:6) = [av.', median.', percentile_low.', percentile_high.'];
            
                cnt = 3; % re-set to calculate with adaptation.
                h=1; % re-set to calculate with adaptation.
            end      

     % Files to SAVE
            if l==1 %climate change only (with and without adaptation)
            fname = sprintf('RD_Results_4_0C_%s_%0.2f_%d_%0.2f.mat', '_CCOnly', RDProbabilityLevel, RDNumberConsecutiveDays, deploymentPercent(deploymentYear,1));
            ResultsAnnualFrequencyEvents_4_0_CC = ResultsAnnualFrequency;
            ResultsExceedDays_4_0_CC = ResultsExceedDays; % reflects exceedance of TMean and role of adaptation on changing thrsholds. SE change has no effect here.
            save(fname, 'ResultsAnnualFrequencyEvents_4_0_CC', 'ResultsExceedDays_4_0_CC')
            
            elseif l==2 %Climate change and population change (with and without adaptation)
            fname = sprintf('RD_Results_4_0C_%s_%s_%d_SSP%d_%0.2f__%d_%0.2f.mat', '_SE', 'Adapt', populationYear4degree, UK_SSP, RDProbabilityLevel, RDNumberConsecutiveDays, deploymentPercent(deploymentYear,1));
            ResultsAnnualFrequencyEvents_4_0_SE = ResultsAnnualFrequency;
            save(fname, 'ResultsAnnualFrequencyEvents_4_0_SE')
            end
            
            if populationYear4degree == 2020
                p=1;
            elseif populationYear4degree == 2030
                p=2;
            elseif populationYear4degree == 2050
                p=3;
            elseif populationYear4degree == 2080 || populationYear4degree == 2100 %2080, or 2100 can be used for end-century. TBC
                p=4;
            end
        end
        h=1; % re-set here if ClimateScenario == 6 and code continues...
    else
        disp ('no input data - 4.0C')
    end
end

disp ('End')