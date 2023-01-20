%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   HARM model (Heat - Adaptation - Risk - Model)
%
%   Estimates heat related mortality, residential discomfort and reduction in
%   labour productivity under future projections of climate change.
%
%   Updated from ARCADIA Impacts model developed in ARCADIA project (ECI, University of Oxford)
%   for OpenCLIM project (Tyndall Centre for Climate Change Research, UEA).
%
%   Uses post-processed UKCP18 data provided from HEAT model (Alan Kennedy-Asser, University of Bristol)
%
%   Author: Katie Jenkins    UEA
%   Date:12/01/2023 [version_05]
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT DATA
%   heat_impact_data.mat
%   1.  Mortality_Threshold_RR (nGridCells * 9)
%       Column 1 and 2: Lon and Lat of grid cell
%       Column 3: NUTS1/GOR name index 1-13 (inc. RoI) (see Region_names)
%       Column 4: Threshold TMean above which heat related mortality will
%       occur by GOR. Source: Vardoulakis et al., (2014) and Hajat et al., (2014)
%       Column 5-9: Exosure-Risk relationship (% change in mortality for
%       each 1degeeC increase in Mean Temp above the MortalityThreshold for
%       all ages and by age group (see Age_group). Source: Vardoulakis et al., 
%       (2014) and Hajat et al., (2014)).
%       
%   2.  Gridded population data {nClimScen}[nGridCells, 7]
%       Column 1&2 lat and long
%       column 3 {1} = Baseline 2011 census based gridded data (1km) 
%       regridded to 12km in ArcGIS. Source: Reis, S. et al.(2017) 
%       https://doi.org/10.5285/0995e94d-6d42-40c1-8ed4-5090d82471e1
%       column 4-7 {1} = gridded poulation per age group based on ONS
%       demographic projections from 2011 census for LAs.
%       Columns 4-7 {2-4} = For 2020, 2030, 2050 UK population projections and demographic data from
%       UK-SSPs SSP5. All ages 1km re-gridded to 12km. Per age group LAs gridded to 12km 
%       Source: https://www.ukclimateresilience.org/products-of-the-uk-ssps-project/
%       
%   3.  population_aggregate [5, 4] - for each climate scenario (columns)
%       gives the aggegrate UK population and split by age groups (rows).
%       Sources as above.
%
%   4.  Gridded dailyDeathRate calculated for 12km grid based on
%       mortality statistics from 2011 for England, Wales, Scotland and NI.
%       
%   5.  Region_names [13,1]: North East, North West, Yorkshire and Humber, East Midlands,
%       West Midlands, East England, London, South East, South West, Wales,
%       Scotland, Northern Ireland, Republic of Ireland.
%
%   6.  Age_group [5,1]: description of the age group classes used.
%      ['All';'0-64';'65-74';'75-84';'85+']
%
%   7. 12km gridded TMean daily data {nGridCells}(12, 14800). .csv files of 
%       daily TMean (one per grid cell) from HEAT model for past (1990-2019); and warming levels 1.5°C
%       2.0°C, 3.0°C and 4.0°C
%
%       These are read in, saved and used for mortality calculations. 
%       Grid coordinates also extracted from file names to provide grid ID for 
%       mapping based on lat_UK_RCM.mat and long_UK_RCM.mat.
%
%   8/9. lat_UK_RCM.mat and long_UK_RCM.mat: [82,112]. From Heat outputs, these are read 
%      in and provide latitude and longitude based on coordinate given in file 
%      names from 7 above. Saved in gridID.
%
%   10. mortality_acclimatisation: {climatescenarios, 1}(nGridCells, 6) Adaptation/acclimatisation values for each warming level.
%       no adaptation; adapt threshold by 1degreeC; adapt threshold by 2 degreeC; 
%       adapt threshold by 95th percentile for 1.5, 2, 3, 4 degree (zero for
%       near past)warming levels.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CURRENT MODEL OUTPUTS:
%   1. Plots average annual heat related deaths (10th & 90th P reflect climate
%   model uncertainty), for CC (Climate change Only Scenario) and CC+SE
%   (Climate Change and future projection of population both accounted for).
%   Scenarios reflect the past, 1.5, 2.0 and 3 degree warming levels.
%   2. Heat related deaths per year per 100,000 population (split
%   by age groups)for CC (Climate change Only Scenario) and CC+SE
%   (Climate Change and future projection of population).
%   3. Average annual heat related deaths per increment of warming above
%   thresholds.
%   4. Plots average annual heat related deaths (10th & 90th P reflect climate
%   model uncertainty), for CC (Climate change Only Scenario) and CC+SE
%   (Climate Change and future projection of population both accounted for).
%   assuming adaptation (acclimatisation adjusted thresholds).
%   5. Daily and annual timeseries as intermediate outputs if needed, split by GCM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('*******************************************************************')
disp('Heat Adaptation Risk Model (HARM)')
disp('version 04')
disp('This version currently estimates UK heat-related mortality')
disp('and results for labour productivity.')
disp('Version for use on DAFNI, June 2022')
disp('*******************************************************************')

tic %stopwatch to measure performance
format short e
dbstop if error %debug if error

load('heat_impact_data') %contains all the 15.m input data required.

%Parameters
socioEcScen = 2; % Socioeconomic scenario options: (Climate change only, Climate change and socio-economic change)
adaptScen = 2; %Used to calculate results WITH NO additional adaptation and WITH additonal adaptation.
nClimScen = 0; %number of climate scenarios - based on input data.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User defined parameters - CHANGE HERE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BiasCorrected = 1; %TRUE = 1, FALSE = 0 %%USER DEFINED AT THIS STAGE IN FUTURE ALWAYS USE BIAS CORRECTED??
impactMetric = 1; %Select impact to run - 1 = mortality; 2 = labour productivity; 3 = residential discomfort.
ClimateScenario = 6; %1 = past; 2 = 1.5 degree; 3 = 2.0 degree; 4 = 3.0 degree; 5 = 4.0 degree; 6 = run ALL.
populationYearbaseline = 2011;
populationYear15degree = 2020; % 2020, 2050, 2080; population year defines when warming level projected to occur. NOTE for residential dicomfort as data from UDM only use 2050 and 2080.
populationYear2degree = 2030; % 2030; 2050; 2080; population year defines when warming level projected to occur. NOTE for residential dicomfort as data from UDM only use 2050 and 2080.
populationYear3degree = 2050; % 2050; 2080; population year defines when warming level projected to occur. NOTE for residential dicomfort as data from UDM only use 2050 and 2080.
populationYear4degree = 2100; % NOTE for residential dicomfort as data from UDM only use 2050 and 2080.
UK_SSP = 5; %data available from UDM for UK-SSPs 2,4 and 5.

%Mortality Related
acclimaScen = 1; %Select approach to behavioural adaptation via natural acclimatisation. 1 = No Adaptation; 2 & 3 = pre-defined thresholds 1 and 2 degrees respectively; 4 = 93rd P of TMean (spatially and temporally explicit based on output from HEAT)

% Residential discomfort related
RDProbabilityLevel = 0.25; % Level at which to define probability of residential thermal discomfort. [See Kingsborough et al: http://dx.doi.org/10.1016/j.crm.2017.01.001. Set values as 0.1, 0.25 and 0.5.
RDNumberConsecutiveDays = 5; % number of consecutive days of overheating above which overheating is calculated.
RDAdaptScenario = 2; % 1 = standard whole house retrofit; 2 = option 1 plus external shutters for shading; 3 = air conditioning uptake only.

%labour productivity related
acclimatised = 1; % 0 = non-acclimatised; 1= acclimatised. Defines which Exposure response Functions (ERFs) are used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Add something to give files more readable names...Mortality
if acclimaScen == 1
    adaptationName = 'No Adaptation';
elseif acclimaScen == 2
    adaptationName = 'Adapt threshold 1degree';
elseif acclimaScen == 3
    adaptationName = 'Adapt threshold 2degree';
elseif acclimaScen == 4
    adaptationName = 'Adapt 93rd Percentile';
end

% Add something to give files more readable names...Residential discomfort
if RDAdaptScenario == 1
    adaptationIntervention = 'Retro';
elseif RDAdaptScenario == 2
    adaptationIntervention = 'Retro_shading';
elseif RDAdaptScenario == 3
    adaptationIntervention = 'a/c';
end

%% Select correct SSP dataset. OFFLINE AT
% MOMENT - IN FUTURE PULL FROM .asc FILES ON DAFNI first then aggregate.
population_aggregate = zeros(5,5);

if UK_SSP == 1
    population = population_SSP1;
    for i=1:5
        population_aggregate(:,i) = sum(population{i}(:,3:7),[]); 
    end
elseif UK_SSP == 2
    population = population_SSP2;
    for i=1:5
        population_aggregate(:,i) = sum(population{i}(:,3:7),[]); 
    end
elseif UK_SSP == 3
    population = population_SSP3;
    for i=1:5
        population_aggregate(:,i) = sum(population{i}(:,3:7),[]); 
    end
elseif UK_SSP == 4
    population = population_SSP4;
    for i=1:5
        population_aggregate(:,i) = sum(population{i}(:,3:7),[]); 
    end
elseif UK_SSP == 5
    population = population_SSP5;
    for i=1:5
        population_aggregate(:,i) = sum(population{i}(:,3:7),[]); 
    end
end

%%
%UKCP18 has 360 days per year (30 days per month)- set to 360 for raw data below %Bias corrected converted to 365 days per year (as leap years removed below)
if BiasCorrected == 0
        daysPerYear = 360;
elseif BiasCorrected == 1
        daysPerYear = 365;
end

if impactMetric == 1
   disp ('Metric: UK heat related mortality')
    
%% 1. Search for available files; Read temperature data from HEAT .csv files; set up arrays; read coordinates
%  Read each gridcell array to a cell Array inc. the cell ID.
%  Each TMean array represents daily data e.g. [1-360] per year [1-30][nCols=10,800] * 12 GCMS [nRows = 12]
%  Only need to do once - if files exist then skip

% RUN EITHER BIAS CORRECTED OR RAW DEPENDENT ON USER DEFINED PARAMETER %
% To account for different number of days and leap years in files

if BiasCorrected == 0 %False, reads raw data which has 30 day months, no leap years
   disp ('Using non-bias corrected data')
 
        %  Baseline / near past TMean
        %  If raw data already read in once and saved then skip initial step
    if (ClimateScenario == 1 || ClimateScenario == 6) %% past or all climate scenarios
        if exist ('input_data_past.mat', 'file') %%can use to speed up
        %if running offline from DAFNI - checks if file already exists.
           load ('input_data_past.mat')
           disp ('loading input data Past TMean')
           nClimScen = nClimScen+1; %counts # climate scenarios being run together

        elseif exist('ARCADIA-Tmean-raw-past', 'file')
            disp ('Read in daily TMean data (Past) for mortality estimates')

            dd = dir('ARCADIA-Tmean-raw-past\*.csv'); % all .csv files in folder
            fileNames = {dd.name}; 
            data_past = cell(numel(fileNames),2);
            data_past(:,1) = regexprep(fileNames, '_UKCP18-12km_T.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)   
                fname = fullfile('ARCADIA-Tmean-raw-past', fileNames(ii)); %allocate each .csv data file to cell array
                data_past{ii,2} = dlmread(fname{1},',',0,0);
            end

            save('input_data_past.mat', 'data_past', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together

        else
            disp ('No Past files exist')
        end
    end

    if (ClimateScenario == 2 || ClimateScenario == 6) %% 1.5degree or all climate scenarios
        %  1.5degreeC
%         if exist ('input_data_1.5C.mat', 'file')
%            load ('input_data_1.5C.mat')
%            disp ('loading input data 1.5C')
%            nClimScen = nClimScen+1; %counts # climate scenarios being run together
% 
%         elseif exist
        if exist('ARCADIA_Tmean_raw_1_5C', 'file')
            disp ('Read in daily TMean data (1.5DegreeC) for mortality estimates')

            dd = dir('ARCADIA_Tmean_raw_1_5C\*.csv'); % all .csv files in folder
            fileNames = {dd.name}; 
            data_1_5C = cell(numel(fileNames),2);
            data_1_5C(:,1) = regexprep(fileNames, '_UKCP18-12km_T.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)
                fname = fullfile('ARCADIA_Tmean_raw_1_5C', fileNames(ii)); %allocate each .csv data file to cell array
                data_1_5C{ii,2} = dlmread(fname{1},',',0,0);
            end

            save('input_data_1.5C.mat', 'data_1_5C', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together

        else
            disp ('No 1.5degreeC files exist')
        end
    end
        %  2degreeC
    if (ClimateScenario == 3 || ClimateScenario == 6) %% 2 degree or all climate scenarios
%         if exist ('input_data_2.0C.mat', 'file')
%            load ('input_data_2.0C.mat')
%            disp ('loading input data 2.0C')
%            nClimScen = nClimScen+1; %counts # climate scenarios being run together
% 
%         elseif
        if exist('ARCADIA_Tmean_raw_2_0C', 'file')
            disp ('Read in daily TMean data (2.0DegreeC) for mortality estimates')

            dd = dir('ARCADIA_Tmean_raw_2_0C\*.csv'); % all .csv files in folder
            fileNames = {dd.name}; 
            data_2_0C = cell(numel(fileNames),2);
            data_2_0C(:,1) = regexprep(fileNames, '_UKCP18-12km_T.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)   
                fname = fullfile('ARCADIA_Tmean_raw_2_0C', fileNames(ii)); %allocate each .csv data file to cell array
                data_2_0C{ii,2} = dlmread(fname{1},',',0,0);
            end

            save('input_data_2.0C.mat', 'data_2_0C', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together

        else
            disp ('No 2.0DegreeC files exist')
        end
    end

        %  3degreeC 
     if (ClimateScenario == 4 || ClimateScenario == 6) %% 3degreeC or all climate scenarios
%         if exist ('input_data_3.0C.mat', 'file')
%            load ('input_data_3.0C.mat')
%            disp ('loading input data 3.0C')
%            nClimScen = nClimScen+1; %counts # climate scenarios being run together
% 
%         elseif exist
        if exist('ARCADIA_Tmean_raw_3_0C', 'file')
            disp ('Read in daily TMean data (3.0DegreeC) for mortality estimates')

            dd = dir('ARCADIA_Tmeanraw_3_0C\*.csv'); % all .csv files in folder
            fileNames = {dd.name}; 
            data_3_0C = cell(numel(fileNames),2);
            data_3_0C(:,1) = regexprep(fileNames, '_UKCP18-12km_T.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)   
                fname = fullfile('ARCADIA_Tmean_raw_3_0C', fileNames(ii)); %allocate each .csv data file to cell array
                data_3_0C{ii,2} = dlmread(fname{1},',',0,0);
            end

            save('input_data_3.0C.mat', 'data_3_0C', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together

        else
            disp ('No 3.0DegreeC files exist')
        end
     end
         %  3degreeC 
     if (ClimateScenario == 5 || ClimateScenario == 6) %% 3degreeC or all climate scenarios
%         if exist ('input_data_4.0C.mat', 'file')
%            load ('input_data_4.0C.mat')
%            disp ('loading input data 4.0C')
%            nClimScen = nClimScen+1; %counts # climate scenarios being run together
% 
%         elseif exist
        if exist('ARCADIA_Tmean_raw_4_0C', 'file')
            disp ('Read in daily TMean data (4.0DegreeC) for mortality estimates')

            dd = dir('ARCADIA_Tmean_raw_4_0C\*.csv'); % all .csv files in folder
            fileNames = {dd.name}; 
            data_4_0C = cell(numel(fileNames),2);
            data_4_0C(:,1) = regexprep(fileNames, '_UKCP18-12km_T.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)   
                fname = fullfile('ARCADIA_Tmean_raw_4_0C', fileNames(ii)); %allocate each .csv data file to cell array
                data_4_0C{ii,2} = dlmread(fname{1},',',0,0);
            end

            save('input_data_4.0C.mat', 'data_4_0C', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together

        else
            disp ('No 4.0DegreeC files exist')
        end
     end
 
elseif BiasCorrected == 1 %reads bias corrected data and removes leap years.
       disp ('Using bias corrected data')
       
    %  Past / baseline
    if (ClimateScenario == 1 || ClimateScenario == 6) %% baseline or all climate scenarios
        if exist ('input_data_past.mat', 'file') %% speed up if running
%         offline from DAFNI for same set of files.
           load ('input_data_past.mat')
           disp ('loading input data Past (Bias corrected)')
           nClimScen = nClimScen+1; %counts # climate scenarios being run together        
        
    elseif exist('ARCADIA-Tmean-BC-past', 'file')
            disp ('Read in daily TMean data (Past) for mortality estimates')

            dd = dir('ARCADIA-Tmean-BC-past\*.csv'); % all .csv files in folder
            fileNames = {dd.name}; 
            data_past = cell(numel(fileNames),2);
            data_past(:,1) = regexprep(fileNames, '_UKCP18-12km_T.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)   
                fname = fullfile('ARCADIA-Tmean-BC-past', fileNames(ii)); %allocate each .csv data file to cell array
                data_past{ii,2} = dlmread(fname{1},',',0,0);
                data_past{ii,2}(:,1461:1461:end) = []; %delete leap year, every 1461 columns. Bias corrected data uses Gregorian Calendar.
                data_past{ii,2}(:,10950) = 0; %Add back in dummy value for 31st December so years are even.
            end

            save('input_data_past.mat', 'data_past', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together

        else
            disp ('No Past (bias corrected) files exist')
        end
    end
    
    if (ClimateScenario == 2 || ClimateScenario == 6) %% 1.5 degree or all climate scenarios
        %  1.5degreeC 
        if exist ('input_data_1.5C.mat', 'file')
           load ('input_data_1.5C.mat')
           disp ('loading input data 1.5C (Bias corrected)')
           nClimScen = nClimScen+1; %counts # climate scenarios being run together

         elseif exist('ARCADIA-Tmean-BC-15', 'file')
            disp ('Read in daily TMean data (1.5DegreeC) for mortality estimates')

            dd = dir('ARCADIA-Tmean-BC-15\*.csv'); % all .csv files in folder
            fileNames = {dd.name}; 
            data_1_5C = cell(numel(fileNames),2);
            data_1_5C(:,1) = regexprep(fileNames, '_UKCP18-12km_T.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)
                fname = fullfile('ARCADIA-Tmean-BC-15', fileNames(ii)); %allocate each .csv data file to cell array
                data_1_5C{ii,2} = dlmread(fname{1},',',0,0);
                data_1_5C{ii,2}(:,1461:1461:end) = []; %delete leap year, every 1461 columns. Bias corrected data uses Gregorian Calendar.
                data_1_5C{ii,2}(:,10950) = 0; %Add back in dummy value for 31st December so years are even.
            end

            save('input_data_1.5C.mat', 'data_1_5C', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together

        else
            disp ('No 1.5degreeC files exist')
        end
    end
    
    %  2degreeC
    if (ClimateScenario == 3 || ClimateScenario == 6) %% 2 degree or all climate scenarios
         if exist ('input_data_2.0C.mat', 'file')
           load ('input_data_2.0C.mat')
           disp ('loading input data 2.0C (Bias corrected)')
           nClimScen = nClimScen+1; %counts # climate scenarios being run together

         elseif exist('ARCADIA-Tmean-BC-20', 'file')
            disp ('Read in daily TMean data (2.0DegreeC) for mortality estimates')

            dd = dir('ARCADIA-Tmean-BC-20\*.csv'); % all .csv files in folder
            fileNames = {dd.name}; 
            data_2_0C = cell(numel(fileNames),2);
            data_2_0C(:,1) = regexprep(fileNames, '_UKCP18-12km_T.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)   
                fname = fullfile('ARCADIA-Tmean-BC-20', fileNames(ii)); %allocate each .csv data file to cell array
                data_2_0C{ii,2} = dlmread(fname{1},',',0,0);
                data_2_0C{ii,2}(:,1461:1461:end) = []; %delete leap year, every 1461 columns. Bias corrected data uses Gregorian Calendar.
                data_2_0C{ii,2}(:,10950) = 0; %Add back in dummy value for 31st December so years are even.
            end

            save('input_data_2.0C.mat', 'data_2_0C', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together

        else
            disp ('No 2.0DegreeC files exist')
        end
    end
    
        %  3degreeC
    if (ClimateScenario == 4 || ClimateScenario == 6) %% 3 degree or all climate scenarios
        if exist ('input_data_3.0C.mat', 'file')
           load ('input_data_3.0C.mat')
           disp ('loading input data 3.0C (Bias corrected)')
           nClimScen = nClimScen+1; %counts # climate scenarios being run together
 
         elseif exist('ARCADIA-Tmean-BC-30', 'file')
            disp ('Read in daily TMean data (3.0DegreeC) for mortality estimates')

            dd = dir('ARCADIA-Tmean-BC-30\*.csv'); % all .csv files in folder
            fileNames = {dd.name}; 
            data_3_0C = cell(numel(fileNames),2);
            data_3_0C(:,1) = regexprep(fileNames, '_UKCP18-12km_T.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)   
                fname = fullfile('ARCADIA-Tmean-BC-30', fileNames(ii)); %allocate each .csv data file to cell array
                data_3_0C{ii,2} = dlmread(fname{1},',',0,0);
                data_3_0C{ii,2}(:,1461:1461:end) = []; %delete leap year, every 1461 columns. Bias corrected data uses Gregorian Calendar.
                data_3_0C{ii,2}(:,10950) = 0; %Add back in dummy value for 31st December so years are even.
            end

            save('input_data_3.0C.mat', 'data_3_0C', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together

        else
            disp ('No 3.0DegreeC files exist')
        end
    end
    
    % 4 degree
    if (ClimateScenario == 5 || ClimateScenario == 6) %% 4 degree or all climate scenarios
        if exist ('input_data_4.0C.mat', 'file')  %%can use to speed up
        %if running offline from DAFNI - checks if file already exists.
           load ('input_data_4.0C.mat')
           disp ('loading input data 4.0C (Bias corrected)')
           nClimScen = nClimScen+1; %counts # climate scenarios being run together

        elseif exist('ARCADIA-Tmean-BC-40', 'file')
            disp ('Read in daily TMean data (4.0DegreeC) for mortality estimates')

             dd = dir('ARCADIA-Tmean-BC-40\*.csv'); % all .csv files in folder
            fileNames = {dd.name}; 
            data_4_0C = cell(numel(fileNames),2);
            data_4_0C(:,1) = regexprep(fileNames, '_UKCP18-12km_T.csv',''); %set name as grid_cell ID (lon_lat)

            for ii = 1:numel(fileNames)   
                fname = fullfile('ARCADIA-Tmean-BC-40', fileNames(ii)); %allocate each .csv data file to cell array
                data_4_0C{ii,2} = dlmread(fname{1},',',0,0);
                data_4_0C{ii,2}(:,1461:1461:end) = []; %delete leap year, every 1461 columns. Bias corrected data uses Gregorian Calendar.
                data_4_0C{ii,2}(:,10950) = 0; %Add back in dummy value for 31st December so years are even.
            end

            save('input_data_4.0C.mat', 'data_4_0C', '-v7.3') %~1GB (cell_IDs and data array)
            nClimScen = nClimScen+1; %counts # climate scenarios being run together

        else
            disp ('No 4.0DegreeC files exist')
        end
    end
end

%% Set lat and long for use in mapping outputs and linking regional mortality thresholds and Relative Risk (RR) to gridded TMean
%(format DS --> currently convert to DMS in ArcGIS for plotting)
% Only need to do once, if file exists skip step.

% if exist ('gridID.mat', 'file')
%     load ('gridID')  
%     disp('loading gridID data')
%     [nGridCells,nVar] = size(gridID);
% else
    disp('calculating gridID data')
    load lat_UK_RCM
    load long_UK_RCM
    

    if(ClimateScenario == 1 || ClimateScenario == 6)
        [nGridCells,nVar] = size(data_past); 
        HEAT_Data = data_past;
    elseif (ClimateScenario == 2)
        [nGridCells,nVar] = size(data_1_5C); 
        HEAT_Data = data_1_5C;
    elseif (ClimateScenario == 3)
        [nGridCells,nVar] = size(data_2_0C); 
        HEAT_Data = data_2_0C;
    elseif (ClimateScenario == 4)
        [nGridCells,nVar] = size(data_3_0C); 
        HEAT_Data = data_3_0C;
    elseif (ClimateScenario == 5)
        [nGridCells,nVar] = size(data_4_0C); 
        HEAT_Data = data_4_0C;
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
   
    lonIDIndex = fix(Mortality_Threshold_RR(:,1)*1e5)/1e5; %truncate decimals to allow match to gridID - different order to IDs in gridID
    latIDIndex = fix(Mortality_Threshold_RR(:,2)*1e5)/1e5; %truncate decimals to allow match to gridID - different order to IDs in gridID
    
    save('gridID', 'gridID', 'lonIDIndex', 'latIDIndex')
    clear cellIndex %not used past this point - use gridID
%end

%[reorderedGridID, ia, ib] = intersect(lonIDIndex (:,1), gridID(:,1),'stable'); %index created to reorder the data from HEAT so rows align with coordinates of the socio-eocnomic data in calculations below.

%% Define size of variables / Cell Arrays for calculations
    [nRows,nCols] = size(HEAT_Data{1,2}); %set up TMean data array size
    years = nCols/daysPerYear; %check years correspond to 30
    ageGroupCnt = size(Age_group,1); %counter to loop through mortality calculations for each age group
    
%% CALCULATE MORTALITY - Average Annual, per 100,000ppl all Ages, and split by age group
%loop through gridded TMean data to calculate daily mortality risk with
%constant and future population for all nClimScens.
   
avMortality = zeros(nRows,nGridCells);
avMortalitySE = zeros(nRows,nGridCells); %SE = Socioeconomic scenario i.e. population and demographic change
avMortPerIncrementUK = zeros(nClimScen*2,30); % results for CC only and CC+SE for plotting
header  = 1:30; % displays the increment exceeded in 1 degree intervals above threshold
avMortPerIncrementUK(1,:) = header; % degrees above threshold in 1 degree intervals
threshold = zeros(nGridCells,1); % re-calculates threshold based on natural acclimatisation scenario.
h=1; %h used to count through nGridCells.

%% PAST
if (ClimateScenario == 1 || ClimateScenario == 6) %% baseline or all climate scenarios
    if exist ('input_data_past.mat', 'file')  %calculate for each time period - based on input files
       disp ('Calculating mortality risk - Near Past')

        %Preallocate cell arrays
        ResultsPast = cell(adaptScen,1); %append all final gridded results for plotting and saving, for each age group, one array for each adaptation scenario
        for n = 1:adaptScen
            ResultsPast{n} = zeros(nGridCells, 7);
        end

        ResultsPastAgeGroup = cell(adaptScen,1); %append all final results aggregated per age group, for plotting and saving, one array for each adaptation scenario
        for n = 1:adaptScen
            ResultsPastAgeGroup{n} = zeros(ageGroupCnt, 5); 
        end
        
        annualMortalityUK_CIs = cell(adaptScen,1);
        for n = 1:adaptScen
            annualMortalityUK_CIs{n} = zeros(1, 2); 
        end

        adapt = 1; %first run assumes no adaptation

        for a = 1:adaptScen

            dailyMort = cell(nGridCells,1);   %store intermediate value of daily additional deaths          
            for n = 1:nGridCells
                dailyMort{n} = zeros(nRows, nCols);
            end

            totMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (summed across all years/GCMs)              
            for n = 1:nGridCells
                totMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
            end

            avMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (average annual)    
            for n = 1:nGridCells
                avMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
            end

            annMortality = cell(nGridCells,1);   %store intermediate value of Mortality risk * population (average annual)      
            for n = 1:nGridCells
                annMortality{n} = zeros(nRows, nCols);
            end

            threshold(:,1) = Mortality_Threshold_RR(:,4) + mortality_acclimatisation{1,1}(:,adapt+2);

            for m = 1:ageGroupCnt %for all ages and per age group
                i=1;
                while (h <= nGridCells)
                    if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) %looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                        for r = 1:nRows
                            for c = 1:nCols
                                if data_past{i,2}(r,c) < threshold(h,1)
                                    % do nothing - speeds up if statement by specifying most likey option first..
                                elseif Mortality_Threshold_RR(h,4) > 0 
                                   dailyMort{h}(r,c) = ((((data_past{i,2}(r,c)- threshold(h,1))*Mortality_Threshold_RR(h,m+4))/100)*dailyDeathRate(h,m+2))*population{1}(h,m+2); % Outputs the additional number of daily deaths %%increments 2 as population + dailyDeathRate includes lat and lon columns first. 
                                   incrementExceeded = round(data_past{i,2}(r,c)- threshold(h,1));
                                   totMortPerIncrement{h}(r, incrementExceeded+1) = totMortPerIncrement{h}(r, incrementExceeded+1) + dailyMort{h}(r,c); %cumulative mortality per each 1 degree increment above the threshold
                                else % do nothing - zero value.
                                end
                            end
                        end
                        h=h+1; %increment through data in Mortality_Threshold_RR
                        i=1; %reset counter
                        r=1; %reset counter
                        c=1; %reset counter

                    else
                            i=i+1;
                    end
                end

    %% average annual mortality per increment %%For plots - UK, All ages (m=1) only, no adaptation (a=1)
               if (m==1) && (a==1)
                    for h = 1:nGridCells
                       for r = 1:nRows
                            avMortPerIncrement{h}(r,:) = totMortPerIncrement{h}(r,:)/years;
                       end
                    end
                       dim = ndims(avMortPerIncrement{1});          % Get the number of dimensions for arrays
                       matrix = cat(dim+1,avMortPerIncrement{:});        % Convert to a (dim+1)-dimensional matrix
                       sumArray = sum(matrix,dim+1);  % Sum mortality across UK grid cells
                       avMortPerIncrementUK(2,:) = mean(sumArray); % average across 12 RCMs
                       
                       %calculating CI across the 12 RCMS per increment
                       CI_increment = zeros(30,2);
                        for q = 1:30
                        stDev = std(sumArray(:,q)); 
                        x = mean(sumArray(:,q));
                        sqRoot = sqrt(12);
                        SEM = stDev/sqRoot;% 
                        ts = tinv([0.025  0.975],12-1);                 % T-Score
                        CI_increment(q,:) = x + ts*SEM;                      %absolute CI
                        %yCI95 = bsxfun(@times, SEM, ts(:));   %value above and below
                        end
                       
                end

     %% Calculate Annual risk (across years and across 12 GCMs) and percentile range
                    annMortalityUK = zeros(nRows, years);                
                   for h = 1:nGridCells
                        annMortality{h} = squeeze(nansum(reshape(dailyMort{h}.',daysPerYear,size(dailyMort{h},2)./daysPerYear,[]))).';
                        annMortalityUK(:,:) = annMortalityUK(:,:) + annMortality{h}(:,:);
                        for r = 1:nRows
                            avMortality(r,h) = mean(annMortality{h}(r,:),[1,2]); %mean across years for each GCM
                        end
                        p_low = prctile(avMortality,(10),1);%10th and 90th percentile across years/GCMs
                        p_high = prctile(avMortality,(90),1);
                        av_ann_rcms = mean(annMortalityUK,2);
                        av = mean(avMortality,1);
                        sum_av = sum(av);
                        sum_p_low = sum(p_low);
                        sum_p_high = sum(p_high);
                    end
                   
                        Column_vector = annMortalityUK(:); %Convert array to column vector                    
                        %calculating CI across the 12 RCMS per day not across
                        %the 2205 gridcells
                        % confidence interval
                        stDev = std(av_ann_rcms); 
                        x = mean(av_ann_rcms);
                        sqRoot = sqrt(12);
                        SEM = stDev/sqRoot;% 
                        ts = tinv([0.025  0.975],12-1);                 % T-Score
                        CI = x + ts*SEM;                      %absolute CI
                        yCI95 = bsxfun(@times, SEM, ts(:));   %value above and below
                        
            
                % Values per 100,000 ppl
                    av_per_population = sum_av/population_aggregate(m,1)*100000; %mean deaths per year/per 100,000 population - as a proportion of age group related population
                    p_low_per_population = sum_p_low/population_aggregate(m,1)*100000; %10th P deaths per year/per 1000 population 
                    p_high_per_population = sum_p_high/population_aggregate(m,1)*100000; %90th P deaths per year/per 1000 population 

                    p_low_per_population_CI = CI(1,1)/population_aggregate(m,1)*100000; %10th P deaths per year/per 1000 population 
                    p_high_per_population_CI = CI(1,2)/population_aggregate(m,1)*100000; %90th P deaths per year/per 1000 population 

                    
                %Save results for plotting:gridded - by age group category (m)
                if m == 1 % all ages
                    ResultsPast{a}(1:h,3:5) = [av.', p_low.', p_high.']; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                    ResultsPast{a}(1:h,1) = lonIDIndex(:,1); %longitude
                    ResultsPast{a}(1:h,2) = latIDIndex(:,1); %latitude
                    annualMortalityUK_CIs{a} = CI;
                %Save results for plotting - UK aggregate
                    ResultsPastAgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                elseif m==2
                    %ResultsPast{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                    ResultsPastAgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                elseif m==3
                    %ResultsPast{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                    ResultsPastAgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                elseif m==4
                    %ResultsPast{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                    ResultsPastAgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                elseif m==5
                    %ResultsPast{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                    ResultsPastAgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                end
                h=1; %set back to 1
            end
            adapt = acclimaScen; %now runs with adaptation.
            
        end 

        % Files to SAVE
        fname = sprintf('ResultsPast_%d.mat', populationYearbaseline);
        save(fname, 'ResultsPast', 'ResultsPastAgeGroup', 'annualMortalityUK_CIs','CI_increment')
        clear('dailyMort','totMortPerIncrement', 'avMortPerIncrement','annMortality', 'p_low', 'p_high', 'av', 'sum_av', 'sum_p_low', 'sum_p_high', 'av_per_population', 'p_low_per_population', 'p_high_per_population')  %daily data not used past this point - tidy up.
    else
        disp ('no input data - past')
    end
    save ('avMortPerIncrementUK.mat', 'avMortPerIncrementUK') % Used to plot deaths per increment for paper.
end

%% Below calculates results for 1.5 degree for both climate change and climate change and SE change (2020 population)scenarios
if (ClimateScenario == 2|| ClimateScenario == 6) %% 1.5degree or all climate scenarios
    if exist ('input_data_1.5C.mat', 'file')
    disp ('Calculating mortality risk - 1.5DegreeC')
       
    p=0; %counter to make sure the correct population cell array is read for years when warming set to occur by user. First run climate change only so p=0.
        for l = 1:socioEcScen % Compute for CC-only and SE+CC

        %set up cell array
        Results_1_5C = cell(adaptScen,1); %append all final results for plotting and saving, for each age group
        for n = 1:adaptScen
            Results_1_5C{n} = zeros(nGridCells, 7); %set automatically in future based on gridcell rows
        end

        Results_1_5C_AgeGroup = cell(adaptScen,1); %append all final results for plotting and saving, for each age group and new array for each adaptation scenario
        for n = 1:adaptScen
            Results_1_5C_AgeGroup{n} = zeros(ageGroupCnt, 5);
        end
        
        annualMortalityUK_CIs_1_5 = cell(adaptScen,1);
        for n = 1:adaptScen
            annualMortalityUK_CIs_1_5{n} = zeros(1, 2); 
        end
        
        adapt = 1; %first run assumes no adaptation

            for a = 1:adaptScen

                dailyMort = cell(nGridCells,1);   %store intermediate value of daily additional deaths          
                for n = 1:nGridCells
                    dailyMort{n} = zeros(nRows, nCols);
                end

                totMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (summed across all years/GCMsCC-only              
                for n = 1:nGridCells
                    totMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                end

                avMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (average annual)   
                for n = 1:nGridCells
                    avMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                end

                annMortality = cell(nGridCells,1);   %store intermediate value of Mortality * population (average annual)        
                for n = 1:nGridCells
                    annMortality{n} = zeros(nRows, nCols);
                end

                threshold(:,1) = Mortality_Threshold_RR(:,4) + mortality_acclimatisation{2,1}(:,adapt+2);  

                for m = 1:ageGroupCnt %for all ages and per age group
                i=1;
                while (h <= nGridCells)
                    if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) %looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                        for r = 1:nRows
                            for c = 1:nCols
                                if data_1_5C{i,2}(r,c) < threshold(h,1)
                                    % do nothing - speeds up if statement by specifying most likey option first..
                                elseif Mortality_Threshold_RR(h,4) > 0 
                                   dailyMort{h}(r,c) = ((((data_1_5C{i,2}(r,c)- threshold(h,1))*Mortality_Threshold_RR(h,m+4))/100)*dailyDeathRate(h,m+2))*population{1+p}(h,m+2); % Outputs the additional number of daily deaths %%increments 2 as population + dailyDeathRate includes lat and lon columns first. 
                                   incrementExceeded = round(data_1_5C{i,2}(r,c)- threshold(h,1));
                                   totMortPerIncrement{h}(r, incrementExceeded+1) = totMortPerIncrement{h}(r, incrementExceeded+1) + dailyMort{h}(r,c); %cumulative mortality per each 1 degree increment above the threshold
                                else % do nothing - zero value.
                                end
                            end
                        end
                        h=h+1; %increment through data in Mortality_Threshold_RR
                        i=1; %reset counter
                        r=1; %reset counter
                        c=1; %reset counter

                    else
                            i=i+1;
                    end
                end          
                 
             %% average annual mortality per increment %%for all ages only and no adaptation     
                   if (m==1) && (a==1)           
                        for h = 1:nGridCells
                           for r = 1:nRows
                                avMortPerIncrement{h}(r,:) = totMortPerIncrement{h}(r,:)/years;
                           end
                        end
                           dim = ndims(avMortPerIncrement{1});          % Get the number of dimensions for arrays
                           matrix = cat(dim+1,avMortPerIncrement{:});        % Convert to a (dim+1)-dimensional matrix
                           sumArray = sum(matrix,dim+1);  % sum mortality across UK grid cells
                           avMortPerIncrementUK(l+2,:) = mean(sumArray);
                           
                       %calculating CI across the 12 RCMS per increment
                       CI_increment = zeros(30,2);
                        for q = 1:30
                        stDev = std(sumArray(:,q)); 
                        x = mean(sumArray(:,q));
                        sqRoot = sqrt(12);
                        SEM = stDev/sqRoot;% 
                        ts = tinv([0.025  0.975],12-1);                 % T-Score
                        CI_increment(q,:) = x + ts*SEM;                      %absolute CI
                        %yCI95 = bsxfun(@times, SEM, ts(:));   %value above and below
                        end
                   end
                    

 %% Calculate Annual risk (across years and across 12 GCMs) and percentile range
                   annMortalityUK = zeros(nRows, years);                
                   for h = 1:nGridCells
                        annMortality{h} = squeeze(nansum(reshape(dailyMort{h}.',daysPerYear,size(dailyMort{h},2)./daysPerYear,[]))).';
                        annMortalityUK(:,:) = annMortalityUK(:,:) + annMortality{h}(:,:);
                        for r = 1:nRows
                            avMortality(r,h) = mean(annMortality{h}(r,:),[1,2]); %mean across years for each GCM
                        end
                        p_low = prctile(avMortality,(10),1);%10th and 90th percentile across years/GCMs
                        p_high = prctile(avMortality,(90),1);
                        av_ann_rcms = mean(annMortalityUK,2);
                        av = mean(avMortality,1);
                        sum_av = sum(av);
                        sum_p_low = sum(p_low);
                        sum_p_high = sum(p_high);
                    end
                   
                        Column_vector = annMortalityUK(:); %Convert array to column vector                    
                        %calculating CI across the 12 RCMS per day not across
                        %the 2205 gridcells
                        % confidence interval
                        stDev = std(av_ann_rcms); 
                        x = mean(av_ann_rcms);
                        sqRoot = sqrt(12);
                        SEM = stDev/sqRoot;% 
                        ts = tinv([0.025  0.975],12-1);                 % T-Score
                        CI = x + ts*SEM;                      %absolute CI
                        yCI95 = bsxfun(@times, SEM, ts(:));   %value above and below

                        %population_aggregate = (sum(population{1,1}(:,3:7))).';
                        av_per_population = sum_av/population_aggregate(m,1+p)*100000; %mean deaths per year/per 100,000 population - as a proportion of total population
                        p_low_per_population = sum_p_low/population_aggregate(m,1+p)*100000; %10th P deaths per year/per 1000 population 
                        p_high_per_population = sum_p_high/population_aggregate(m,1+p)*100000; %10th P deaths per year/per 1000 population 
                        
                        p_low_per_population_CI = CI(1,1)/population_aggregate(m,1+p)*100000; %10th P deaths per year/per 1000 population 
                        p_high_per_population_CI = CI(1,2)/population_aggregate(m,1+p)*100000; %90th P deaths per year/per 1000 population 

                    %Save results for plotting -
                    if m == 1 % all ages
                        Results_1_5C{a}(1:h,3:5) = [av.', p_low.', p_high.']; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                        Results_1_5C{a}(1:h,1) = lonIDIndex(:,1); %longitude
                        Results_1_5C{a}(1:h,2) = latIDIndex(:,1); %latitude
                        annualMortalityUK_CIs_1_5{a} = CI;
                        Results_1_5C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.']; % UK aggregate per 100,000
                    elseif m==2
                       %Results_1_5C{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                        Results_1_5C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                    elseif m==3
                       %Results_1_5C{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                        Results_1_5C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                    elseif m==4
                       %Results_1_5C{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                        Results_1_5C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                    elseif m==5
                       %Results_1_5C{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                        Results_1_5C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                    end
                    h=1; %set back to 1
                end
                adapt = acclimaScen; %now runs with adaptation.
            end

            % Files to SAVE %%IN future update to include index on UDM run and
            % other types of adaptation....
            if l==1 %climate change only
            fname = sprintf('Results_1_5C_%s_%d_SSP%d.mat', adaptationName, populationYearbaseline, UK_SSP); 
            Results_1_5C_CC = Results_1_5C;
            annualMortalityUK_CIs_1_5_CC = annualMortalityUK_CIs_1_5;
            incrementUK_CIs_1_5_CC = CI_increment;
            Results_1_5C_AgeGroup_CC = Results_1_5C_AgeGroup;
            save(fname, 'Results_1_5C_CC', 'Results_1_5C_AgeGroup_CC', 'annualMortalityUK_CIs_1_5_CC', 'incrementUK_CIs_1_5_CC')

            elseif l==2 %Climate change and population change
            fname = sprintf('Results_1_5C%s_%s_%d_SSP%d.mat', '_SE', adaptationName, populationYear15degree, UK_SSP);
            Results_1_5C_SE = Results_1_5C;
            Results_1_5C_AgeGroup_SE = Results_1_5C_AgeGroup;
            annualMortalityUK_CIs_1_5_SE = annualMortalityUK_CIs_1_5;
            incrementUK_CIs_1_5_SE = CI_increment;
            save(fname,'Results_1_5C_SE', 'Results_1_5C_AgeGroup_SE', 'annualMortalityUK_CIs_1_5_SE', 'incrementUK_CIs_1_5_SE') 
            %clear('dailyMort','totMortPerIncrement', 'avMortPerIncrement','annMortality', 'p_low', 'p_high', 'av', 'sum_av', 'sum_p_low', 'sum_p_high', 'av_per_population', 'p_low_per_population', 'p_high_per_population')  %data not used past this point
            end
            
            if populationYear15degree == 2020
                p=1;
            elseif populationYear15degree == 2030
                p=2;
            elseif populationYear15degree == 2050
                p=3;
            elseif populationYear15degree == 2080
                p=4;
            end
        end
        
    else
        disp ('no input data - 1.5degree')
    end
    save ('avMortPerIncrementUK.mat', 'avMortPerIncrementUK') % Used to plot deaths per increment for paper.
end

%%  Calculate results for 2 degree
if (ClimateScenario == 3|| ClimateScenario == 6) %% 2degree or all climate scenarios
    if exist ('input_data_2.0C.mat', 'file')  %calculate for each time period - based on input files
        disp ('Calculating mortality risk - 2.0DegreeC')

        %sorted_data_2_0C = data_2_0C(ib,:,:); %sorts so that rows match the coordinate order of socio-economic data stored in heat_impact_data.mat
        p=0; %counter to make sure the correct population cell array is read for years when warming set to occur by user.

        for l = 1:socioEcScen % Compute for CC-only and SE+CC      

            Results_2C = cell(adaptScen,1); %append all final results for plotting and saving, for each age group
            for n = 1:adaptScen
                Results_2C{n} = zeros(nGridCells, 7); %set automatically in future based on gridcell rows
            end

            Results_2C_AgeGroup = cell(adaptScen,1); %append all final results for plotting and saving, for each age group and new array for each adaptation scenario
            for n = 1:adaptScen
                Results_2C_AgeGroup{n} = zeros(ageGroupCnt, 5); %set automatically in future based on gridcell rows
            end

            annualMortalityUK_CIs_2 = cell(adaptScen,1);
            for n = 1:adaptScen
            annualMortalityUK_CIs_2{n} = zeros(1, 2); 
            end
            
            adapt = 1; %first run assumes no adaptation

            for a = 1:adaptScen

                % Pre-allocate Cell Arrays to store intermediate results 1:adaptScen
                dailyMort = cell(nGridCells,1);   %store intermediate value of daily additional deaths         
                for n = 1:nGridCells
                    dailyMort{n} = zeros(nRows, nCols);
                end

                totMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (summed across all years/GCMs)             
                for n = 1:nGridCells
                    totMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                end

                annMortality = cell(nGridCells,1);   %store intermediate value of Mortality risk * population (average annual)      
                for n = 1:nGridCells
                    annMortality{n} = zeros(nRows, nCols);
                end

                avMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (average annual)    
                for n = 1:nGridCells
                avMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                end

                threshold(:,1) = Mortality_Threshold_RR(:,4) + mortality_acclimatisation{3,1}(:,adapt+2);  

                for m = 1:ageGroupCnt %for all ages and per age group
                i=1;
                while (h <= nGridCells)
                    if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) %looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                        for r = 1:nRows
                            for c = 1:nCols
                                if data_2_0C{i,2}(r,c) < threshold(h,1)
                                    % do nothing - speeds up if statement by specifying most likey option first..
                                elseif Mortality_Threshold_RR(h,4) > 0 
                                   dailyMort{h}(r,c) = ((((data_2_0C{i,2}(r,c)- threshold(h,1))*Mortality_Threshold_RR(h,m+4))/100)*dailyDeathRate(h,m+2))*population{1+p}(h,m+2); % Outputs the additional number of daily deaths %%increments 2 as population + dailyDeathRate includes lat and lon columns first. 
                                   incrementExceeded = round(data_2_0C{i,2}(r,c)- threshold(h,1));
                                   totMortPerIncrement{h}(r, incrementExceeded+1) = totMortPerIncrement{h}(r, incrementExceeded+1) + dailyMort{h}(r,c); %cumulative mortality per each 1 degree increment above the threshold
                                else % do nothing - zero value.
                                end
                            end
                        end
                        h=h+1; %increment through data in Mortality_Threshold_RR
                        i=1; %reset counter
                        r=1; %reset counter
                        c=1; %reset counter

                    else
                            i=i+1;
                    end
                end          

    %% average annual mortality per increment %%for all ages only

                    if (m==1) && (a==1)        
                        for h = 1:nGridCells
                           for r = 1:nRows
                                avMortPerIncrement{h}(r,:) = totMortPerIncrement{h}(r,:)/years;
                           end
                        end
                           dim = ndims(avMortPerIncrement{1});          % Get the number of dimensions for arrays
                           matrix = cat(dim+1,avMortPerIncrement{:});        % Convert to a (dim+1)-dimensional matrix
                           sumArray = sum(matrix,dim+1);  % Get the mean across 3d array (summing mortality across UK grid cells)
                           avMortPerIncrementUK(l+4,:) = mean(sumArray);
                           
                       %calculating CI across the 12 RCMS per increment
                       CI_increment = zeros(30,2);
                        for q = 1:30
                        stDev = std(sumArray(:,q)); 
                        x = mean(sumArray(:,q));
                        sqRoot = sqrt(12);
                        SEM = stDev/sqRoot;% 
                        ts = tinv([0.025  0.975],12-1);                 % T-Score
                        CI_increment(q,:) = x + ts*SEM;                      %absolute CI
                        %yCI95 = bsxfun(@times, SEM, ts(:));   %value above and below
                        end
                    end

 %% Calculate Annual risk (across years and across 12 GCMs) and percentile range
                    annMortalityUK = zeros(nRows, years);                
                   for h = 1:nGridCells
                        annMortality{h} = squeeze(nansum(reshape(dailyMort{h}.',daysPerYear,size(dailyMort{h},2)./daysPerYear,[]))).';
                        annMortalityUK(:,:) = annMortalityUK(:,:) + annMortality{h}(:,:);
                        for r = 1:nRows
                            avMortality(r,h) = mean(annMortality{h}(r,:),[1,2]); %mean across years for each GCM
                        end
                        p_low = prctile(avMortality,(10),1);%10th and 90th percentile across years/GCMs
                        p_high = prctile(avMortality,(90),1);
                        av_ann_rcms = mean(annMortalityUK,2);
                        av = mean(avMortality,1);
                        sum_av = sum(av);
                        sum_p_low = sum(p_low);
                        sum_p_high = sum(p_high);
                    end
                   
                        Column_vector = annMortalityUK(:); %Convert array to column vector                    
                        %calculating CI across the 12 RCMS per day not across
                        %the 2205 gridcells
                        % confidence interval
                        stDev = std(av_ann_rcms); 
                        x = mean(av_ann_rcms);
                        sqRoot = sqrt(12);
                        SEM = stDev/sqRoot;% 
                        ts = tinv([0.025  0.975],12-1);                 % T-Score
                        CI = x + ts*SEM;                      %absolute CI
                        yCI95 = bsxfun(@times, SEM, ts(:));   %value above and below
                    
                    %population_aggregate = (sum(population{1,1}(:,3:7))).';
                    av_per_population = sum_av/population_aggregate(m,1+p)*100000; %mean deaths per year/per 100,000 population - as a proportion of total population
                    p_low_per_population = sum_p_low/population_aggregate(m,1+p)*100000; %10th P deaths per year/per 1000 population 
                    p_high_per_population = sum_p_high/population_aggregate(m,1+p)*100000; %10th P deaths per year/per 1000 population 

                    p_low_per_population_CI = CI(1,1)/population_aggregate(m,1+p)*100000; %10th P deaths per year/per 1000 population 
                    p_high_per_population_CI = CI(1,2)/population_aggregate(m,1+p)*100000; %90th P deaths per year/per 1000 population 

                    %Save results for plotting -
                    if m == 1
                        Results_2C{a}(1:h,3:5) = [av.', p_low.', p_high.']; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                        Results_2C{a}(1:h,1) = lonIDIndex(:,1); %longitude
                        Results_2C{a}(1:h,2) = latIDIndex(:,1); %latitude
                        annualMortalityUK_CIs_2{a} = CI;
                        Results_2C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.']; % UK aggregate per 100,000
                    elseif m==2
                        %Results_2C{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                        Results_2C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                    elseif m==3
                        %Results_2C{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                        Results_2C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                    elseif m==4
                        %Results_2C{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                        Results_2C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                    elseif m==5
                        %Results_2C{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                        Results_2C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                    end
                    h=1; %set back to 1
                end
                adapt = acclimaScen; %now runs with adaptation.
            end
                % Files to SAVE
            if l==1
            fname = sprintf('Results_2C_%s_%d_SSP%d.mat', adaptationName,populationYearbaseline,UK_SSP);
            Results_2C_CC = Results_2C;
            Results_2C_AgeGroup_CC = Results_2C_AgeGroup;
            annualMortalityUK_CIs_2_CC = annualMortalityUK_CIs_2;
            incrementUK_CIs_2_CC = CI_increment;
            save(fname, 'Results_2C_CC', 'Results_2C_AgeGroup_CC', 'annualMortalityUK_CIs_2_CC', 'incrementUK_CIs_2_CC')
            elseif l==2
            fname = sprintf('Results_2C%s_%s_%d_SSP%d.mat', '_SE', adaptationName, populationYear2degree,UK_SSP);
            Results_2C_SE = Results_2C;
            Results_2C_AgeGroup_SE = Results_2C_AgeGroup;
            annualMortalityUK_CIs_2_SE = annualMortalityUK_CIs_2;
            incrementUK_CIs_2_SE = CI_increment;
            save(fname, 'Results_2C_SE', 'Results_2C_AgeGroup_SE', 'annualMortalityUK_CIs_2_SE', 'incrementUK_CIs_2_SE') 
            clear('dailyMort','totMortPerIncrement', 'avMortPerIncrement','annMortality', 'p_low', 'p_high', 'av', 'sum_av', 'sum_p_low', 'sum_p_high', 'av_per_population', 'p_low_per_population', 'p_high_per_population')  %data not used past this point
            end

            if populationYear2degree == 2030
                p=2;
            elseif populationYear2degree == 2050
                p=3;
            elseif populationYear2degree == 2080
                p=4;
            end
        end
     else
        disp ('no input data - 2degree')
     end
    save ('avMortPerIncrementUK.mat', 'avMortPerIncrementUK') % Used to plot deaths per increment for paper.
end

%% Calculate results for 3 degree
if (ClimateScenario == 4|| ClimateScenario == 6) %% 3degree or all climate scenarios
    if exist ('input_data_3.0C.mat', 'file')  %calculate for each time period - based on input files
        disp ('Calculating mortality risk - 3.0DegreeC')

       %sorted_data_3_0C = data_3_0C(ib,:,:); %sorts so that rows match the coordinate order of socio-economic data stored in heat_impact_data.mat
       p=0; %counter to make sure the correct population cell array is read for 2050s.

        for l = 1:socioEcScen % Compute for CC-only and SE+CC

        %set up cell array
        Results_3C = cell(adaptScen,1); %append all final results for plotting and saving, for each age group
        for n = 1:adaptScen
            Results_3C{n} = zeros(nGridCells, 7); %set automatically in future based on gridcell rows
        end

        Results_3C_AgeGroup = cell(adaptScen,1); %append all final results for plotting and saving, for each age group and new array for each adaptation scenario
        for n = 1:adaptScen
            Results_3C_AgeGroup{n} = zeros(ageGroupCnt, 5); %set automatically in future based on gridcell rows
        end
        
        annualMortalityUK_CIs_3 = cell(adaptScen,1);
        for n = 1:adaptScen
        annualMortalityUK_CIs_3{n} = zeros(1, 2); 
        end

        adapt = 1; %first run assumes no adaptation
            for a = 1:adaptScen

                dailyMort = cell(nGridCells,1);   %store intermediate value of daily additional deaths         
                for n = 1:nGridCells
                    dailyMort{n} = zeros(nRows, nCols);
                end

                totMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (summed across all years/GCMs)            
                for n = 1:nGridCells
                    totMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                end

                avMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (average annual) 
                for n = 1:nGridCells
                    avMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                end

                annMortality = cell(nGridCells,1);   %store intermediate value of Mortality risk * population (average annual)       
                for n = 1:nGridCells
                    annMortality{n} = zeros(nRows, nCols);
                end

                threshold(:,1) = Mortality_Threshold_RR(:,4) + mortality_acclimatisation{4,1}(:,adapt+2);  

                for m = 1:ageGroupCnt %for all ages and per age group
                i=1;
                while (h <= nGridCells)
                    if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) %looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                        for r = 1:nRows
                            for c = 1:nCols
                                if data_3_0C{i,2}(r,c) < threshold(h,1)
                                    % do nothing - speeds up if statement by specifying most likey option first..
                                elseif Mortality_Threshold_RR(h,4) > 0 
                                   dailyMort{h}(r,c) = ((((data_3_0C{i,2}(r,c)- threshold(h,1))*Mortality_Threshold_RR(h,m+4))/100)*dailyDeathRate(h,m+2))*population{1+p}(h,m+2); % Outputs the additional number of daily deaths %%increments 2 as population + dailyDeathRate includes lat and lon columns first. 
                                   incrementExceeded = round(data_3_0C{i,2}(r,c)- threshold(h,1));
                                   totMortPerIncrement{h}(r, incrementExceeded+1) = totMortPerIncrement{h}(r, incrementExceeded+1) + dailyMort{h}(r,c); %cumulative mortality per each 1 degree increment above the threshold
                                else % do nothing - zero value.
                                end
                            end
                        end
                        h=h+1; %increment through data in Mortality_Threshold_RR
                        i=1; %reset counter
                        r=1; %reset counter
                        c=1; %reset counter

                    else
                            i=i+1;
                    end
                end  

    %% average annual mortality per increment %%for all ages only 
                    if (m==1) && (a==1)            
                        for h = 1:nGridCells
                           for r = 1:nRows
                                avMortPerIncrement{h}(r,:) = totMortPerIncrement{h}(r,:)/years;
                           end   
                        end
                           dim = ndims(avMortPerIncrement{1});          % Get the number of dimensions for arrays
                           matrix = cat(dim+1,avMortPerIncrement{:});        % Convert to a (dim+1)-dimensional matrix
                           sumArray = sum(matrix,dim+1);  % sum mortality across UK grid cells
                           avMortPerIncrementUK(l+6,:) = mean(sumArray);
                           
                       %calculating CI across the 12 RCMS per increment
                       CI_increment = zeros(30,2);
                        for q = 1:30
                        stDev = std(sumArray(:,q)); 
                        x = mean(sumArray(:,q));
                        sqRoot = sqrt(12);
                        SEM = stDev/sqRoot;% 
                        ts = tinv([0.025  0.975],12-1);                 % T-Score
                        CI_increment(q,:) = x + ts*SEM;                      %absolute CI
                        %yCI95 = bsxfun(@times, SEM, ts(:));   %value above and below
                        end
                    end

 %% Calculate Annual risk (across years and across 12 GCMs) and percentile range
                    annMortalityUK = zeros(nRows, years);                
                   for h = 1:nGridCells
                        annMortality{h} = squeeze(nansum(reshape(dailyMort{h}.',daysPerYear,size(dailyMort{h},2)./daysPerYear,[]))).';
                        annMortalityUK(:,:) = annMortalityUK(:,:) + annMortality{h}(:,:);
                        for r = 1:nRows
                            avMortality(r,h) = mean(annMortality{h}(r,:),[1,2]); %mean across years for each GCM
                        end
                        p_low = prctile(avMortality,(10),1);%10th and 90th percentile across years/GCMs
                        p_high = prctile(avMortality,(90),1);
                        av_ann_rcms = mean(annMortalityUK,2);
                        av = mean(avMortality,1);
                        sum_av = sum(av);
                        sum_p_low = sum(p_low);
                        sum_p_high = sum(p_high);
                    end
                   
                        Column_vector = annMortalityUK(:); %Convert array to column vector                    
                        %calculating CI across the 12 RCMS per day not across
                        %the 2205 gridcells
                        % confidence interval
                        stDev = std(av_ann_rcms); 
                        x = mean(av_ann_rcms);
                        sqRoot = sqrt(12);
                        SEM = stDev/sqRoot;% 
                        ts = tinv([0.025  0.975],12-1);                 % T-Score
                        CI = x + ts*SEM;                      %absolute CI
                        yCI95 = bsxfun(@times, SEM, ts(:));   %value above and below
       
                    
                    %population_aggregate = (sum(population{1,1}(:,3:7))).';
                    av_per_population = sum_av/population_aggregate(m,1+p)*100000; %mean deaths per year/per 100,000 population - as a proportion of total population
                    p_low_per_population = sum_p_low/population_aggregate(m,1+p)*100000; %10th P deaths per year/per 1000 population 
                    p_high_per_population = sum_p_high/population_aggregate(m,1+p)*100000; %10th P deaths per year/per 1000 population 

                    p_low_per_population_CI = CI(1,1)/population_aggregate(m,1+p)*100000; %10th P deaths per year/per 1000 population 
                    p_high_per_population_CI = CI(1,2)/population_aggregate(m,1+p)*100000; %90th P deaths per year/per 1000 population 

                    %Save results for plotting -
                    if m == 1
                        Results_3C{a}(1:h,3:5) = [av.', p_low.', p_high.']; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                        Results_3C{a}(1:h,1) = lonIDIndex(:,1); %longitude
                        Results_3C{a}(1:h,2) = latIDIndex(:,1); %latitude
                        annualMortalityUK_CIs_3{a} = CI;
                        Results_3C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.']; % UK aggregate per 100,000
                    elseif m==2
                       %Results_3C{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                        Results_3C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                    elseif m==3
                        %Results_3C{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                        Results_3C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                    elseif m==4
                        %Results_3C{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                        Results_3C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                    elseif m==5
                        %Results_3C{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                        Results_3C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                    end
                    h=1; %set back to 1
                end
                adapt = acclimaScen; %now runs with adaptation.
            end

            % Files to SAVE
            if l==1
            fname = sprintf('Results_3C_%s_%d_SSP%d.mat', adaptationName, populationYearbaseline, UK_SSP);
            Results_3C_CC = Results_3C;
            Results_3C_AgeGroup_CC = Results_3C_AgeGroup;
            annualMortalityUK_CIs_3_CC = annualMortalityUK_CIs_3;
            incrementUK_CIs_3_CC = CI_increment;
            save(fname, 'Results_3C_CC', 'Results_3C_AgeGroup_CC', 'annualMortalityUK_CIs_3_CC', 'incrementUK_CIs_3_CC')
            clear('dailyMort','totMortPerIncrement', 'avMortPerIncrement','annMortality', 'p_low', 'p_high', 'av', 'sum_av', 'sum_p_low', 'sum_p_high', 'av_per_population', 'p_low_per_population', 'p_high_per_population')  %daily data not used past this point
            elseif l==2
            fname = sprintf('Results_3C%s_%s_%d_SSP%d.mat', '_SE', adaptationName, populationYear3degree, UK_SSP);
            Results_3C_SE = Results_3C;
            Results_3C_AgeGroup_SE = Results_3C_AgeGroup;
            annualMortalityUK_CIs_3_SE = annualMortalityUK_CIs_3;
            incrementUK_CIs_3_SE = CI_increment;
            save(fname,'Results_3C_SE', 'Results_3C_AgeGroup_SE', 'annualMortalityUK_CIs_3_SE', 'incrementUK_CIs_3_SE') 
            clear('dailyMort','totMortPerIncrement', 'avMortPerIncrement','annMortality', 'p_low', 'p_high', 'av', 'sum_av', 'sum_p_low', 'sum_p_high', 'av_per_population', 'p_low_per_population', 'p_high_per_population')  %data not used past this point
            end
            
            if  populationYear3degree == 2050
                p=3;
            elseif populationYear3degree == 2080
                p=4;
            end
        end
    else
        disp ('no input data - 3degree')
    end
    save ('avMortPerIncrementUK.mat', 'avMortPerIncrementUK') % Used to plot deaths per increment for paper.
end

    
%% Calculate results for 4 degree
if (ClimateScenario == 5|| ClimateScenario == 6) %% 4degree or all climate scenarios
    if exist ('input_data_4.0C.mat', 'file')  %calculate for each time period - based on input files
        disp ('Calculating mortality risk - 4.0DegreeC')

       %sorted_data_4_0C = data_4_0C(ib,:,:); %sorts so that rows match the coordinate order of socio-economic data stored in heat_impact_data.mat
       p=0; %counter to make sure the correct population cell array is read for 2050s.

        for l = 1:socioEcScen % Compute for CC-only and SE+CC

        %set up cell array
        Results_4C = cell(adaptScen,1); %append all final results for plotting and saving, for each age group
        for n = 1:adaptScen
            Results_4C{n} = zeros(nGridCells, 7); %set automatically in future based on gridcell rows
        end

        Results_4C_AgeGroup = cell(adaptScen,1); %append all final results for plotting and saving, for each age group and new array for each adaptation scenario
        for n = 1:adaptScen
            Results_4C_AgeGroup{n} = zeros(ageGroupCnt, 5); %set automatically in future based on gridcell rows
        end

        annualMortalityUK_CIs_4 = cell(adaptScen,1);
        for n = 1:adaptScen
        annualMortalityUK_CIs_4{n} = zeros(1, 2); 
        end
        
        adapt = 1; %first run assumes no adaptation
            for a = 1:adaptScen

                dailyMort = cell(nGridCells,1);   %store intermediate value of daily additional deaths         
                for n = 1:nGridCells
                    dailyMort{n} = zeros(nRows, nCols);
                end

                totMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (summed across all years/GCMs)            
                for n = 1:nGridCells
                    totMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                end

                avMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (average annual) 
                for n = 1:nGridCells
                    avMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                end

                annMortality = cell(nGridCells,1);   %store intermediate value of Mortality risk * population (average annual)       
                for n = 1:nGridCells
                    annMortality{n} = zeros(nRows, nCols);
                end

                threshold(:,1) = Mortality_Threshold_RR(:,4) + mortality_acclimatisation{5,1}(:,adapt+2);  

                for m = 1:ageGroupCnt %for all ages and per age group
                i=1;
                while (h <= nGridCells)
                    if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) %looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                        for r = 1:nRows
                            for c = 1:nCols
                                if data_4_0C{i,2}(r,c) < threshold(h,1)
                                    % do nothing - speeds up if statement by specifying most likey option first..
                                elseif Mortality_Threshold_RR(h,4) > 0 
                                   dailyMort{h}(r,c) = ((((data_4_0C{i,2}(r,c)- threshold(h,1))*Mortality_Threshold_RR(h,m+4))/100)*dailyDeathRate(h,m+2))*population{1+p}(h,m+2); % Outputs the additional number of daily deaths %%increments 2 as population + dailyDeathRate includes lat and lon columns first. 
                                   incrementExceeded = round(data_4_0C{i,2}(r,c)- threshold(h,1));
                                   totMortPerIncrement{h}(r, incrementExceeded+1) = totMortPerIncrement{h}(r, incrementExceeded+1) + dailyMort{h}(r,c); %cumulative mortality per each 1 degree increment above the threshold
                                else % do nothing - zero value.
                                end
                            end
                        end
                        h=h+1; %increment through data in Mortality_Threshold_RR
                        i=1; %reset counter
                        r=1; %reset counter
                        c=1; %reset counter

                    else
                            i=i+1;
                    end
                end  

    %% average annual mortality per increment %%for all ages only
    annMortalityUK = zeros(nRows, years); 
                    if (m==1) && (a==1)            
                        for h = 1:nGridCells
                           for r = 1:nRows
                                avMortPerIncrement{h}(r,:) = totMortPerIncrement{h}(r,:)/years;
                           end   
                        end
                           dim = ndims(avMortPerIncrement{1});          % Get the number of dimensions for arrays
                           matrix = cat(dim+1,avMortPerIncrement{:});        % Convert to a (dim+1)-dimensional matrix
                           sumArray = sum(matrix,dim+1);  % sum mortality across UK grid cells
                           avMortPerIncrementUK(l+8,:) = mean(sumArray);
                           
                       %calculating CI across the 12 RCMS per increment
                       CI_increment = zeros(30,2);
                        for q = 1:30
                        stDev = std(sumArray(:,q)); 
                        x = mean(sumArray(:,q));
                        sqRoot = sqrt(12);
                        SEM = stDev/sqRoot;% 
                        ts = tinv([0.025  0.975],12-1);                 % T-Score
                        CI_increment(q,:) = x + ts*SEM;                      %absolute CI
                        %yCI95 = bsxfun(@times, SEM, ts(:));   %value above and below
                        end
                    end

    %% Calculate Annual risk (across years and across 12 GCMs) and percentile range
                    annMortalityUK = zeros(nRows, years);
                    CI_spatial = zeros(nGridCells,2);
                   for h = 1:nGridCells
                        annMortality{h} = squeeze(nansum(reshape(dailyMort{h}.',daysPerYear,size(dailyMort{h},2)./daysPerYear,[]))).';
                        annMortalityUK(:,:) = annMortalityUK(:,:) + annMortality{h}(:,:);
                        for r = 1:nRows
                            avMortality(r,h) = mean(annMortality{h}(r,:),[1,2]); %mean across years for each GCM
                        end
                        p_low = prctile(avMortality,(10),1);%10th and 90th percentile across years/GCMs
                        p_high = prctile(avMortality,(90),1);
                        av_ann_rcms = mean(annMortalityUK,2);
                        av = mean(avMortality,1);
                        sum_av = sum(av);
                        sum_p_low = sum(p_low);
                        sum_p_high = sum(p_high);
                        
                        stDev = std(avMortality(:, h)); 
                        x = mean(avMortality(:, h));
                        sqRoot = sqrt(12);
                        SEM = stDev/sqRoot;% 
                        ts = tinv([0.025  0.975],12-1);                 % T-Score
                        CI_spatial(h,:) = x + ts*SEM;                      %absolute CI
                        %yCI95 = bsxfun(@times, SEM, ts(:));   %value above and below
                        
                    end
                   
                        %Column_vector = annMortalityUK(:); %Convert array to column vector                    
                        %calculating CI across the 12 RCMS per day not across
                        %the 2205 gridcells
                        % confidence interval
                        stDev = std(av_ann_rcms); 
                        x = mean(av_ann_rcms);
                        sqRoot = sqrt(12);
                        SEM = stDev/sqRoot;% 
                        ts = tinv([0.025  0.975],12-1);                 % T-Score
                        CI = x + ts*SEM;                      %absolute CI
                        yCI95 = bsxfun(@times, SEM, ts(:));   %value above and below
                        
                    %population_aggregate = (sum(population{1,1}(:,3:7))).';
                    av_per_population = sum_av/population_aggregate(m,1+p)*100000; %mean deaths per year/per 100,000 population - as a proportion of total population
                    p_low_per_population = sum_p_low/population_aggregate(m,1+p)*100000; %10th P deaths per year/per 1000 population 
                    p_high_per_population = sum_p_high/population_aggregate(m,1+p)*100000; %10th P deaths per year/per 1000 population 

                    p_low_per_population_CI = CI(1,1)/population_aggregate(m,1+p)*100000; %10th P deaths per year/per 1000 population 
                    p_high_per_population_CI = CI(1,2)/population_aggregate(m,1+p)*100000; %90th P deaths per year/per 1000 population 

                    %Save results for plotting -
                    if m == 1
                        Results_4C{a}(1:h,3:7) = [av.', p_low.', p_high.', CI_spatial]; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                        Results_4C{a}(1:h,1) = lonIDIndex(:,1); %longitude
                        Results_4C{a}(1:h,2) = latIDIndex(:,1); %latitude
                        annualMortalityUK_CIs_4{a} = CI;
                        Results_4C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.']; % UK aggregate per 100,000
                    elseif m==2
                       %Results_3C{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                        Results_4C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                    elseif m==3
                        %Results_3C{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                        Results_4C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                    elseif m==4
                        %Results_4C{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                        Results_4C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                    elseif m==5
                        %Results_4C{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                        Results_4C_AgeGroup{a}(m, 1:5) = [av_per_population.', p_low_per_population.', p_high_per_population.', p_low_per_population_CI.', p_high_per_population_CI.'];
                    end
                    h=1; %set back to 1
                end
                adapt = acclimaScen; %now runs with adaptation.
            end

            % Files to SAVE
            if l==1
            fname = sprintf('Results_4C_%s_%d_SSP%d.mat', adaptationName, populationYearbaseline, UK_SSP);
            Results_4C_CC = Results_4C;
            Results_4C_AgeGroup_CC = Results_4C_AgeGroup;
            annualMortalityUK_CIs_4_CC = annualMortalityUK_CIs_4;
            incrementUK_CIs_4_CC = CI_increment;
            save(fname, 'Results_4C_CC', 'Results_4C_AgeGroup_CC', 'annualMortalityUK_CIs_4_CC', 'incrementUK_CIs_4_CC')
            clear('dailyMort','totMortPerIncrement', 'avMortPerIncrement','annMortality', 'p_low', 'p_high', 'av', 'sum_av', 'sum_p_low', 'sum_p_high', 'av_per_population', 'p_low_per_population', 'p_high_per_population')  %daily data not used past this point
            elseif l==2
            fname = sprintf('Results_4C%s_%s_%d_SSP%d.mat', '_SE', adaptationName, populationYear4degree, UK_SSP);
            Results_4C_SE = Results_4C;
            Results_4C_AgeGroup_SE = Results_4C_AgeGroup;
            annualMortalityUK_CIs_4_SE = annualMortalityUK_CIs_4;
            incrementUK_CIs_4_SE = CI_increment;
            save(fname,'Results_4C_SE', 'Results_4C_AgeGroup_SE', 'annualMortalityUK_CIs_4_SE', 'incrementUK_CIs_4_SE') 
            clear('dailyMort','totMortPerIncrement', 'avMortPerIncrement','annMortality', 'p_low', 'p_high', 'av', 'sum_av', 'sum_p_low', 'sum_p_high', 'av_per_population', 'p_low_per_population', 'p_high_per_population')  %data not used past this point
            end
            
            p=4; 
        end
    else
        disp ('no input data - 4degree')
    end
    save ('avMortPerIncrementUK.mat', 'avMortPerIncrementUK') % Used to plot deaths per increment for paper.
end

%% Amalgamate results for plotting - UK data based on sum of gridded data (All ages and by age group)Average, 10th and 90th P summed.
%Absolute values also split by age group - differs to mortality per 100,000
%population per age group.
% Only run when all climate scenarios run.
    if ClimateScenario == 6

        AvAnnualDeathsUK = zeros(2,15);
        AvAnnualDeathsUK_SE = zeros(2,15);
        AvAnnualDeathsAdapt = zeros(2,5);
        AvAnnualDeathsAdapt_SE = zeros(2,5);
        CIs_CC = zeros(2,10);
        CIs_SE = zeros(2,10);

        for a = 1:2
            % For all ages
            AvAnnualDeathsUK(a,1:3) = [sum(ResultsPast{a}(:,3)), sum(ResultsPast{a}(:,6)), sum(ResultsPast{a}(:,7))]; % Near Past
            AvAnnualDeathsUK(a,4:6) = [sum(Results_1_5C_CC{a}(:,3)), sum(Results_1_5C_CC{a}(:,6)), sum(Results_1_5C_CC{a}(:,7))]; % 1.5 degree
            AvAnnualDeathsUK(a,7:9) = [sum(Results_2C_CC{a}(:,3)), sum(Results_2C_CC{a}(:,6)), sum(Results_2C_CC{a}(:,7))]; % 2 degree C
            AvAnnualDeathsUK(a,10:12) = [sum(Results_3C_CC{a}(:,3)), sum(Results_3C_CC{a}(:,6)), sum(Results_3C_CC{a}(:,7))]; % 3 degree C
            AvAnnualDeathsUK(a,13:15) = [sum(Results_4C_CC{a}(:,3)), sum(Results_4C_CC{a}(:,6)), sum(Results_4C_CC{a}(:,7))]; % 4 degree C
            
            AvAnnualDeathsUK_SE(a,1:3) = [sum(ResultsPast{a}(:,3)), sum(ResultsPast{a}(:,6)), sum(ResultsPast{a}(:,7))];% Near Past (population same in both scenarios)
            AvAnnualDeathsUK_SE(a,4:6) = [sum(Results_1_5C_SE{a}(:,3)), sum(Results_1_5C_SE{a}(:,6)), sum(Results_1_5C_SE{a}(:,7))]; % 1.5 degree CC+SE
            AvAnnualDeathsUK_SE(a,7:9) = [sum(Results_2C_SE{a}(:,3)), sum(Results_2C_SE{a}(:,6)), sum(Results_2C_SE{a}(:,7))]; %2 degree C CC+SE
            AvAnnualDeathsUK_SE(a,10:12) = [sum(Results_3C_SE{a}(:,3)), sum(Results_3C_SE{a}(:,6)), sum(Results_3C_SE{a}(:,7))]; % 3 degree C CC+SE
            AvAnnualDeathsUK_SE(a,13:15) = [sum(Results_4C_SE{a}(:,3)), sum(Results_4C_SE{a}(:,6)), sum(Results_4C_SE{a}(:,7))]; % 4 degree C CC+SE

            CIs_CC(a, 1:2) = [sum(ResultsPast{a}(:,3))- annualMortalityUK_CIs{a}(1,1),annualMortalityUK_CIs{1}(1,2)- sum(ResultsPast{a}(:,3))]; 
            CIs_CC(a, 3:4) = [sum(Results_1_5C_CC{a}(:,3))- annualMortalityUK_CIs_1_5_CC{a}(1,1),annualMortalityUK_CIs_1_5_CC{a}(1,2)- sum(Results_1_5C_CC{a}(:,3))]; 
            CIs_CC(a, 5:6) = [sum(Results_2C_CC{a}(:,3))- annualMortalityUK_CIs_2_CC{a}(1,1),annualMortalityUK_CIs_2_CC{a}(1,2)- sum(Results_2C_CC{a}(:,3))]; 
            CIs_CC(a, 7:8) = [sum(Results_3C_CC{a}(:,3))- annualMortalityUK_CIs_3_CC{a}(1,1),annualMortalityUK_CIs_3_CC{a}(1,2)- sum(Results_3C_CC{a}(:,3))]; 
            CIs_CC(a, 9:10) = [sum(Results_4C_CC{a}(:,3))- annualMortalityUK_CIs_4_CC{a}(1,1),annualMortalityUK_CIs_4_CC{a}(1,2)- sum(Results_4C_CC{a}(:,3))]; 
            
            CIs_SE(a, 1:2) = [sum(ResultsPast{a}(:,3))- annualMortalityUK_CIs{a}(1,1),annualMortalityUK_CIs{1}(1,2)- sum(ResultsPast{a}(:,3))]; 
            CIs_SE(a, 3:4) = [sum(Results_1_5C_SE{a}(:,3))- annualMortalityUK_CIs_1_5_SE{a}(1,1),annualMortalityUK_CIs_1_5_SE{a}(1,2)- sum(Results_1_5C_SE{a}(:,3))]; 
            CIs_SE(a, 5:6) = [sum(Results_2C_SE{a}(:,3))- annualMortalityUK_CIs_2_SE{a}(1,1),annualMortalityUK_CIs_2_SE{a}(1,2)- sum(Results_2C_SE{a}(:,3))]; 
            CIs_SE(a, 7:8) = [sum(Results_3C_SE{a}(:,3))- annualMortalityUK_CIs_3_SE{a}(1,1),annualMortalityUK_CIs_3_SE{a}(1,2)- sum(Results_3C_SE{a}(:,3))]; 
            CIs_SE(a, 9:10) = [sum(Results_4C_SE{a}(:,3))- annualMortalityUK_CIs_4_SE{a}(1,1),annualMortalityUK_CIs_4_SE{a}(1,2)- sum(Results_4C_SE{a}(:,3))]; 
        end

            %difference between NO adaptation and WITH adaptation strategies for stacked
            %plot (rows = adapt scenarios, columns = time-periods).
            AvAnnualDeathsAdapt(1,1:5) = [sum(ResultsPast{2}(:,3)), sum(Results_1_5C_CC{2}(:,3)), sum(Results_2C_CC{2}(:,3)), sum(Results_3C_CC{2}(:,3)), sum(Results_4C_CC{2}(:,3))];
            AvAnnualDeathsAdapt(2,1:5) = [sum(ResultsPast{1}(:,3))- sum(ResultsPast{2}(:,3)), sum(Results_1_5C_CC{1}(:,3))-sum(Results_1_5C_CC{2}(:,3)), sum(Results_2C_CC{1}(:,3))-sum(Results_2C_CC{2}(:,3)), sum(Results_3C_CC{1}(:,3)) - sum(Results_3C_CC{2}(:,3)), sum(Results_4C_CC{1}(:,3)) - sum(Results_4C_CC{2}(:,3))]; 
            AvAnnualDeathsAdapt_SE(1,1:5) = [sum(ResultsPast{2}(:,3)), sum(Results_1_5C_SE{2}(:,3)), sum(Results_2C_SE{2}(:,3)), sum(Results_3C_SE{2}(:,3)),sum(Results_4C_SE{2}(:,3))];
            AvAnnualDeathsAdapt_SE(2,1:5) = [sum(ResultsPast{1}(:,3))-sum(ResultsPast{2}(:,3)), sum(Results_1_5C_SE{1}(:,3))-sum(Results_1_5C_SE{2}(:,3)), sum(Results_2C_SE{1}(:,3))-sum(Results_2C_SE{2}(:,3)), sum(Results_3C_SE{1}(:,3))-sum(Results_3C_SE{2}(:,3)), sum(Results_4C_SE{1}(:,3))-sum(Results_4C_SE{2}(:,3))];

            %% PLOTS
            % Create some initial output figures

            %1. Average annual deaths for each scenario (All ages) CC-only and CC+SE

            figure;
            M = max(AvAnnualDeathsUK_SE(1,:)*1.1); %set y axis max limit across sub-plots
            x=1:5;
            subplot(1,2,1);

            % Create bar and set individual colours
            b = bar(AvAnnualDeathsUK(1,[1,4,7,10,13]), 'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);
            b.FaceColor = 'flat';
            b.CData(1,:) = [1 0.843137264251709 0];
            b.CData(2,:) = [0.87058824300766 0.490196079015732 0];
            b.CData(3,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
            b.CData(4,:) = [0.600000023841858 0.200000002980232 0];
            b.CData(5,:) = [0.40,0.00,0.00];

            hold on
            p_errhigh =  CIs_CC(1,[2,4,6,8,10]); %90th p - average
            p_errlow =  CIs_CC(1,[1,3,5,7,9]); % average - 10th p)
            %p_errhigh = AvAnnualDeathsUK(1,[3,6,9,12,15])- AvAnnualDeathsUK(1,[1,4,7,10,13]); %90th p - average
            %p_errlow = AvAnnualDeathsUK(1,[1,4,7,10,13])-AvAnnualDeathsUK(1,[2,5,8,11,14]); % average - 10th p)

            er = errorbar(x, AvAnnualDeathsUK(1,[1,4,7,10,13]), p_errlow, p_errhigh, 'LineWidth',1); %CI95%
            er.Color = [0 0 0];                            
            er.LineStyle = 'none';

            % Create ylabel
            ylabel('Average Annual Heat Related Deaths');

            % Create xlabel
            xlabel('Climate Scenario');

            % Create title
            title('Climate Change Only');

            % Set the remaining axes properties
            set (gca, 'ylim', [0 M]);
            set(gca,'XTick',[1 2 3 4 5],'XTickLabel',{'Past','1.5','2.0','3.0','4.0'},'YGrid',...
                'on');
            ax = gca;
            ax.YAxis.Exponent = 0;

            subplot(1,2,2);
            % Create bar
            b= bar(AvAnnualDeathsUK_SE(1,[1,4,7,10,13]), 'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);
            b.FaceColor = 'flat';
            b.CData(1,:) = [1 0.843137264251709 0];
            b.CData(2,:) = [0.87058824300766 0.490196079015732 0];
            b.CData(3,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
            b.CData(4,:) = [0.600000023841858 0.200000002980232 0];
            b.CData(5,:) = [0.40,0.00,0.00];

            hold on
            p_errhighSE =  CIs_SE(1,[2,4,6,8,10]); %90th p - average
            p_errlowSE =  CIs_SE(1,[1,3,5,7,9]); % average - 10th p)
%             p_errhighSE = AvAnnualDeathsUK_SE(1,[3,6,9,12,15])- AvAnnualDeathsUK_SE(1,[1,4,7,10,13]); %90th p - average
%             p_errlowSE = AvAnnualDeathsUK_SE(1,[1,4,7,10,13])-AvAnnualDeathsUK_SE(1,[2,5,8,11,14]); % average - 10th p)

            er = errorbar(x, AvAnnualDeathsUK_SE(1,[1,4,7,10,13]), p_errlowSE, p_errhighSE, 'LineWidth',1);    
            er.Color = [0 0 0];                            
            er.LineStyle = 'none';

            % Create ylabel
            %ylabel('Average Annual Heat Related Deaths');

            % Create xlabel
            xlabel('Climate Scenario');

            % Create title
            title('Climate and Population Change');

            % Set the remaining axes properties
            set (gca, 'ylim', [0 M])
            set(gca,'XTick',[1 2 3 4 5],'XTickLabel',{'Past','1.5','2.0','3.0','4.0'},'YGrid',...
                'on');
            ax = gca;
            ax.YAxis.Exponent = 0;

           
%             %% 2 versus 4 degree for paper
%             
%             %twovsfour
%             
%             figure;
%             x=1:3;
% 
%             % Create bar and set individual colours
%             b = bar(twovsfour(1,1:3), 'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);
%             b.FaceColor = 'flat';
%             b.CData(1,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
%             b.CData(2,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
%             b.CData(3,:) = [0.40,0.00,0.00];
% 
%             hold on
%             p_errhigh = twovsfour(3,1:3)- twovsfour(1,1:3); %90th p - average
%             p_errlow = twovsfour(1,1:3)-twovsfour(2,1:3); % average - 10th p)
% 
%             er = errorbar(x, twovsfour(1,1:3), p_errlow, p_errhigh, 'LineWidth',1);    
%             er.Color = [0 0 0];                            
%             er.LineStyle = 'none';
% 
%             % Create ylabel
%             ylabel('Average Annual Heat Related Deaths');
% 
%             % Create xlabel
%             xlabel('Climate Scenario');
% 
%             % Create title
%             title('Climate and Population Change');
% 
%             % Set the remaining axes properties
%             %set (gca, 'ylim', [0 M]);
%             set(gca,'XTick',[1 2 3],'XTickLabel',{'2.0 2050','2.0 2080','4.0 2080'},'YGrid',...
%                 'on');
%             ax = gca;
%             ax.YAxis.Exponent = 0;
            
                       
            %% Clustered bars per age group, climate scenario (cc only)
            % Per 100,000 population
            %%NOT USED IN PAPER

            AvAgeResults = zeros(nClimScen, ageGroupCnt); %ignore all ages.
            %AvAgeResults_10percentile = zeros(nClimScen, ageGroupCnt);
            %AvAgeResults_90percentile = zeros(nClimScen, ageGroupCnt);
            AvAgeResults_95CI_low = zeros(nClimScen, ageGroupCnt);
            AvAgeResults_95CI_high = zeros(nClimScen, ageGroupCnt);

            AvAgeResults(1,:) = ResultsPastAgeGroup{1}(1:5,1);
            AvAgeResults(2,:) = Results_1_5C_AgeGroup_CC{1}(1:5,1);
            AvAgeResults(3,:) = Results_2C_AgeGroup_CC{1}(1:5,1);
            AvAgeResults(4,:) = Results_3C_AgeGroup_CC{1}(1:5,1);
            AvAgeResults(5,:) = Results_3C_AgeGroup_CC{1}(1:5,1);

            AvAgeResults_95CI_low(1,:) = ResultsPastAgeGroup{1}(1:5,4);
            AvAgeResults_95CI_low(2,:) = Results_1_5C_AgeGroup_CC{1}(1:5,4); 
            AvAgeResults_95CI_low(3,:) = Results_2C_AgeGroup_CC{1}(1:5,4); 
            AvAgeResults_95CI_low(4,:) = Results_3C_AgeGroup_CC{1}(1:5,4);
            AvAgeResults_95CI_low(5,:) = Results_3C_AgeGroup_CC{1}(1:5,4);

            AvAgeResults_95CI_high(1,:) = ResultsPastAgeGroup{1}(1:5,5);
            AvAgeResults_95CI_high(2,:) = Results_1_5C_AgeGroup_CC{1}(1:5,5); 
            AvAgeResults_95CI_high(3,:) = Results_2C_AgeGroup_CC{1}(1:5,5); 
            AvAgeResults_95CI_high(4,:) = Results_3C_AgeGroup_CC{1}(1:5,5); 
            AvAgeResults_95CI_high(5,:) = Results_3C_AgeGroup_CC{1}(1:5,5); 

            AvAgeResults = AvAgeResults.'; %Transpose for clustered bar
            AvAgeResults_95CI_low = AvAgeResults_95CI_low.';
            AvAgeResults_95CI_high = AvAgeResults_95CI_high.';

            M = max(AvAgeResults_95CI_high(:,5)*1.1); %set y axis max limit across sub-plots
            figure

            b = bar(AvAgeResults, 'FaceColor', 'flat');
            b(1).CData(1,:) = [1 0.843137264251709 0];
            b(2).CData(1,:) = [0.87058824300766 0.490196079015732 0];
            b(3).CData(1,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
            b(4).CData(1,:) = [0.600000023841858 0.200000002980232 0];
            b(1).CData(2,:) = [1 0.843137264251709 0];
            b(2).CData(2,:) = [0.87058824300766 0.490196079015732 0];
            b(3).CData(2,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
            b(4).CData(2,:) = [0.600000023841858 0.200000002980232 0];
            b(1).CData(3,:) = [1 0.843137264251709 0];
            b(2).CData(3,:) = [0.87058824300766 0.490196079015732 0];
            b(3).CData(3,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
            b(4).CData(3,:) = [0.600000023841858 0.200000002980232 0];
            b(1).CData(4,:) = [1 0.843137264251709 0];
            b(2).CData(4,:) = [0.87058824300766 0.490196079015732 0];
            b(3).CData(4,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
            b(4).CData(4,:) = [0.600000023841858 0.200000002980232 0];
            

            hold on
            p_errhigh = AvAgeResults_95CI_high - AvAgeResults; %90th p - average
            p_errlow = AvAgeResults - AvAgeResults_95CI_low; % average - 10th p)

            % Calculate the number of groups and number of bars in each group
            [ngroups,nbars] = size(AvAgeResults);
            % Get the x coordinate of the bars
            x = nan(nbars, ngroups);
            for i = 1:nbars
                x(i,:) = b(i).XEndPoints;
            end

            er = errorbar(x', AvAgeResults, p_errlow, p_errhigh, 'k', 'LineWidth',1,'linestyle','none');                                

            % Create ylabel
            ylabel('Average Annual Heat Related Deaths per 100,000 population');

            % Create xlabel
            xlabel('Age Group');
            %create title
            title ('Climate Change Only')

            set (gca, 'ylim', [0 M])
            set(gca,'XTick',[1 2 3 4],'XTickLabel',{'0-64', '65-74', '75-84','85+'},'YGrid',...
                'on');

            %%%
            %% Clustered bars per age group, climate scenario (cc and SE)
            % Per 100,000 population

            AvAgeResultsSE = zeros(nClimScen, ageGroupCnt); %ignore all ages.
            %AvAgeResults_10percetileSE = zeros(nClimScen, ageGroupCnt-1);
            %AvAgeResults_90percetileSE = zeros(nClimScen, ageGroupCnt-1);
            AvAgeResults_95CI_lowSE = zeros(nClimScen, ageGroupCnt);
            AvAgeResults_95CI_highSE = zeros(nClimScen, ageGroupCnt);

            AvAgeResultsSE(1,:) = ResultsPastAgeGroup{1}(:,1);
            AvAgeResultsSE(2,:) = Results_1_5C_AgeGroup_SE{1}(:,1);
            AvAgeResultsSE(3,:) = Results_2C_AgeGroup_SE{1}(:,1);
            AvAgeResultsSE(4,:) = Results_3C_AgeGroup_SE{1}(:,1);
            AvAgeResultsSE(5,:) = Results_4C_AgeGroup_SE{1}(:,1);

            AvAgeResults_95CI_lowSE(1,:) = ResultsPastAgeGroup{1}(:,4); 
            AvAgeResults_95CI_lowSE(2,:) = Results_1_5C_AgeGroup_SE{1}(:,4); 
            AvAgeResults_95CI_lowSE(3,:) = Results_2C_AgeGroup_SE{1}(:,4); 
            AvAgeResults_95CI_lowSE(4,:) = Results_3C_AgeGroup_SE{1}(:,4);
            AvAgeResults_95CI_lowSE(5,:) = Results_4C_AgeGroup_SE{1}(:,4);

            AvAgeResults_95CI_highSE(1,:) = ResultsPastAgeGroup{1}(:,5); 
            AvAgeResults_95CI_highSE(2,:) = Results_1_5C_AgeGroup_SE{1}(:,5); 
            AvAgeResults_95CI_highSE(3,:) = Results_2C_AgeGroup_SE{1}(:,5); 
            AvAgeResults_95CI_highSE(4,:) = Results_3C_AgeGroup_SE{1}(:,5);
            AvAgeResults_95CI_highSE(5,:) = Results_4C_AgeGroup_SE{1}(:,5);
            
            AvAgeResultsSE = AvAgeResultsSE.'; %Transpose for clustered bar
            AvAgeResults_95CI_lowSE = AvAgeResults_95CI_lowSE.';
            AvAgeResults_95CI_highSE = AvAgeResults_95CI_highSE.';

            M = max(AvAgeResults_95CI_highSE(:,5)*1.1); %set y axis max limit across sub-plots

            figure
            b = bar(AvAgeResultsSE, 'FaceColor', 'flat');
            b(1).CData(1,:) = [1 0.843137264251709 0];
            b(2).CData(1,:) = [0.87058824300766 0.490196079015732 0];
            b(3).CData(1,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
            b(4).CData(1,:) = [0.600000023841858 0.200000002980232 0];
            b(5).CData(1,:) = [0.40,0.00,0.00];
            b(1).CData(2,:) = [1 0.843137264251709 0];
            b(2).CData(2,:) = [0.87058824300766 0.490196079015732 0];
            b(3).CData(2,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
            b(4).CData(2,:) = [0.600000023841858 0.200000002980232 0];
            b(5).CData(2,:) = [0.40,0.00,0.00];
            b(1).CData(3,:) = [1 0.843137264251709 0];
            b(2).CData(3,:) = [0.87058824300766 0.490196079015732 0];
            b(3).CData(3,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
            b(4).CData(3,:) = [0.600000023841858 0.200000002980232 0];
            b(5).CData(3,:) = [0.40,0.00,0.00];
            b(1).CData(4,:) = [1 0.843137264251709 0];
            b(2).CData(4,:) = [0.87058824300766 0.490196079015732 0];
            b(3).CData(4,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
            b(4).CData(4,:) = [0.600000023841858 0.200000002980232 0];
            b(5).CData(4,:) = [0.40,0.00,0.00];
            b(1).CData(5,:) = [1 0.843137264251709 0];
            b(2).CData(5,:) = [0.87058824300766 0.490196079015732 0];
            b(3).CData(5,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
            b(4).CData(5,:) = [0.600000023841858 0.200000002980232 0];
            b(5).CData(5,:) = [0.40,0.00,0.00];
            
            hold on
            
            p_errhigh = AvAgeResults_95CI_highSE - AvAgeResultsSE; %90th p - average
            p_errlow = AvAgeResultsSE - AvAgeResults_95CI_lowSE; % average - 10th p)

            % Calculate the number of groups and number of bars in each group
            [ngroups,nbars] = size(AvAgeResultsSE);
            % Get the x coordinate of the bars
            x = nan(nbars, ngroups);
            for i = 1:nbars
                x(i,:) = b(i).XEndPoints;
            end

            er = errorbar(x', AvAgeResultsSE, p_errlow, p_errhigh, 'k', 'LineWidth',1,'linestyle','none');                                

            % Create ylabel
            ylabel('Average Annual Heat Related Deaths per 100,000 population');


            % Create xlabel
            xlabel('Age Group');
            %create title
            title ('Climate and Population Change')


            set (gca, 'ylim', [0 M])
            set(gca,'XTick',[1 2 3 4 5],'XTickLabel',{'All Ages', '0-64', '65-74', '75-84','85+'},'YGrid',...
                'on');

            %% clustered bar - mortality per each 1 degree increment for each scenario
            load('avMortPerIncrementUK')

            figure
            subplot(2,1,1);

            avMortPerIncrementUK_0dp = round(avMortPerIncrementUK); % round to nearest integer
            M = 3250; %max(avMortPerIncrementUK_CI_high(10,:)*1.1); %set y axis max limit across sub-plots
            p_errlow = round(avMortPerIncrementUK_CI_low([2,7,8,9,10], 1:16).');
            p_errhigh = round(avMortPerIncrementUK_CI_high([2,7,8,9,10], 1:16).');

            b = bar(avMortPerIncrementUK_0dp([2:3, 5, 7, 9],1:16).');
            set(b(1),'FaceColor',[1 0.843137264251709 0]);
            set(b(2),'FaceColor',[0.87058824300766 0.490196079015732 0]);
            set(b(3),'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);
            set(b(4),'FaceColor',[0.600000023841858 0.200000002980232 0]);
            set(b(5),'FaceColor',[0.40,0.00,0.00]);


            % Create ylabel
            ylabel('Average Annual Heat Related Deaths');

            % Create xlabel
            xlabel('Degrees warming above regional thresholds');
            %create title
            title ('Climate Change Only')
            set (gca, 'ylim', [0 M])


            subplot(2,1,2);
            b = bar(avMortPerIncrementUK_0dp([2,4, 6, 8, 10],1:16).');
            set(b(1),'FaceColor',[1 0.843137264251709 0]);
            set(b(2),'FaceColor',[0.87058824300766 0.490196079015732 0]);
            set(b(3),'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);
            set(b(4),'FaceColor',[0.600000023841858 0.200000002980232 0]);
            set(b(5),'FaceColor',[0.40,0.00,0.00]);
            hold on
            labels = {'Past','1.5','2.0','3.0', '4.0'};
            legend(labels);
                    
            % Calculate the number of groups and number of bars in each group
            [ngroups,nbars] = size(p_errlow);
            % Get the x coordinate of the bars
            x = nan(nbars, ngroups);
            for i = 1:nbars
                x(i,:) = b(i).XEndPoints;
            end
            
            er = errorbar(x.', avMortPerIncrementUK_0dp([2,4, 6, 8, 10],1:16).', p_errlow, p_errhigh, 'k', 'LineWidth',1,'linestyle','none'); 

            % Create ylabel
            ylabel('Average Annual Heat Related Deaths');
            xlabel('Degrees warming above regional thresholds');
            %create title
            title ('Climate and Population Change')
            set (gca, 'ylim', [0 M])

            %% Plots for Adaptation (scenarios for natural acclimatisation - 93rd Percentile or 1 & 2 degree threshold)
            figure;
            subplot(1,2,1);
            b= bar (AvAnnualDeathsAdapt.', 'stacked', 'FaceColor','flat');

            %colours for adapt thresholds no adaptation and adapScen.
            b(2).CData(1,:) = [0.87058824300766 0.490196079015732 0];
            b(1).CData(1,:) = [0.800000011920929 0.800000011920929 0.800000011920929];
            b(1).CData(2,:) = [0.87058824300766 0.490196079015732 0];
            b(2).CData(2,:) = [0.800000011920929 0.800000011920929 0.800000011920929];
            b(1).CData(3,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
            b(2).CData(3,:) = [0.800000011920929 0.800000011920929 0.800000011920929];
            b(1).CData(4,:) = [0.600000023841858 0.200000002980232 0];
            b(2).CData(4,:) = [0.800000011920929 0.800000011920929 0.800000011920929];
            b(1).CData(5,:) = [0.600000023841858 0.200000002980232 0];
            b(2).CData(5,:) = [0.800000011920929 0.800000011920929 0.800000011920929];
            hold on

            M = AvAnnualDeathsUK_SE(1,13)*1.3; %set y axis max limit across sub-plots
            % Create ylabel
            ylabel('Average Annual Heat Related Deaths');

            % Create xlabel
            xlabel('Climate Scenario');

            % Create title
            title('Climate Change Only');

             labels = {'No Adaptation','Adaptation'};
             legend(labels);
             set(legend,...
                 'Position',[0.157187955560562 0.762734641133266 0.276785708750997 0.13589211247274]);

            % Set the remaining axes properties
            set (gca, 'ylim', [0 M])
            set(gca,'XTick',[1 2 3 4 5],'XTickLabel',{'Past','1.5','2.0','3.0', '4.0'},'YGrid',...
                'on');

            subplot(1,2,2);
            b=bar (AvAnnualDeathsAdapt_SE.', 'stacked', 'FaceColor','flat');     

            b(2).CData(1,:) = [0.87058824300766 0.490196079015732 0];
            b(1).CData(1,:) = [0.800000011920929 0.800000011920929 0.800000011920929];
            b(1).CData(2,:) = [0.87058824300766 0.490196079015732 0];
            b(2).CData(2,:) = [0.800000011920929 0.800000011920929 0.800000011920929];
            b(1).CData(3,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
            b(2).CData(3,:) = [0.800000011920929 0.800000011920929 0.800000011920929];
            b(1).CData(4,:) = [0.600000023841858 0.200000002980232 0];
            b(2).CData(4,:) = [0.800000011920929 0.800000011920929 0.800000011920929];
            b(1).CData(5,:) = [0.600000023841858 0.200000002980232 0];
            b(2).CData(5,:) = [0.800000011920929 0.800000011920929 0.800000011920929];
            hold on

            % Create ylabel
            ylabel('Average Annual Heat Related Deaths');

            % Create xlabel
            xlabel('Climate Scenario');

            % Create title
            title('Climate and Population Change');

            % Set the remaining axes properties
            set (gca, 'ylim', [0 M])
            set(gca,'XTick',[1 2 3 4 5],'XTickLabel',{'Past','1.5','2.0','3.0', '4.0'},'YGrid',...
                'on');
% 
%             % adaptation comparison
%             AvAdaptResultsSE = AvAdaptResultsSE.'; %Transpose for clustered bar
%             AvAdaptResults_10percentileSE = AvAdaptResults_10percentileSE.';
%             AvAdaptResults_90percentileSE = AvAdaptResults_90percentileSE.';
% 
%             M = 22000; %set y axis max limit across sub-plots
% 
%             figure
%             b = bar(AvAdaptResultsSE, 'FaceColor', 'flat');
%             b(1).CData(1,:) = [1 0.843137264251709 0];
%             b(2).CData(1,:) = [1 0.843137264251709 0];
%             b(3).CData(1,:) = [1 0.843137264251709 0];
%             b(4).CData(1,:) = [1 0.843137264251709 0];
%             b(1).CData(2,:) = [0.87058824300766 0.490196079015732 0];
%             b(2).CData(2,:) = [0.87058824300766 0.490196079015732 0];
%             b(3).CData(2,:) = [0.87058824300766 0.490196079015732 0];
%             b(4).CData(2,:) = [0.87058824300766 0.490196079015732 0];
%             b(1).CData(3,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
%             b(2).CData(3,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
%             b(3).CData(3,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
%             b(4).CData(3,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
%             b(1).CData(4,:) = [0.600000023841858 0.200000002980232 0];
%             b(2).CData(4,:) = [0.600000023841858 0.200000002980232 0];
%             b(3).CData(4,:) = [0.600000023841858 0.200000002980232 0];
%             b(4).CData(4,:) = [0.600000023841858 0.200000002980232 0];
%             b(1).CData(5,:) = [0.40,0.00,0.00];
%             b(2).CData(5,:) = [0.40,0.00,0.00];
%             b(3).CData(5,:) = [0.40,0.00,0.00];
%             b(4).CData(5,:) = [0.40,0.00,0.00];
% 
%             hold on
%            % p_errhigh = AvAdaptResults_90percentileSE - AvAdaptResultsSE; %90th p - average
%            % p_errlow = AvAdaptResultsSE - AvAdaptResults_10percentileSE; % average - 10th p)
% 
%             % Calculate the number of groups and number of bars in each group
%             [ngroups,nbars] = size(AvAdaptResultsSE);
%             % Get the x coordinate of the bars
%             x = nan(nbars, ngroups);
%             for i = 1:nbars
%                 x(i,:) = b(i).XEndPoints;
%             end
% 
%             er = errorbar(x', AvAdaptResultsSE,  AvAdaptResults_10percentileSE,  AvAdaptResults_90percentileSE, 'k', 'LineWidth',1,'linestyle','none');                                
% 
%             % Create ylabel
%             ylabel('Average Annual Heat Related Deaths');
% 
% 
%             % Create xlabel
%             xlabel('Climate Scenario');
%             %create title
%             title ('Climate and Population Change')
% 
%             labels = {'No Adaptation', '1°C', '2°C','93rd Percentile'};
%             legend(labels);
% 
%             set (gca, 'ylim', [0 M])
%             set(gca,'XTick',[1 2 3 4 5],'XTickLabel',{'Past','1.5','2.0','3.0', '4.0'},'YGrid',...
%                 'on');
%             ax = gca;
%             ax.YAxis.Exponent = 0;
    end
elseif impactMetric == 2
    run HARM_labour_productivity.m
elseif impactMetric == 3
     run HARM_residential_discomfort.m
end

toc  %stopwatch to measure performance
disp ('End')