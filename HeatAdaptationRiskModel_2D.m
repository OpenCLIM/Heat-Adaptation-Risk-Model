function [] = HeatAdaptationRiskModel_2D(InputDir,varname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   HEAT-HARM model (Heat - Adaptation - Risk - Model)
%
%   Estimates heat related mortality under future projections of climate change.
%   TO ADD: Residential dicomfort and reduction in labour productivity
%
%   Updated from ARCADIA Impacts model developed in ARCADIA project (ECI, University of Oxford)
%   for OpenCLIM project (Tyndall Centre for Climate Change Research, UEA).
%
%   Uses post-processed UKCP18 data provided from HEAT model (Alan Kennedy-Asser, University of Bristol)
%
%   Author: Katie Jenkins    UEA    (Edits by Alan Kennedy-Asser, UoB)
%   Date:12/06/2023 [version_3.1]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('*******************************************************************')
disp('Heat Adaptation Risk Model (HARM)')
disp('version 3.1')
disp('This version currently estimates UK heat-related mortality')
disp('Version for use on DAFNI, March 2023')
disp('*******************************************************************')

% Start timer
startt = now;

format short e

%% Set some defaults
% If the model has been passed no variable mapping to input directory,
% figure out if it is on DAFNI or running locally and update accordingly
if ~exist('InputDir','var')
%     % For testing purposes, to show the correct data has copied to the Docker
%     % container (remove this later):
%     disp(' ')
%     disp('You are here:')
%     pwd
     if strcmp(pwd,'/code')
%         disp('Running in Docker container with these files:')
%         ls
%         disp(' ')
%         disp('-----')
%         
%         disp('These are data files:')
%         ls('/data/')
%         disp('These are input files:')
%         ls('/data/inputs/')
%         disp('These are climate files:')
%         ls('/data/inputs/Climate/')        
%         
%         disp(' ')
%         disp('-----')
        InputDir = '/data/inputs/Climate/';
        OutputDir = '/data/outputs/';
        mkdir(OutputDir)
        
    else
        % Update these as necessary to a sensible local default
        InputDir = 'inputs/Climate/';
        OutputDir = 'outputs/';
        
    end
end

if ~exist('OutputDir','var')
    OutputDir = 'outputs/';
end


% Set default variable name if one has not been provided
if ~exist('varname','var')
    % Take the environment variable defined name if possible
    varname = char(getenv('VARNAME'));
    if isempty(varname)
        % Otherwise set tas as default
        varname = 'tas';
    end
end


%% Load input data in 2D
% Check for input files
% Check if T data has been passed to the model
files = dir([InputDir '*.nc']); % Then check if any files exist with this root
if isempty(files)
    disp('No climate input data. Stopping.') % Cancel if not
    return
end

% Load first file as a template
file = [files(1).folder,'/',files(1).name];
T = ncread(file,varname);
S = size(T);


% Pass useful input data
if ~exist([InputDir,'LSM.mat'],'file')
    if S(1) == 82 && S(2) == 112
        load('PreProcessedData/LSM12.mat') % Default = UKCP18 RCM LSM
        LSM = LSM12;
        load('PreProcessedData/lat_UK_RCM.mat')
        lats = lat_UK_RCM;
        load('PreProcessedData/long_UK_RCM.mat')
        lons = long_UK_RCM;
        CPM = 0;
    elseif S(1) == 484 && S(2) == 606
        load('PreProcessedData/LSM2.mat') % Default = UKCP18 RCM LSM
        LSM = LSM2;
        load('PreProcessedData/lat_UK_CPM.mat')
        lats = lat_UK_CPM;
        load('PreProcessedData/long_UK_CPM.mat')
        lons = long_UK_CPM;
        CPM = 1;
    end
else
    load([InputDir,'LSM.mat']);
    load([InputDir,'lats.mat']);
    load([InputDir,'lons.mat']);
end

if exist([InputDir,'grid_idx.mat'],'file')
    load([InputDir,'grid_idx.mat'])
    load([InputDir,'grid_idy.mat'])
else
    disp('Grid data is missing. Assuming UKCP18 RCM domain.')
    grid_idx = 1:82;
    grid_idy = 1:112;
end


try
    projection_x_coordinate = ncread(file,'projection_x_coordinate');
    projection_y_coordinate = ncread(file,'projection_y_coordinate');
catch
    disp('Projection of input data unknown, setting generic values')
    projection_x_coordinate = grid_idx;
    projection_y_coordinate = grid_idy;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User defined parameters from environment variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Heat/cold related mortality

% Read in environment variable parameters
impactMetric = getenv('IMPMET');

% impactMetric = 'ColdMort' % For testing
if strcmp(impactMetric,'HeatMort') || strcmp(impactMetric,'ColdMort')
    
    % Load defaults
    if CPM == 0
        load('PreProcessedData/population_SSP2_2D.mat')
        load('PreProcessedData/population_SSP4_2D.mat')
        load('PreProcessedData/population_SSP5_2D.mat')
        load('PreProcessedData/dailyDeathRate_2D.mat')
        
        if strcmp(impactMetric,'HeatMort')
            disp('Calculating heat-related mortality')
            disp(' ')
            load('PreProcessedData/Mortality_RR_2D.mat')
            load('PreProcessedData/Mortality_MMT_2D.mat')
        elseif strcmp(impactMetric,'ColdMort')
            disp('Calculating cold-related mortality')
            disp(' ')
            load('PreProcessedData/Mortality_RR_2D_cold.mat')
            load('PreProcessedData/Mortality_MMT_2D_cold.mat')
            Mortality_RR_2D = Mortality_RR_2D_cold;
            Mortality_MMT_2D = Mortality_MMT_2D_cold;
        end
    elseif CPM == 1
        load('PreProcessedData/population_SSP2_2D_CPM.mat')
        load('PreProcessedData/population_SSP4_2D_CPM.mat')
        load('PreProcessedData/population_SSP5_2D_CPM.mat')
        load('PreProcessedData/dailyDeathRate_2D_CPM.mat')
        population_SSP2_2D = population_SSP2_2D_CPM;
        population_SSP4_2D = population_SSP4_2D_CPM;
        population_SSP5_2D = population_SSP5_2D_CPM;
        dailyDeathRate_2D = dailyDeathRate_2D_CPM;
        
        if strcmp(impactMetric,'HeatMort')
            disp('Calculating heat-related mortality')
            disp(' ')
            load('PreProcessedData/Mortality_RR_2D_CPM.mat')
            load('PreProcessedData/Mortality_MMT_2D_CPM.mat')
            Mortality_RR_2D = Mortality_RR_2D_CPM;
            Mortality_MMT_2D = Mortality_MMT_2D_CPM;
        elseif strcmp(impactMetric,'ColdMort')
            disp('Calculating cold-related mortality')
            disp(' ')
            load('PreProcessedData/Mortality_RR_2D_CPM_cold.mat')
            load('PreProcessedData/Mortality_MMT_2D_CPM_cold.mat')
            Mortality_RR_2D = Mortality_RR_2D_CPM_cold;
            Mortality_MMT_2D = Mortality_MMT_2D_CPM_cold;
        end
    end
    
    % If user has provided a MMT and RR field, load these instead
    if exist('/data/inputs/Mortality/MMT.mat','file') && exist('/data/inputs/Mortality/RR.mat','file')
        disp('Using alternative MMT and RR data: Loading')
        disp(' ')
        load('/data/inputs/Mortality/MMT.csv')
        Mortality_MMT_2D = MMT;
        load('/data/inputs/Mortality/RR.csv')
        Mortality_RR_2D = RR;
    end
    % Rename
    Mortality_RR = Mortality_RR_2D;
    
    % Choose whether to run climate change only (1) or climate + socioeconomic
    % changes (2): (default = climate change only)
    socioEcScen = getenv('SOCIOECON');

    
    % If interested in socioeconomics, choose scenario and time period:
    if strcmp(socioEcScen,'CConly')
        disp('Calculating effect of climate change only (population remains unchanged from 2011 baseline)')
        disp(' ')
        population = squeeze(population_SSP5_2D(:,:,:,1));
        % Find total population for each age category and time period
        population_aggregate_SSP5 = squeeze(nansum(nansum(population_SSP5_2D,1),2));
        population_aggregate = population_aggregate_SSP5(:,1);
        
    elseif strcmp(socioEcScen,'SSP2-2050')
        disp('Using default population from UK-SSPs for 2050 (SSP2)')
        disp(' ')
        popID = 4;
        population = squeeze(population_SSP2_2D(:,:,:,popID));
        % Find total population for each age category and time period
        population_aggregate_SSP2 = squeeze(nansum(nansum(population_SSP2_2D,1),2));
        population_aggregate = population_aggregate_SSP2(:,popID);
        
    elseif strcmp(socioEcScen,'SSP2-2080')
        disp('Using default population from UK-SSPs for 2080 (SSP2)')
        disp(' ')
        popID = 5;
        population = squeeze(population_SSP2_2D(:,:,:,popID));
        % Find total population for each age category and time period
        population_aggregate_SSP2 = squeeze(nansum(nansum(population_SSP2_2D,1),2));
        population_aggregate = population_aggregate_SSP2(:,popID);
                
    elseif strcmp(socioEcScen,'SSP4-2050')
        disp('Using default population from UK-SSPs for 2050 (SSP4)')
        disp(' ')
        popID = 4;
        population = squeeze(population_SSP4_2D(:,:,:,popID));
        % Find total population for each age category and time period
        population_aggregate_SSP4 = squeeze(nansum(nansum(population_SSP4_2D,1),2));
        population_aggregate = population_aggregate_SSP4(:,popID);
        
    elseif strcmp(socioEcScen,'SSP4-2080')
        disp('Using default population from UK-SSPs for 2080 (SSP4)')
        disp(' ')
        popID = 5;
        population = squeeze(population_SSP4_2D(:,:,:,popID));
        % Find total population for each age category and time period
        population_aggregate_SSP4 = squeeze(nansum(nansum(population_SSP4_2D,1),2));
        population_aggregate = population_aggregate_SSP4(:,popID);
                
    elseif strcmp(socioEcScen,'SSP5-2050')
        disp('Using default population from UK-SSPs for 2050 (SSP5)')
        disp(' ')
        popID = 4;
        population = squeeze(population_SSP5_2D(:,:,:,popID));
        % Find total population for each age category and time period
        population_aggregate_SSP5 = squeeze(nansum(nansum(population_SSP5_2D,1),2));
        population_aggregate = population_aggregate_SSP5(:,popID);
        
    elseif strcmp(socioEcScen,'SSP5-2080')
        disp('Using default population from UK-SSPs for 2080 (SSP5)')
        disp(' ')
        popID = 5;
        population = squeeze(population_SSP5_2D(:,:,:,popID));
        % Find total population for each age category and time period
        population_aggregate_SSP5 = squeeze(nansum(nansum(population_SSP5_2D,1),2));
        population_aggregate = population_aggregate_SSP5(:,popID);        
        
    elseif strcmp(socioEcScen,'UDM')
        disp('Using population data from UDM: Loading')
        disp(' ')
        filespop = dir(['/data/inputs/Population/population_total_uk_*.tif']); % Then check if any files exist with this root
        
        if exist('/data/inputs/Population/population_total_uk-12km.asc','file')
            population = rot90(arcgridread('population_total_uk-12km.asc'),3);
            population = cat(3,population,rot90(arcgridread('/data/inputs/Population/population_demographics_0-64.asc'),3));
            population = cat(3,population,rot90(arcgridread('/data/inputs/Population/population_demographics_65-74.asc'),3));
            population = cat(3,population,rot90(arcgridread('/data/inputs/Population/population_demographics_75-84.asc'),3));
            population = cat(3,population,rot90(arcgridread('/data/inputs/Population/population_demographics_85.asc'),3));
            population_aggregate = squeeze(nansum(nansum(population,1),2));
            
        elseif ~isempty(filespop)        
            % Load first file as a template
            filepop = [filespop(1).folder,'/',filespop(1).name];
            
            population = rot90(geotiffread(filepop),3);
            population = cat(3,population,rot90(geotiffread('/data/inputs/Population/population_total_demographic_0_64.tif'),3));
            population = cat(3,population,rot90(geotiffread('/data/inputs/Population/population_total_demographic_65_74.tif'),3));
            population = cat(3,population,rot90(geotiffread('/data/inputs/Population/population_total_demographic_75_84.tif'),3));
            population = cat(3,population,rot90(geotiffread('/data/inputs/Population/population_total_demographic_85.tif'),3));
            population_aggregate = squeeze(nansum(nansum(population,1),2));
            
        else
            disp('Invlaid UDM data: cancelling')
            disp('-----')
            return
        end
        
        if length(population(:,1,1)) == 55 && length(population(1,:,1)) == 100
            disp('UDM data provided on standard UDM grid (55 x 100) -> convert to 112 x 82')
            populationtemp = population;
            population = nan(82,112,5);
            population(19:73,11:110,:) = populationtemp;
            disp(' ')
        end
        
    elseif strcmp(socioEcScen,'Other')
        disp('Using alternative population data: Loading')
        % If user has provided a population and daily death rate field, load these instead
        if exist('/data/inputs/Population/population.mat','file') && exist('/data/inputs/Population/deathrate.mat','file')
            load('/data/inputs/Population/deathrate.mat');
            load('/data/inputs/Population/population.mat');
        else
            disp('Invlaid population data (ensure files are named population.csv and deathrate.csv): cancelling')
            disp('-----')
        end
    end
    
    % Choose how much adaptation is possible, either by a numeric value which
    % is added to the MMT or by the string 'acclim', which loads the variable
    % reg_acclim.mat computed by HEAT: (default = 0 / no acclimatisation)
    if strcmp(impactMetric,'HeatMort') % Acclimatisation only applicable for heat-related mortality
        
        acclimaScen = getenv('ACCLIMSCEN');
        if isempty(acclimaScen)
            acclimaScen = 'None'; %Select approach to behavioural adaptation via natural acclimatisation. 1 = No Adaptation; 2 & 3 = pre-defined thresholds 1 and 2 degrees respectively; 4 = 93rd P of TMean (spatially and temporally explicit based on output from HEAT)
        end
        
        if strcmp(acclimaScen,'None')
            disp('Assuming no acclimatisation or adaptation')
            disp(' ')
            adaptationName = 'No Adaptation';
            acc = 0;
            RRchange = 0;
            
        elseif strcmp(acclimaScen,'93rdMMTRR') || strcmp(acclimaScen,'93rdMMT')
            % Load the updated acclimatisation
            if exist([InputDir,'reg_acclim_2D.mat'],'file')
                adaptationName = 'Acclimatisation to percentile shift';
                disp('Acclimatisation file available: loading this')
                load([InputDir,'reg_acclim_2D.mat']);
                
                disp('Updating MMT')
                acc = reg_acclim_2D;
                
                if strcmp(acclimaScen,'93rdMMTRR')
                    disp('Updating RR')
                    RRchange = acc*1.07; % Value 1.07 from regression of regional RR vs MMT
                else
                    RRchange = 0;
                end
                disp(' ')
            else
                disp('Acclimatisation file is missing. Defaulting to 0.')
                disp(' ')
                acc = 0;
                RRchange = 0;
            end
            
        elseif strcmp(acclimaScen,'Custom')
            adaptationName = 'Custom acclimatisation of MMT and RR';
            disp('Custom acclimatisation defined')
            disp(' ')
            acc = str2double(string(getenv('MMT')));
            if isempty(acc)
                acc = 0;
            end
            RRchange = str2double(string(getenv('RR')));
            if isempty(RRchange)
                RRchange = 0;
            end
        end
        
        % Update MMT and RR
        Mortality_MMT = Mortality_MMT_2D + acc;
        Mortality_RR = Mortality_RR + RRchange;
        
    elseif strcmp(impactMetric,'ColdMort')
        
        acclimaScen = getenv('ACCLIMSCEN');
        if isempty(acclimaScen)
            acclimaScen = 'None'; %Select approach to behavioural adaptation via natural acclimatisation. 1 = No Adaptation; 2 & 3 = pre-defined thresholds 1 and 2 degrees respectively; 4 = 93rd P of TMean (spatially and temporally explicit based on output from HEAT)
        end
        
        if strcmp(acclimaScen,'Custom')
            adaptationName = 'Custom acclimatisation of MMT and RR';
            disp('Custom acclimatisation defined')
            disp(' ')
            acc = str2double(string(getenv('MMT')));
            if isempty(acc)
                acc = 0;
            end
            
            % Update MMT and RR
            Mortality_MMT = Mortality_MMT_2D + acc;
            
        else
            Mortality_MMT = Mortality_MMT_2D;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% Go through each dataset (model simulation) and calculate mortality
    % Based on tempalte file, create template output arrays
    totMortPerIncrement = nan(S(1),S(2),5,30); % long x lat x ages x increments
    avMortPerIncrement = nan(S(1),S(2),5,30); % long x lat x ages x increments
    avMortality = nan(S(1),S(2),5); % long x lat x ages 
    
    disp(['Loading dataset ',files(1).name])
    T = ncread(file,varname);
    
    % Find length of record in years
    dates = ncread(file,'yyyymmdd');
    
    nyears = str2double(dates(end,1:4)) - str2double(dates(1,1:4)) +1;
    
    % Find when T exceeds MMT
    if strcmp(impactMetric,'HeatMort')
        T_ex = T - Mortality_MMT;
    elseif strcmp(impactMetric,'ColdMort')
        % T should be an average over a 28 day period - first 27 days need
        % a shorter averaging window
        T_temp = T;
        for Td = 2:length(T(1,1,:))
            if Td > 27
                T(:,:,Td) = nanmean(T_temp(:,:,Td-27:Td),3);
            else
                T(:,:,Td) = nanmean(T_temp(:,:,1:Td),3);
            end
        end
        T_ex = Mortality_MMT - T;
    end
    
    % Otherwise make 0
    T_ex(T_ex<0) = 0;
    
    % Go through each age category
    ages = {'all ages','0-64','65-74','75-84','85+'};
    for a = 1:5
        disp(['Calculating for age category: ',char(ages(a))])
        % Do the mortality calculation
        dailyMort = (((T_ex) .* Mortality_RR(:,:,a))/100).*dailyDeathRate_2D(:,:,a).*population(:,:,a);
        
        % Round T exceeding threshold to nearest degree
        incrementExceeded = round(T_ex);
        
        % Find how many increments are needed
        increments = nanmax(nanmax(nanmax(nanmax(incrementExceeded))));
        
        % Do calculation of mortality per increment
        for i = 1:increments
            totMortPerIncrement(:,:,a,i) = nansum(dailyMort.*(incrementExceeded == i),3);
        end
        
        % Calculate averages and percentiles
        avMortPerIncrement(:,:,a,:) = totMortPerIncrement(:,:,a,:)/nyears;
        avMortality(:,:,a) = (sum(dailyMort,3)/nyears);
        
    end
    
    disp(' ')
    
    
    %% Saving outputs
    disp('Saving outputs')
    disp(' ')
    
    % For reproducing Katie's ERL figures
    Fig3vals = squeeze(nansum(nansum(avMortPerIncrement(:,:,1,:))));
    Fig1_2vals = squeeze(nansum(nansum(avMortality(:,:,1))));
    Fig4vals = squeeze(nansum(nansum(avMortality(:,:,:),1),2))./population_aggregate*100000;
    Fig5vals = squeeze(avMortality(:,:,1));
    
    dlmwrite([OutputDir,'F1_2-',impactMetric,'.csv'],Fig1_2vals, 'delimiter', ',', 'precision', '%i')
    dlmwrite([OutputDir,'F3-',impactMetric,'.csv'],Fig3vals, 'delimiter', ',', 'precision', '%i')
    dlmwrite([OutputDir,'F4-',impactMetric,'.csv'],Fig4vals, 'delimiter', ',', 'precision', '%i')
    dlmwrite([OutputDir,'F5-',impactMetric,'.csv'],Fig5vals, 'delimiter', ',', 'precision', '%i')
    
    
    % Set netCDF details
    units = 'deaths';
    label_units = 'deaths';
    plot_label = 'Annual average deaths';
    if strcmp(impactMetric,'HeatMort')
        standard_name = 'average heat mortality per increment';
        long_name = 'Average heat-related mortality per year, for each age category and degree increment above MMT';
        description = 'Average heat-related mortality per year, for each age category and degree increment above MMT';
    elseif strcmp(impactMetric,'ColdMort')
        standard_name = 'average cold mortality per increment';
        long_name = 'Average cold-related mortality per year, for each age category and degree increment above MMT';
        description = 'Average cold-related mortality per year, for each age category and degree increment above MMT';
        
    end
    
    ages_char = ['00-99';'00-64';'65-74';'75-84';'85-99']';
    
    % Write netCDF
    fname_long = [OutputDir,'HARM_',impactMetric,'_',dates(1,1:8),'-',dates(end,1:8),'.nc'];
    if exist(fname_long,'file')
        delete(fname_long)
    end
    
    nccreate(fname_long,standard_name,'Dimensions',{'projection_x_coordinate',length(projection_x_coordinate),'projection_y_coordinate',length(projection_y_coordinate),'ages',5,'increments',30},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,standard_name,avMortPerIncrement);
    ncwriteatt(fname_long,standard_name,'standard_name',standard_name);
    ncwriteatt(fname_long,standard_name,'long_name',long_name);
    ncwriteatt(fname_long,standard_name,'units',units);
    ncwriteatt(fname_long,standard_name,'description',description);
    ncwriteatt(fname_long,standard_name,'label_units',label_units);
    ncwriteatt(fname_long,standard_name,'plot_label',plot_label);
    ncwriteatt(fname_long,standard_name,'missing_value',-9999);
    
    % Add lat and long data
    nccreate(fname_long,'projection_x_coordinate','Dimensions',{'projection_x_coordinate',length(projection_x_coordinate)},'Datatype','single','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,'projection_x_coordinate',projection_x_coordinate);
    ncwriteatt(fname_long,'projection_x_coordinate','axis','X');
    
    nccreate(fname_long,'projection_y_coordinate','Dimensions',{'projection_y_coordinate',length(projection_y_coordinate)},'Datatype','single','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,'projection_y_coordinate',projection_y_coordinate);
    ncwriteatt(fname_long,'projection_y_coordinate','axis','Y');
    
    nccreate(fname_long,'ages','Dimensions',{'ages',5,'string64',64},'Datatype','char')
    ncwrite(fname_long,'ages',ages_char);
    ncwriteatt(fname_long,'ages','standard_name','age_bands');
    ncwriteatt(fname_long,'ages','axis','Z');
    
    nccreate(fname_long,'increments','Dimensions',{'increments',30},'Datatype','single','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,'increments',1:30);
    ncwriteatt(fname_long,'increments','long_name','degrees C above MMT');
    ncwriteatt(fname_long,'increments','units','\B0C');
    
    % Write some general attributes
    ncwriteatt(fname_long,'/','collection','HARM derived variable')
    ncwriteatt(fname_long,'/','creation_date',datestr(now))
    ncwriteatt(fname_long,'/','title','Heat impact metric')
    ncwriteatt(fname_long,'/','version','HARM v2.0')
    
    % Add meta data if same format as UKCP18 data
    if length(totMortPerIncrement(:,1,1,1)) == 82 && length(totMortPerIncrement(1,:,1,1)) == 112
        ncwriteatt(fname_long,standard_name,'grid_mapping','transverse_mercator');
        
        ncwriteatt(fname_long,'projection_x_coordinate','standard_name','easting');
        ncwriteatt(fname_long,'projection_x_coordinate','long_name','easting');
        ncwriteatt(fname_long,'projection_x_coordinate','units','m');
        
        ncwriteatt(fname_long,'projection_y_coordinate','standard_name','northing');
        ncwriteatt(fname_long,'projection_y_coordinate','long_name','northing');
        ncwriteatt(fname_long,'projection_y_coordinate','units','m');
    end
    
    
    
    if strcmp(impactMetric,'HeatMort')
        standard_name = 'average heat mortality';
        long_name = 'Average annual heat-related mortality per age category';
        description = 'Average annual heat-related mortality per age category';
    elseif strcmp(impactMetric,'ColdMort')
        standard_name = 'average cold mortality';
        long_name = 'Average annual cold-related mortality per age category';
        description = 'Average annual cold-related mortality per age category';
    end
    
    nccreate(fname_long,standard_name,'Dimensions',{'projection_x_coordinate',length(projection_x_coordinate),'projection_y_coordinate',length(projection_y_coordinate),'ages',5},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,standard_name,avMortality);
    ncwriteatt(fname_long,standard_name,'standard_name',standard_name);
    ncwriteatt(fname_long,standard_name,'long_name',long_name);
    ncwriteatt(fname_long,standard_name,'units',units);
    ncwriteatt(fname_long,standard_name,'description',description);
    ncwriteatt(fname_long,standard_name,'label_units',label_units);
    ncwriteatt(fname_long,standard_name,'plot_label',plot_label);
    ncwriteatt(fname_long,standard_name,'missing_value',-9999);
    %
    nccreate(fname_long,'Population','Dimensions',{'projection_x_coordinate',length(projection_x_coordinate),'projection_y_coordinate',length(projection_y_coordinate),'ages',5},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,'Population',population);
    ncwriteatt(fname_long,'Population','standard_name','Population');
    ncwriteatt(fname_long,'Population','long_name','Population in each age category');
    ncwriteatt(fname_long,'Population','units','People');
    ncwriteatt(fname_long,'Population','description','Population in each age category');
    ncwriteatt(fname_long,'Population','label_units','Population in age category');
    ncwriteatt(fname_long,'Population','plot_label','Population');
    ncwriteatt(fname_long,'Population','missing_value',-9999);
    %
    nccreate(fname_long,'Relative risk','Dimensions',{'projection_x_coordinate',length(projection_x_coordinate),'projection_y_coordinate',length(projection_y_coordinate),'ages',5},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,'Relative risk',Mortality_RR);
    ncwriteatt(fname_long,'Relative risk','standard_name','Relative risk');
    ncwriteatt(fname_long,'Relative risk','long_name','Relative risk of increased mortality when MMT is exceeded');
    ncwriteatt(fname_long,'Relative risk','description','Relative risk of increased mortality when MMT is exceeded');
    ncwriteatt(fname_long,'Relative risk','plot_label','Relative risk');
    ncwriteatt(fname_long,'Relative risk','missing_value',-9999);
    %
    nccreate(fname_long,'Minimum mortality temperature','Dimensions',{'projection_x_coordinate',length(projection_x_coordinate),'projection_y_coordinate',length(projection_y_coordinate)},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,'Minimum mortality temperature',Mortality_MMT);
    ncwriteatt(fname_long,'Minimum mortality temperature','standard_name','Minimum mortality temperature');
    ncwriteatt(fname_long,'Minimum mortality temperature','long_name','Minimum mortality temperature');
    ncwriteatt(fname_long,'Minimum mortality temperature','units','\B0C');
    ncwriteatt(fname_long,'Minimum mortality temperature','description','Minimum mortality temperature');
    ncwriteatt(fname_long,'Minimum mortality temperature','label_units','Minimum mortality temperature');
    ncwriteatt(fname_long,'Minimum mortality temperature','plot_label','Minimum mortality temperature');
    ncwriteatt(fname_long,'Minimum mortality temperature','missing_value',-9999);
    
    % Create plot
    if length(grid_idx) == 82 && length(grid_idy) == 112
        figure
        UK_subplot(avMortality(:,:,1),'Mean annual mortality, all population',OutputDir,lats,lons)
        set(gcf,'position',[0 0 350 350],'units','normalized')
        
        % Save output figure
        filename = [OutputDir,'MortalityMap.png'];
        warning('off','all')
        export_fig(filename,  '-png', '-nocrop', '-m5', '-zbuffer');
        warning('on','all')
        close all
    end

    
    
elseif strcmp(impactMetric,'Res')
    disp('Calculating residential discomfort')
    
    % Load defaults
    load('PreProcessedData/building_numbers_SSP2_2D.mat')
    load('PreProcessedData/building_numbers_SSP4_2D.mat')
    load('PreProcessedData/building_numbers_SSP5_2D.mat')
    load('PreProcessedData/population_SSP2_2D.mat')
    load('PreProcessedData/population_SSP4_2D.mat')
    load('PreProcessedData/population_SSP5_2D.mat')
    load('PreProcessedData/building_threshold_2D.mat')
    
    % Choose whether to run climate change only (1) or climate + socioeconomic
    % changes (2): (default = climate change only)
    socioEcScen = getenv('SOCIOECON');
    
    % If interested in socioeconomics, choose scenario and time period:
    if strcmp(socioEcScen,'CConly')
        disp('Calculating effect of climate change only (population remains unchanged from 2011 baseline)')
        disp(' ')
        population = squeeze(population_SSP5_2D(:,:,:,1));
        building_numbers = squeeze(building_numbers_SSP5_2D(:,:,:,1));
        building_threshold = squeeze(building_threshold_2D(:,:,:,:,1));
        % Find total population for each age category and time period
        population_aggregate_SSP5 = squeeze(nansum(nansum(population_SSP5_2D,1),2));
        population_aggregate = population_aggregate_SSP5(:,1);
        
    elseif strcmp(socioEcScen,'SSP2-2050')
        disp('Using default population from UK-SSPs for 2050 (SSP2)')
        disp(' ')
        popID = 4;
        population = squeeze(population_SSP2_2D(:,:,:,popID));
        building_numbers = squeeze(building_numbers_SSP2_2D(:,:,:,popID));
        building_threshold = squeeze(building_threshold_2D(:,:,:,:,popID));
        population_aggregate = population_aggregate_SSP2(:,popID);
        
    elseif strcmp(socioEcScen,'SSP2-2080')
        disp('Using default population from UK-SSPs for 2080 (SSP2)')
        disp(' ')
        popID = 5;
        population = squeeze(population_SSP2_2D(:,:,:,popID));
        building_numbers = squeeze(building_numbers_SSP2_2D(:,:,:,popID));
        building_threshold = squeeze(building_threshold_2D(:,:,:,:,popID));
        population_aggregate = population_aggregate_SSP2(:,popID);
        
    elseif strcmp(socioEcScen,'SSP4-2050')
        disp('Using default population from UK-SSPs for 2050 (SSP4)')
        disp(' ')
        popID = 4;
        population = squeeze(population_SSP4_2D(:,:,:,popID));
        building_numbers = squeeze(building_numbers_SSP4_2D(:,:,:,popID));
        building_threshold = squeeze(building_threshold_2D(:,:,:,:,popID));
        population_aggregate = population_aggregate_SSP4(:,popID);
        
    elseif strcmp(socioEcScen,'SSP4-2080')
        disp('Using default population from UK-SSPs for 2080 (SSP4)')
        disp(' ')
        popID = 5;
        population = squeeze(population_SSP4_2D(:,:,:,popID));
        building_numbers = squeeze(building_numbers_SSP4_2D(:,:,:,popID));
        building_threshold = squeeze(building_threshold_2D(:,:,:,:,popID));
        population_aggregate = population_aggregate_SSP4(:,popID);
        
    elseif strcmp(socioEcScen,'SSP5-2050')
        disp('Using default population from UK-SSPs for 2050 (SSP5)')
        disp(' ')
        popID = 4;
        population = squeeze(population_SSP5_2D(:,:,:,popID));
        building_numbers = squeeze(building_numbers_SSP5_2D(:,:,:,popID));
        building_threshold = squeeze(building_threshold_2D(:,:,:,:,popID));
        population_aggregate = population_aggregate_SSP5(:,popID);
        
    elseif strcmp(socioEcScen,'SSP2-2080')
        disp('Using default population from UK-SSPs for 2080 (SSP2)')
        disp(' ')
        popID = 5;
        population = squeeze(population_SSP5_2D(:,:,:,popID));
        building_numbers = squeeze(building_numbers_SSP5_2D(:,:,:,popID));
        building_threshold = squeeze(building_threshold_2D(:,:,:,:,popID));
        population_aggregate = population_aggregate_SSP5(:,popID);
        
    elseif strcmp(socioEcScen,'UDM')
        disp('Using population data from UDM: Loading')
        disp(' ')
        if exist('/data/inputs/Population/population_total-12km.asc','file')
            [population,RefMat] = arcgridread('population_total-12km.asc');
            population_aggregate = squeeze(nansum(nansum(population,1),2));
            
        else
            disp('Invlaid UDM data: cancelling')
            disp('-----')
            return
        end
        if exist('/data/inputs/UrbanData/dwelling_detached_total-12km.asc','file')
            [building_numbers,RefMat] = arcgridread('dwelling_detached_total-12km.asc');
            building_numbers = cat(3,building_numbers,arcgridread('dwelling_semi-detached_total-12km.asc'));
            building_numbers = cat(3,building_numbers,arcgridread('dwelling_flat_total-12km.asc'));
            building_numbers = cat(3,building_numbers,arcgridread('dwelling_terraced_total-12km.asc'));
            building_numbers_aggregate = squeeze(nansum(nansum(building_numbers,1),2));
            
        else
            disp('Invlaid UDM data: cancelling')
            disp('-----')
            return
        end
    end
    
    % Choose how much adaptation is possible, either by a numeric value which
    % is added to the MMT or by the string 'acclim', which loads the variable
    % reg_acclim.mat computed by HEAT: (default = 0 / no acclimatisation)
    
    acclimaUptake = getenv('UPTAKE');
    if isempty(acclimaUptake)
        disp('Assuming no acclimatisation or adaptation')
        disp(' ')
        adaptationName = 'No Adaptation';
        building_threshold = squeeze(building_threshold(:,:,:,1));
        
    else
        disp(['Acclimatisation or adaptation assumed in ',char(acclimaUptake),'% of housing stock'])
        disp(' ')
        deploymentPercent = double(acclimaUptake)/100;
        building_threshold_a = squeeze(building_threshold(:,:,:,2));
        building_threshold = squeeze(building_threshold(:,:,:,1));
    end
    
    
    
    %% Check for input files
    % Check if T data has been passed to the model
    files = dir([InputDir '*.nc']); % Then check if any files exist with this root
    if isempty(files)
        disp('No climate input data. Stopping.') % Cancel if not
        return
    end
    
    % Load first file as a template
    file = [files(1).folder,'/',files(1).name];
    T = ncread(file,varname);
    S = size(T);
    try
        projection_x_coordinate = ncread(file,'projection_x_coordinate');
        projection_y_coordinate = ncread(file,'projection_y_coordinate');
    catch
        disp('Projection of input data unknown, setting generic values')
        projection_x_coordinate = grid_idx;
        projection_y_coordinate = grid_idy;
    end
    
    % Find length of record in years
    dates = ncread(file,'yyyymmdd');
    
    nyears = str2double(dates(end,1:4)) - str2double(dates(1,1:4)) +1;
    
    % Based on this, create template output arrays
    population = squeeze(population(:,:,1));
    building_numbers_nan = sum(building_numbers,3);
    building_numbers_nan(building_numbers_nan == 0) = nan;
    ppl_per_dwellings = population./building_numbers_nan;
    
    days_exceeded = zeros(S(1),S(2),4);
    ppl_per_event = zeros(S(1),S(2),S(3),4);
    events = zeros(S(1),S(2),S(3));
    
    % Calculate the number of days exceeding temperature threshold
    disp('Finding days exceeding building comfort thresholds')
    for b = 1:4
        thresh = repmat(building_threshold(:,:,b),1,1,S(3));
        % Seperately, calculate for adapted houses if required
        if exist('deploymentPercent','var')
            thresh_a = repmat(building_threshold_a(:,:,b),1,1,S(3));
        end
        
        days_exceeded(:,:,b) =  sum(T >= thresh,3);
        % Calculate how many people that is affecting each day
        disp('Calculating number of people affected')
        ppl_per_event(:,:,:,b) = repmat(building_numbers(:,:,b) .* ppl_per_dwellings,1,1,S(3));
        if ~exist('deploymentPercent','var')
            ppl_per_event(:,:,:,b) = ppl_per_event(:,:,:,b).*(T >= thresh);
        else
            ppl_per_event(:,:,:,b) = ppl_per_event(:,:,:,b).*(T >= thresh)*(1-deploymentPercent) + ppl_per_event(:,:,:,b).*(T >= thresh_a)*deploymentPercent;
        end
    end
        
    disp('Identifying if population affected exceed thresholds to qualify as an "event"')
    % Set threshold for classifying an event (e.g. if 25% of the population experience discomfort)
    RDProbabilityLevel = getenv('RDTHRESH');
    if isempty(RDProbabilityLevel)
        disp('No residential discomfort threshold set: defaulting to 25%')
        disp(' ')
        RDProbabilityLevel = 0.25;
    else
        RDProbabilityLevel = str2double(string(RDProbabilityLevel))/100;
    end
    % Set threshold for classifying an event length (e.g. if 5 or more days count as 'events')
    RDNumberConsecutiveDays = getenv('RDCONSEC');
    if isempty(RDNumberConsecutiveDays)
        disp('No event length threshold set: defaulting to 5 days')
        disp(' ')
        RDNumberConsecutiveDays = 5;
    else
        RDNumberConsecutiveDays = str2double(string(RDNumberConsecutiveDays));
    end
    
    ppl_per_event_thres = sum(ppl_per_event,4)./repmat(population,1,1,S(3)) >= RDProbabilityLevel;
    
    % Check the number of consecutive days
    disp('Identifying if population affected exceed thresholds to qualify as an "extreme event"')
    consecutive_days = 0;
    for i = 1:S(1)
        for j = 1:S(2)
            for k = 1:S(3)
                if ppl_per_event_thres(i,j,k) == 0
                    % do nothing - speeds up if statement by specifying most likely option first..
                    consecutive_days = 0;
                elseif ppl_per_event_thres(i,j,k) == 1
                    consecutive_days = consecutive_days+1;
                    if consecutive_days == RDNumberConsecutiveDays
                        events(i,j,k) = 1;
                        consecutive_days = 0;
                    end
                end
            end
            consecutive_days = 0;
        end
    end
    
    disp(' ')
    
    
    %% Saving outputs
    disp('Saving outputs')
    disp(' ')
    
    % Some basic summary values for plotting
    MaxPplAffected = nan(nyears,1);
    for y = 1:nyears
        ids = 1+(y-1)*floor(S(3)/nyears):y*floor(S(3)/nyears);
        MaxPplAffected(y) = nanmax(nansum(nansum(ppl_per_event(:,:,ids),1),2),[],3);
    end
    AnnAveMaxPplAffected = nanmean(MaxPplAffected);
    AveNumEvents = nansum(nansum(nansum(events,1),2)>0,3)/nyears;
    
    dlmwrite([OutputDir,'Ann_average_people_affected.csv'],AnnAveMaxPplAffected, 'delimiter', ',', 'precision', '%i')
    dlmwrite([OutputDir,'Number_of_extreme_res_discomfort_events.csv'],AveNumEvents, 'delimiter', ',', 'precision', '%i')    
    
    % Set netCDF details
    units = 'people';
    label_units = 'people';
    plot_label = 'Number of people experiencing residential discomfort';
    standard_name = 'Number of people experiencing residential discomfort';
    long_name = 'Number of people experiencing residential discomfort';
    description = 'Number of people experiencing residential discomfort';
    
    building_char = ['Flats        ';'Detached     ';'semi-Detached';'Terraced     ']';
    
    % Write netCDF
    fname_long = [OutputDir,'HARM_',impactMetric,'_',dates(1,1:8),'-',dates(end,1:8),'.nc'];
    if exist(fname_long,'file')
        delete(fname_long)
    end
    
    nccreate(fname_long,standard_name,'Dimensions',{'projection_x_coordinate',length(projection_x_coordinate),'projection_y_coordinate',length(projection_y_coordinate),'yyyymmdd',length(dates(:,1)),'building_type',4},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,standard_name,ppl_per_event);
    ncwriteatt(fname_long,standard_name,'standard_name',standard_name);
    ncwriteatt(fname_long,standard_name,'long_name',long_name);
    ncwriteatt(fname_long,standard_name,'units',units);
    ncwriteatt(fname_long,standard_name,'description',description);
    ncwriteatt(fname_long,standard_name,'label_units',label_units);
    ncwriteatt(fname_long,standard_name,'plot_label',plot_label);
    ncwriteatt(fname_long,standard_name,'missing_value',-9999);
    
    nccreate(fname_long,'events','Dimensions',{'projection_x_coordinate',length(projection_x_coordinate),'projection_y_coordinate',length(projection_y_coordinate),'yyyymmdd',length(dates(:,1))},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,'events',events);
    ncwriteatt(fname_long,'events','standard_name','Major discomfort events');
    ncwriteatt(fname_long,'events','long_name','Days when a major discomfort event occurs');
    ncwriteatt(fname_long,'events','units','days');
    ncwriteatt(fname_long,'events','description','Days when a major discomfort event occurs exceeding population % and duration thresholds');
    ncwriteatt(fname_long,'events','label_units','days');
    ncwriteatt(fname_long,'events','plot_label','Days when a major discomfort event occurs');
    ncwriteatt(fname_long,'events','missing_value',-9999);
   
    % Add lat and long data
    nccreate(fname_long,'projection_x_coordinate','Dimensions',{'projection_x_coordinate',length(projection_x_coordinate)},'Datatype','single','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,'projection_x_coordinate',projection_x_coordinate);
    ncwriteatt(fname_long,'projection_x_coordinate','axis','X');
    
    nccreate(fname_long,'projection_y_coordinate','Dimensions',{'projection_y_coordinate',length(projection_y_coordinate)},'Datatype','single','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,'projection_y_coordinate',projection_y_coordinate);
    ncwriteatt(fname_long,'projection_y_coordinate','axis','Y');
    
    nccreate(fname_long,'yyyymmdd','Dimensions',{'yyyymmdd',length(dates(:,1)),'string64',64},'Datatype','char')
    ncwrite(fname_long,'yyyymmdd',dates);
    ncwriteatt(fname_long,'yyyymmdd','standard_name','dates');
    ncwriteatt(fname_long,'yyyymmdd','axis','Z');
    
    nccreate(fname_long,'building_type','Dimensions',{'building_type',4,'string64',64},'Datatype','char')
    ncwrite(fname_long,'building_type',building_char');
    ncwriteatt(fname_long,'building_type','standard_name','building_type');
    ncwriteatt(fname_long,'building_type','axis','ZZ');
    
    
    % Write some general attributes
    ncwriteatt(fname_long,'/','collection','HARM derived variable')
    ncwriteatt(fname_long,'/','creation_date',datestr(now))
    ncwriteatt(fname_long,'/','title','Heat impact metric')
    ncwriteatt(fname_long,'/','version','HARM v2.0')
    
    % Add meta data if same format as UKCP18 data
    if length(ppl_per_event(:,1,1,1)) == 82 && length(ppl_per_event(1,:,1,1)) == 112
        ncwriteatt(fname_long,standard_name,'grid_mapping','transverse_mercator');
        
        ncwriteatt(fname_long,'projection_x_coordinate','standard_name','easting');
        ncwriteatt(fname_long,'projection_x_coordinate','long_name','easting');
        ncwriteatt(fname_long,'projection_x_coordinate','units','m');
        
        ncwriteatt(fname_long,'projection_y_coordinate','standard_name','northing');
        ncwriteatt(fname_long,'projection_y_coordinate','long_name','northing');
        ncwriteatt(fname_long,'projection_y_coordinate','units','m');
    end
    
    % Create plot
    if length(grid_idx) == 82 && length(grid_idy) == 112
        figure
        
        UK_subplot(sum(events,3)/nyears,'Mean annual no. of extreme residential discomfort events',OutputDir,lats,lons)
        set(gcf,'position',[0 0 350 350],'units','normalized')
        
        % Save output figure
        filename = [OutputDir,'ResidentialDiscomfortMap.png'];
        warning('off','all')
        export_fig(filename,  '-png', '-nocrop', '-m5', '-zbuffer');
        warning('on','all')
        close all
    end
    
end

%% Finish up
endt = now;
fprintf('Total time taken to run: %s\n', datestr(endt-startt,'HH:MM:SS'))

disp ('HARM finished')
