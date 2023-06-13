#!/bin/bash

echo "Unzipping data if required"
tar xfv /data/inputs/Climate/*.tar.gz -C /data/inputs/Climate
tar xfv /data/inputs/UrbanData/*.tar.gz -C /data/inputs/UrbanData
tar xfv /data/inputs/Mortality/*.tar.gz -C /data/inputs/Mortality
tar xfv /data/inputs/Population/*.tar.gz -C /data/inputs/Population
unzip /data/inputs/Climate/*.zip -d /data/inputs/Climate
unzip /data/inputs/UrbanData/*.zip -d /data/inputs/UrbanData
unzip /data/inputs/Mortality/*.zip -d /data/inputs/Mortality
unzip /data/inputs/Population/*.zip -d /data/inputs/Population

echo "Copying and unzipping PreProcessedData to working directory"
cp -r /data/inputs/PreProcessedData /code/
cp /data/inputs/PreProcessedData/PreProcessedData.tar.gz /code/PreProcessedData.tar.gz
tar -xf /code/PreProcessedData.tar.gz -C /code/PreProcessedData/

echo "Removing climate data from nested directory, if necessary"
cp /data/inputs/Climate/Climate/* /data/inputs/Climate/ 

echo "Running containerised model"
/usr/bin/mlrtapp/HeatAdaptationRiskModel_2D
