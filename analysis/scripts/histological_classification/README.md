## Prediction of smoking status through histological images

This folder contains all the code used to predict the smoking status of individuals using Neural Networks. The required data is also provided in the /data folder, with exception of the smoking status data. 

The analysis was performed for 4 tissues: Lung, Esophagus Mucosa, Pancreas and thyroid

Before running the code, install the all the depedencies: 

``` 
pip install -r requirements.txt
``` 

The pipeline main code in the notebook main.ipynb. This pipeline is split into x parts


1 - Download the images/tiles
2 - Data Preprocessing
3 - Modeling (xception_ft was used, but other networks are also implemented here. See src.setup_model function).
4 - Prediction 
5- Evaluation
