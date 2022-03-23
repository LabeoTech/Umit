### Extracting data from ROIs

Here we show how to extract the imaging data from regions of interest (ROIs) created using the *ROImanager* app.   
There are two ways to extract the ROI data: using the main gui (*umIToolbox*) or through the standalone version of *DataViewer* app. The sections below describe the steps to retrieve the ROI data using an event-triggered calcium imaging data as example. The data consists on visual responses to drifting gratings. The data was preprocessed using a band-pass temporal filter, split by trials and normalized to obtain the calcium signal values as &#916;F/F.

##### Extracting ROI data using the *umIToolbox* main GUI
___

The following steps show how to extract the ROI data using the main GUI. If you already created the  ROI file (```ROImasks_xxxx.mat```), skip steps 2 to 4.

1. Open the main GUI (```umIToolbox```) and load the the project file.
2. In the *Visualization* tab, choose the file to be used to create the ROIs and click on *ROImanager* button.
3. In the *ROImanager* interface, create the ROIs and save the ROI file. For more information, click [here](../../docs/devDocs/roimanager.md).
  >Note: You can also launch the *ROImanager* from the *DataViewer* app (for details, click [here](../../docs/devDocs/dataviewer.md)). In this case, note that you will have to manually save the ROI file to the subject's folder. This is important because the function used to extract the ROI data will look for the ROI file inside the subject's folder only.   

  <p align="center">
    <img src="../../assets/gifs/extractdatafromroi_steps2to3.gif" alt="ROImanager_ROIShapeCreationGif" width = "1200"/>
    <em>Steps 2 to 3.</em>
  </p>
4. Close the *ROImanager* app.
5. In the *Pipeline control panel* tab, select the recordings to extract the data launch the *Pipeline Builder* interface.
6. In the *Pipeline Builder*, add the function ```getDataFromROI``` to the pipeline.
7. Click on *Set Opts* button.
8. Type the name of the ROI file to be used by the function (default : ```ROImasks_data.mat```).
  >Note : By default, the function will calculate the spatial average of the ROIs. To change this, type the name of the aggregation function in the  ```SpatialAgg``` field. The available options are: *none* (no spatial aggregation is performed), *median*, *mean*, *min*, *max*,*std* and *mode*.   

9. Click on *Save Config* button and select the input file that will be used in the ```getDataFromROI``` function.
10. In the *Pipeline control panel* tab, click on *Run Pipeline* to execute the function over the selected elements. The function will create a .mat file (default : ```ROI_data.mat```) in the object's save folder. This file contains the data extracted from each ROI.

<p align="center">
  <img src="../../assets/gifs/extractdatafromroi_steps5to10.gif" alt="ROImanager_ROIShapeCreationGif" width = "1200"/>
  <em>Steps 5 to 10.</em>
</p>
11. In the *Analysis* tab,  select one or multiple objects containing the *ROI_data* files. Then, in the file list box, select the *ROI_data* file created in the previous step. Select the observations (i.e. ROIs) exported and click on *Select*.
12. In the *Apps* panel, the *Utilities* section one has the option to visualize the data as a table in *TableView* and has the option to save the data to a .mat file (*Save to MAT*) or to export the data to a .csv file (*Export to CSV*).
>Note: The *TableView* and *Export to CSV* options are available only for unidimensional data (with size 1 x N)!   

<p align="center">
  <img src="../../assets/gifs/extractdatafromroi_steps11to12.gif" alt="ROImanager_ROIShapeCreationGif" width = "1200"/>
  <em>Steps 11 to 12.</em>
</p>

##### Extracting ROI data using the standalone version of *DataViewer* app.
___

Info available soon.

##### How the ROI data is organized in .mat files
___

Info available soon.

##### Links to functions' documentation
___
[umIToolbox (main GUI)](../../docs/userDocs/maingui.md) | [DataViewer](../../docs/devDocs/dataviewer.md) | [ROImanager](../../docs/devDocs/roimanager.md)




[**<< Home**](../../index.md)
