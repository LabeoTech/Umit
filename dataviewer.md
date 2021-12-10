### DataViewer
Visualization of imaging data
##### Description
___
```DataViewer``` is a graphical interface to visualize functional imaging recordings. The app provides basic tools to explore distinct types of recordings such as image time series, event-triggered time series and correlation maps. ```DataViewer``` can be called from the [Visualization](/visualization_tab.md) tab or as **standalone**. 
The app opens the following data (dimension names):
* Image time series (Y,X,T)
* Event-triggered image time series (E,Y,X,T)
* Single frames (Y,X)
* Seed-pixel correlation maps (Y,X,S)

##### Syntax
___

```DataViewer()```: Opens the main GUI.   

```DataViewer(datFile)```: Opens the app with the imaging data from ```datFile```. The input ```datFile``` is a string with the name of a .DAT file located in the current directory or the full path of a .DAT file created by one of the [analysis functions](/index.md/#analysisfunctions) of ***umIT***.

##### The interface
___   

The graphical interface is composed of a main window, containing file and data control options, an image window and a temporal signal window (if the time dimension exists).
![DataViewer_MainFigs](/assets/img/dataviewer_mainFigs.png)
##### 1. Menu Bar
* **File >>**
    * ***Open:***  Opens dialog to load previously saved .DAT files.
    * ***Import Raw data:*** Import raw data from **LabeoTech Imaging Systems**.
    * ***Export images:*** Opens dialog to export Image and temporal signal windows (PNG (default), JPEG, PDF or EPS formats).
* **Settings >>**
    * ***Image Options:*** Opens dialog to set image colormap and clipping values.
    * ***Playback speed:*** Sets movie speed for playback button (3. Controls). If the selected speed is too high, the app will reduce it to avoid lag.
* ***Utilities >>***
    * ***Data processing:*** **(Standalone only)** Tool to setup image processing pipeline on the data using the [analysis functions](/index.md/#analysisfunctions) from the toolbox. 
    *  ***Save as...:*** Open dialog to save current data into a .DAT file.
* ***Add-ons >>***
    * ***ROImanager:*** Launches the [ROImanager](/ROImanager.md) app as add-on.

##### 2. File info
Displays the current file path and the experiment type.
##### 3. Controls
The main window has some control options for image time series: 
* ***Slider:*** interactive selection of frames.
* ***Play button:*** plays a movie.   

In addition, one can move the white crosshair in the *Image Window* using the X and Y coordinates.

##### 4. Image window
Click on the image to select a pixel. For image time series, the signal amplitude of the selected pixel will be updated in the temporal profile figure. Use the figure's *Data Tips* button to add data tips containing the pixel's location and it's value.
<img src="https://s-belanger.github.io/Umit/assets/gifs/dataviewer_imagFig_clicking.gif" alt="DataViewer_Anim1" />

##### 4. Time profile window
This window will be visible only for data containing **T**ime dimension. Double-click over the plotted line to shift the current frame to the point in time. As for the *Image Window*, use the figure's *Data Tips* button to add data tips  showing the signal amplitude and time delay values.
<img src="https://s-belanger.github.io/Umit/assets/gifs/dataviewer_timeFig_clicking.gif" alt="DataViewer_Anim2" />











