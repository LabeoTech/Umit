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

```DataViewer()``` Opens the main GUI.
```DataViewer(datFile)``` Opens the app with the imaging data from ```datFile```. The input ```datFile``` is a string with the name of a .DAT file located in the current directory or the full path of a .DAT file created by one of the [analysis functions](/index.md/#analysisfunctions) of ***umIT***.

##### Exploring Image time series (Y,X,T)
If the data was not provided as input, first open the .DAT file.
![DataViewerMenu1](/assets/img/dataviewer_menu1.png)
Once the data is loaded, two figures are displayed: the *Image Window*  and the *Temporal signal*. The *Image Window* shows each frame of the movies while the *Temporal signal* window shows the amplitude of the selected pixels across time.










