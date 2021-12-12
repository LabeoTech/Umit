### ROImanager   

This app creates and edits Regions of Interest (ROIs).

##### Description
___
```ROImanager``` is a graphical interface to create, edit and save regions of interest (ROIs) associated to an image. In addition, it provides preset ROIs of the top projection of mouse cortical areas from the Mouse Allen Brain Atlas. The ```ROImanager``` app can be called from the [Visualization](/visualization_tab.md) tab, from the *Add-ons* tab in the [DataViewer](/dataviewer.md) app or as **standalone** (see syntax below). 

##### Syntax
___

```ROImanager()```: Opens only the main GUI.   

```ROImanager(datFile)```: Opens the *first frame* of the imaging data from ```datFile```. The input ```datFile``` is a string with the name of a .DAT file located in the current directory or the full path of a .DAT file created by one of the [analysis functions](/index.md/#analysisfunctions) of ***umIT***.

```ROImanager(data)```: Opens the image contained in the variable ```data``` from the Matlab's workspace. ```data``` can be either a 2-D numeric matrix or a 3-D matrix encoding an RGB image. In the latter case, the app will automatically transform the data in a grayscale image using Matlab's built-in function ```rgb2gray```.

```ROImanager(ROIfile)```: Opens the image and ROIs contained in the ```ROIfile```. The ```ROIfile``` is a string with the name of a .mat file with prefix *"ROI_"* located in the current directory or its full path. The ```ROIfile``` is created by ```ROImanager``` to store ROI information.

##### Sections
___   

* [The interface](#theinterface)
* [Image settings](#imagesettings)
* [Drawing and editing polygonal ROIs](#drawingandeditingpolygonalrois)
* [Creating point ROIs](#creatingpointrois)
* [Using the Mouse Brain Atlas preset ROIs](#usingthemousebrainatlaspresetrois)
* [Rules for ROI fitting](#rulesforroifitting)

##### The interface
___   

The graphical interface is composed of a main window containing file and data control options as well as a table listing all ROIs over the image.

![ROImanager_MainFig](/assets/img/roimanager_mainFig.png)

##### 1. Menu Bar

* **File >>**
    * ***Open File:***  Opens dialog to load an ```ROIfile```.
    * ***New:*** Creates a new ROI table. Existing ROIs will be erased.
    * ***Save File:*** Saves the ROI data to the current ```ROIfile```.
    * ***Save as...:*** Opens dialog to save ROI data to a new ```ROIfile```.
    * ***Import from CSV:*** **(Not implemented yet)** Creates new table wit ROIs from a CSV file.
    * ***Export to CSV:*** **(Not implemented yet)** Exports ROI info to a CSV file.
    * ***Open Image:*** Opens dialog to load a new image. Images can be obtained from a ```datFile```, or from a variable ```data```. The ```data``` can be a variable stored in a *.mat* file or a RGB image in the formats *PNG*, *JPEG* or *TIF*. RGB images will be transformed to grayscale 2-D images.
* **Draw >>**
    * ***Polygon:*** Creates a new customizable polygonal ROI.
    * ***Point:*** Creates a new point ROI. A point ROI consists of a single pixel from the image.
    * **Mouse Allen Brain Atlas>>** 
    Preset ROIs based on the top projection of cortical areas areas obtained from the *Mouse Allen Brain Atlas*. The ROIs are available in two formats: **Areas** or **Centroids**. Click on *>>(Areas/Centroids)>>Select (areas/centroids)* to select the one or more ROIs from the list. Access *>>(Areas/Centroids)>>Options...* to set optional parameters.
        *  ***Areas:*** Uses the cortical area boundaries as ROIs. **Option**: shrink the ROI by a given amount of pixels.
        *  ***Centroids*** Draws a circle around each area's centroid. **Option**: select the radius of the circle in pixels or milimeters.

* ***Selection >>*** ROIs are selected by checking the boxes in the *Selection* column of the ROI table.
    * **Edit>>**
        *  ***Color:*** Change color of the selected ROI. Single selection only.
        *  **Shape>>** Edit ROI shape (*works only with polygonal ROIs!*). 
            * ***Constrained Edit:*** Performs scaling, rotation and translation of selected ROIs. Can edit more than one ROI at once. In this case, the central point of rotation and scaling is the centroid of the ROI set.
            * ***Unconstrained Edit:*** Free edit of the shape and position of an ROI  by changing/adding/deleting its vertices. Single selection only.
    * ***Delete:*** Erases selected ROI(s) from the table. Multi-selection allowed.
* **Image >>**
    * **Set Origin>>** Sets the axis origin (0,0) of the image
        *  ***New:*** Allows the selection of a new origin.
        *  ***Import from file:*** Import the origin coordinates stored in a variable inside a .mat file. The coordinates must be inside the image’s limits.
        *  ***Align image to origin:*** Allows the rotation of the around the origin. This option is disabled when ```ROImanager``` is used as *Add-on* in [DataViewer](/dataviewer.md) app.
    * ***Set pixel size:*** Sets the pixel ratio in pixel per millimeter.
    * ***Set colormap:*** Sets image colormap and clipping values.
    * **Mask>>** Creates a logical mask to help with ROI drawing.
        * ***Draw new:*** Creates a new mask by drawing a polygon. Existing masks will be overwritten.
        * ***Import from file...*** Loads a logical mask from a variable inside a .mat file. 
        * ***Show/Hide:*** Toggles mask display over the image.
    * ***Export Image as Reference:*** Exports the image parameters as an *ImageReferenceFrame.mat* file which is used in the automatic and manual alignment functions of ***umIT***.

##### 2. ROI display modes
There are 3 options to display the ROIs over the image:
1. ***normal:*** ROIs are displayed with a semi-transparent face color with black contours.
2. ***contours only:*** ROIs are displayed as coloured contours and no face color.
3. ***highlight pixels:*** ROI shapes are shown as *contours only* while the pixels inside the ROIs are highlighted.

##### 3. ROI table.
A table containing the ROI name, type (Shape or point), date of creation, centroid coordinates as well as descriptive statistics of the shape and and selected pixels. 

![ROImanager_Table](/assets/img/roimanager_roitable.png)

#####  Image settings
___

#####  Drawing and editing polygonal ROIs
___

To create a new polygonal ROI:
1. Click on Draw>>Polygon.
2. Enter the name of the new ROI
3. Draw the polygon. 
4. Once the drawing is finished, one can edit the shape by adding/deleting vertices or translating the polygon across the image.
5. Finally, double-click inside the ROI to save.

<img src="https://s-belanger.github.io/Umit/assets/gifs/roimanager_creatingShapeROI.gif" alt="ROImanager_ROIShapeCreationGif"/>

Now, to edit the ROI: 
1. Select the ROI by checking the correspondent box in the *Selection* column of the ROI table.
2. Click on Selection>>Edit>> (Constrained Edit/ Unconstrained Edit). 
3. Edit the ROI shape.
4. Double-click inside the ROI to save.

###### Example of constrained edit:
<img src="https://s-belanger.github.io/Umit/assets/gifs/roimanager_constEditShapeROI.gif" alt="ROImanager_ROIShapeEdit1Gif"/>

###### Example of unconstrained edit:
<img src="https://s-belanger.github.io/Umit/assets/gifs/roimanager_UnconstEditShapeROI.gif" alt="ROImanager_ROIShapeEdit2Gif"/>

#####  Creating point ROIs
___

Point ROIs are simply ROIs consisted of a single pixel from the image. This type of ROI does not allow any type of editing. Thus, if you want to change the location of the ROI, you can simply create a new point with the same name of the one that you want to change.

To create a point ROI:
1. Click on Draw>>Point
2. Enter the name of the new ROI
3. Select the point on the image.
4. Double-click to save.

<img src="https://s-belanger.github.io/Umit/assets/gifs/roimanager_creatingPointROI.gif" alt="ROImanager_ROIShapeEdit2Gif"/>

[**<< Home**](/index.md)

#####  Using the Mouse Brain Atlas preset ROIs
___

#####  Rules for ROI fitting
___






