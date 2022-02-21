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

* [The interface](#the-interface)
* [Drawing and editing polygonal ROIs](#drawing-and-editing-polygonal-rois)
* [Creating point ROIs](#creating-point-rois)
* [Using the Mouse Brain Atlas preset ROIs](#using-the-mouse-brain-atlas-preset-rois)
* [Rules for ROI fitting](#rules-for-roi-fitting)

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


#####  Using the Mouse Brain Atlas preset ROIs
___

The preset ROIs were created from the top projection of the mouse cortical areas (see image below) obtained from the *Mouse Allen Brain Atlas*.

> Tip: For a more accurate result, it is advisable to set **Bregma** as the image's **origin** and to set the image's **pixel size** before drawing the ROIs from the Mouse Brain Atlas. Once these parameters are set, the ROI mask will be automatically place the mask's Bregma over the origin point and resize it to approximate the mask's real size.

![ROImanager_MouseAreas](/assets/img/roimanager_mouse-areas.png)

> Note: The position of the **Bregma** is provided here as a rough estimate. The data from the Mouse Brain atlas does not provide any anatomical landmark coordinates from the mouse skull. For more info on this, see this [discussion](https://community.brain-map.org/t/why-doesnt-the-3d-mouse-brain-atlas-have-bregma-coordinates/158) from the Allen Brain Map Community Forum.

There are two modes for the Mouse Brain Atlas preset ROIs: *Areas* and *Centroids*. The *Areas* option draws the full surface of each cortical area while the *Centroids* option creates a circle around the area's centroid.
Optional parameters for each mode can be set by clicking on *Draw >> Mouse Allen Brain Atlas >> (Areas/Centroids) >> Options...*.
#####  Area options
One can choose to use the full surface of each ROI (default) or to shrink it by N pixels. This only affects the selected pixels inside the ROIs and not the shape themselves. The shrinking algorithm removes pixels from the ROI border.

<img src="https://s-belanger.github.io/Umit/assets/img/roimanager_area-options.png" alt="ROImanager_AreaOptions" width = "400"/>

######  Example of shrinking of the Left primary visual cortex ROI by removing 15 pixels from the ROI border:

<img src="https://s-belanger.github.io/Umit/assets/img/roimanager_area-shrink-example.png" alt="ROImanager_AreaShrinkEx" width = "400"/>

#####  Centroid options
One can choose the radius of the circles in pixels (default = 1 px) and in millimeters. The latter is only available if the image's *pixel size* is already set.

<img src="https://s-belanger.github.io/Umit/assets/img/roimanager_centroid-options.png" alt="ROImanager_CentroidOptions" width = "400"/>

> Note: If the ROIs are larger than the cortical areas themselves, the circles will automatically intersect with the areas' borders to avoid selecting pixels outside the areas.


##### Area selection
Once the options are set (either for Areas or Centroids), click on Select (Areas/Centroids) button.
A table containing the columns acronyms, names, functional modality and Selection will appear.

![ROImanager_MouseAreasSelectionBox](/assets/img/roimanager_area-selection-fig.png)

> Tip: For checking/uncheking multiple chekboxes from the *Selected* column, first highlight the cells, then click outside the table and press *Enter*.

After checking the areas, close the window to save.
##### Fitting the mask
After areas selection, a mask is plotted over the image and can be fitted using a constrained edit method.
<img src="https://s-belanger.github.io/Umit/assets/gifs/roimanager_fittingMouseAreas.gif" alt="ROImanager_MouseAreasFitGif"/>

#####  Rules for ROI fitting
___

The ```ROImanager``` app will automatically detect if the ROIs are over an invalid region of the image.
Invalid regions consist of :
* Regions outside the image boundaries.
* Areas of the image containing ```NaN```s.
* Areas outside the selected area using a logical mask.

If the ROI is completely inside an invalid region, it will be automatically deleted. However, if the ROI is partially inside an invalid region, a window will appear giving a choice to delete, try to fix or ignore.
If *try to fix* is selected, ```ROImanager``` will redraw the ROI shape based on the valid pixels inside the image.
> Note: Depending on the size of the remaining ROIs, the fixing algorithm may create empty ROIs (without any pixels inside) or split the ROI. Check the descriptive statistics on the table to spot potentially broken ROIs!

If *ignore* is selected, some of the ROI shape statistics, such as the centroid coordinates or area, may be wrong. Regarding the selected pixels inside the ROIs, the pixels located outside the image will be removed, however the ones in regions with ```NaN```s or outside the logical mask areas will be kept.

###### Example of ROI fix partially outside the logical mask
<img src="https://s-belanger.github.io/Umit/assets/gifs/roimanager_fixingROIs.gif" alt="ROImanager_FixingROIGif"/>





















[**<< Home**](../../index.md)
