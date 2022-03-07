### run_HemoCorrection
Performs hemodynamic correction of the fluorescence signal.

##### Description
___
This function removes the hemodynamic fluctuations from a fluorescence signal using data acquired with the LabeoTech's optical imaging systems.\
The function uses multi-channel recordings to measure the light absorption of the hemoglobin in order to remove it from the fluorescence channel. The algorithm used in this function is described by [Valley et al., 2020](#references).

##### Syntax
___

`outData = run_HemoCorrection(SaveFolder)`: Performs the hemodynamic correction of the fluorescence data (`fluo.dat` or `fluo_475.dat`) using 3 channels: red (`red.dat`), green (`green.dat`) and amber (`yellow.dat`) located in `SaveFolder`. It returns the corrected data  as `outData`.

`outData = run_HemoCorrection(SaveFolder, opts)`: Performs the hemodynamic correction of the fluorescence data (`fluo.dat` or `fluo_475.dat`) using the channels indicated in the [opts](#opts---optional-parameters) structure.

`[outData, metaData] = run_HemoCorrection(__)`: Outputs the meta data `metaData` associated with the data stored in `outData`.

>Important!\
This function assumes that all `.dat` files containing the fluorescence and hemodynamic data already exist in the Save Folder. Please refer to the documentation on [*run_ImagesClassification*](../../docs/devDocs/run_imagesclassification.md) for details on how to create the `.dat` files.

##### Inputs
___
###### SaveFolder - save folder path
*character vector*   
Full path to the folder where the imaging data `.dat` files are stored.

###### opts - optional parameters
*structure*   
Structure containing the optional parameters of the function.   
The `opts` structure contains the list of possible channels to be used in the hemodynamic correction:
* *Red* (bool) default = `true`: If **true**, uses the data from the *Red* channel stored in the `red.dat` file.
* *Green* (bool) default = `true`: If **true**, uses the data from the *Green* channel stored in the `green.dat` file.
* *Amber* (bool) default = `true`: If **true**, uses the data from the *Amber* channel stored in the `yellow.dat` file.

##### Output(s)
___

###### outData - output array
*3D numerical matrix*   
Matrix with dimensions *Y,X,T* containing the corrected fluorescence data.

###### metaData - data's meta data
*structure*   
Structure containing meta data associated with `outData`. For details, read the documentation of [loadDatFile](../../docs/devDocs/loaddatfile.md) function.

##### Examples
___

###### Perform hemodynamic correction using all 3 channels
To perform the hemodynamic correction using the *Red*, *Green* and *Amber* channels:\
`outData = run_HemoCorrection('C:\ROOT\SAVEFOLDER');`   

###### Perform hemodynamic correction using the Red and Green channels only:
First, create the *opts* structure:\
`opts = struct('Red',true, 'Green',true, 'Amber',false);`   
Now, the `opts` structure is passed as an input argument to the *run_HemoCorrection* function:\
`outData = run_HemoCorrection('C:\ROOT\SAVEFOLDER',opts);`   

Here, the function will use only the data stored in `red.dat` and `green.dat` to perform the hemodynamic correction in the fluorescence data.


##### References
___

1. Valley, M. T., M. G. Moore, J. Zhuang, N. Mesa, D. Castelli, D. Sullivan, M. Reimers, and J. Waters. ‘Separation of Hemodynamic Signals from GCaMP Fluorescence Measured with Wide-Field Imaging’. Journal of Neurophysiology 123, no. 1 (1 January 2020): 356–66. https://doi.org/10.1152/jn.00304.2019.

##### Links to related functions
___

[run_ImagesClassification](../../docs/devDocs/run_imagesclassification.md)









[**<< Home**](../../index.md)
