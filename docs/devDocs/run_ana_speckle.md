### run_Ana_Speckle

Calculates the blood flow from the Laser Contrast Speckle Imaging (LSCI) data.

##### Description
___

This function reads the raw speckle data stored in a `speckle.dat` file and calculates flow of imaged blood vessels using the algorithm described in [[1]](https://doi.org/10.1117/1.3285504).

##### Syntax
___

`outFile = run_Ana_Speckle(SaveFolder)`: Calculates the blood flow (in arbitrary units) from the LSCI data from the file `speckle.dat` located in the folder `SaveFolder`.

##### Inputs
___
###### SaveFolder - save folder path
*character vector*   
Full path to the folder where the LSCI raw data `speckle.dat` file is stored and where the blood flow data will be saved.

##### Output(s)
___

###### outFile - blood flow file name
*cell*   
Cell containing the default name (`Flow.dat`) of the file containing blood flow data.\
> Note: This output is used by *umIT* to create and manage analysis pipelines.

The file `Flow.dat`, and it's associated meta data file `Flow.mat`, are created during the execution of the function `run_Ana_Speckle` and automatically saved in the folder `SaveFolder`.

##### Examples
___

###### Calculate the blood flow from raw LSCI data:
To calculate the blood flow of raw LSCI data stored in `speckle.dat` file:\
`outFile = run_Ana_Speckle('C:\ROOT\SAVEFOLDER');`\

##### References
___

1. Boas, David A., and Andrew K. Dunn. ‘Laser Speckle Contrast Imaging in Biomedical Optics’. Journal of Biomedical Optics 15, no. 1 (February 2010): 011109. https://doi.org/10.1117/1.3285504.

##### Links to related functions
___

[run_ImagesClassification]((../../docs/devDocs/run_imagesclassification.md)) \| [run_SpeckleMapping]((../../docs/devDocs/run_SpeckleMapping.md))


[**<< Home**](../../index.md)
