### loadDatFile
Loads the data from a `.dat` file to the workspace.

##### Description
___

This function opens the data stored in `.dat` files created by *umIT* functions in the current workspace. Optionally, opens the meta data stored in the `.mat` associated with the `.dat` file.

##### Syntax
___

`data = loadDatFile(DatFileName)`: loads the data stored in `DatFileName` in the variable `data` in Matlab's workspace.

`data = loadDatFile(DatFileName,metaDatFileName)`: loads the data stored in `DatFileName` in `data` using the meta data stored in `metaDatFileName`.

`[data, metaData] = loadDatFile(__)`: also returns the content of the meta data file (`.mat` file) associated with `DatFileName` and stores it in `metaData`.
##### Inputs
___
###### DatFileName - .dat file path
*character vector*   
Full path for the `.dat` file containing the data to be opened.

###### metaDatFileName - .mat file path (optional)
*character vector*   
Full path for the `.mat` file associated with the `.dat` file stored in `DatFileName`.\
Optional parameter. If not provided, the function will look for the `.mat` file in the same folder as `DatFileName`. Otherwise, it will use the information in `metaDatFileName` to read the data in `DatFileName`.

##### Outputs
___

###### data - output array
*numerical array*   
Multi-dimensional numerical array with dimensions set by meta data's variables `datSize` and `datLength` and data type set by meta data's `Datatype` variable.

###### metaData - data's meta data
*structure*   
Structure containing meta data associated with `data`. The `metaData` structure contains variables that describe the organization of the `DatFileName` file and other variables that may be added by an analysis function from *umIT*. The minimal variables required to open a `.dat` files are:
* *Datatype*: type of numerical data (e.g. 'single' or 'double')
* *datSize*: array with the first two dimensions of `data`. For example, `datSize = [512,512]` for an image time series of dimensions 512 by 512 by 30000.
* *datLength*: 1xN array with the size 3rd to the Nth dimension of `data`. For instance, given the example above, `datLength = 30000`.

##### Examples
___

###### Load data to workspace
To load the data from a `fluo.dat` file in Matlab's current workspace:
`data = loadDatFile('C:\ROOT\SAVEDIR\fluo.dat')`;
Here, the function looks for the meta data file `fluo.mat` in the folder `C:\ROOT\SAVEDIR\` to load `data`.

###### Load data and meta data to worskpace
Following the example above, to load the data array and meta data structure of `fluo.dat`:
`[data, metaData] = loadDatFile('C:\ROOT\SAVEDIR\fluo.dat')`;

###### Specify meta data file to load data:
To load the data from a `.dat` file with an specific meta data file:
`data = loadDatFile('C:\ROOT\SAVEDIR\fluo.dat', 'C:\ROOT\SAVEDIR\metaDataFile.mat')`;

##### Links to related functions
___

[mapDatFile](../../docs/devDocs/mapdatfile.md)


[**<< Home**](../../index.md)
