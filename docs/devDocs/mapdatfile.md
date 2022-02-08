### mapDatFile
Creates memory map of data to a `.dat` file.

##### Description
___

This function creates a memory map to a `.dat` file created by *umIT* functions. Optionally, does the same to the `.mat` meta data file associated with the `.dat` file. In *umIT*, memory mapping is used mainly by visualization apps, where a fast and recurrent access to the data is needed. For more information on memory mapping, please consult [Matlab's documentation](https://www.mathworks.com/help/matlab/import_export/overview-of-memory-mapping.html).

##### Syntax
___

`mmFile = mapDatFile(DatFileName)`: maps the data to `DatFileName`.

`mmFIle = mapDatFile(DatFileName,metaDatFileName)`: maps the data to `DatFileName` using the meta data stored in `metaDatFileName`.

`[mmFile, metaData] = mapDatFile(__)`: also returns the content of the meta data file (`.mat` file) associated with `DatFileName` as a Matlab's [*MAT-file object*](https://www.mathworks.com/help/matlab/ref/matlab.io.matfile.html?searchHighlight=MAT-file&s_tid=srchtitle_MAT-file_2) `metaData`.
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

###### mmFile - Memory map
*memmapfile object*   
Memmory map file that maps to `DatFileName`. The data inside `mmFile` is accessed through the `mmFile.Data` property. For more information on memmapfile objects, please consult Matlab's documentation of the [*memmapfile*](https://www.mathworks.com/help/matlab/ref/memmapfile.html?searchHighlight=memmapfile&s_tid=srchtitle_memmapfile_1) built-in function.

###### metaData - MAT-file object
*matlab.io.MatFile*   
MAT-file containing meta data associated with the data mapped in `mmFile`. For more information on MAT-file objects, please consult Matlab's documentation of the [*matfile*](https://www.mathworks.com/help/matlab/ref/matlab.io.matfile.html?searchHighlight=matfile&s_tid=srchtitle_matfile_1) built-in function.
##### Examples
___

###### Map data
To map the data from a `fluo.dat` file:\
`mmFile = mapDatFile('C:\ROOT\SAVEDIR\fluo.dat')`\
Here, the function looks for the meta data file `fluo.mat` in the folder `C:\ROOT\SAVEDIR\` to memory map the data.

###### Map data and meta data
Following the example above, to obtain a memory map and a MAT-file objects pointing to data array and meta data structure of `fluo.dat`, respectively:\
`[mmFile, metaData] = mapDatFile('C:\ROOT\SAVEDIR\fluo.dat')`

###### Specify meta data file to load data:
To map to `.dat` file with an specific meta data file:\
`mmFile = mapDatFile('C:\ROOT\SAVEDIR\fluo.dat', 'C:\ROOT\SAVEDIR\metaDatFile.mat')`

##### Links to related functions
___

[loadDatFile](../../docs/devDocs/loaddatfile.md)


[**<< Home**](../../index.md)
