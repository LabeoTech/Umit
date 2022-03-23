## umIToolbox installation

The current version of the umIToolbox was written in Matlab and is compatible with versions ***R2019a*** and above.
The following Matlab toolboxes are necessary to run all modules of ***umIT*** :
* Imaging processing toolbox
* Signal processing toolbox
* Statistics and Machine Learning toolbox (\*not necessary for standalone versions of *DataViewer* and *ROImanager*)

To install ***umIT*** , first clone or download the source code from our [GitHub repository](https://github.com/S-Belanger/Umit).
Some imaging analysis functions are stored in another GitHub repository ([IOI_ana](https://github.com/flesage/ioi_ana)). Please, clone or download the *master* repository if you wish to use the built-in analysis functions.

Once everything is installed, create an environment variable named **Umitoolbox** with the path for the toolbox folder.
Here are some guides to create environment variables in [Linux](https://phoenixnap.com/kb/linux-set-environment-variable) and [MacOS](https://phoenixnap.com/kb/set-environment-variable-mac).
However, if your are running Matlab in Windows, execute the function *Umitoolbox_setupEnv.m* to automatically set the environment variable to your Windows account.

Now the toolbox is ready to use!   



[**<< Home**](/index.md)                                