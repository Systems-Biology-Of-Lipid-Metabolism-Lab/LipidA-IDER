# LipidAIDER

LipidA-IDER is an automated structure annotation tool for systems-level scale identification of lipid A from high resolution tandem mass spectrometry for .ms2 files.

# How to Install and Run

Step 0: Install git (if you have not done so before). Instructions [learn how to set up git here](https://github.com/git-guides/install-git)

Step 1: Open the terminal in the directory where LipidA-IDER is in. You can use the following command to do so:
~~~
$ cd \PATH_TO_LipidA-IDER
~~~

Step 2: To clone the repository enter the command into the terminal:
~~~
$ git clone https://github.com/Systems-Biology-Of-Lipid-Metabolism-Lab/LipidA-IDER
~~~

Step 3a: For users who prefer to use a Graphical User Interface, ensure that the terminal is opened in the directory where LipidA-IDER is found, then enter the following command into the terminal :
~~~
$ python LipidAIDER_GUI.py
~~~

Step 3b: For users who prefer using the Command Line, ensure that the terminal is opened in the directory where LipidA-IDER is found, and then refer to `For generic users (CLI/Terminal inputs)` for usage guide.

As of current, LipidA-IDER has only been tested on Windows. Use of LipidA-IDER on MacOS and Linux has *not* been tested. 

# Testing LipidA_IDER using demo files
In the event that you do not have initial .ms2 data, we have provided the some demo input files for your use here: [\Source\demo_input_files](/Source/demo_input_files) 

# Preparing data files
LipidAIDER currently accepts data in .ms2 format. To convert to .ms2, use MSConvert (ProteoWizard) (download from: https://proteowizard.sourceforge.io/download.html)
Settings: Output format: .ms2, Filters: Peak Picking, Algorithm: CWT, MS Levels: 2 - 

# Preparing parameter settings
For users who will be using default settings, you may skip this section.

The default parameters are provided in [\Settings\LipidAIDER_AnalysisParam.csv]. 
Users may edit this file to according to the nature of their data. The default settings for Lipid A has been recorded in [\Settings\README.md](/Settings/README.md)
It is also possible to create a duplicate settings file to save preset settings, and then selecting a saved preset it in the LipidA-IDER GUI/CLI. More info on how to select the settings as mentioned in the `For generic users (GUI option)` and `For generic users (CLI/Terminal inputs)` sections before running the analysis. 

For more info on the definition of each parameter, refer to: [\Settings\README.md](/Settings/README.md) 

# Running LipidA-IDER
## For generic users (GUI option):
Upon running [LipidAIDER_GUI.py](LipidAIDER_GUI.py) as seen from Step 3a, there will be a pop-up window that appears like this:
<br><img src='/images/init_ss.png' width="800" height="150">


To load files, simply click the 'select' button here to open the file directory to select the input .ms2 file(s) for analysis:
<br><img src='/images/load_ms2.png' width="800" height="150"><br>
Once the files are loaded, LipidA-IDER should show that the files are loaded:
<br><img src='/images/selected_files.png' width="800" height="150"><br>

If you would like to modify parameters to better suit your data, you may do so by changing the relevant setting values in the csv file found in [\Settings\LipidAIDER_AnalysisParam.csv](/Settings/LipidAIDER_AnalysisParam.csv). 
You may refer to the previous section on `Preparing parameter settings` for more information. 

Once both the desired input .ms2 files and settings files have been selected, clicking the 'Run Analysis' will prompt LipidA-IDER to start the analysis of the input files.

Upon completion of the analysis. the generated output will be found in the '\Logger\<YYYYMMDD_HHMMSS>\Batch Output' subfolder. 
i.e. if the code was run on 10 Jun 2022 at 1823 hours, the 'Batch Output' folder would be found here: [\Logger\20220610_182340](/Logger/20220610_182340).

## For generic users (CLI/Terminal inputs):
To [LIPIDAIDER_main.py](LIPIDAIDER_main.py), prepare the batch file using templates [\Settings\LipidAIDER_BatchAnalysisFiles.csv] to input the analysis parameter file and data files, and [\Settings\LipidAIDER_AnalysisParam.csv] to edit analysis parameters, based on the MS-data nature. 
Using the command python LIPIDAIDER_main.py the settings used in these files will be used. Alternatively, refer to python LIPIDAIDER_main.py -h for command functions if using user-defined file names

Upon completion of the analysis. the generated output will be found in the '\Logger\<YYYYMMDD_HHMMSS>\Batch Output' subfolder. 
i.e. if the code was run on 10 Jun 2022 at 1823 hours, the 'Batch Output' folder would be found here: [\Logger\20220610_182340](/Logger/20220610_182340).

## For developers:
The functions and scripts for further development can be found in the [\LipidA-IDER](/LipidA-IDER) subfolder. 

# Credits
When using code(s) in this repository, or any part-thereof, do help us by citing the authors of this paper:
>X.L. Guan,
>J.Y-X Loh,
>M. Lizwan,
>S.C.M. Chan,
>J. M C Kwan,
>T.P. Lim,
>T.H. Koh,
>L-Y. Hsu,
>B.T.K. Lee,
doi: TBC

## Contact
microbialipidb@gmail.com

# License
LipidA-IDER is released under the [NTUItive Dual License](LICENSE.txt).