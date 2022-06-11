# LipidAIDER

## Project Description
LipidAIDER is an automated and robust structure annotation tool for systems-level scale identification of lipid A from high resolution tandem mass spectrometry for .ms2 files.

## How to Install and Run

Step 0: Install git (if you have not done so before). Instructions [learn how to set up git here](https://github.com/git-guides/install-git)

Step 1: Open the terminal in the directory where LipidAIDER is in. You can use the following command to do so:
~~~
$ cd \PATH_TO_LIPIDAIDER
~~~

Step 2: To clone the repository enter the command into the terminal:
~~~
$ git clone https://github.com/Systems-Biology-Of-Lipid-Metabolism-Lab/LipidAIDER
~~~

Step 3: To run/start LipidAIDER, simply enter the following into the terminal:
~~~
$ python LipidAIDER_GUI.py
~~~

As of current, LipidAIDER has only been tested on Windows. Use of LipidAIDER on MacOS and Linux has *not* been tested. 

## How to use LipidAIDER
For explanation, the following was performed using the demo files found here: [\Source\demo_input_files](\Source\demo_input_files)

### For generic users:
Upon running [LipidAIDER_GUI.py](LipidAIDER_GUI.py), there will be a pop-up window that appears like this:
<br><img src='/images/init_ss.png' width="800" height="100">


To load files, simply click the 'select' button here to open the file directory to select the input .ms2 file(s) for analysis:
<br><img src='/images/load_ms2.png' width="800" height="100"><br>
Once the files are loaded, LipidAIDER should show that the files are loaded:
<br><img src='/images/selected_files.png' width="800" height="100"><br>


If in any case, you would like to use non-default parameters, you may refer to the csv file in [\Settings\AnalysisParam_CGH.csv](\Settings\AnalysisParam_CGH.csv). 
We suggest that you create a new copy following the same format as the original file, and then select it in the LipidAIDER GUI before running the analysis. 

Once both the desired input .ms2 files and settings files have been selected, clicking the 'Run Analysis' will prompt LipidAIDER to start the analysis of the input files.

Upon completion of the analysis. the generated output will be found in the '\Logger\<YYYYMMDD_HHMMSS>\Batch Output' subfolder. 
i.e. if the code was run on 10 Jun 2022 at 1823 hours, then the 'Batch Output' folder would be found here: [\Logger\20220610_182340](\Logger\20220610_182340).

### For developers:
The functions and scripts for development can be found in the [\LipidAIDER](\LipidAIDER) subfolder. 

## Credits
When using this code, or any part-thereof, please cite:

### Contact

# License
LipidAIDER is released under the [NTUItive Dual License](LICENSE.txt).