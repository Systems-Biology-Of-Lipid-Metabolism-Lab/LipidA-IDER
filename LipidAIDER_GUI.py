#!/usr/bin/env python3

import multiprocessing
import sys
import datetime
import os
import logging
from LipidAIDER.Setup import *
from LipidAIDER.Ion_Data import *
from LipidAIDER.Misc import *
import shutil
from multiprocessing import Pool
from LipidAIDER.Core import *
import tkinter as tk
from tkinter import END, HORIZONTAL, filedialog as fd
from tkinter import messagebox

def selectMs2File():
    # initialise variables
    ms2_files_list = []
    filenames = []
    filestr = ''

    if ms2_filename_var.get():
        ms2_filename.delete('1.0', END)
        ms2_filename_var.set(filestr)
        ms2_filename.insert(float(1),'')

    # open file selection dialog
    filenames = fd.askopenfiles()
    
    for i,file in enumerate(filenames):
        ms2_filename.insert(float(i),file.name + '; \n')
        ms2_files_list.append(file.name)
        if filestr == '':
            filestr = filestr + file.name
        else:
            filestr = filestr + '; \n' + file.name
    ms2_filename_var.set(filestr)
    #ms2_filename.config(text= width=len(ms2_filename_var.get()))

    # if ms2_filename_var.get():
    #     ms2_files_list = ms2_filename_var.get()
    #     #files_list = list(ms2_files_list.split('; \n'))
    #     files_list = list(ms2_files_list)
    #     for i, text_file in enumerate(files_list):
    #         #position = f'{i}.0'
    #         ms2_filename.insert(float(i),text_file)

def selectSettingsFile():
    filename = fd.askopenfile()
    settings_filename_var.set(filename.name)
    settings_filename.config(width=len(settings_filename_var.get()))

class myDialog:
    def __init__(self, parent):
        self.top = tk.Toplevel(parent)
        self.top.title("Running")
        tk.Label(self.top, text="Code is now running.").pack()
        self.top.wait_visibility()
        self.top.transient(parent)
        self.top.grab_set()
    def destroy(self):
        self.top.destroy()

def runAnalysis():
    ms2_files_list = ms2_filename_var.get()
    ms2_files_list = list(ms2_files_list.split('; \n'))
    print('output here')
    for file in ms2_files_list:
        print(file)
    if ms2_files_list:
        for file in ms2_files_list:
                if not os.path.exists(file):
                    messagebox.showerror("Error", "Unable to find MS2 file.")
                    return

    # check settings file
    if not os.path.exists(settings_filename_var.get()):
        messagebox.showerror("Error", "Unable to find settings file.")
        return

    # change current working to where this script is located
    os.chdir(os.path.dirname(__file__))

    # run analysis
    time_start = datetime.datetime.now()
    no_of_files = 0
    for file in ms2_files_list:
        analysis = AutoAnalyser(file=settings_filename_var.get(), ms2file=file, cores=multiprocessing.cpu_count())
        analysis.set_experiment_name("Expt1")
        analysis.set_batch_name("Batch1")
        analysis.scans = readMs2File(file, threshold_ppm=analysis.settings["fragments_threshold"])
        analysis.ions = [ [f"scan_{0}_{i}_{analysis.scans[i]['parent_mz']:.3f}_{analysis.scans[i]['retention_time']:.3f}", i] for i in range(len(analysis.scans)) ]
        analysis.log_dir = log_dir
        no_of_files += analysis.return_no_of_ion_files()
        analysis.print_names()
        analysis.analyse_all()    
    time_end = datetime.datetime.now()
    logger.info(f"Started at {time_start.strftime('%H:%M:%S')} Ended at {time_end.strftime('%H:%M:%S')}")
    seconds = (time_end-time_start).total_seconds()
    minutes = round(seconds/ 60,2)
    logger.info(f"Time elapsed {minutes} minutes")
    logger.info(f"No of Files {no_of_files}")
    logger.info(f"Time taken per file: {seconds/no_of_files:.3f} sec/file")
    shutil.rmtree(os.path.join(analysis.log_dir, "Libraries"))

logger = logging.getLogger(__name__)

if __name__ == "__main__":

    # set up logging
    log_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Logger", datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
    os.makedirs(log_dir)
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(os.path.join(log_dir, "main.log"), "w"),
            logging.StreamHandler(sys.stdout)
        ]
    )

    root = tk.Tk()
    root.title(string='LipidA-IDER')
    root.resizable(False, False)
    tk.Label(root, text="MS2 file:").grid(row=0, column=0, sticky=tk.W, pady=5, padx=5)
    tk.Label(root, text="Settings file:").grid(row=2, column=0, sticky=tk.W, pady=5, padx=5)
    # field for .ms2 file
    filename_scrollbar = tk.Scrollbar(root,orient="vertical")
    filename_scrollbar.grid(row=0, column=2, sticky=tk.NS)
    ms2_filename_var = tk.StringVar()
    #ms2_filename = tk.Entry(root, textvariable=ms2_filename_var)
    ms2_filename = tk.Text(root,height=10,yscrollcommand=filename_scrollbar)
    ms2_filename.grid(row=0, column=1, sticky=tk.W, pady=5, padx=5)
    #ms2_filename = tk.Entry(root, textvariable=ms2_filename_var,xscrollcommand=filename_scrollbar.set)
    #ms2_filename.grid(row=0, column=1, sticky=tk.W, pady=5, padx=5)
    #ms2_filename.config(width=100)
    #scrollbar = tk.Scrollbar(root, orient='vertical', command=text.yview)
    #ms2_filename.focus()
    #ms2_filename.pack(side='bottom')

    filename_scrollbar.config(command=ms2_filename.yview)
    ms2_filename.config()

    # field for the settings file, use the default
    settings_scrollbar = tk.Scrollbar(root,orient="horizontal")
    settings_scrollbar.grid(row=3, column=1, sticky=tk.EW)
    settings_filename_var = tk.StringVar(value=os.path.abspath(os.path.join(os.path.dirname(__file__), "Settings", "LipidAIDER_AnalysisParam.csv")))
    settings_filename = tk.Entry(root, textvariable=settings_filename_var,xscrollcommand=settings_scrollbar.set)
    settings_filename.config(width=100)
    settings_filename.grid(row=2, column=1, sticky=tk.W, pady=5, padx=5)

    tk.Button(root, text="Select .ms2 File(s)", command=selectMs2File).grid(row=0, column=3, sticky=tk.W, pady=5, padx=5)
    tk.Button(root, text="Select Settings File", command=selectSettingsFile).grid(row=2, column=3, sticky=tk.W, pady=5, padx=5)

    tk.Button(root, text="Run analysis", command=runAnalysis).grid(row=4, column=0, sticky=tk.W, pady=5, padx=5)
    tk.Button(root, text="Quit", command=root.destroy).grid(row=4, column=1, sticky=tk.W, pady=5, padx=5)

    root.mainloop()