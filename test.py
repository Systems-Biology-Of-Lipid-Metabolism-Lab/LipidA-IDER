#!/usr/bin/env python3

from random import random, seed
from LipidAIDER.Misc import *
import logging
import sys
import tkinter as tk
from tkinter import filedialog as fd
from tkinter import messagebox
import os

def selectMs2File():
    filename = fd.askopenfile()
    ms2_filename_var.set(filename.name)
    ms2_filename.config(width=len(ms2_filename_var.get()))

def selectSettingsFile():
    filename = fd.askopenfile()
    settings_filename_var.set(filename.name)
    settings_filename.config(width=len(settings_filename_var.get()))

def runAnalysis():
    # check the ms2 files
    if not os.path.exists(ms2_filename_var.get()):
        messagebox.showerror("Error", "Unable to find MS2 file.")
        return
    # check the settings file
    if not os.path.exists(settings_filename_var.get()):
        messagebox.showerror("Error", "Unable to find settings file.")
        return

root = tk.Tk()
root.resizable(False, False)
tk.Label(root, text="MS2 file:").grid(row=0, column=0, sticky=tk.W, pady=5, padx=5)
tk.Label(root, text="Settings file:").grid(row=1, column=0, sticky=tk.W, pady=5, padx=5)
# field for the ms2 file
ms2_filename_var = tk.StringVar()
#ms2_filename = tk.Entry(root, textvariable=ms2_filename_var, state=tk.DISABLED)
ms2_filename = tk.Entry(root, textvariable=ms2_filename_var)
ms2_filename.grid(row=0, column=1, sticky=tk.W, pady=5, padx=5)
# field for the settings file, use the default
settings_filename_var = tk.StringVar(value=os.path.abspath(os.path.join(os.path.dirname(__file__), "Settings", "AnalysisParam_CGH.csv")))
#settings_filename = tk.Entry(root, textvariable=settings_filename_var, state=tk.DISABLED)
settings_filename = tk.Entry(root, textvariable=settings_filename_var)
settings_filename.config(width=len(settings_filename.get()))
settings_filename.grid(row=1, column=1, sticky=tk.W, pady=5, padx=5)

tk.Button(root, text="Select", command=selectMs2File).grid(row=0, column=2, sticky=tk.W, pady=5, padx=5)
tk.Button(root, text="Select", command=selectSettingsFile).grid(row=1, column=2, sticky=tk.W, pady=5, padx=5)

tk.Button(root, text="Run analysis", command=runAnalysis).grid(row=2, column=0, sticky=tk.W, pady=5, padx=5)
tk.Button(root, text="Quit", command=root.destroy).grid(row=2, column=1, sticky=tk.W, pady=5, padx=5)

root.mainloop()