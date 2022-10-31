#!/usr/bin/env python3

import multiprocessing
import sys
import datetime
import os
import logging
from LipidAIDER.Setup import *
from LipidAIDER.Ion_Data import *
from LipidAIDER.Misc import *
import argparse
import shutil
from multiprocessing import Pool
from LipidAIDER.Core import *

def main(analysis_files, verbose=True):    
    time_start = datetime.datetime.now()
    no_of_files = 0
    analysis_line_counter = 0
    for analysis_file in analysis_files:
        expt_name = analysis_file["expt_name"]
        batch_name = analysis_file["batch_name"]
        ms2file = analysis_file["ms2_filename"]
        paramfile = analysis_file["settings_filename"]
        pause = 0.1
        test = AutoAnalyser(file=paramfile, ms2file=ms2file, cores=analysis_file["cores"], verbose=verbose)
        test.set_experiment_name(expt_name)
        test.set_batch_name(batch_name)
        try:
            pause = float(pause)
        except:
            pause = 0
        # read in the ms2file
        logger.info(f"Reading {ms2file} MS2 file for analysis")
        test.scans = readMs2File(ms2file, threshold_ppm=test.settings["fragments_threshold"]) # save the scan data to the test AutoAnalyser object
        logger.info(f"Found {len(test.scans)} ms1 ions for analysis")
        test.ions = [ [f"scan_{analysis_line_counter}_{i}_{test.scans[i]['parent_mz']:.3f}_{test.scans[i]['retention_time']:.3f}", i] for i in range(len(test.scans)) ]
        test.log_dir = log_dir
        batch_no_of_files = test.return_no_of_ion_files()
        no_of_files += batch_no_of_files
        test.print_names()
        test.analyse_all()
        analysis_line_counter += 1 # increment the analysis line counter by 1
    time_end = datetime.datetime.now()
    logger.info(f"Started at {time_start.strftime('%H:%M:%S')} Ended at {time_end.strftime('%H:%M:%S')}")
    seconds = (time_end-time_start).total_seconds()
    minutes = round(seconds/ 60,2)
    logger.info(f"Time elapsed {minutes} minutes")
    try:
        logger.info(f"No of Files {no_of_files}")
        logger.info(f"Time taken per file: {seconds/no_of_files:.3f} sec/file")
    except:
        logger.error("Error")    
    # remove the libraries folder unless told not to
    if not verbose:
        logger.info("Removing libraries sub directory")
        shutil.rmtree(os.path.join(log_dir, "Libraries"))
    else:
        logger.info("Keeping libraries sub directory")

logger = logging.getLogger(__name__)

if __name__ == '__main__':
    # set up the parser and parse the commandline input
    parser = argparse.ArgumentParser("LipidAIDER")
    parser.add_argument("-i", "--input", action="store", dest="input", type=str, default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "Settings", "LipidAIDER_BatchAnalysisFiles.csv"), help="Input analysis CSV file. Overwites all settings.")
    parser.add_argument("-m", "--ms2", action="store", dest="ms2", type=str, help="Input MS2 file.")
    parser.add_argument("-s", "--settings", action="store", dest="settings", type=str, default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "Settings", "LipidAIDER_AnalysisParam.csv"), help="Input MS2 file analysis setting file. Only used when the ms2 option is set.")
    parser.add_argument("-e", "--experiment", action="store", dest="experiment", type=str, default="Experiment1", help="Input MS2 file experiment name. Only used when the ms2 option is set.")
    parser.add_argument("-b", "--batch", action="store", dest="batch", type=str, default="Batch1", help="Input MS2 file batch name. Only used when the ms2 option is set.")
    parser.add_argument("-c", "--cores", action="store", dest="cores", type=int, default=multiprocessing.cpu_count(), help="Number of computational cores to use. Defaults to the number of cores on the system.")
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", help="Verbose flag. Keep intermediate output for debugging and checks.")
    parser.add_argument("-o", "--output", action="store", dest="output", type=str, default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "Logger", datetime.datetime.now().strftime("%Y%m%d_%H%M%S")), help="Output directory. Will default to the Logger directory where this script is located with a date time stamp")
    parser.add_argument("-f", "--force", action="store_true", dest="force", help="Overwrite the output directory.")
    opt = parser.parse_args()
    # process the input
    # perform checks on the output directory before creating it
    if os.path.exists(opt.output):
        if opt.force:
            print(f"Removing output directory due to -f option.")
            shutil.rmtree(opt.output)
        else:
            print(f"Unable to create a new output directory {opt.output} as it already exists.")
            sys.exit(1)
    # set up the logging
    log_dir = os.path.abspath(opt.output)
    os.makedirs(log_dir)
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(os.path.join(log_dir, "main.log"), "w"),
            logging.StreamHandler(sys.stdout)
        ]
    )

    # use specified ms2 file
    analysis_files = [] # define the analysis files
    if opt.ms2 != None:
        analysis_files.append({
            "expt_name" : opt.experiment,
            "batch_name" : opt.batch,
            "ms2_filename" : opt.ms2,
            "settings_filename" : opt.settings,
            "cores" : opt.cores,
        })        
    # use the input setting files
    else:
        for line in read_csv(opt.input)[1:]:
            analysis_files.append({
                "expt_name" : line[0],
                "batch_name" : line[1],
                "ms2_filename" : line[2],
                "settings_filename" : line[4],
                "cores" : opt.cores,
            })

    # check files to make sure they exist
    for analysis_file in analysis_files:
        # check the ms2 file
        if os.path.exists(analysis_file["ms2_filename"]):
            analysis_file["ms2_filename"] = os.path.abspath(analysis_file["ms2_filename"]) # change to abs path
        else:
            # try to find the file in the Source directory
            npath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Source", analysis_file["ms2_filename"])
            if os.path.exists(npath):
                analysis_file["ms2_filename"] = os.path.abspath(npath)
            else:
                logger.error(f"Unable to locate the MS2 file: {analysis_file['ms2_filename']}")
                sys.exit(1)
        # check the analysis settings file
        if os.path.exists(analysis_file["settings_filename"]):
            analysis_file["settings_filename"] = os.path.abspath(analysis_file["settings_filename"]) # change to abs path
        else:
            # try to find file in the Source directory
            npath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Settings", analysis_file["settings_filename"])
            if os.path.exists(npath):
                analysis_file["settings_filename"] = os.path.abspath(npath)
            else:
                logger.error(f"Unable to locate the analysis settings file: {analysis_file['settings_filename']}")
                sys.exit(1)

    # change the current working to where this script is located
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    logger.info(f"Processing with {len(analysis_files)} analysis files.")
    print(opt.verbose)
    main(analysis_files, verbose=opt.verbose)