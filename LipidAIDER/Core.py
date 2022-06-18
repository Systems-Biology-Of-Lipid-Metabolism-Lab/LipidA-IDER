import sys
import datetime
import os
from LipidAIDER.Setup import *
from LipidAIDER.Ion_Data import *
from LipidAIDER.Misc import *
from multiprocessing import Pool

class AutoAnalyser:
    def __init__ (self, file, ms2file, cores, verbose=False):
        self.ions = []
        #Load Param
        logger.info("Auto-Analysis")
        
        self.filename = file
        self.ms2file = ms2file
        self.settings= {}
        Param_file = read_csv(file)
        mapping = {
            "neutral_loss_calc_threshold":float,
            "neutral_loss_FA_lib_filename":str,
            "neutral_loss_HG_lib_filename":str,
            "anion_HG_lib_filename":str,
            "level_2_lib_filename":str,
            "fragment_rules_filename":str,
            "neutral_loss_mode":str,
            "neutral_loss_threshold":float,       
            "fragments_type": split_to_list,
            "fragments_mode":str,
            "fragments_threshold": float,
            "fragments_lib":str,
            "filter_ion_data_threshold":int,
            "filter_ion_data_mode":str,
            "filter_NL_threshold": float,
            "filter_lipid_components_threshold":float,
            "filter_top_lipid_hits":int,
            "filter_lipid_hit_min_score":float,
            "filter_lipid_hit_min_peaks": int,
            "LpxA_threshold":int,
            "LpxD_threshold":int,
            "carbon_count_FA_lower_threshold":int,
            "carbon_count_FA_upper_threshold":int,
            "non_labile_secondary_FA":split_to_list,
            "missing_fragment_score_penalty":float,
            "backbone_fragment_penalty_multiplier":float,
            "frag_types_with_penalty":split_to_list,
            "parallel_analysis":int,
            "clear_previous_data":bool,
            "clear_libraries_after_analysis":bool,
            "RT_bin_window":float,
            "filter_penalty":float,
            "filter_FA_NL_Da_threshold":float,
            "score_cut_off_filename":str,
            "HG_chemical_composition_filename":str,
            "total_intensity_minimum_cutoff":float,
            "peaks_considered_top_n_ratio":float,
        }

        for setting,value in Param_file[1:]:
            self.settings[setting] = mapping[setting](value)
        self.settings["parallel_analysis"] = cores
        self.settings["dummy ion mz"] = 123456
        self.settings["abort_if_no_FA"] = True
        self.settings["abort_if_no_fragment"] = True
        self.verbose = verbose

        ###Continue Init###   
        self.preload_libraries()
        self.reset_batch_files()
        self.create_lipid_frag_rules_lib()

        self.bin_error_list = ("No fatty acids found", "No fragments found")

        # load the score cut off
        self.score_cut_offs = {}
        for fields in self.readSettingFile(self.settings["score_cut_off_filename"], main_setting_folder="Libraries"):
            if fields[0] != "L2HG":
                self.score_cut_offs[fields[0]] = float(fields[1])
        
        # load the element mass
        self.settings["element_mass"] = {}
        file = self.readSettingFile("Elements.csv", "Libraries")
        for i in range(len(file[0])):
            if file[0][i] != "Element":
                self.settings["element_mass"][file[0][i]] = float(file[1][i])
        logger.info(f"Loaded element mass: {self.settings['element_mass']}")
        
    def set_clear_previous(self,boo):
        if boo == True or boo == False:
            self.settings["clear_previous_data"] = boo
            
    def set_clear_library_after_analysis(self,boo):
         if boo == True or boo == False:
            self.settings["clear_libraries_after_analysis"] = boo
            
    def is_hit_logical_check(self,hit_ion,hit_dict): #check against LpxA,LpxD,parent,NL...
        NL_list = [x for x in self.ion_data.get_neutral_loss_from_parent_fragment(hit_ion) if x in self.FA_list]
        if len(NL_list) > 0:
            return neutral_loss_check(hit_dict,NL_list)
        else:
            return True
            
    def preload_libraries(self): #We preload and clean some library files for permutation generation later to save time
        #Sugars#
        logger.info(f"Loading fragments_lib")
        file = self.readSettingFile(self.settings["fragments_lib"], main_setting_folder="Libraries")[1:]
        headers = file[0]
        self.Frag_lib ={}
        for row in file[1:]:
            zipped_row = list(zip(headers,row))
            name = zipped_row[1][1]
            self.Frag_lib[name] ={}
            self.Frag_lib[name]["Mass"] = convert_string_to_float_smart(zipped_row[0][1])
            components = []
            for header, data in zipped_row[2:]:
                if "Yes" in data:
                    components.append(header)
                data = convert_string_to_float_smart(data)
                if type(data) == float:
                    self.Frag_lib[name][header]=data/100
            self.Frag_lib[name]["Components"]=components
        
        # Fatty Acids
        logger.info(f"Loading neutral_loss_FA_lib_filename")
        self.FA_lib_raw = self.readSettingFile(self.settings["neutral_loss_FA_lib_filename"], main_setting_folder="Libraries")[2:]

        # HG
        logger.info(f"Loading neutral_loss_HG_lib_filename")
        self.HG_NL_lib_raw = self.readSettingFile(self.settings["neutral_loss_HG_lib_filename"], main_setting_folder="Libraries")[2:]

        self.mz_FA_lib,self.mz_HG_NL_lib = {},{} ### in the format: {mz: name}
        self.FA_lib,self.HG_NL_lib = {},{} ### in the format: {name: mz}
        
        for lib in [self.FA_lib_raw,self.HG_NL_lib_raw]:
            for line in lib:
                mass = float(line[0])
                name = line[1]
                if lib == self.FA_lib_raw:
                    self.mz_FA_lib[mass] = name
                    self.FA_lib[name] = mass
                elif lib == self. HG_NL_lib_raw:
                    self.mz_HG_NL_lib[mass] = name
                    self.HG_NL_lib[name] = mass                   
        
        # list of all hydroxy ketenes in lib
        self.FA_hydroxy_ketenes_list = list(filter(lambda x: is_hydroxy_ketene(x) and carbon_count_between(x, self.settings["carbon_count_FA_lower_threshold"], self.settings["carbon_count_FA_upper_threshold"]), self.FA_lib.keys()))
        # list of all unsat ketenes in lib
        self.FA_unsat_ketenes_list = list(filter(lambda x: is_unsat_ketene(x) and carbon_count_between(x,self.settings["carbon_count_FA_lower_threshold"],self.settings["carbon_count_FA_upper_threshold"]), self.FA_lib.keys()))
        # list of all sat ketenes in lib
        self.FA_sat_ketenes_list = list(filter(lambda x: is_sat_ketene(x) and carbon_count_between(x,self.settings["carbon_count_FA_lower_threshold"],self.settings["carbon_count_FA_upper_threshold"]), self.FA_lib.keys()))
        
        # anion HG lib
        logger.info(f"Loading anion_HG_lib_filename")
        self.anion_HG_lib = self.readSettingFile(self.settings["anion_HG_lib_filename"], main_setting_folder="Libraries")[2:]

        # HG chemical composition
        logger.info(f"Loading HG_chemical_composition_filename")
        file = self.readSettingFile(self.settings["HG_chemical_composition_filename"], main_setting_folder="Libraries")
        header = [ x.replace("HG.", "") for x in file[0] ]
        self.HG_chemical_composition = {}
        for line in file[1:]:            
            line = { header[i] : line[i] for i in range(len(header)) }
            self.HG_chemical_composition[line["HG"]] = { k : int(v) for (k, v) in line.items() if k != "HG" }

        #Level 2 IDs#
        self.reset_L2D()
        
    def reset_L2D(self):
        # global L2D
        self.L2D = L2D_Manager()
        
    def reset_batch_files(self):
        self.batch_name = None
        self.experiment_name = None
        self.batch_results = [[
            "AnalysisParam File",
            "MS2 File",
            "Feature Scan ID",
            "RT (s)",
            "Total Height",
            "Top N Peak Total Height",
            "% Fragment Area",
            "Frag. Types Found",
            "FA found",
            "HG found",
            "Level 2 ID",
            "Chemical Formula",
            "Exp. Parent Ion",
            "Theoretical Parent Ion",
            "Difference",
            "Rank"] + self.components_in_frag_type("Di-glucosamine")+ ["Rel. Penalty","Frag","Peaks","Final Score", "Penalty number", "Total Hit FA", "Penalty", "Backbone Fragment Penalty number", "Backbone Fragment Penalty multiplier", "Peaks considered", "Score Dist.", "Found HG", "Found FA", "L2 Name","L3 Name"]]
            
    def set_batch_name(self,name):
        self.batch_name = str(name).replace(" ","_")

    def set_experiment_name(self,name):
        self.experiment_name = str(name).replace(" ","_")

    def print_names(self):
        logger.info(f"Experiment name: {self.experiment_name} batch: {self.batch_name}")
        
    def announce_settings(self):
        self.batch_analysis_parameters = [self.filename, self.ms2file]
        logger.info("Calculating neutral losses to a maximum of {neutral_loss_calc_threshold} and comparing against {neutral_loss_FA_lib_filename} and {neutral_loss_HG_lib_filename} at {neutral_loss_threshold}".format(**self.settings))
        logger.info("Building {} fragment libraries with {} and comparing at {} {}. Fragments: {}".format(len(self.settings["fragments_type"]), self.settings["fragments_lib"], self.settings["fragments_threshold"], self.settings["fragments_mode"], self.settings["fragments_type"]))
        logger.info("Peak Area filter: {} {}; \nNeutral Loss filter: {}%; \nLipid Components filter: {}% except for very labile components".format(self.settings["filter_ion_data_threshold"], self.settings["filter_ion_data_mode"], self.settings["filter_NL_threshold"], self.settings["filter_lipid_components_threshold"]))
        logger.info("LpxA filter: {}; LpxD filter: {}; Searching for fatty acids between C{} and C{}. Non-labile fatty acids include: {}".format(self.settings["LpxA_threshold"], self.settings["LpxD_threshold"], self.settings["carbon_count_FA_lower_threshold"], self.settings["carbon_count_FA_upper_threshold"], self.settings["non_labile_secondary_FA"]))
        logger.info("Fragment penalty: {}%; Backbone multiplier: {}".format(self.settings["missing_fragment_score_penalty"], self.settings["backbone_fragment_penalty_multiplier"]))
        logger.info("Saving and displaying top {} lipid A candidates".format(self.settings["filter_top_lipid_hits"]))

    def get_analysis_name(self,x):
        if ".csv" in x:
            return x.replace(".csv","")
        elif ".txt" in x:
            return x.replace(".txt","")
        else:
            return False
        
    def add_ion_files (self,list_of_files):               
         #ions should be in the form (filename,analysis name)
        if type(list_of_files) != list:
            list_of_files = [list_of_files]
        
        self.ions += list_of_files
        logger.info("Added {} ions to Auto-Analyser".format(len(list_of_files)))

    def add_ion_files_from_directory(self,directory,filter_dir = lambda x: True):
        directory = os.path.join(os_source_path,directory)
        ion_files = os.listdir(directory)
        list_of_files = []
        for ion_file in filter(filter_dir,ion_files):            
            analysis_name = self.get_analysis_name(ion_file)
            if analysis_name:
                list_of_files.append((analysis_name,os.path.join(directory,ion_file)))
        logger.info("Loaded files from directory {} to Auto-Analyser".format(directory))
        self.add_ion_files(list_of_files)

    def return_no_of_ion_files(self):
        return len(self.ions)
    
    def guess_no_of_FA(self,parent_ion_mass):
        no = (parent_ion_mass - 300)/250
        return math.ceil(no)

    def readSettingFile(self, filename, main_setting_folder=None):
        """
        Method to read a settings file. The method to attempt to use the filename as an absolute path and then search the location of the settings file (self.filename) and then lastly the mains script folders as specified by main_setting_folder.        
        """
        # initialize the possible filenames as teh absolute filename and a filename in the same directory as the settings file.
        filenames = [filename, os.path.join(os.path.dirname(self.filename), filename)]
        # add more if the main_setting_folder is not null
        if main_setting_folder != None:
            filenames.append(os.path.join(os.getcwd(), main_setting_folder, filename))
        # search through the possibilities
        for filename in filenames:
            if os.path.exists(filename):
                output = read_csv(filename)
                logger.info(f"Found setting file {filename}")
                return output
        # unable to locate the setting file
        logger.error(f"Unable to locate setting file {filenames[0]}")
        sys.exit(1)

    def pFunc(self, x):
        return self.individual_analysis(x[0], x[1])

    def pll_analysis(self):#TODO 
        with Pool(self.settings["parallel_analysis"]) as p:
            parallel_results = p.map(self.pFunc, self.ions)
        total_results = [x[0] for x in parallel_results]
        total_id = [x[1] for x in parallel_results if x[1] != None]
        list(map(lambda line: self.batch_results.extend(line), total_results))
        list(map(lambda line: self.L2D.set_library_with_score_check(line),total_id)) ###Un-adjusted score check
            
    def individual_analysis(self,ion_data_name,ion_file):
        self.current_no_ions += 1
        try:
            ion_data_name = truncate_filename(ion_data_name)
            logger.info(ion_data_name)            
            logger.info("File {} out of {} for batch {}...".format(self.current_no_ions,self.total_no_ions,self.batch_name))
            self.analyse(ion_file,ion_data_name)
        except Exception as exc:    
            error_string = sys.exc_info()[1]
            traceback.print_exc()
            logger.info("An error occurred for {}, file: {}, {}".format(ion_data_name,ion_file,error_string))
            temp = self.batch_analysis_parameters.copy()
            temp += [
                self.ion_data_name,
                self.ion_data.get_retention_time(),
                round(self.ion_data.get_total_area()),
                round(self.ion_data.get_total_area_of_fragments()),
                round(self.ion_data.get_fragment_area_ratio(),2),
                "N/A",
                self.FA_list,
                self.HG_list,
                "Unknown",
                "N/A",
                self.ion_data.get_parent_ion(),
                error_string,
            ]
            self.batch_results.append(temp)
        if self.settings["parallel_analysis"]!= 1:
            return self.batch_results[1:],self.hit_level_2_id ###Un-adjusted score check
    
    def bin_RT(self):
        ## Define position index of RT, area, RT_bin_id, bin_ref_flag
        RT_index=3
        area_index=4
        parent_mz_index=12
        RT_bin_id_index=5
        bin_ref_flag_index=6
        RT_bin_window=self.settings['RT_bin_window']
        MZ_bin_window=self.settings["fragments_threshold"]
        
        #### Pop header 
        header_line=self.batch_results.pop(0)
        header_line.insert(RT_bin_id_index,"RT_Bin_SN")
        header_line.insert(bin_ref_flag_index,"RT_Bin_TopPeak")
        
        #### Cast RT and Area as float
        self.batch_results=[[float(sublist[i]) if i in [RT_index,area_index] else sublist[i] for i in range(len(sublist))] for sublist in self.batch_results]
        
        ### Sort by RT
        self.batch_results=sorted(self.batch_results, key=lambda x: x[RT_index])

        # get a table of rt and mz for clustering
        rts = [[self.batch_results[i][RT_index], self.batch_results[i][parent_mz_index], self.batch_results[i][area_index], i] for i in range(len(self.batch_results)) ]

        # use hclust to cluster the retention time by a threshold
        clusters = getHclustByDistanceClustersForSortedNumericalList([ x[0] for x in rts ], thresholdFunc=lambda x: x <= RT_bin_window, distanceFunc=calcDist)        
        cluster_ids = []
        for i in range(len(clusters)):
            if isinstance(clusters[i], float):
                cluster_ids.append(i)
            else:
                cluster_ids += [i] * len(clusters[i])
        if len(cluster_ids) != len(rts):
            raise Exception("Hclust algorithm did not produce the same length for retention time")
        rts = [ rts[i] + [cluster_ids[i]] for i in range(len(cluster_ids)) ]

        # key the rts using the cluster id
        rts_dict = {}
        for i in range(len(rts)):
            if not rts[i][4] in rts_dict:
                rts_dict[rts[i][4]] = []
            rts_dict[rts[i][4]].append(rts[i])

        # within each cluster, check for the mz
        rts = []
        for cid in rts_dict:
            rts_dict[cid].sort(key=lambda x: x[1]) # sort by the mz
            clusters = getHclustByDistanceClustersForSortedNumericalList([ x[1] for x in rts_dict[cid] ], thresholdFunc=lambda x: x <= MZ_bin_window, distanceFunc=calcPpmDist)
            cluster_ids = []
            for i in range(len(clusters)):
                if isinstance(clusters[i], float):
                    cluster_ids.append(i)
                else:
                    cluster_ids += [i] * len(clusters[i])
            if len(cluster_ids) != len(rts_dict[cid]):
                raise Exception("Hclust algorithm did not produce the same length for m/z")
            rts_dict[cid] = [ rts_dict[cid][i] + [cluster_ids[i]] for i in range(len(cluster_ids)) ]
            rts += rts_dict[cid]

        # add the rt bin id
        rts.sort(key=lambda x: (x[4], x[5]))
        rt_id = 1
        rts[0].append(rt_id)
        for i in range(1, len(rts)):
            if rts[i][4] != rts[i-1][4] or rts[i][5] != rts[i-1][5]:
                rt_id += 1
            rts[i].append(rt_id)
        rts.sort(key=lambda x: x[3]) # sort back to the original order

        # find the max area for each of the rt bin id
        max_area = {}
        for row in rts:
            if row[6] not in max_area:
                max_area[row[6]] = row[2]
            if row[2] > max_area[row[6]]:
                max_area[row[6]] = row[2]

        # populate the ref_flag for the max area
        for row in rts:
            if max_area[row[6]] == row[2]:
                row.append(1)
            else:
                row.append(0)

        # write out the binning results for checking later
        if self.verbose:
            os_binning_path = os.path.join(self.log_dir, "Binning")
            os.makedirs(os_binning_path, exist_ok=True)
            filename = os.path.join(os_binning_path, f"{self.experiment_name}_{self.batch_name}_binning.csv")
            write_csv(filename, [["retention_time", "parent_ion_mz", "total_area", "original_order", "rt_bin", "mz_bin", "final_bin", "bin_ref_flag"]] + rts)

        # add the rt bin to the data
        for i in range(len(self.batch_results)):
            self.batch_results[i].insert(RT_bin_id_index, rts[i][6]) # add the RT bin id
            self.batch_results[i].insert(bin_ref_flag_index, rts[i][7])
        
        ### Put back header line
        self.batch_results.insert(0,header_line)
                
    def analyse_all(self):
        self.announce_settings()
        self.total_no_ions = len(self.ions)
        logger.info("Total no of files: {} for experiment {} and batch {}".format(self.total_no_ions,self.experiment_name,self.batch_name))
        self.current_no_ions = 0
        if self.settings["parallel_analysis"] != 1:
            self.pll_analysis()
        else:
            results = []
            for (ion_data_name, ion_file) in self.ions:
                self.individual_analysis(ion_data_name, ion_file)
                results += self.batch_results[1:]
            self.batch_results = [self.batch_results[0]] + results
        self.bin_RT()
        self.print_batch_results()
        self.ions = []
        self.reset_batch_files()
        
    def print_batch_results(self):
        date = datetime.date.today().strftime("%Y%m%d")
        time = datetime.datetime.now().strftime("%H_%M_%S")
        filename = date
        if self.experiment_name != None:
            filename+="_"+self.experiment_name
        filename += "_Batch_Results_"+time
        if self.batch_name != None:
            filename += "_"+self.batch_name
        os_batch_path = os.path.join(self.log_dir, "Batch Output")
        os.makedirs(os_batch_path, exist_ok=True)
        filename = os.path.join(os_batch_path, filename+".csv")

        # old full version output (not used)
        if self.verbose:
            write_csv(filename, self.batch_results)
        # new version output (currently in use) 
        else:
            selected = [i for i in range(len(self.batch_results[0])) if self.batch_results[0][i] in ("AnalysisParam File", "MS2 File", "Feature Scan ID","RT (s)", "Total Height", "RT_Bin_SN", "RT_Bin_TopPeak", "Top N Peak Total Height", "Level 2 ID", "Chemical Formula", "Exp. Parent Ion", "Theoretical Parent Ion", "Difference", "Rank", "Sugar", "C4' Headgroup", "C1 Headgroup", "C3' FA", "C3' Sec FA", "C2' ketene", "C2' Sec FA", "C3 FA", "C2 ketene", "C2 Sec FA", "Peaks", "Final Score", "Score Dist.", "L2 Name", "L3 Name")]
            o = list()
            for row in self.batch_results:
                output_result = [i < len(row) and row[i] or "" for i in selected]
                output_result.append(self.determine_if_HG_position_is_resolved(str(output_result[8])))
                o.append(output_result)
            write_csv(filename, o)
        return

    def determine_if_HG_position_is_resolved(self, L2_ID):
        # return True or False
        result = L2_ID
        non_resolved_HG_types = ['PP', 'PPP', 'PPPEtN','PPAra4N']
        if '/' not in L2_ID:
            if L2_ID == "Level 2 ID":
                return 'HG Position Resolved?'
            else:
                return ''
        L2_ID_split = L2_ID.split('/')
        HG_type = L2_ID_split[0]
        if any(non_resolved_HG_type in HG_type for non_resolved_HG_type in non_resolved_HG_types):
            return 'Type3B'
        else:
            return 'Type3A'

    def analyse(self,ion_file,ion_data_name):
        self.batch_results = [[
            "AnalysisParam File",
            "MS2 File",
            "Feature Scan ID",
            "RT (s)",
            "Total Height",
            "Top N Peak Total Height",
            "% Fragment Area",
            "Frag. Types Found",
            "FA found",
            "HG found",
            "Level 2 ID",
            "Chemical Formula",
            "Exp. Parent Ion",
            "Theoretical Parent Ion",
            "Difference",
            "Rank"] + self.components_in_frag_type("Di-glucosamine")+ ["Rel. Penalty","Frag","Peaks","Final Score", "Penalty number", "Total Hit FA", "Penalty", "Backbone Fragment Penalty number", "Backbone Fragment Penalty multiplier", "Peaks considered", "Score Dist.", "Found HG", "Found FA", "L2 Name","L3 Name"]]
        #Setting Up
        logger.info("Opening {}".format(ion_file))
        self.ion_data_name = ion_data_name
        self.ion_data = Ion_Data(self.ion_data_name, self.scans[ion_file], data_path=os.path.join(self.log_dir, "Intermediate Data"), dummy_ion_mz = self.settings["dummy ion mz"])
        
        parent_ion_mass = self.ion_data.get_parent_ion()
        self.hit_level_2_id = None
        self.msp_dict = {}
            
        #Ion Filtering and Lib Directory Creation
        self.ion_data.filter_ion_data(mode = self.settings["filter_ion_data_mode"],threshold = self.settings["filter_ion_data_threshold"])

        self.ion_lib_path = os.path.join(self.log_dir, "Libraries", ion_data_name)       
        create_directory(self.ion_lib_path)
        
        #Neutral Loss
        self.ion_data.write_ion_list()
        self.ion_data.calc_neutral_loss(self.settings["neutral_loss_calc_threshold"])

        #Build a Lipid Dict
        components = self.components_in_frag_type("Di-glucosamine")
        #Each component has its own subdictionary; with possible component:frequency
        self.Lipid_dict = {component:{} for component in components}
        self.Lipid_dict["Sugar"] = {"Di-glucosamine":1}

        #(1) Neutral Loss Analysis
        '''Compare ion data with fatty acid NL and headgroup NL libraries.
        Hits are saved in an internal library for permutation in fragment analysis'''
        
        #Comparison with HG
        self.ion_data.compare_neutral_loss(self.settings["neutral_loss_mode"],self.settings["neutral_loss_threshold"],self.settings["neutral_loss_HG_lib_filename"],self.mz_HG_NL_lib)
        HG_NL_report = self.ion_data.write_neutral_loss()[1:]
        self.HG_list = [] #list of all headgroups, no frequency filtering required
        for line in HG_NL_report:
            HG_hit = line[3]
            HG = line[4]
            parent,daughter = float(line[0]),float(line[1])

            if is_modified_headgroup(HG) == False:
                HG = convert_dehyd_to_hyd(HG)                     
                if HG not in self.HG_list and is_headgroup(HG):
                    self.HG_list.append(HG)
                    logger.info("Added {} as headgroup".format(HG))
                
        #Comparison with HG Anions
        self.ion_data.compare_ions(self.settings["neutral_loss_mode"], self.settings["neutral_loss_threshold"], self.anion_HG_lib, "anion_HG_lib")
        HG_Ion_Report = self.ion_data.write_fragment_ions()[1:]
        for line in HG_Ion_Report:
            mz = line[0]
            HG = line[2]
            threshold = float(line[3])
            HG = convert_dehyd_to_hyd(HG)          
            if HG not in self.HG_list and is_headgroup(HG):
                self.HG_list.append(HG)
                logger.info("Added {} as headgroup @ {} @ {} {}".format(HG,mz,round(threshold,3),self.settings["neutral_loss_mode"]))
        
        self.ion_data.compare_neutral_loss(self.settings["neutral_loss_mode"],self.settings["neutral_loss_threshold"],self.settings["neutral_loss_FA_lib_filename"],self.mz_FA_lib)
        logger.info(f"FA NL row counts before filtering: {len(self.ion_data.neutral_loss_report)}")
        self.ion_data.neutral_loss_report = [ x for x in self.ion_data.neutral_loss_report if isinstance(x[1], str) or x[1] >= self.settings["filter_FA_NL_Da_threshold"] ]
        logger.info(f"FA NL row counts after filtering: {len(self.ion_data.neutral_loss_report)}")

        FA_NL_report = [line[4] for line in self.ion_data.write_neutral_loss()[1:]]
        FA_CC_dict = {}
        self.FA_list = []
    
        if len(FA_NL_report) == 0 and self.settings["abort_if_no_FA"]:
            raise Exception ("No fatty acids found")

        for FA in FA_NL_report:
            if carbon_count_between(FA,self.settings["carbon_count_FA_lower_threshold"],self.settings["carbon_count_FA_upper_threshold"]):
                add_key_to_dict(FA_CC_dict,get_carbon_count(FA))
        
        #(2) Neutral Loss Frequency Filter
        # Clean the data for Fatty Acids; Note down the carbon count and add it to an internal library.       
        # filter out low frequency FA
        FA_total = sum(FA_CC_dict.values())
        for FA in FA_CC_dict:
            percent = get_percentage(FA_CC_dict[FA],FA_total,dp = 1)
            if percent > self.settings["filter_NL_threshold"]: #5% cutoff
                self.FA_list.append(FA)
                logger.info("Added {}, {}%".format(FA,percent))
            else:
                logger.info("Filtered {}, {}%".format(FA,percent))

        # Add some default headgroups
        for default_HG in ["Water","Phosphate"]: #"P+GalN"
            if default_HG not in self.HG_list:
                self.HG_list.append(default_HG)

        #sort lists...
        self.HG_list.sort();self.FA_list.sort()
        #Generate two separate FA lists: (1) ketene loss (2) acid losses
        self.FA_ketene_list_from_NL_analysis,self.FA_acid_list_from_NL_analysis = [],[]
        for count in self.FA_list: #carbon_count
            ketene = create_ketene_from_carbon_count(count)
            acid = create_acid_from_carbon_count(count)
            self.FA_ketene_list_from_NL_analysis += [ketene,convert_unhydroxy_to_hydroxy(ketene)]#secondary acids
            self.FA_acid_list_from_NL_analysis += [convert_unhydroxy_to_hydroxy(acid),convert_sat_to_unsat(acid)]

        found_hg_loss = False
        found_fa_loss = False
        # check for HG loss
        hits = [ x for x in HG_NL_report if x[0] >= self.ion_data.get_parent_ion() - self.settings["neutral_loss_threshold"]/1000 and x[0] <= self.ion_data.get_parent_ion() + self.settings["neutral_loss_threshold"]/1000 ]
        if len(hits) > 0:
            logger.info("Found a HG hit")
            logger.info(hits)
            found_hg_loss = True
        # compute the possible first neutral loss
        possible_first_neutral_loss = [ self.ion_data.get_parent_ion() - x for x in self.HG_NL_lib.values() ]
        possible_first_neutral_loss = [ [x, x - self.settings["neutral_loss_threshold"]/1000, x + self.settings["neutral_loss_threshold"]/1000] for x in possible_first_neutral_loss ]        
        # check whether it is found in the parent column of the neutral loss report
        for row in possible_first_neutral_loss:
            # check against the FA loss
            hits = [ x for x in self.ion_data.neutral_loss_report[1:] if x[0] >= row[1] and x[0] <= row[2]]
            if len(hits) > 0:
                logger.info("Found a FA hit")
                logger.info(row)
                logger.info(hits)
                found_fa_loss = True
        
        #(3) Fragment Analysis
        self.all_fragment_types_with_hits = {}#All fragment types with hits will be used for fragmentation analysis
        self.all_fragment_types_with_hits["Di-glucosamine"] = {}#we save for searching later (fragment score)
        self.all_fragment_types_with_hits["Di-glucosamine"]["No of Hit Ions"] = 999
        self.all_fragment_types_with_hits["Di-glucosamine"]["Max Mass"]= round(parent_ion_mass,3)+1
        self.all_fragment_types_with_hits["Di-glucosamine"]["Min Mass"]= 0

        self.comparison_func = lambda hit_dict: LpxA_check(hit_dict,self.settings["LpxA_threshold"]) and LpxD_check(hit_dict,self.settings["LpxD_threshold"]) and parent_secondary_check(hit_dict)

        for frag_type in self.settings["fragments_type"]:
            
            logger.info("Analysing {} fragments...".format(frag_type))
            #Generate a permutation library
            perm_lib_gen = self.create_perm_lib_gen_csv(frag_type)
            perm_filename = "{}_Fragments_{}".format(self.ion_data_name,frag_type)

            max_mass_threshold = self.Frag_lib[frag_type]["Max Mass %"]*parent_ion_mass
            min_mass_threshold = self.Frag_lib[frag_type]["Min Mass %"]*parent_ion_mass          
            perm_lib = generate_permutation_as_csv(perm_lib_gen,perm_filename,sub_folder = self.ion_lib_path,\
                                                   max_mass = max_mass_threshold, min_mass = min_mass_threshold,#exclude any mases outside the thresholds
                                                   frag_comp = self.components_in_frag_type(frag_type),#components in fragment, used for making hit-dicts for the next 2 checks
                                                   comparison_func = self.comparison_func)            
            self.ion_data.compare_ions(self.settings["fragments_mode"], self.settings["fragments_threshold"], self.readSettingFile(os.path.join(self.ion_lib_path, perm_lib))[2:], perm_lib)
            
            #Analyse the fragment ions that match and add their components to self.Lipid_dict
            Fragment_Ion_Report =self.ion_data.write_fragment_ions()[1:]
            fragment_components = self.components_in_frag_type(frag_type)
            logger.info("Verifying Hits for {} fragments...".format(frag_type))
            filtered_fragment_ion_report = {} #organised by hit ion
            
            #Removing Illogical Fragments
            for line in Fragment_Ion_Report:
                hit_ion = line[0]
                hit_name = line[2]
                hit_components = line[2].split(".")           
                hit_dict = dict(zip(fragment_components,hit_components))
                #Logical Checks
                if self.is_hit_logical_check(hit_ion,hit_dict):# checking against neutral losses, previous checks done during permutation
                    if hit_ion not in filtered_fragment_ion_report:
                        filtered_fragment_ion_report[hit_ion]={"No of Hits": 1, "Hits":[(hit_name,hit_dict)]}
                    else:
                        filtered_fragment_ion_report[hit_ion]["Hits"].append((hit_name,hit_dict))
                        filtered_fragment_ion_report[hit_ion]["No of Hits"]+=1
                else:
                    logger.info("Filtered off {} for absent fatty acids".format(hit_name,frag_type))
                    self.ion_data.del_fragment_from_ion(hit_ion,hit_name)
                    
            logger.info("Verified Hits for {}".format(frag_type))
                
            if len(filtered_fragment_ion_report) > 0:#dictionary of all fragments that have a hit
                list_of_keys = filtered_fragment_ion_report.keys()
                self.all_fragment_types_with_hits[frag_type] = {"No of Hit Ions": len(list_of_keys),\
                                                                "Max Mass": round(max(list_of_keys),3)*1.1,\
                                                                "Min Mass": round(min(list_of_keys),3)*0.9 }#Save for searching later (fragment score)
            else:
               continue
            
            for hit_ion in filtered_fragment_ion_report:
                #Weighting
                hit_weight = 100/filtered_fragment_ion_report[hit_ion]["No of Hits"]#Weight result by how many other hits are found
                ###hit_weight *= (hit_ion/max_mass_threshold)**2 #Weight result by mass
                hit_weight = round(hit_weight,3)
                logger.info("Ion Mass: {}, {} Verified Hits, Weight: {}".format(round(hit_ion,3), filtered_fragment_ion_report[hit_ion]["No of Hits"], hit_weight))

                for hit_name,hit_dict in filtered_fragment_ion_report[hit_ion]["Hits"]:
                    logger.info("Hit {}".format(hit_name))
                    
                    #Break up Components
                    for component,component_hit in hit_dict.items():
                        ###Secondary FA###
                        if is_secondary_FA(component):
                            parent_hit = hit_dict[get_parent(component)]
                            if is_component_hit_blank(component_hit):
                                exceptions = []
                                if is_unsaturated(parent_hit) or is_component_hit_blank(parent_hit):
                                    supplementary = self.get_supplementary_hits_for_component(component) #there must be a secondary FA
                                    
                                    for hit in supplementary:
                                        add_key_to_dict(self.Lipid_dict[component],hit,exception_list = exceptions,freq = round(hit_weight/len(supplementary),3))
                                else:
                                    add_key_to_dict(self.Lipid_dict[component],component_hit,exception_list = exceptions,freq = hit_weight)
                            else:
                                exceptions = []
                                add_key_to_dict(self.Lipid_dict[component],component_hit,exception_list = exceptions,freq = hit_weight)
                            
                        ###Headgroups##    
                        elif is_headgroup_component(component):
                            add_key_to_dict(self.Lipid_dict[component],component_hit,exception_list = ["NIL"],freq = hit_weight)
                            if component_hit == "Phosphate":#accounting for headgroup losses
                                supplementary = list(filter(contains_P,self.get_supplementary_hits_for_component(component)))
                                for hit in supplementary:
                                    add_key_to_dict(self.Lipid_dict[component],hit,exception_list = ["NIL"],freq = round(hit_weight/len(supplementary),3))       
                            
                        ###Amide-Linked###       
                        elif is_amide_linked_FA(component):#can be blank                            
                           #if unsat_dehydroxy, add sat_hydroxy version instead due to the possibility of dehydration and sec FA loss
                            if is_unsaturated(component_hit) and not is_hydroxy(component_hit):
                                component_hit_hydroxy = convert_unsat_to_hydroxy(component_hit)
                                add_key_to_dict(self.Lipid_dict[component],component_hit_hydroxy,exception_list = [],freq = hit_weight)
                                add_key_to_dict(self.Lipid_dict[component],component_hit,exception_list = [],freq = hit_weight)
                            else:
                                add_key_to_dict(self.Lipid_dict[component],component_hit,exception_list = ["Water"],freq = hit_weight)

                        ###Ester-Linked###
                        elif is_ester_linked_FA(component):#C3 and C3' FA
                        #if unsat_dehydroxy, add sat_hydroxy version instead due to the possibility of dehydration and sec FA loss
                        #Also, possible for #C3 and C3' FA to be H2O
                            if is_unsaturated(component_hit) and not is_hydroxy(component_hit): 
                                component_hit_hydroxy = convert_unsat_to_hydroxy(component_hit)
                                add_key_to_dict(self.Lipid_dict[component],component_hit_hydroxy,exception_list = [],freq = hit_weight)
                                add_key_to_dict(self.Lipid_dict[component],component_hit,exception_list = [],freq = hit_weight)
                            elif is_hydroxy(component_hit):
                                add_key_to_dict(self.Lipid_dict[component],component_hit,exception_list = ["NIL"],freq = hit_weight)
                            elif component_hit == "Water": #Water is added as a possibility
                                add_key_to_dict(self.Lipid_dict[component],component_hit,exception_list = ["NIL"],freq = hit_weight)
                            elif component_hit == "NIL": #there must be a FA
                                exceptions = ["NIL"]                               
                                supplementary = self.get_supplementary_hits_for_component(component)
                                for hit in supplementary:
                                    add_key_to_dict(self.Lipid_dict[component],hit,exception_list = exceptions,freq = round(hit_weight/len(supplementary),3))                                
                            else:
                                raise Exception ("Unknown hit {} for component {}".format(component_hit,component))
                            
                        elif not is_sugar(component):
                            raise Exception ("Unknown component {}".format(component))
                    
            logger.info("###self.Lipid_dict Modifications for {}###".format(frag_type))
            for key,val in self.Lipid_dict.items():
                if len(val) >0:
                    logger.info("{} \n {}".format(key,val))
                 
        ###End of Fragment###
        if len(self.all_fragment_types_with_hits) == 1 and self.settings["abort_if_no_fragment"]: #No fragments found; abort analysis
            raise Exception ("No fragments found for fragment search: {}".format(self.settings["fragments_type"]))
        
        #(4) Cleaning up of self.Lipid_dict, removing low frequency hits
        new = {}
        empty_components = []#note down

        for component in self.Lipid_dict:
            new[component] = {}
            logger.info("###Filtering {} at {}%###".format(component,self.settings["filter_lipid_components_threshold"]))
            
            #Check for empty components...
            if len(self.Lipid_dict[component]) == 0:
                component_empty = True
                empty_components.append(component)
                supplementary = self.get_supplementary_hits_for_component(component)                  
                logger.info("No hits found for {}, added following items".format(component))
                for item in supplementary:
                   add_key_to_dict(self.Lipid_dict[component],item,freq = 1/len(supplementary))
            else:
                component_empty = False
 
            total = sum(self.Lipid_dict[component].values())
            max_percent = get_percentage(max(self.Lipid_dict[component].values()),total)
            if max_percent < self.settings["filter_lipid_components_threshold"]:
                threshold = max_percent
                logger.info("Reduced threshold for component {} from {} to max. freq. {}".format(component,self.settings["filter_lipid_components_threshold"],max_percent))
            else:
                threshold = self.settings["filter_lipid_components_threshold"]
            
            for hit,count in self.Lipid_dict[component].items():
                percent = get_percentage(count,total)
                    
                if percent >= threshold or is_very_labile(component) or component_empty:
                    new[component][hit] = percent
                    logger.info("Added in {} , {}% in {}".format(hit,percent,component))

                else:
                    logger.info("Filtered off {} , {}% in {}".format(hit,percent,component))
                    
        self.Lipid_dict = new
        
        #Generate the permutations
        Lipid_perm_gen = self.create_lipid_perm_lib_gen_csv()
        Lipid_filename = "{}_Lipid_Di-glucosamine".format(self.ion_data_name)
        
        max_mass_threshold = self.Frag_lib["Di-glucosamine"]["Max Mass %"]*parent_ion_mass
        min_mass_threshold = self.Frag_lib["Di-glucosamine"]["Min Mass %"]*parent_ion_mass
        Lipid_lib = generate_permutation_as_csv(Lipid_perm_gen,Lipid_filename,sub_folder = self.ion_lib_path,\
                                                max_mass = max_mass_threshold,min_mass = min_mass_threshold,
                                                frag_comp = self.components_in_frag_type("Di-glucosamine"),
                                                comparison_func = self.comparison_func)
        self.ion_data.compare_ions(self.settings["fragments_mode"], self.settings["fragments_threshold"], self.readSettingFile(os.path.join(self.ion_lib_path, Lipid_lib))[2:], Lipid_lib, parent_only = True)

        #Find and rank the hits
        Lipid_final_data = self.ion_data.write_fragment_ions()[1:]
        Lipid_final_hits = len(Lipid_final_data)
        if Lipid_final_hits == 0:
            raise Exception("No lipid final hits found for library")
        logger.info("{} hits found for Auto-Analysis of {} from file {}".format(Lipid_final_hits,self.ion_data_name,ion_file))
        logger.info("Filtering hits...")
        fragment_components = self.components_in_frag_type("Di-glucosamine")
        contents = []

        # troubleshooting files if verbose is true
        if self.verbose:
            penalty_file = open(os.path.join(self.log_dir, "Intermediate Data", self.ion_data_name, "hit_penalty_details.tsv"), "w")
            penalty_header = ["hit_name", "fragment_name", "lost_FA", "loss_type", "rounded_mass", "hit_no_of_FA", "penalty", "backbone_fragment_penalty_multipler"]
            penalty_file.write("\t".join(penalty_header) + "\n")
            hit_file = open(os.path.join(self.log_dir, "Intermediate Data", self.ion_data_name, "hit_details.tsv"), "w")
            hit_file.write("\t".join(["hit_name", "level_2_id_string", "hit_score_frag_mod", "lipid_frag_total_penalty", "found_hg_loss", "found_fa_loss", "total_area", "total_area_ok", "peaks_considered", "top_n", "peaks_considered_top_n_ratio", "hit_peaks_considered_top_n_ratio_ok", "hit_score"]) + "\n")        

        self.penalty_data = {} # dict to store all the penalty data
        for line in Lipid_final_data:
            hit_mass = line[0]#experimental
            hit_ion = round(hit_mass,3) 
            hit = line[1] #theoretical
            hit_name = line[2]          
            hit_components = hit_name.split(".")
            hit_dict = dict(zip(fragment_components,hit_components))
            hit_threshold = line[3]+" "+line[4]

            #Level 2 ID#
            hit_level_2_id = self.L2D.get_l2d_from_hit_dict(float(hit),hit_dict,ion_file,self.ion_data_name,0,self.HG_NL_lib)
            self.level_2_id_string = self.L2D.get_l2d_string(hit_level_2_id)
            
            #Lipid A Categorisation#
            lipid_type = get_lipid_type(hit_dict)
            self.hit_no_of_FA = len([component_hit for component,component_hit in hit_dict.items() if is_FA(component) and not is_component_hit_blank(component_hit)])

            self.lipid_frag_counted_area = 0 #accounted area
            self.lipid_frag_counted_peaks = 0 #accounted peaks
            self.lipid_frag_counted_peaks_list = []#avoid double-counting
            self.lipid_frag_total_area = self.ion_data.get_total_area_of_fragments()#total area to account for
            self.lipid_frag_total_penalty = 0
            self.penalty_data[hit_name] = [] # add a new hit to the penalty data
            self.current_penalty_data = self.penalty_data[hit_name] # set the current penalty data to the current one
            
            if lipid_type not in self.lipid_frag_rules_lib:
                logger.info("Rules for this type ({}) are not confirmed yet, it will be ignored".format(lipid_type))
                continue
            
            logger.info("###Accepted {} at mass {}###".format(hit_name,hit_ion))
            self.msp_dict[hit_name] = MSPfile(hit_name,hit, "[M-H]-", rt=self.ion_data.get_retention_time(), comment=f"Scan={self.ion_data_name}; filename={os.path.basename(self.ms2file)}", L2D=self.level_2_id_string)
            self.current_msp = self.msp_dict[hit_name]
            
            for frag_type in self.lipid_frag_rules_lib[lipid_type]["Fragments Analysed"]:#Fragment Score
                hit_dict_copy = copy.deepcopy(hit_dict)
                if frag_type == "Di-glucosamine":
                    self.check_all_lipid_fragments (frag_type,hit_dict_copy,lipid_type,"[M-H]")
                else:
                    frag_dict = {i: hit_dict_copy[i] for i in self.components_in_frag_type(frag_type)}
                    frag_dict["Sugar"] = frag_type
                    self.check_all_lipid_fragments (frag_type,frag_dict,lipid_type,frag_type)
                logger.info("{} Fragment Score: {}%".format(frag_type,100*round(self.lipid_frag_counted_area/self.lipid_frag_total_area,5)))
                
            if self.lipid_frag_total_area > 0:
                hit_score_frag_mod = round(self.lipid_frag_counted_area/self.lipid_frag_total_area,5)
            else:
                hit_score_frag_mod = 0
                
            logger.info("###Fragment Score: {} % for {} peaks###".format(100*hit_score_frag_mod,self.lipid_frag_counted_peaks))
            # ensure that the penalty does not exceed 100%, set at 100% when it exceed 100%
            if self.lipid_frag_total_penalty > 100:
                logger.info(f"Total penalty {self.lipid_frag_total_penalty} exceeds 100%, setting it to 100%")
                self.lipid_frag_total_penalty = 100.0
            logger.info("Total penalty: {}".format(round(self.lipid_frag_total_penalty,2)))
            hit_score = round(100*hit_score_frag_mod * (1-self.lipid_frag_total_penalty/100), 1)

            # subtract the penalty for the 2 filters
            if found_hg_loss == False and found_fa_loss == False:
                logger.info("HG and FA not found filter activated, subtracting score")
                hit_score -= self.settings["filter_penalty"]

            # set the score to 0 if it is negative
            if hit_score < 0:
                logger.info("Score is less than 0, setting it to 0.")
                hit_score = 0

            logger.info("###Tentative Score: {}###".format(hit_score)) 

            # check for the total area
            total_area_ok = True
            if self.ion_data.get_total_area_of_fragments() < self.settings["total_intensity_minimum_cutoff"]:
                logger.info(f"Total area {self.ion_data.get_total_area_of_fragments():.6f} is less than the minimum {self.settings['total_intensity_minimum_cutoff']:.2f}")
                total_area_ok = False

            # check for the peaks considered
            hit_peaks_considered_top_n_ratio = float(len(self.ion_data.ion_list)) / float(self.settings["filter_ion_data_threshold"])
            hit_peaks_considered_top_n_ratio_ok = True
            if hit_peaks_considered_top_n_ratio < self.settings["peaks_considered_top_n_ratio"]:
                logger.info(f"Peaks considered top n ratio {hit_peaks_considered_top_n_ratio:.6f} is less than the minimum {self.settings['peaks_considered_top_n_ratio']:.2f}")
                hit_peaks_considered_top_n_ratio_ok = False

            # if verbose is true, write out hit details
            if self.verbose:
                # write out the penalty details
                for row in self.penalty_data[hit_name]:
                    o = [hit_name]
                    o += [ x in row and str(row[x]) or "" for x in penalty_header[1:] ]
                    penalty_file.write("\t".join(o) + "\n")
                # write out the hit details
                hit_file.write(f"{hit_name}\t{self.level_2_id_string}\t{hit_score_frag_mod}\t{self.lipid_frag_total_penalty}\t{found_hg_loss}\t{found_fa_loss}\t{self.ion_data.get_total_area_of_fragments():.6f}\t{total_area_ok}\t{len(self.ion_data.ion_list)}\t{self.settings['filter_ion_data_threshold']}\t{hit_peaks_considered_top_n_ratio}\t{hit_peaks_considered_top_n_ratio_ok}\t{hit_score}\n")

            #Top 5 only      
            if hit_score > self.settings["filter_lipid_hit_min_score"] and self.lipid_frag_counted_peaks >= self.settings["filter_lipid_hit_min_peaks"] + how_many_labile_components(hit_dict) and len(self.current_msp.get_fragment_types()) > 0 and (found_hg_loss == True or found_fa_loss == True) and total_area_ok and hit_peaks_considered_top_n_ratio_ok:
                if len(contents) < self.settings["filter_top_lipid_hits"] or hit_score > contents[-1][-1]:
                    if len(contents) == self.settings["filter_top_lipid_hits"]:
                        contents.pop(-1)
                    hit_dict = dict(zip(fragment_components,hit_components))  
                    contents.append([hit_name, hit_dict, hit_mass, hit, hit_threshold]+\
                                    hit_components+\
                                    [round(self.lipid_frag_total_penalty,1),round(hit_score_frag_mod,3),self.lipid_frag_counted_peaks,hit_score])
                    contents.sort(key = lambda x: x[-1],reverse = True)

        # close the troubleshooting files
        if self.verbose:
            penalty_file.close()
            hit_file.close()

        # write out the penalty data to the data directory for troubleshooting if verbose is true
        if len(contents) > 0 and self.verbose:
            with open(os.path.join(self.log_dir, "Intermediate Data", self.ion_data_name, "penalty_details.tsv"), "w") as outfile:
                hit_names = [ x[0] for x in contents ]
                headers = ["hit_name"] + list(self.penalty_data[hit_names[0]][0].keys())
                # write out headers
                outfile.write("\t".join(headers) + "\n")
                # write out the penalties
                for hit_name in hit_names:
                    for row in self.penalty_data[hit_name]:
                        o = [hit_name] + [ str(row[x]) for x in headers[1:] ]
                        outfile.write("\t".join(o) + "\n")

        #Results#
        if len(contents)>0:
            logger.info("Top Hits for {}".format(self.ion_data_name))
            min_penalty = min(line[-4] for line in contents)
            logger.info("Minimum Penalty: {}".format(min_penalty)) ###Adjusted Score###
            # compute the score ranking
            hit_scores = [ x[-1] for x in contents ]
            if len(hit_scores) > 1:
                hit_ranks = rankByRevStdCompetition(hit_scores)
            else:
                hit_ranks = [1]
            # iterate through the hits
            for rank,line in enumerate(contents):
                hit_name = line[0]
                hit_dict = line[1]
                hit_mass = line[2] #experimental
                hit = line[3] #theoretical
                hit_threshold = line[4]
                line[-1] = round(line[-1] ,ndigits=1)
                hit_score = line[-1]
                level_3_id_name = self.L2D.get_l3d_name_from_hit_dict(hit_dict)
                self.current_penalty_data = self.penalty_data[hit_name]
                
                #Level 2 ID#
                self.hit_level_2_id = self.L2D.get_l2d_from_hit_dict(float(hit),hit_dict,ion_file,self.ion_data_name,float(hit_score),self.HG_NL_lib)
                #Rank 1#
                if rank == 0:
                     #Level 2 ID Addition#
                    self.L2D.set_library_with_score_check(self.hit_level_2_id) ###Un-adjusted score check
                    self.level_2_id_string = self.L2D.get_l2d_string(self.L2D.get_library(float(hit)))
                    self.level_2_id_name = self.hit_level_2_id["L2 Name"]
                    
                #MSPfile#
                msp_obj = self.msp_dict[hit_name]
                for HG_Report_line in HG_Ion_Report: ### Addition of HG to msp
                    HG_mz,HG_hit,HG_name = HG_Report_line[:3]
                    if self.hit_level_2_id["HG mz"] >= HG_mz:
                        if is_ion_from_HG(hit_dict["C4' Headgroup"],HG_name) or is_ion_from_HG(hit_dict["C1 Headgroup"],HG_name):
                            HG_percent = round(100*self.ion_data.get_area(HG_mz)/self.lipid_frag_total_area,2)
                            msp_obj.add_peak(float(HG_hit),HG_percent,HG_name)

                msp_obj.add_comment("L2D={}".format(msp_obj.get_L2D()))                
                msp_obj.add_L3_name(level_3_id_name)
                
                ###Addition of Parent Ion###
                msp_parent_ion = msp_obj.get_precursor_mz()
                msp_parent_ion_area = self.ion_data.compare_ions_get_total_area(self.settings["fragments_mode"],self.settings["fragments_threshold"],"[M-H]-",msp_parent_ion,save=False)[1]
                if msp_parent_ion_area >= self.settings["dummy ion mz"]:
                    msp_parent_ion_area -= self.settings["dummy ion mz"] #adjust for dummy parent ion
                if msp_parent_ion_area > 0:
                    msp_parent_ion_percent = round(100*msp_parent_ion_area/self.lipid_frag_total_area,2)
                    msp_obj.add_peak(msp_parent_ion,msp_parent_ion_percent,"[M-H]-")

                ###Print & Save##
                msp_contents = msp_obj.print_msp()
                msp_mz = round(msp_obj.get_precursor_mz(),3)
                msp_mz_path = os.path.join(self.log_dir,"msp",str(msp_mz))
                os.makedirs(msp_mz_path, exist_ok=True)
                msp_path = os.path.join(msp_mz_path, self.ion_data_name)
                os.makedirs(msp_path, exist_ok=True)
                msp_filename = os.path.join(msp_path, self.ion_data_name + "_rank_" + str(rank+1) + ".txt")
                write_txt(msp_filename, msp_contents)
                msp_fragment_types = msp_obj.get_fragment_types()
                
                # compute the hit_score_dist
                if len(contents) != rank+1:
                    logger.info(f"Computing score_dist {hit_score} {contents[rank+1][-1]}")
                    hit_score_dist = round(hit_score - (contents[rank+1][-1]), 1)
                    str_diff = compare_hit_dict(hit_dict,contents[rank+1][1])
                else:
                    level_2_HG = self.level_2_id_string[:self.level_2_id_string.find("/")]
                    if level_2_HG in self.score_cut_offs:
                        hit_score_dist = round(hit_score - self.score_cut_offs[level_2_HG],1)
                    else:
                        hit_score_dist = round(hit_score - self.settings["filter_lipid_hit_min_score"],1)
                
                if len(contents) > 0 and rank != 0:
                    str_diff = compare_hit_dict(hit_dict,contents[0][1])
                else:
                    str_diff = ""
                
                hit_name = clean_name(hit_name)
                logger.info("Hit Name: {}, Hit Score: {}, Score Dist: {}".format(hit_name,hit_score,hit_score_dist))
                ###Batch Results##
                (chem_comp, hit) = computeL2ChemicalFormula(self.level_2_id_string, self.HG_chemical_composition, self.settings["element_mass"])
                # recompute the difference
                hit_threshold = f"{abs(hit_mass - hit) / hit * 1E6:0.3f} ppm"
                temp = self.batch_analysis_parameters.copy()
                temp += [
                    self.ion_data_name,
                    self.ion_data.get_retention_time(),
                    round(self.ion_data.get_total_area()),
                    round(self.ion_data.get_total_area_of_fragments()),
                    round(self.ion_data.get_fragment_area_ratio(),2),
                    msp_fragment_types,
                    self.FA_list,
                    self.HG_list,
                    self.level_2_id_string,
                    chem_comp,
                    hit_mass,
                    f"{hit:0.6f}",
                    hit_threshold,
                    hit_ranks[rank],
                ]
                temp += list(map(clean_name,line[5:]))
                temp += [
                    len([ x for x in self.current_penalty_data if x["penalty"] != None ]),
                    len(self.current_penalty_data) > 0 and self.current_penalty_data[0]["hit_no_of_FA"] or 0,
                    self.settings["missing_fragment_score_penalty"],
                    len([ x for x in self.current_penalty_data if x["backbone_fragment_penalty_multipler"] != None ]),
                    self.settings["backbone_fragment_penalty_multiplier"],
                    len(self.ion_data.ion_list),
                    hit_score_dist,
                    found_hg_loss,
                    found_fa_loss,
                    self.level_2_id_name,
                    level_3_id_name
                ]
                temp[28] = msp_obj.num_peaks
                self.batch_results.append(temp)
            contents
            filename = "###Final Results_{}###.csv".format(self.ion_data_name)
            filename = os.path.join(self.ion_data.get_path(),filename)
            write_csv(filename,contents)
        else:
            (chem_comp, hit) = computeL2ChemicalFormula(self.level_2_id_string, self.HG_chemical_composition, self.settings["element_mass"])
            #if got hit but score all very low
            temp = self.batch_analysis_parameters.copy()
            temp += [
                self.ion_data_name,
                self.ion_data.get_retention_time(),
                round(self.ion_data.get_total_area()),
                round(self.ion_data.get_total_area_of_fragments()),
                round(self.ion_data.get_fragment_area_ratio(),2),
                "N/A",
                self.FA_list,
                self.HG_list,
                self.level_2_id_string,
                chem_comp,
                self.ion_data.get_parent_ion(),
                "Putative L2 proposed, insufficient data for Level 3",
            ]
            self.batch_results.append(temp)             
        logger.info("###End of Analysis for {}###".format(self.ion_data_name))

    def get_supplementary_hits_for_component(self,component):
        if is_headgroup_component(component):
            supplementary = self.HG_list
        elif is_ester_linked_FA(component):
            supplementary = list(filter(lambda x: is_hydroxy(x) and is_saturated(x) and is_acid(x),self.FA_acid_list_from_NL_analysis))+["Water"]
        elif is_amide_linked_FA(component):
            supplementary = ["NIL"]+list(filter(is_even_carbon_count,\
                                            self.FA_hydroxy_ketenes_list+self.FA_unsat_ketenes_list))
        elif is_secondary_FA(component):
            if component in self.settings["non_labile_secondary_FA"]:
                supplementary = ["NIL"]+self.FA_sat_ketenes_list
            else:
                supplementary = ["NIL"]+self.FA_ketene_list_from_NL_analysis

        return supplementary
    
    def components_in_frag_type(self,frag_type):
        if frag_type in self.Frag_lib:
            return self.Frag_lib[frag_type]["Components"]
        else:
            raise Exception("Fragment type {} not supported".format(frag_type))
    
    def create_perm_lib_gen_csv(self,frag_type):
        components = self.components_in_frag_type(frag_type)

        contents = []
        contents.append(["Created by Auto-Analyser for file {} for analysis of {} fragments".format(self.ion_data_name,frag_type)])
        contents.append(["m/z","Name"])

        for component in components:
            contents.append(["Component",component])
            
            if is_sugar(component):
                mass =self.Frag_lib[frag_type]["Mass"]
                contents.append([mass,frag_type])

            elif is_headgroup_component(component):
                if "Water" not in self.HG_list:
                    self.HG_list.append("Water")
                for HG in self.HG_list:
                    mass = self.HG_NL_lib[HG]
                    contents.append([mass,HG])
            
            elif is_ester_linked_FA(component): #acid loss, limited range
                contents.append([0,"NIL"]) #account for acid loss of C3 FA
                contents.append([self.HG_NL_lib["Water"],"Water"]) #account for ketene loss of C3 FA
                for FA in self.FA_acid_list_from_NL_analysis:
                    mass =self.FA_lib[FA]
                    contents.append([mass,FA])
                        
            elif is_secondary_FA(component): #ketene loss, limited range
                contents.append([0,"NIL"]) #account for acid loss of Sec FA or its absence
                if component == "C2' Sec FA":
                    ketene_list = self.FA_sat_ketenes_list+self.FA_unsat_ketenes_list#any C# is possible for sec FA2
                else:
                    ketene_list = self.FA_ketene_list_from_NL_analysis #else look at neutral loss
                for FA in ketene_list:
                    mass =self.FA_lib[FA]
                    contents.append([mass,FA])
                    
            elif is_amide_linked_FA(component): #ketene loss, any C# is possible
                contents.append([0,"NIL"]) #unlikely but possible
                for FA in filter(is_even_carbon_count,\
                                 self.FA_hydroxy_ketenes_list+self.FA_unsat_ketenes_list):
                    mass =self.FA_lib[FA]
                    contents.append([mass,FA])
        filename = "PermLibGen_{}.csv".format(frag_type)
        pathname = os.path.join(self.ion_lib_path,filename)
        write_csv(pathname,contents)
        return os.path.join(self.ion_lib_path,filename)

    def create_lipid_perm_lib_gen_csv(self,frag_type = "Di-glucosamine"):
        contents = []
        contents.append(["Created by Auto-Analyser for file {} for analysis of {} fragments".format(self.ion_data_name,frag_type)])
        contents.append(["m/z","Name"])
        
        for component in self.Lipid_dict:
            contents.append(["Component",component])
            
            if is_sugar(component):
                mass =self.Frag_lib[frag_type]["Mass"]
                contents.append([mass,frag_type])
                    
            elif is_FA(component) or is_headgroup_component(component):
                if is_headgroup_component(component) and "Water" not in self.Lipid_dict[component]:
                    self.Lipid_dict[component]["Water"] = "N/A"
                    logger.info("Added Water to {}".format(component))
                for item in self.Lipid_dict[component]:
                    if item in self.HG_NL_lib:
                        mass = self.HG_NL_lib[item]
                    else:
                        mass = self.FA_lib[item]
                    contents.append([mass,item])     
                         
        filename = "PermLibGen_Lipid_{}.csv".format(frag_type)
        pathname = os.path.join(self.ion_lib_path,filename)
        write_csv(pathname,contents)
        return os.path.join(self.ion_lib_path,filename)

    ###Fragmentation Analysis###
    def calc_mass_from_dict(self,hit_dict):
        mass = 0
        for item in hit_dict:
            if item in self.FA_lib:
                mass += self.FA_lib[item]
            elif item in self.HG_NL_lib:
                mass += self.HG_NL_lib[item]
            elif item in self.Frag_lib:
                mass += float(self.Frag_lib[item]["Mass"])
            else:
                raise Exception ("Unknown item {}".format(item))
        return mass

    def create_lipid_frag_rules_lib(self):
        self.HG_frag_rules_lib = load_fragmentHG_from_settings()###Headgroup Fragmentations###
        self.lipid_frag_rules_type_filter = load_fragmentexcl_from_settings()###Special Exclusions###
        ###Edit Fragmentation Rules Here###
        self.lipid_frag_rules_lib = {}
        for line in self.readSettingFile(self.settings["fragment_rules_filename"], main_setting_folder="Settings"):
            comment, HG1, HG2, analysed, penalised, *LipidA_frag = line
            lst = []
            for frag in LipidA_frag:
                if frag != "":
                    lst.append(frag.split(","))
            dct = {"Fragments Analysed": analysed.split(","),"Fragments Penalised": penalised.split(","),"Component Order":lst}
            self.lipid_frag_rules_lib[HG1,HG2] = dct
        
    def check_all_lipid_fragments (self,frag_type,hit_dict,lipid_type,fragment_name = None):
        '''Checks the rule library, remembers the fragment name and initiates fragment analysis through group/sequence'''
        self.lipid_frag_dict = {} #None --> there's no component True --> there's a component and its loss is found; False --> not found
        logger.info("Analysing fragmentation order for {} and {} fragments...".format(lipid_type,frag_type))
        
        self.fragment_lipid(frag_type,fragment_name,hit_dict,None)#check for fragments without any losses
        self.fragment_lipid(frag_type,fragment_name,hit_dict,"H2O")#check for a single degree of dehydration always
        if lipid_type in self.lipid_frag_rules_lib:
            rules = self.lipid_frag_rules_lib[lipid_type]["Component Order"]
            self.settings["frag_types_with_penalty"] = self.lipid_frag_rules_lib[lipid_type]["Fragments Penalised"]
                         
            for line in rules:
                fragment_name,hit_dict = self.process_lipid_fragmentation_rules_by_line(line,frag_type,hit_dict,lipid_type,fragment_name)
            lost_FA_sequence = self.lipid_frag_rules_type_filter.get(frag_type,[])#check the exclusions last
            if lost_FA_sequence != []:
                fragment_name,hit_dict = self.check_all_lipid_fragments_consecutive(frag_type,hit_dict,lipid_type,lost_FA_sequence,fragment_name)
        else:
            logger.info("Rules for this type ({}) are not confirmed yet, it will be ignored".format(lipid_type))
            self.lipid_frag_counted_area = 0
            
    def process_lipid_fragmentation_rules_by_line(self,line,frag_type,hit_dict,lipid_type,fragment_name):
        if line[0] == "group":
            lost_FA_group = [cpt for cpt in line[1:] if cpt not in self.lipid_frag_rules_type_filter.get(frag_type,[])]
            fragment_name,hit_dict = self.check_all_lipid_fragments_competitive(frag_type,hit_dict,lipid_type,lost_FA_group,fragment_name)
        elif line[0] == "seq":
            lost_FA_sequence = [cpt for cpt in line[1:] if cpt not in self.lipid_frag_rules_type_filter.get(frag_type,[])]
            fragment_name,hit_dict = self.check_all_lipid_fragments_consecutive(frag_type,hit_dict,lipid_type,lost_FA_sequence,fragment_name)
        elif line[0] == "branch":
            self.process_lipid_fragmentation_rules_by_line(line[1:],frag_type,hit_dict,lipid_type,fragment_name) #recurse function, do not modify fragment_name,hit_dict
        else:
            raise Exception("Incorrect fragment rules for lipid_type {}: {} and {}".format(lipid_type,line[0],line[1]))

        return fragment_name,hit_dict
        
    def check_all_lipid_fragments_consecutive (self,frag_type,hit_dict,lipid_type,lost_FA_sequence,fragment_name = None):
        logger.info("Consecutive losses of {}".format(lost_FA_sequence))
        for lost_FA in lost_FA_sequence:
            self.fragment_lipid(frag_type,fragment_name,hit_dict,lost_FA)
            status = self.lipid_frag_dict.get(lost_FA,None)
            if status == None: #if component absent, continue
                continue
            else: #remove component and continue
                fragment_name,hit_dict = self.remove_component_from_lipid_fragment(hit_dict,lost_FA,fragment_name)[0]###only final pathway given   
        return fragment_name,hit_dict

    def check_all_lipid_fragments_competitive(self,frag_type,hit_dict,lipid_type,lost_FA_group,fragment_name = None):
        #check individually
        logger.info("Competitive individual losses of {}".format(lost_FA_group))
        lost_FA_group_missing = []
        
        for lost_FA in lost_FA_group:
            self.fragment_lipid(frag_type,fragment_name,hit_dict,lost_FA)
            FA_bool = self.lipid_frag_dict.get(lost_FA,None)
            if FA_bool== None:#components not present
                lost_FA_group_missing.append(lost_FA)
        
        for lost_FA in lost_FA_group_missing:
            lost_FA_group.remove(lost_FA)#remove them so all combinations later are legitimate
    
        #combo
        logger.info("Combination losses of {}...".format(lost_FA_group))
        fragment_combinations_evaluated= []
        
        for no_of_combo_components in range(1,len(lost_FA_group)):#iterate through all possible combinations
            comb = list(itertools.combinations(lost_FA_group,no_of_combo_components+1))
            for combination in comb:
                hit_dict_mod  = copy.deepcopy(hit_dict) #create copies
                fragment_name_mod = fragment_name
                
                ###Re-arrangement###
                #we arrange components with multiple fragmentations to be the last if it is present
                if "C3' FA" in combination:
                    re_arrangement_index = combination.index("C3' FA")
                else:
                     re_arrangement_index = None                   
                    
                if re_arrangement_index != None:
                    combination = combination[:re_arrangement_index]+combination[re_arrangement_index+1:]+(combination[re_arrangement_index],)

                combination_list = [[fragment_name_mod,hit_dict_mod]]
                combination_list_mod = []
                for lost_FA in combination[:-1]: #remove each lost_FA in combination except for the last one
                    for fragment_name_mod,hit_dict_mod in combination_list:
                        combination_list_mod += self.remove_component_from_lipid_fragment(hit_dict_mod,lost_FA,fragment_name_mod,remove = False)
                    combination_list = copy.deepcopy(combination_list_mod)
                
                for fragment_name_mod,hit_dict_mod in combination_list:
                    lost_FA = combination[-1]
                    if [fragment_name_mod,lost_FA] not in fragment_combinations_evaluated:
                        self.fragment_lipid(frag_type,fragment_name_mod,hit_dict_mod,lost_FA)#submit last one as component to be removed
                        fragment_combinations_evaluated.append([fragment_name_mod,lost_FA])
                    else:
                        pass
                         
        #remove all
        for lost_FA in lost_FA_group:
            fragment_name,hit_dict = self.remove_component_from_lipid_fragment(hit_dict,lost_FA,fragment_name)[0] ###only final pathway given
        return fragment_name,hit_dict

    def remove_component_from_lipid_fragment(self,hit_dict,lost_FA,fragment_name,remove = True):
        ###remove = False: give intermediate pathways for headgroup; remove = True: give final pathway for headgroup
        output = [] #list of name,dict
        
        ###remove components without checking mass and score###
        if lost_FA == None:
            output.append([fragment_name,hit_dict])
        elif is_headgroup_component(lost_FA):
            
            for HG_fragment,loss_type in fragment_HG(hit_dict[lost_FA],self.HG_frag_rules_lib,remove):
                hit_dict_mod = copy.deepcopy(hit_dict)
                hit_dict_mod[lost_FA] = HG_fragment
                fragment_name_mod = "{} -{}({})".format(fragment_name,lost_FA,loss_type)
                output.append([fragment_name_mod,hit_dict_mod])
                
        elif is_primary_FA(lost_FA):
            hit_dict[lost_FA] =  "NIL"
            fragment_name += " -{}(acid)".format(lost_FA)
            output.append([fragment_name,hit_dict])
        elif is_secondary_FA(lost_FA):
             hit_dict[lost_FA] =  "NIL"
             parent_FA = get_parent(lost_FA)
             hit_dict[parent_FA] = convert_hydroxy_to_unsat(hit_dict[parent_FA]) #convert to blank + dehydrate parent
             fragment_name += " -{}(acid)".format(lost_FA)
             output.append([fragment_name,hit_dict])
        else:
            raise Exception ("Not a valid component {}".format(lost_FA))

        return output
        
    def fragment_lipid(self,frag_type,fragment_name,hit_dict,lost_FA):
        ###removes, calc mass and score
        if frag_type == "Di-glucosamine":
            save = True
        else:
            save = False
            
        hit_dict_mod = copy.deepcopy(hit_dict)
        if lost_FA == None:
            ###No Loss###
            loss_type = "fragment"
            hit_mass = self.calc_mass_from_dict(hit_dict.values())
            hit_name = ".".join(hit_dict.values())
            if frag_type != "Di-glucosamine":
                self.score_lipid_fragment(lost_FA,hit_name,hit_mass,hit_dict,save,loss_type,fragment_name,frag_type)

        elif lost_FA == "H2O":
            ###Dehydration###
            loss_type = "dehydration"
            hit_mass = self.calc_mass_from_dict(hit_dict.values())
            hit_mass -= float(self.HG_NL_lib["Water"])
            hit_name = ".".join(hit_dict.values())
            self.score_lipid_fragment("H2O",hit_name,hit_mass,hit_dict,save,loss_type,fragment_name,frag_type)

        elif "adduct" in lost_FA:
            ###Adducts###
            loss_type = "adduct"
            lost_FA =  lost_FA.replace("adduct","")
            hit_mass = self.calc_mass_from_dict(hit_dict.values())
            hit_mass += float(self.HG_NL_lib[lost_FA])
            hit_name = ".".join(hit_dict.values())
            hit_name += "+{}".format(lost_FA)

            self.score_lipid_fragment(lost_FA,hit_name,hit_mass,hit_dict,save,loss_type,fragment_name,frag_type,penalty = False)
            
        elif "loss" in lost_FA:
            ###Unusual Losses###
            loss_type = "alkane"
            lost_FA =  lost_FA.replace("loss","")
            hit_mass = self.calc_mass_from_dict(hit_dict.values())
            hit_mass -= float(self.HG_NL_lib[lost_FA])
            hit_name = ".".join(hit_dict.values())
            #hit_name += "-{}".format(lost_FA)

            self.score_lipid_fragment(lost_FA,hit_name,hit_mass,hit_dict,save,loss_type,fragment_name,frag_type,penalty = False)
            
        elif lost_FA not in hit_dict:
            self.lipid_frag_dict[lost_FA] = None
            logger.info("Component {} not present in <{}>".format(lost_FA,fragment_name))
            return

        ###Secondary Fatty Acids###   
        elif is_secondary_FA(lost_FA) and hit_dict[lost_FA] != "NIL":
            
            ###Acid Loss###
            hit_dict_mod[lost_FA] = "NIL"
            parent_FA = get_parent(lost_FA)
            hit_dict_mod[parent_FA] = convert_hydroxy_to_unsat(hit_dict[parent_FA]) #convert to blank + dehydrate parent
            hit_name = ".".join(hit_dict_mod.values())
            hit_mass_mod = self.calc_mass_from_dict(hit_dict_mod.values())
            self.score_lipid_fragment(lost_FA,hit_name,hit_mass_mod,hit_dict_mod,save,"acid",fragment_name,frag_type)

        ###Ester-Linked###       
        elif is_ester_linked_FA(lost_FA) and hit_dict[lost_FA] not in ["NIL","Water"]: #check for acid and ketene loss
            
            ###Ketene Loss###
            if lost_FA == "C3' FA":
                hit_dict_mod[lost_FA] = "Water" #convert to Water
                hit_mass_mod = self.calc_mass_from_dict(hit_dict_mod.values())
                hit_name = ".".join(hit_dict_mod.values())         
                self.score_lipid_fragment(lost_FA,hit_name,hit_mass_mod,hit_dict_mod,save,"ketene",fragment_name,frag_type)
                
            ###Acid Loss###
            hit_dict_mod[lost_FA] = "NIL" #convert to blank
            hit_mass_mod = self.calc_mass_from_dict(hit_dict_mod.values())
            hit_name = ".".join(hit_dict_mod.values())
            self.score_lipid_fragment(lost_FA,hit_name,hit_mass_mod,hit_dict_mod,save,"acid",fragment_name,frag_type)
            
        ###Amide-Linked###       
        elif is_amide_linked_FA(lost_FA) and hit_dict[lost_FA] not in ["NIL","Water"]: 
            
            ###Ketene Loss###
            hit_dict_mod[lost_FA] = "NIL" #convert to blank
            hit_mass_mod = self.calc_mass_from_dict(hit_dict_mod.values())
            hit_name = ".".join(hit_dict_mod.values())         
            self.score_lipid_fragment(lost_FA,hit_name,hit_mass_mod,hit_dict_mod,save,"ketene",fragment_name,frag_type)
                
        ###Headgroups###        
        elif is_headgroup_component(lost_FA) and hit_dict[lost_FA] not in ["NIL","Water"]:
            hit_dict_mod = copy.deepcopy(hit_dict)
            ###Headgroup Loss###
            for HG_fragment,loss_type in fragment_HG(hit_dict[lost_FA],self.HG_frag_rules_lib,remove = False):
                hit_dict_mod[lost_FA] = HG_fragment
                hit_mass_mod = self.calc_mass_from_dict(hit_dict_mod.values())
                hit_name = ".".join(hit_dict_mod.values())
                self.score_lipid_fragment(lost_FA,hit_name,hit_mass_mod,hit_dict_mod,save,loss_type,fragment_name,frag_type)
        
        else:
            self.lipid_frag_dict[lost_FA] = None
            logger.info("Component {} not present in <{}>".format(lost_FA,fragment_name))
            
    def score_lipid_fragment(self,lost_FA,hit_name,hit_mass,hit_dict,save,loss_type,fragment_name,frag_type,penalty = True):
        if lost_FA not in self.lipid_frag_dict:
            self.lipid_frag_dict[lost_FA] =  False
          
        matched_ions,total_matching_area = self.ion_data.compare_ions_get_total_area(self.settings["fragments_mode"],self.settings["fragments_threshold"],hit_name,hit_mass,save)
        rounded_mass = round(hit_mass,5)
        self.lipid_frag_added_area = 0
        
        if total_matching_area > 0:
            logger.info(hit_name)
            self.lipid_frag_added_area = total_matching_area
            self.lipid_frag_dict[lost_FA] = True
            double_count = False
            
            for ion in matched_ions:#avoid double counting
                if ion not in self.lipid_frag_counted_peaks_list:
                    self.lipid_frag_counted_peaks_list.append(ion)
                else:
                    double_count = True
                    
            if double_count == False: 
                self.lipid_frag_counted_peaks += 1
                self.lipid_frag_counted_area += self.lipid_frag_added_area
                percent = round(100*self.lipid_frag_added_area/self.lipid_frag_total_area,2)
                logger.info( "FOUND: #{}-<{}># as #{}# loss at mass {} for {} %".format(fragment_name, lost_FA, loss_type, rounded_mass, percent))
                self.current_msp.add_peak(rounded_mass,percent,fragment_name,lost_FA,loss_type)           
            else:
                logger.info( "FOUND: #{}-<{}># as #{}# loss at mass {} for 0 %".format(fragment_name,lost_FA,loss_type,rounded_mass))
                self.current_msp.add_peak(rounded_mass,0,fragment_name,lost_FA,loss_type)
        else:
            logger.info(hit_name)
            logger.info( "NOT FOUND: {}-<{}> as {} loss at mass {}".format(fragment_name, lost_FA, loss_type, rounded_mass))
            if penalty:
                self.current_penalty_data.append({
                    "fragment_name" : fragment_name,
                    "lost_FA" : lost_FA,
                    "loss_type" : loss_type,
                    "rounded_mass" : rounded_mass,
                    "hit_no_of_FA" : self.hit_no_of_FA,
                    "penalty" : None,
                    "backbone_fragment_penalty_multipler" : None,
                })
                self.penalise_lipid_fragment(lost_FA, frag_type)

    def penalise_lipid_fragment(self,lost_FA,frag_type):
        if frag_type in self.settings["frag_types_with_penalty"] and is_FA(lost_FA) and lost_FA != "C2' Sec FA": #no penalty for headgroups/C2' Sec FA
            self.lipid_frag_total_penalty_increment = self.settings["missing_fragment_score_penalty"]/self.hit_no_of_FA
            self.current_penalty_data[-1]["penalty"] = self.settings["missing_fragment_score_penalty"]
            if frag_type != "Di-glucosamine":
                self.lipid_frag_total_penalty_increment *= self.settings["backbone_fragment_penalty_multiplier"] #penalty bonus for cross-ring fragments
                self.current_penalty_data[-1]["backbone_fragment_penalty_multipler"] = self.settings["backbone_fragment_penalty_multiplier"]
            self.lipid_frag_total_penalty += self.lipid_frag_total_penalty_increment
            logger.info("Penalty incurred: {} %; current penalty: {} %".format(round(self.lipid_frag_total_penalty_increment,2), round(self.lipid_frag_total_penalty, 2)))
            
    def find_unusual_neutral_losses_for_matched_ions(self,matched_ions,frag_type,fragment_name,hit_dict,lost_FA):
        ###Update Fragment Name##
        fragment_name,hit_dict = self.remove_component_from_lipid_fragment(hit_dict,lost_FA,fragment_name)[0]
        daughter_neutral_losses =[]
        for ion in matched_ions:
            daughter_neutral_losses+= self.ion_data.get_daughter_neutral_losses_for_ion(ion)
            
        daughter_neutral_losses = set(filter(is_alkane,daughter_neutral_losses))
        for alkane in daughter_neutral_losses:
            self.fragment_lipid(frag_type,fragment_name,hit_dict,alkane+"loss")

def filter_many(test_string,str_lst):
    if len(str_lst) > 1 and str_lst[0] == "OR":
        return any(map(lambda x: filter_string(test_string,x),str_lst[1:]))
    else:
        return all(map(lambda x: filter_string(test_string,x),str_lst))

def filter_string(test_string,string):
    if string[:2] != "--" and string in test_string:
        return True
    elif string[:2] == "--" and string[2:] not in test_string:
        return True
    else:
        return False
