from LipidAIDER.Setup import *
import logging

logger = logging.getLogger(__name__)

###Ion Data###
#Set-up
Ion_Report = {"Parents": {}, "Daughters": {}, "Fragment Ions": {"No of Hits":0} }

class Ion_Data:
    def __init__ (self, name, scan, data_path, dummy_ion_mz = 123456):
        self.name = name
        self.ion_dict = {}
        self.ion_raw = {}
        self.readme_file = {} #Saving of parameters for writing reports
        self.neutral_loss = [] #This is a list of all neutral losses... misses and hits
        self.reset_neutral_loss() #This is a list of hits only, resetting it brings it back to default headers
        self.reset_fragment_ions() #Setting the headers
        self.os_path = os.path.join(data_path, self.name)
        os.makedirs(self.os_path)
        self.dummy_ion_mz = dummy_ion_mz #dummy peak area value
        for fragment in scan["fragments"]:
            self.ion_dict[fragment[0]] = copy.deepcopy(Ion_Report)
            self.ion_dict[fragment[0]]["Area"] = fragment[1]
            self.ion_raw[fragment[0]] = copy.deepcopy(Ion_Report)
            self.ion_raw[fragment[0]]["Area"] = fragment[1]
        self.retention_time = scan["retention_time"] * 60
        parent = pick_parent_ion(scan["fragments"], scan["parent_mz"])
        # if the parent fragment is not found, create a dummy one
        if parent not in self.ion_dict:
                self.ion_dict[parent] = copy.deepcopy(Ion_Report)
                self.ion_dict[parent]["Area"] = self.dummy_ion_mz #dummy value
                self.ion_dict[parent]["Nature"] = "Substitute"
        self.parent_ion_mass = parent
        self.ion_list = list(self.ion_dict)
        self.ion_list.sort(reverse = True) # descending order
        self.ion_raw_list = list(self.ion_raw)
        self.ion_raw_list.sort(reverse = True) # descending order

    def adjust_threshold(self,threshold):
        #Adjust threshold by a max of 10 if there are a lot of isotopic ions
        prev_mass = None
        threshold_new = threshold
        for ion in self.ion_list:
            if prev_mass == None:
                prev_mass = ion
            else:
                if calc_mDa(prev_mass,ion) < 200:
                    threshold_new += 1
                prev_mass = ion
        return min(threshold_new,threshold+10)
          
    def get_parent_ion_top_flag(self,threshold):
        parent_ion = self.get_parent_ion()
        self.ion_raw_list.sort(reverse = True, key = self.get_raw_area)
        adjusted_threshold = self.adjust_threshold(threshold)
        approved = len(self.ion_raw_list[:adjusted_threshold])
        removed = len(self.ion_raw_list[adjusted_threshold:])
        removed_list = self.ion_raw_list[adjusted_threshold:]
        logger.info("log_type='XUELI' /".join(map(str,self.ion_raw_list[:adjusted_threshold])))
        if parent_ion in self.ion_raw_list[:adjusted_threshold]:
            return(1)
        else:
            return(0)
             
        ###XUELI COMMENT 0: IF POSSIBLE INSTEAD OF CHOOSING TOP IONS, WE NEED TO CHANGE THIS TO FIND MINIMAL AREA/INTENSITY, THEN OBTAIN LIST WITH ALL SIGNALS THAT IS MINIMAL AREA/ INTENSITY +1
    def filter_ion_data(self, mode, threshold):
        if mode == "Top Ions":
            parent_ion = self.get_parent_ion()
            self.ion_list.sort(reverse = True, key = self.get_area)
            adjusted_threshold = self.adjust_threshold(threshold)
            approved = len(self.ion_list[:adjusted_threshold])
            removed = len(self.ion_list[adjusted_threshold:])
            removed_list = self.ion_list[adjusted_threshold:]

            #We increase sensitivity for the 0 -260 m/z range to pick up HG anions
            supplement_low_mass = [x for x in self.ion_list if x< 260]
            supplement_low_mass.sort(reverse = True, key = self.get_area)
            supplement_no = 0

            for ion in supplement_low_mass:
                if ion in removed_list:
                    removed_list.remove(ion)
                    supplement_no +=1
                    approved +=1
                    removed -= 1                 
                if supplement_no == 10:#max out at 10
                    break
                
            for ion in removed_list:
                if abs(ion-parent_ion) > 1:
                    del self.ion_dict[ion]
                else:
                    approved +=1
                    removed -= 1
        logger.info("Filtered peak area data at {} {} threshold; removed {} peaks and approved {} peaks".format(threshold, mode, removed, approved))
        self.ion_list = list(self.ion_dict)
        self.ion_list.sort(reverse = True)#sorted in descending order
        
    ###Getters###
    def get_parent_ion(self):
        if self.parent_ion_mass == None:
            precursor = self.ion_list[0]
            peaks = [(key,self.ion_dict[key]["Area"]) for key in self.ion_dict if calc_mDa(precursor,key)<1000]
            parent = pick_parent_ion(peaks,precursor)
            logger.info("Parent ion m/z: {}".format(round(parent,3)))
            if parent not in self.ion_dict: #TODO Jason when would this happen
                self.ion_dict[parent] = copy.deepcopy(Ion_Report)
                self.ion_dict[parent]["Area"] = self.dummy_ion_mz #dummy value
                self.ion_dict[parent]["Nature"] = "Substitute"
                logger.info("Precursor ion used as substitute for parent ion, peak area set to {}".format(self.dummy_ion_mz))
            self.parent_ion_mass = parent
            return self.parent_ion_mass
        else:
            return self.parent_ion_mass
    ##XUELI COMMENT 1: THESE PARTS ARE WHERE NEUTRAL LOSS ARE CALCULATED. TWO FUNCTIONS NEEDED. 1) ONLY NEED TO CALULATE FOR IONS BETWEEN PARENT M/Z and M/Z 250
      #2) QUALITY CHECK. 
 #if parent mz loss of FA present? (ie a hit for FA loss from the parent)
 #If parent mz loss of HG present? (ie a hit needed for HG loss from the parent)
 #if 1) absent but 2) present > ok
 #if 1) present but 2) absent > ok
 #But if both absent > quality of spectra is not perfect. Will have to penalize cos it likely is a lipid A with L2 ID, but can't be a L3 ID. 

    def get_neutral_loss_from_parent_fragment(self,parent_mass):
        neutral_losses = []
        for daughter_data in self.ion_dict[parent_mass]["Daughters"].values():
            NL = daughter_data[2]
            if is_acid(NL) or is_ketene (NL):
                FA_CC = get_carbon_count(NL)
                if FA_CC not in neutral_losses: neutral_losses.append(FA_CC)
        return neutral_losses
    
    def get_daughter_neutral_losses_for_ion(self,ion):
        ##Index = ["Daughters,"Neutral Loss", "Hit","Hit Name", "Diff.","Mode"]
        ion_dict_daughters = self.ion_dict[ion]["Daughters"]
        if ion_dict_daughters != {}:
            return [ion_dict_daughters[key][2] for key in ion_dict_daughters] ###returns list of neutral losses by name
        else:
            return []
        
    def get_parent_neutral_losses_for_ion(self,ion):
        ##Index = ["Daughters,"Neutral Loss", "Hit","Hit Name", "Diff.","Mode"]
        ion_dict_parents = self.ion_dict[ion]["Parents"]
        if  ion_dict_parents != {}:
            return  [ion_dict_parents[key][2] for key in  ion_dict_parents] ###returns list of neutral losses by name
        else:
            return []
    ##XUELI COMMENT 1: THESE PARTS ARE WHERE NEUTRAL LOSS ARE CALCULATED.END. 
    def get_path(self):
        return self.os_path

    def get_parent_ion_area(self):
        if self.ion_dict[self.parent_ion_mass]["Nature"] == "Substitute":
            return 0
        else:
            return self.get_area(self.parent_ion_mass)
        
    def get_area(self,mass):
        return self.ion_dict[mass]["Area"]
        
    def get_raw_area(self,mass):
        return self.ion_raw[mass]["Area"]
     
    def get_total_area_of_fragments(self):
        return sum(map(lambda mass: self.get_area(mass),
                    filter(lambda mass: mass > 260 and abs(mass-self.get_parent_ion())>17.9,
                           self.ion_dict.keys())))
    
    def get_total_area(self):
        return sum(map(lambda mass: 0 if self.get_area(mass) == self.dummy_ion_mz else self.get_area(mass),self.ion_dict.keys()))

    def get_fragment_area_ratio(self):
        if self.get_total_area() > 0:
            return self.get_total_area_of_fragments()/self.get_total_area()
        else:
            return 0
    
    def get_retention_time(self):
        return self.retention_time
    
    ###Neutral Loss###
    def reset_neutral_loss(self):
        self.neutral_loss_report = [["Parent Ion","Daughter Ion","Neutral Loss", "Hit","Hit Name", "Diff."]]
        
    def calc_neutral_loss (self,cap):
        if len(self.neutral_loss) == 0:
            length = len(self.ion_list)
            parent_daughter_combinations = list(itertools.combinations(self.ion_list,2))
            for combination in parent_daughter_combinations:
                parent = max(combination)
                daughter = min(combination)
                diff = parent-daughter
                if diff< cap and diff> 15: 
                    self.neutral_loss.append ([diff,parent,daughter])
        else:
            logger.info("Neutral Loss Calculation Complete")
        
    def compare_neutral_loss (self,mode,threshold,library_name,library_dict): #library_dict is in format {name:mz_string}
        self.reset_neutral_loss() #remove prior data
        comparison_string = str(library_name) + " at " + str(threshold) + mode
        self.readme_file["Neutral Loss"] = comparison_string
        logger.info("Comparing with {}".format(comparison_string))
        library = list(library_dict.keys())
        
        #Comparing with Library and Writing Data Internally
        for neutral_loss,parent,daughter in self.neutral_loss:
            hit_list = binary_search_shallow_array(library,neutral_loss,threshold,mode)
            #format of hit_list => [ ( [lib mz,lib name,lib formula,""], ppm/mDa ) ]
            if len(hit_list)>0:
                mass,difference= hit_list[0]
                name = library_dict[mass]
                self.ion_dict[parent]["Daughters"][daughter] = [neutral_loss,mass,name,str(round(difference,3)),mode]
                self.ion_dict[daughter] ["Parents"][parent] = [neutral_loss,mass,name,str(round(difference,3)),mode]
                self.neutral_loss_report.append([parent,daughter,neutral_loss,mass,name,str(round(difference,3)),mode])
        
    ###Ion###
    def reset_fragment_ions(self):
        self.fragment_ion_report = [["Ion", "Hit","Hit Name", "Diff.","Mode"]]

    def compare_ions_get_total_area(self,mode,threshold,fragment_name,fragment_mass,save):
        total_area = 0
        matched_ions = []
        hit_list = binary_search_shallow_array(self.ion_list,fragment_mass,threshold,mode)
        if len(hit_list) > 0:
            for ion,diff in hit_list:
                total_area+= self.get_area(ion)
                matched_ions.append(ion)
                if save:
                    self.ion_dict[ion]["Fragment Ions"]["No of Hits"] += 1 #we have multiple hits, we must use a different key for each
                    self.ion_dict[ion]["Fragment Ions"][fragment_name] = [round(fragment_mass,6),round(diff,3),mode]
                    self.fragment_ion_report.append([ion,fragment_mass,fragment_name,round(diff,3),mode])
            
        return matched_ions,total_area
        
    def compare_ions(self, mode, threshold, library, library_name, parent_only=False):
        self.reset_fragment_ions() # remove prior data
        comparison_string = "{} at {} {}".format(library_name,threshold,mode)
        self.readme_file["Fragment Ions"] = comparison_string
        logger.info("Comparing with {}".format(comparison_string))
        if parent_only :
            selection = [self.get_parent_ion()]
        else:
            selection = self.ion_list
        for ion in selection:
            hit_list = binary_search_deep_array(library,ion,threshold,mode)
            if len(hit_list) >0:
                for hit,diff in hit_list:
                    fragment_mass,fragment_name = hit[0],hit[1]
                    self.ion_dict[ion]["Fragment Ions"]["No of Hits"] += 1 #we have multiple hits, we must use a different key for each
                    self.ion_dict[ion]["Fragment Ions"][fragment_name] = [str(round(fragment_mass,6)),str(round(diff,3)),mode]
                    self.fragment_ion_report.append([ion,fragment_mass,fragment_name,str(round(diff,3)),mode])
        
    #Removing Data
    def del_fragment_from_ion(self,ion,fragment_name):
        del self.ion_dict[ion]["Fragment Ions"][fragment_name]
        
    #Exporting Data
    ###Machine Readable###
    def write_ion_list(self):
        contents = [["m/z","Area"]]
        for ion in self.ion_list:
            if "Area" in self.ion_dict[ion]:
                contents.append([ion,self.ion_dict[ion]["Area"]])
            else: contents.append([ion])
                
        filename = self.name+"_Ion_List.csv"
        filename = os.path.join(self.os_path,filename)        
        write_csv(filename,contents)
        return contents
    
    #Jason:     
    def write_neutral_loss(self):
        if "Neutral Loss" not in self.readme_file:
            raise Exception ("Neutral Loss incomplete")      
        filename = self.name+"_NL_"+ self.readme_file["Neutral Loss"] + ".csv"
        filename = os.path.join(self.os_path,filename)
        write_csv(filename,self.neutral_loss_report)
        return self.neutral_loss_report

    def write_fragment_ions(self):
        if "Fragment Ions" not in self.readme_file:
            raise Exception ("Fragment Ions incomplete")
        filename = self.name+"_FI_"+ self.readme_file["Fragment Ions"]  + ".csv"
        filename = os.path.join(self.os_path,filename)
        write_csv(filename,self.fragment_ion_report)
        return self.fragment_ion_report
        
     #Ion_Report Styles   
    def write_ion_report_neutral_loss (self):
        if "Neutral Loss" not in self.readme_file:
            raise Exception ("Neutral Loss incomplete")
        
        contents = []
        parameters = self.readme_file["Neutral Loss"]
        contents.append (["Ion Report - Neutral Loss", parameters])
        for ion in sorted(self.ion_dict.keys(), reverse=True):
            ion_dict = self.ion_dict[ion]
            contents.append([])
            contents.append(["Ion:",ion])
            if "Area" in ion_dict:
                area = ion_dict["Area"]
                contents.append(["Area",area])
            
            if ion_dict["Parents"] != {}:
                contents.append(["Parents:","Neutral Loss", "Hit","Hit Name", "Diff.","Mode"])
                for key,val in ion_dict["Parents"].items():
                    contents.append([key]+val)
                    
            if ion_dict["Daughters"] != {}:
                contents.append(["Daughters:","Neutral Loss", "Hit","Hit Name", "Diff.","Mode"])   
                for key,val in ion_dict["Daughters"].items():
                    contents.append([key]+val)
        
        filename = self.name+"_IonReport_NL_"+ parameters + ".csv"
        filename = os.path.join(self.os_path,filename)
        write_csv(filename,contents)
        return filename

    def write_ion_report_fragment_ions (self):
        if "Fragment Ions" not in self.readme_file:
            raise Exception ("Fragment Ions incomplete")
        
        contents = []
        parameters = self.readme_file["Fragment Ions"]
        contents.append (["Ion Report - Fragment Ions", parameters])
        for ion in sorted(self.ion_dict.keys(), reverse=True):
            ion_dict = self.ion_dict[ion]
            contents.append([])
            
            if ion_dict["Fragment Ions"] != {}:
                contents.append(["Ion:",ion])
                if "Area" in ion_dict:
                    area = ion_dict["Area"]
                    contents.append(["Area",area])
                contents.append(["Fragment:"])
                for key,val in ion_dict["Fragment Ions"].items():
                   contents.append([val[0]]+[key]+val[1:])
                        
        filename = self.name+"_IonReport_Frag_Ions_"+ parameters + ".csv"
        filename = os.path.join(self.os_path,filename)
        write_csv(filename,contents)
#XUELI COMMENT1: THE NEUTRAL LOSS REPORT IS WHERE WE CAN TARGET TO FIND IF PARENT OR PARENT minus HEADGROUP ARE PRESENT.THIS IS WHERE THE REPORT IS BEING GENERATED. 
    def write_ion_report(self):
        if "Fragment Ions" not in self.readme_file:
            raise Exception ("Fragment Ions incomplete")
        if "Neutral Loss" not in self.readme_file:
            raise Exception ("Neutral Loss incomplete")

        contents = []
        NL_parameters = self.readme_file["Neutral Loss"]
        contents.append (["Ion Report - Neutral Loss"])
        contents.append ([NL_parameters])
        F_parameters = self.readme_file["Fragment Ions"]
        contents.append (["Ion Report - Fragment Ions"])
        contents.append ([F_parameters])

        for ion in sorted(self.ion_dict.keys(), reverse=True):
            ion_dict = self.ion_dict[ion]
            contents.append([])
            contents.append(["Ion:",ion])
            if "Area" in ion_dict:
                area = ion_dict["Area"]
                contents.append(["Area",area])
            
            if ion_dict["Parents"] != {}:
                contents.append(["Parents:","Neutral Loss", "Hit","Hit Name", "Diff.","Mode"])
                for key,val in ion_dict["Parents"].items():
                    contents.append([key]+val)
                    
            if ion_dict["Daughters"] != {}:
                contents.append(["Daughters:","Neutral Loss", "Hit","Hit Name", "Diff.","Mode"])   
                for key,val in ion_dict["Daughters"].items():
                   contents.append([key]+val)

            if ion_dict["Fragment Ions"] ["No of Hits"] != 0:
              contents.append(["Fragments:","","Hits:","Hit Name","Diff.","Mode"])
              for key,val in ion_dict["Fragment Ions"].items():
                  if key != "No of Hits":
                    contents.append(["","",]+[val[0]]+[key]+val[1:])
                    
        filename = self.name+"_Ion_Report_Complete.csv"
        filename = os.path.join(self.os_path,filename)
        write_csv(filename,contents)

    def write_readme(self):
        contents = self.readme_file.items()
        filename = self.name+"_readme.csv"
        filename = os.path.join(self.os_path,filename)
        write_csv(filename,contents)
