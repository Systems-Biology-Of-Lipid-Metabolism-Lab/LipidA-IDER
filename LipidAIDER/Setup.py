import math
import bisect
import traceback
import itertools
import heapq
import functools
import os
import datetime
import copy
import logging

###Import Own Modules###
from LipidAIDER.Hit_Dict import *
from LipidAIDER.File_Handling import *

logger = logging.getLogger(__name__)

###Math###
def sigmoid(x):
  intermediate = -1000*(x-0.001)
  return 1 / (1 + math.exp(intermediate))

def logit(x):
  return math.log(x/(1.000001-x))

def get_percentage(n,total,dp = 1):
  return round((100*n/total),dp)

def calc_ppm(value1,value2):
    difference = abs(value1-value2)*1000000
    ppm = difference / max(value1,value2)
    return ppm

def upper_limit_ppm(num,ppm):
    return (1+(ppm/1000000))*num

def lower_limit_ppm(num,ppm):
    return (1-(ppm/1000000))*num 

def binary_search_deep_array(array,num,threshold,mode):#array[][0] should be mz
    upper_limit = return_upper_limit_function_from_mode(mode)(num,threshold)
    lower_limit = return_lower_limit_function_from_mode(mode)(num,threshold)
    array = [[float(k[0])]+k[1:] for k in array]
    array.sort(key = lambda k:k[0])
    array_bisect = [k[0] for k in array]
    
    upper_index = bisect.bisect_right(array_bisect,upper_limit)
    lower_index = bisect.bisect_left(array_bisect,lower_limit)

    if lower_index == upper_index:
        return []
    else:
      hit_list = [(hit,return_comparison_function_from_mode(mode)(num,hit[0])) for hit in array[lower_index:upper_index] \
                  if return_comparison_function_from_mode(mode)(num,hit[0])<= threshold]
      hit_list.sort(key = lambda x: x[1])
      return hit_list

def binary_search_shallow_array(array,num,threshold,mode): #mzlist
    upper_limit = return_upper_limit_function_from_mode(mode)(num,threshold)
    lower_limit = return_lower_limit_function_from_mode(mode)(num,threshold)
    array.sort()
    upper_index = bisect.bisect_right(array,upper_limit)
    lower_index = bisect.bisect_left(array,lower_limit)

    if lower_index == upper_index:
        return []
    else:
      hit_list = [(hit,return_comparison_function_from_mode(mode)(num,hit)) for hit in array[lower_index:upper_index] \
                  if return_comparison_function_from_mode(mode)(num,hit)<= threshold]
      hit_list.sort(key = lambda x: x[1])
      return hit_list

def calc_mDa(value1,value2):
    return abs(value1-value2)*1000

def upper_limit_mDa(num,mDa):
  return num + mDa/1000

def lower_limit_mDa(num,mDa):
  return num - mDa/1000

def return_comparison_function_from_mode(mode):
      if mode == "ppm":
          return calc_ppm
      elif mode == "mDa":
          return calc_mDa
      else:
          raise ValueError("Incorrect mode used.")
        
def return_upper_limit_function_from_mode(mode):
      if mode == "ppm":
          return upper_limit_ppm
      elif mode == "mDa":
          return upper_limit_mDa
      else:
          raise ValueError("Incorrect mode used.")

def return_lower_limit_function_from_mode(mode):
      if mode == "ppm":
          return lower_limit_ppm
      elif mode == "mDa":
          return lower_limit_mDa
      else:
          raise ValueError("Incorrect mode used.")
        
def convert_string_to_float_smart(string):
  try:
    return float(string)
  except:
    return string

###File Directory###
os_path = os.getcwd()
os_data_path = os.path.join(os_path,"Data")
os_lib_path = os.path.join(os_path,"Libraries")
os_logger_path = os.path.join(os_path,"Logger")
os_source_path = os.path.join(os_path,"Source")
os_batch_path = os.path.join(os_logger_path, "Batch Output")

def create_directory(path):
  if not os.path.exists(path):
    os.makedirs(path)
    logger.info("Successfully created the directory {}".format(path))    
  else:
    logger.info("Directory {} already exists".format(path))    

#TODO Jason: this is where we want to comment out so all intermediate files are preserved
# def clear_directory(path,emphasis = False,delete_directory = False):
#   try:
#     for root, dirs, files in os.walk(path):
#       no_files = len(files)
#       for file in files:
#           os.remove(os.path.join(root, file))
#     if delete_directory:
#       os.rmdir(path)
#   except OSError:
#       # logger.print ("Unable to clear directory {}".format(path),emphasis)
#       logging.error("Unable to clear directory {}".format(path))
#   else:
#       # logger.print ("Successfully cleared the directory {} of {} files".format(path,no_files),emphasis)
#       logging.error("Successfully cleared the directory {} of {} files".format(path,no_files))

###Logger###
class Logger():
    def __init__(self,log_type_exceptions = ["Dev Notes","Ion Data","XUELI","Neutral Loss","Anions","Neutral Loss Filter","Frag. Analysis","Lipid Constr.","Lipid Analysis"]):
        self.date = datetime.date.today().strftime("%Y%m%d")
        self.path = os.path.join(os_logger_path,self.date)
        self.batch_path = os.path.join(self.path,"Batch Output")
        self.msp_path = os.path.join(self.path,"msp")
        self.mute = False
        self.hold = False
        self.name = ""
        self.exceptions = []
        self.log_categories = ["General"]
        self.reset_log()
        self.held_prints = {}
        self.held_name = None
        time = datetime.datetime.now().strftime("%H:%M:%S")
        self.print("Logger initialised on {} at time {}".format(self.date,time),emphasis = True)
        self.create_logger_directory()
        self.create_log_directory()
        self.create_batch_directory()
        self.create_msp_directory()
        self.exclude_log_types(log_type_exceptions)
        
    def create_logger_directory(self):#Directory Creation
      try:
        os.mkdir(self.path)
      except:
        pass

    def create_log_directory(self):
      try:
        os.mkdir(os.path.join(self.path,"Log"))
      except:
        pass
      
    def create_batch_directory(self):
      try:
        os.mkdir(self.batch_path)
      except:
        pass

    def create_msp_directory(self):
      try:
        os.mkdir(self.msp_path)
      except:
        pass
      
    def set_name(self,name):
        self.name = name
        self.print("Logger name set to {}".format(name))

    def reset_log(self):
        self.log = [["Log"]]
   
    def print(self,string,log_type="General",emphasis=False,hold=False):
      if type(log_type)!= type("str"):
        log_type = "General"
      if hold or self.hold:
        self.held_prints[self.held_name].append((string,log_type,emphasis))
      else:
        time = datetime.datetime.now().strftime("%H:%M:%S")
        dict_string = str(string).replace("\n"," ")
        self.log.append([time,log_type,dict_string])
        if emphasis:
            string = "\n###"+string+"###"       
        if log_type not in self.exceptions and self.mute == False:
          print(string)
          
    def keep_upcoming_prints(self,name):
      self.hold = True
      self.held_name = name
      self.held_prints[name]=[]

    def release_upcoming_prints(self):
      self.hold = False
      self.held_name = None
      
    def release_held_prints(self,name):
      if name in self.held_prints:
        for string,log_type,emphasis in self.held_prints[name]:
          self.print(string,log_type,emphasis)
        del self.held_prints[name]

    def release_held_prints_as_list(self):
      prints = copy.deepcopy(self.held_prints)
      self.held_prints = []
      return prints
      
    def exclude_log_type(self,log_type):
        self.exceptions.append(log_type)

    def exclude_log_types(self,types_lst):
      for log_type in types_lst:
        self.exceptions.append(log_type)

    def remove_all_exceptions(self,log_type):
        self.exceptions = []

    def mute_logger(self):
        self.mute = True

    def unmute_logger(self):
        self.mute = False

    def write_log(self):
        time = datetime.datetime.now().strftime("%H_%M_%S_")
        self.print("Log file exported at {} Name: {}".format(time,self.name),emphasis = True)
        filename = os.path.join(self.path,"Log",time+self.name+ ".csv")
        contents = self.log
        write_csv(filename,contents)
        self.reset_log()

    def get_logger_path(self):
      return self.path

    def get_logger_batch_path(self):
      return self.batch_path

#TODO Jason not sure about initializing the log here..    
# logger = Logger()

###Frequency Dictionary###
def add_key_to_dict(dictionary,output_key,exception_list = None,freq = 1):
    '''Add input key to dictionary as output key if input key is not in the exception list.
    If output key is already present, increase the count by one'''
    if exception_list == None or output_key not in exception_list:
        if output_key not in dictionary:
            dictionary[output_key] = freq
        elif type(dictionary[output_key]) == int or type(dictionary[output_key]) == float  :
            dictionary[output_key] += freq
            dictionary[output_key] = round(dictionary[output_key],5)
        else:
            exception_string = "Tried to add {} key but {} is not an integer or float".format(output_key,freq)
            raise Exception (exception_string)

        if dictionary[output_key] == 0: ###edge case where freq == 0 after rounding
          dictionary.pop(output_key,None) 

###Mass Modifications###
def get_anion_mass(mass):
  return float(mass) - H_atom.get_mass()
def convert_mass(mass,charge): #converts to monocharged mass
  if charge == 1:
    return mass
  elif charge == 2:
    return mass*2+H_atom.get_mass()

def convert_2charge_mass(mass):
  return (mass-H_atom.get_mass())/2

###Level 2 ID### 
class L2D_Manager:
    def __init__(self, file_number=None):
        self.path = os.path.join(os.getcwd(), "Libraries", "Level 2 ID")
        self.library = {}
        if file_number == None:
            file_number = max(map(self.get_file_number_from_filename,(filter(lambda x: "Lib_Level_2_ID" in x,os.listdir(self.path)))))
        self.file_number = file_number  
        self.filename= os.path.join(self.path,self.create_filename_from_file_number(self.file_number))
        self.initial_library = read_csv(self.filename)
        self.introduction = self.initial_library[0]
        self.headers = self.initial_library[1]
        convert_to_float = ["HG mz","mz","Score"]
        convert_to_int = ["Acyl","C#","Sat#","DBE#","Hydroxyl#"]
        for item in self.initial_library[2:]:
            item_dict = dict(zip(self.headers,item))
            for key in convert_to_float:
                item_dict[key] = float(item_dict[key])
            for key in convert_to_int:
                item_dict[key] = int(item_dict[key])
            if item_dict["Score"] > 0: self.set_library (item_dict)
            else: continue
        self.string_conversion_dict = {
          "1": "mono", "2": "di", "3": "tri", "4": "tetra", "5": "penta", "6": "hexa", "7": "hepta", "8":"octo", #Number conversions
          "Water": "", "Phosphate": "phosphoryl", "Pyrophosphate": "phosphoryl,phosphoryl", "Triphosphate": "phosphoryl,phosphoryl,phosphoryl", #Phosphates
          "P": "phosphoryl","Pyr": "pyrophosphoryl", #Phosphorylation
          "EtN": "ethanolamine", "GalN": "aminogalactose", "Ara4N":"aminoarabinose" #saccharides + EtN
        }
        self.L2HG_conversion_dict = {
          "18.01": "NIL",
          "36.02": "FREE_OH",
          "115.99": "P", "195.95": "PP", "275.92": "PPP", "355.89": "PPPP",
          "240.00": "PPEtN","362.00": "PPEtNPEtN","318.96":"PPPEtN","441.97":"PPEtNPPEtN",
          "247.05":"PAra4N","378.10":"PAra4NAra4N","327.01":"PPAra4N","458.07":"PAra4NPAra4N",
          "518.09":"PGalNPGalN"
        }
        
    def get_file_number_from_filename(self,filename):
      try:
        return int(filename.replace("Lib_Level_2_ID_","").replace(".csv",""))
      except Exception:
        return 0
      
    def create_filename_from_file_number(self,file_number):
        return "Lib_Level_2_ID_{}.csv".format(file_number)

    def set_library(self,item_dict):
        mz = item_dict["mz"]
        self.library[mz] = item_dict
        
    def set_library_with_score_check(self,item_dict):
        mz = item_dict["mz"]
        new_lib_score = item_dict.get("Score",0)
        if new_lib_score == 0: ###Score = 0 can be misleading...
          return
        
        existing_lib = self.get_library(mz)
        if existing_lib == None or existing_lib["Score"] < new_lib_score:
            self.update_item_dict_with_datetime(item_dict)
            self.set_library(item_dict)
            logger.info("New Level 2 ID accepted for mz {} with score {}; shorthand: {}".format(mz,new_lib_score, self.get_l2d_string(item_dict)))
        else:
          pass

    def update_item_dict_with_datetime(self,item_dict):
        item_dict["Last Accessed Date"] = datetime.date.today().strftime("%Y%m%d")
        item_dict["Last Accessed Time"] = datetime.datetime.now().strftime("%H_%M_%S")

    def get_library(self,mz,ppm = 30, mode = "ppm"):
        hit_list = binary_search_shallow_array(list(self.library.keys()),mz,ppm,mode)
        if len(hit_list) == 0:
          return None
        else:
          lib_mz = hit_list[0][0]
          return self.library.get(lib_mz,None)

    def get_l2d_from_hit_dict(self,mz,hit_dict,source_filename,source_name,score,HG_NL_lib):
        acyl,hgmz,c,sat,dbe,hydroxyl = 0,0,0,0,0,0
        headgroups = []
        sugar = None
        for component in hit_dict:
            component_hit = hit_dict[component]
            if is_FA(component) and not is_component_hit_blank(component_hit):
                acyl += 1
                c += get_carbon_count_number(component_hit)
                if is_hydroxy(component_hit):
                    hydroxyl +=1
                elif is_saturated(component_hit):
                    sat += 1
                elif is_unsaturated(component_hit):
                    dbe += 1
                else:
                    raise Exception ("Unknown component hit {}".format(component_hit))
            elif is_headgroup_component(component):
                hgmz += HG_NL_lib[component_hit]
                headgroups.append(component_hit)
            elif is_sugar(component):
              sugar = component_hit

        hgmz = round(hgmz,2)
        l2d_name = self.get_l2d_name(headgroups,acyl,sugar,hgmz,c,dbe,hydroxyl)
        level_2_id = {"L2 Name": l2d_name,"mz": round(mz,3), "Acyl": acyl, "HG mz": hgmz, "C#":c, "Sat#":sat,\
                     "DBE#": dbe, "Hydroxyl#": hydroxyl, "Headgroups": headgroups, "Sugar": sugar,\
                      "Source Filename": source_filename, "Source Name": source_name, "Score": score}
        return level_2_id

    def get_l2d_string(self,level_2_id):
        return "{}/{}/{}/{}/{}".format(self.convert_hgmz_to_string(level_2_id["HG mz"]),level_2_id["Acyl"],level_2_id["C#"],\
                                      level_2_id["DBE#"],level_2_id["Hydroxyl#"])
        
    def get_l2d_name(self,headgroups,acyl,sugar,hgmz,c,dbe,hydroxyl):
        name = "Lipid A {}-acyl {} {} ({}/{}/{}/{})".format(
          self.convert_string(acyl),
          self.convert_headgroups_to_L3namestrings(headgroups),
          self.uncapitalise_first(sugar),
          self.convert_hgmz_to_string(hgmz),c,dbe,hydroxyl)
        return name

    def get_l3d_name_from_hit_dict(self,hit_dict):
      acyl = 0
      sugar = None
      for component in hit_dict:
          component_hit = hit_dict[component]
          if is_FA(component) and not is_component_hit_blank(component_hit):
              acyl += 1
          elif is_sugar(component):
            sugar = component_hit
      
      hg1_string = self.convert_headgroup_to_L3namestring("C4'",hit_dict["C4' Headgroup"])
      hg2_string = self.convert_headgroup_to_L3namestring("C1",hit_dict["C1 Headgroup"])

      name = "Lipid A {}-acyl {} {} {} [{}-{}/{}-{}/{}-{}-{}/{}-{}]".format(
        self.convert_string(acyl),
        hg1_string,
        hg2_string,
        self.uncapitalise_first(sugar),
        self.convert_headgroup_to_string(hit_dict["C4' Headgroup"]),
        *[self.convert_fatty_acid_to_string(hit_dict[x]) for x in
        ["C3' FA","C3' Sec FA","C2' ketene","C2' Sec FA","C3 FA","C2 ketene","C2 Sec FA"]],
        self.convert_headgroup_to_string(hit_dict["C1 Headgroup"]))
      return name       

    def convert_string(self,string):
      try:
        return self.string_conversion_dict[str(string)]
      except KeyError:
        return ""

    def convert_headgroup_to_L3namestring(self,C,headgroup): #single HG
      HG_string = self.convert_headgroups_to_L3namestrings([headgroup])
      if HG_string != "":
        return "{}-{}".format(C,HG_string)
      else:
        return ""
      
    def convert_headgroups_to_L3namestrings(self,headgroups):
        headgroups_string = ""
        temp_dict = {}
        for headgroup in headgroups:
          components = headgroup.split("+")
          headgroup_name = "-".join(map(self.convert_string,components))
          for part in headgroup_name.split(","):
            if part in temp_dict:
              temp_dict[part]+=1
            elif part != "":
              temp_dict[part]=1

        for part in sorted(temp_dict.keys()):
          headgroups_string += "{}-{}".format(self.convert_string(temp_dict[part]),part)
        return headgroups_string

    def convert_hgmz_to_string(self,hgmz):
      try:
        return self.L2HG_conversion_dict[str(hgmz)]
      except KeyError:
        return hgmz

    def convert_headgroup_to_string(self,HG):
      HG = HG.replace("+","")
      HG = HG.replace("Pyrophosphate","PP")
      HG = HG.replace("Pyr","PP")
      HG = HG.replace("Phosphate","P")
      HG = HG.replace("Water","NIL")
      return HG
    
    def convert_fatty_acid_to_string(self,FA):
      if is_component_hit_blank(FA):
        return "NIL"
      else:
        return remove_ketene_acid(FA)
      
    def print_l2d_lib(self,experiment_name = None):
        self.file_number +=1
        if experiment_name == None:
          self.filename = os.path.join(self.path,self.create_filename_from_file_number(self.file_number))
        else:
          self.filename = os.path.join(self.path,"Level_2_ID_{}.csv".format(experiment_name))

        contents = []
        contents.append(self.introduction)
        contents.append(self.headers)
        for mz in sorted(self.library.keys()):
            line = [self.library[mz][key] for key in self.headers]
            contents.append(line)
        write_csv(self.filename,contents)

    def uncapitalise_first(self,s):
      if len(s) == 0:
        return s
      else:
        return s[0].lower() + s[1:]
      
###Permutation Generator###
'''Takes in a .csv of known components (name and m/z) and spits out all possible permutations for use by Ion_Data object'''

def generate_catalogue(library_name,sub_folder):    
    lib_path = os.path.join(os_lib_path,sub_folder,library_name)
    library = read_csv(lib_path)[2:]
    #First we scan the list and catalogue the components as a dictionary of lists
    catalogue = {}
    for mass,name in library:
        if mass == "Component":
            current_key = name         
            catalogue[current_key] = []#add an empty list
        else:
            catalogue[current_key].append([float(mass),name])
    return catalogue
        
def generate_permutation_as_csv(library_name,filename,sub_folder = None, min_mass = 0, max_mass = 10000,frag_comp = None, comparison_func = None):
    '''Generation of permutations from the base file'''
    logger.info("Generating permutations from base file "+library_name)
    catalogue = generate_catalogue(library_name,sub_folder)
    components = list(catalogue.keys())
    components_no = len(components)
    components_string = ".".join(components)
    no = max(itertools.accumulate([len(catalogue[k]) for k in catalogue],lambda x,y:x*y))
    time_start = datetime.datetime.now()
  
    headers = [["Generated from " + library_name],["m/z","Name",components_string]]
    results = []

    def check_permutation(permutation):
      permutation = list(permutation)
      permutation_mass = sum(map(lambda x: x[0],permutation))
      if min_mass <= permutation_mass <= max_mass:
        hit_components = list(map(lambda x: x[1],permutation))
        hit_dict = dict(zip(frag_comp,hit_components))
        if comparison_func == None or comparison_func(hit_dict):
          permutation_name = ".".join(hit_components)
          return [permutation_mass,permutation_name]

    try:
      results = list(filter(
        lambda x: type(x) == list,
        map(check_permutation,
            itertools.product(*catalogue.values()))))
      
    except Exception as x:
        raise Exception ("Too many combinations {} {}".format(no,x))  

    time_end = datetime.datetime.now()
    time_taken = (time_end-time_start).total_seconds()
    logger.info("Evaluated {} combinations in {} seconds".format(no,time_taken))
    
    results.sort()
    filename = "Lib_{}_{}_{}.csv".format(truncate_filename(filename),components_no,len(results)) #Write .csv file for future use
    if sub_folder != None:
        file_path = os.path.join(os_lib_path,sub_folder,filename)
    else:
        file_path = os.path.join(os_lib_path,filename)
    write_csv(file_path,headers+results)
    return filename

###Parent Ion Picking Algorithm###
def pick_parent_ion(peak_list,precursor_peak):
###Each peak is represented as a tuple (m/z,peak area)
  peak_list_org = []
  peak_list.sort(reverse = True, key = lambda x:x[1])
  for i,(peak,area) in enumerate(peak_list):
    ppm = calc_ppm(peak,precursor_peak)
    if i != len(peak_list)-1:
      area_dist = area/peak_list[i+1][1]
    else:
      area_dist = 0
    peak_list_org.append((peak,area,ppm,area_dist))#sorted by descending area

  for peak,area,ppm,area_dist in peak_list_org:
    #if ppm <= 30 and area_dist > 5:
    if ppm <= 30:
      return peak
    else: break #only look at the largest
      
  # for peak,area,ppm,area_dist in peak_list_org:
  #   if area_dist > 5:
  #     return peak
  #   else: break #only look at the largest

  return precursor_peak

###msp###
class MSPfile:
    def __init__ (self,name,precursor_mz,precursor_type,rt = 0, comment = "", L2D = ""):
        self.name = name + " " + precursor_type
        self.precursor_mz = float(precursor_mz)
        self.precursor_type = precursor_type
        self.rt = rt
        self.comment = comment
        self.L2D = L2D

        self.peak_list= {} #mass, area, fragment_name
        self.num_peaks = 0
        self.fragment_types = set()

    def get_precursor_mz(self):
      return self.precursor_mz

    def get_fragment_types(self):
      return list(self.fragment_types)

    def get_name(self):
      return self.name

    def get_L2D(self):
      return self.L2D
    
    def add_peak (self,mass,area,fragment_name,lost_FA = None,loss_type = None):
        if lost_FA != None and loss_type != None: ### Handle Lipid A fragmentation output
            fragment_name = fragment_name + " -{}({})".format(lost_FA,loss_type)
            if "[M-H]" not in fragment_name:
              if "-" in fragment_name:
                fragment_type = fragment_name.split(" -")[0]
                self.fragment_types.add(fragment_type)
              else:
                 self.fragment_types.add(fragment_name)
                
        if mass in self.peak_list:
          if fragment_name not in self.peak_list[mass]["name"]:
            self.peak_list[mass]["name"].append(fragment_name)
        else:
            self.peak_list[mass] = {"area":area,"name":[fragment_name]}
            self.num_peaks += 1

    def add_comment (self,comment):
      self.comment += "||{}".format(comment)

    def add_L3_name (self,name):
      self.L3_name = name
      
    def print_msp (self):
        contents = []
        contents.append('NAME: {}\n'.format(self.L3_name))
        contents.append('PRECURSORMZ: {}\n'.format(self.precursor_mz))
        contents.append('PRECURSORTYPE: {}\n'.format(self.precursor_type))
        contents.append('RETENTIONTIME: {}\n'.format(self.rt))
        contents.append('Comment: {}\n'.format(self.comment))
        contents.append('Num Peaks: {}\n'.format(self.num_peaks))
        for mass in sorted(self.peak_list.keys(),reverse=True):
            area = round(self.peak_list[mass]["area"]*100) #scale by 100
            fragment_name = ' || '.join(self.peak_list[mass]["name"]).replace(" ","")
            contents.append("{} {} \"{}\"\n".format(
              mass,
              area,
              fragment_name
              ))
        #contents.append("\n")
        return contents
