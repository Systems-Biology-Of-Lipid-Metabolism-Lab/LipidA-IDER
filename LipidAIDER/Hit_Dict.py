###Hit_dict Operators for Components & Component Hits###
from LipidAIDER.File_Handling import *
import os
import logging

logger = logging.getLogger(__name__)

###Operators for Fatty Acids###
def is_unsaturated (FA_string):
    return type(FA_string) == str and ":1" in FA_string
def is_saturated (FA_string):
  return type(FA_string) == str and ":0" in FA_string
def is_ketene (FA_string):
  return type(FA_string) == str and "ketene" in FA_string
def is_acid (FA_string):
  return type(FA_string) == str and "acid" in FA_string
def is_hydroxy(FA_string):
    return type(FA_string) == str and "_(OH)" in FA_string
def is_saturated_unhydroxy(FA_string):
  return type(FA_string) == str and not is_hydroxy(FA_string) and is_saturated(FA_string)
def is_hydroxy_ketene(FA_string):
    return is_hydroxy(FA_string) and is_ketene(FA_string)
def is_unsat_ketene(FA_string):
     return is_unsaturated(FA_string) and is_ketene(FA_string)   
def is_sat_ketene(FA_string):
     return is_saturated(FA_string) and is_ketene(FA_string)
def is_blank(FA_string):
      return FA_string == "NIL"
    
def create_acid_from_carbon_count (count):
    return ("{}:0 acid".format(count))
def create_ketene_from_carbon_count (count):
    return convert_acid_to_ketene(create_acid_from_carbon_count (count))
def clean_name(string): #remove ketene, acid
    if type(string) == str:
        rems = [" ketene"," acid","Water","NIL"]
        reps = ["","","NIL","NIL"]
        for rem, rep in zip(rems,reps):
            string = string.replace(rem,rep)
    return string


def get_carbon_count (FA_string):
    if is_acid(FA_string) or is_ketene(FA_string):
        return FA_string[:3]#returns "CXX"
    elif is_alkane(FA_string):
        return FA_string[2:5]
    else:
        raise Exception("Carbon count error for {}".format(FA_string))
  
def get_carbon_count_number(FA_string):
    try:
        return int(get_carbon_count (FA_string)[1:])
    except Exception:
        logger.error("Unknown carbon count for {}, returned 0. log_type='Dev Notes'".format(FA_string))
        return 0

def is_even_carbon_count(FA_string):
    return float(get_carbon_count_number(FA_string))%2 == 0
def carbon_count_between(FA_string,min_count = 0,max_count = 999, ):
  carbon_count = int(get_carbon_count(FA_string)[1:])
  return carbon_count <= max_count and carbon_count >= min_count
def difference_in_carbon_count(FA_string1,FA_string2):
  return abs(get_carbon_count_number(FA_string1) - get_carbon_count_number(FA_string2))
    
def convert_unsat_to_sat (FA_string):
    return FA_string.replace(":1",":0")
def convert_sat_to_unsat (FA_string):
    return FA_string.replace(":0",":1")
def convert_unhydroxy_to_hydroxy(FA_string):
    if is_ketene(FA_string):
        return FA_string.replace(" ketene", "_(OH) ketene")
    else:
        return FA_string.replace(" acid", "_(OH) acid")
def convert_hydroxy_to_unhydroxy(FA_string):
    if is_ketene(FA_string):
        return FA_string.replace("_(OH) ketene"," ketene")
    else:
        return FA_string.replace( "_(OH) acid"," acid")   
def convert_unsat_to_hydroxy (FA_string):
    return convert_unhydroxy_to_hydroxy(convert_unsat_to_sat(FA_string))
def convert_hydroxy_to_unsat (FA_string):
    return convert_hydroxy_to_unhydroxy(convert_sat_to_unsat(FA_string))
def convert_ketene_to_acid (FA_string):
    return FA_string.replace("ketene","acid")
def convert_acid_to_ketene (FA_string):
    return FA_string.replace("acid","ketene")
def convert_ketene_acid_vv(FA_string):
    if is_ketene(FA_string):
        return convert_ketene_to_acid(FA_string)
    else:
        return convert_acid_to_ketene(FA_string)
def remove_ketene_acid(FA_string):
    return FA_string.replace(" acid","").replace(" ketene","")

###Operators for Lipid_dict / Hit_dict Components###
def get_parent(component):
    secondary = ("C3' Sec FA", "C2' Sec FA", "C2 Sec FA")
    primary = ("C3' FA", "C2' ketene", "C2 ketene")
    if component in secondary:
        return primary[secondary.index(component)]
    else:
        return None
def get_lipid_type(hit_dict):
  return (hit_dict["C4' Headgroup"],hit_dict["C1 Headgroup"])
def is_FA (component):
    return is_amide_linked_FA(component) or is_ester_linked_FA(component) \
           or is_secondary_FA(component)
def is_secondary_FA(component):
    return type(component) == str and "Sec" in component
def is_primary_FA(component):
  return type(component) == str and not is_secondary_FA(component) and not is_headgroup_component(component)
def is_sugar(component):
    return component == "Sugar"
def is_headgroup_component(component):
    return "Headgroup" in component
def is_amide_linked_FA(component):
    return type(component) == str and component in ("C2' ketene", "C2 ketene")
def is_ester_linked_FA(component):
    return type(component) == str and component in ("C3' FA", "C3 FA")
def is_very_labile(component):
    return component in ("C3' Sec FA", "C3 FA")
def is_component_hit_blank(component_hit):
    return component_hit in ("NIL","Water")
def is_component_hit_water(component_hit):
  return component_hit == "Water"

###Hit_dict checks###
def how_many_labile_components(hit_dict):
  n = 0    
  for component in hit_dict:
    if is_ester_linked_FA(component) or is_headgroup_component(component) or is_secondary_FA(component):
      if not is_component_hit_blank(hit_dict[component]):
        n+=1
  return n
    
def LpxA_check(hit_dict,threshold): #C3' and C3 FAs
  if "C3' FA" in hit_dict and "C3 FA" in hit_dict and \
     not is_component_hit_blank(hit_dict["C3' FA"]) and not is_component_hit_blank(hit_dict["C3 FA"]):
    boolean = difference_in_carbon_count(hit_dict["C3' FA"],hit_dict["C3 FA"]) <= threshold
    return boolean
  else:
    return True
  
def LpxD_check(hit_dict,threshold): #C2' and C2 ketenes
  if "C2' ketene" in hit_dict and "C2 ketene" in hit_dict and \
     not is_component_hit_blank(hit_dict["C2' ketene"]) and not is_component_hit_blank(hit_dict["C2 ketene"]):
    boolean = difference_in_carbon_count(hit_dict["C2' ketene"],hit_dict["C2 ketene"]) <= threshold
    return boolean
  else:
    return True

def parent_secondary_check(hit_dict):
  for component in hit_dict:
    if is_secondary_FA(component):
        parent_component = get_parent(component)
        parent_component_hit = hit_dict[parent_component]
        component_hit = hit_dict[component]
        if not is_component_hit_blank(component_hit): #secondary is not empty
          if is_component_hit_blank(parent_component_hit) or is_unsaturated(parent_component_hit):
            return False
  return True

def neutral_loss_check(hit_dict,NL_list):
  for NL_CC in NL_list: #Iterate through all neutral losses
    for component,component_hit in hit_dict.items(): #iterate through all components
      if NL_CC in component_hit and not is_amide_linked_FA(component):
        return True #stop when any NL is found; True
  return False

def check_hit_dict_for_FA(boolean,hit_dict,FA_lst):
  if boolean == False: ##don't check
    return True
  
  for component,*allowed_no in FA_lst:
    component_hit = hit_dict.get(component,"NIL")
    if is_component_hit_blank (component_hit) or get_carbon_count_number(component_hit) in allowed_no:
      continue
    else:
      return False
  return True

def compare_hit_dict(hit_dict1,hit_dict2):
    FA_list1,FA_list2 = [],[]
    for FA_list,hit_dict in zip([FA_list1,FA_list2],[hit_dict1,hit_dict2]):
        for component in hit_dict:
            if is_FA(component):
                try:
                    CC = get_carbon_count(hit_dict[component])
                    CC = remove_ketene_acid(hit_dict[component])
                except Exception as x:
                    CC = "nil"
                FA_list.append(CC)
        FA_list.sort()

    if FA_list1 == FA_list2:
        change = "Acyl Shuffle"
    else:
        acyl_common = []
        for item in FA_list1:
            if item in FA_list2:
                FA_list2.remove(item)
                acyl_common.append(item)
        acyl_change = len(FA_list1)-len(acyl_common)
        change = "{} Acyl Change".format(acyl_change)

    return change
            
###Operators for Headgroups###
def is_headgroup(HG_string):
  return not(is_acid(HG_string) or is_ketene(HG_string) or "$$" in HG_string)
def contains_P(HG_string):
  return is_monoP(HG_string) or is_pyrP(HG_string) or HG_string in ["Phosphate","Pyrophosphate","Triphosphate"]
def is_monoP (HG_string):
    return "P+" in HG_string
def is_pyrP (HG_string):
    return "Pyr+" in HG_string
def is_noP (HG_string):
    return is_monoP (HG_string) == False and is_pyrP (HG_string) == False \
       and HG_string != "Water"
def is_onlyP (HG_string):
    return convert_dehyd_to_hyd(HG_string) in ["Phosphate","Pyrophosphate","Triphosphate"]
def is_dehydrated (HG_string):
    return "*d" in HG_string
def is_modified_headgroup(HG_string):
  for i in ("EtN","GalN","Ara4N"):
    if i in HG_string:
      return True
  return False
def get_modification(HG_string):
  for i in ("EtN","GalN","Ara4N"):
    if i in HG_string:
      return i
  return False    

def add_P(P,HG_string):
  return P+"+"+HG_string
def remove_P(HG_string):
  return HG_string.replace("Pyr+","").replace("P+","")
  
def convert_P(end,HG_string):
  HG_string = remove_P(HG_string)
  return add_P(end,HG_string)
def convert_dehyd_to_hyd(HG_string):
    return HG_string.replace("*d","")
def convert_hyd_to_dehyd(HG_string):
    return HG_string+"*d"

def is_ion_from_HG(HG_string,ion_name):
    ion_name = convert_dehyd_to_hyd(ion_name)
    return any([ion_name == HG_string,#exact match
                is_onlyP(ion_name) == is_onlyP(HG_string) == True, #phosphates
                get_modification(ion_name) == get_modification (HG_string) != False, #share the same modification
                is_onlyP(ion_name) and is_modified_headgroup(HG_string)]) #phosphate from modified headgroup

###Get Fragment HG Settings###

def load_fragmentHG_from_settings(file = "Headgroup Fragmentation.csv"):
    HG_Frag_file = read_csv(os.path.join(os.getcwd(),"Settings",file))
    glob_HG_Frag_dict = {}

    for line in HG_Frag_file[1:]:
        commment, HG_string, *HG_frag = line
        lst = []
        for frag in HG_frag:
            if frag != "": lst.append(frag.split(","))
        glob_HG_Frag_dict[HG_string]=lst
        

    return glob_HG_Frag_dict

def fragment_HG(HG_string,dct,remove = True):
    #FALSE gives ALL fragmentation possibilities; #TRUE gives final fragment (last element)
    fragment = dct.get(HG_string,[["NIL","unknown"]])
    return [fragment[-1]] if remove else fragment

def load_fragmentexcl_from_settings(file="FragmentExceptions.csv"):
    LipidA_Frag_file = read_csv(os.path.join(os.getcwd(),"Settings",file))
    LipidA_Frag_dict = {}

    for line in LipidA_Frag_file[1:]:
        fragment,excl = line
        LipidA_Frag_dict[fragment] = excl.split(",")
    return LipidA_Frag_dict
