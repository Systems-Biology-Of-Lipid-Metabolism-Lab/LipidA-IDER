###Imports###
import csv
import sys

csv.field_size_limit(min(sys.maxsize, 2147483646))

###Parse .txt file###
def extract_mass_from_txt(string):
    i = string.find("\t")
    return string[:i]

def extract_area_from_txt(string):
    i = string.find("\t")
    string2 = string[i+1:]
    j = string2.find("\t")
    return string2[:j]
  
###CSV Import and Writing###
def read_csv(csvfilename):
    """
    Reads a csv file and returns a list of list
    containing rows in the csv file and its entries.
    """
    rows = []

    with open(csvfilename) as csvfile:
        file_reader = csv.reader(csvfile)
        for row in file_reader:
            rows.append(row)
    return rows

def write_csv(csvfilename,contents):
    with open(csvfilename,'w', newline="") as file:
        writer = csv.writer(file)
        writer.writerows(contents)

def truncate_filename(filename):
  if len(filename) > 40:
    return filename[:10]+"..."+filename[-30:]
  else:
    return filename

def write_txt (txtfilename, contents):
    with open(txtfilename,"w") as file:
        file.writelines(contents)
        file.close()

###Mapping CSV to var###
def split_to_list(string):
    return list(filter(lambda x: x!="",string.split(",")))




