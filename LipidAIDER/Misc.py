import logging
import sys

logger = logging.getLogger(__name__)

##############################################
# New functions added by Bernett
##############################################

# Method to sort a descending series of floats by the standard competition ranking (https://en.wikipedia.org/wiki/Ranking).
def rankByRevStdCompetition(scores):
    score_ranks = [1]
    score_offset = 0
    for i in range(1, len(scores)):
        if scores[i] < scores[i-1]:
            score_ranks.append(score_ranks[i-1] + 1 + score_offset)
            score_offset = 0
        elif scores[i] == scores[i-1]:
            score_ranks.append(score_ranks[i-1])
            score_offset += 1
        else:
            logging.error("Scores are not in descending order")
            sys.exit(1)
    return score_ranks

# function to compute the chemical formula from the L2 text
def computeL2ChemicalFormula(l2_text, HG_chem_comp, element_mass):
    # counts of the atoms
    counts = {
        "C" : 0,
        "H" : 0,
        "O" : 0,
        "N" : 0,
        "P" : 0,
    }
    fields = l2_text.split("/")
    if len(fields) == 5:
        # convert the counts to integers
        fields[1:] = [ int(x) for x in fields[1:] ]
        (HG, Acyl, C, DB, OH) = fields 
        # HG
        if HG in HG_chem_comp:
            for (k, v) in HG_chem_comp[HG].items():
                counts[k] += v
        else:
            logger.error(f"Unable to locate HG {HG} in HG chemical composition data")
            return ("Unable to compute chemical composition", 0.0)
        counts["C"] += C + 12
        counts["H"] += C * 2 - DB * 2 - Acyl * 2 + 24
        counts["O"] += Acyl * 2 + OH - Acyl + 9
        counts["N"] += 2
        o = ""
        mass = 0.0
        for x in ("C", "H", "O", "N", "P"):
            if counts[x] != 0:
                o += f"{x}{counts[x]}"
                mass += counts[x] * element_mass[x]
        mass -= element_mass["H"]
        return (o, mass)
    else:
        raise Exception("Unable to parse L2 string")

# function to read ms2 file
def readMs2File(path, threshold_ppm=30):
    scans = []
    with open(path) as infile:
        for line in infile:
            line = line.strip()
            if line.startswith("S\t"):
                fields = line.split("\t")
                scans.append({
                    "parent_mz" : float(fields[3]),
                    "retention_time" : None,
                    "parent_intensity" : None,
                    "fragments" : [],
                })
            elif line.startswith("I"):
                fields = line.split("\t")
                if fields[1] == "TIC":
                    scans[-1]["parent_intensity"] = float(fields[2])
                elif fields[1] == "RTime":
                    scans[-1]["retention_time"] = float(fields[2])
            elif line.startswith("D"):
                continue
            elif line.startswith("Z"):
                continue
            elif line.startswith("H"):
                continue
            else:
                fields = line.split(" ")
                scans[-1]["fragments"].append((float(fields[0]), float(fields[1])))
    # filter away for the fragments which are greater than the parent ion + a threshold (threshold_ppm)
    for scan in scans:
        threshold = (threshold_ppm/1E6 * scan["parent_mz"]) + scan["parent_mz"]
        nfragments = []
        for fragment in scan["fragments"]:
            if fragment[0] > threshold:
                logger.info(f"Fragment {fragment[0]:.3f} greater than threshold {threshold:.3f} for scan {scan['parent_mz']:.3f} m/z {scan['retention_time']:.3f} RT")
            else:
                nfragments.append(fragment)
        scan["fragments"] = nfragments
    return scans

# function to compute the distance where a < b
def calcDist(a, b):
    if isinstance(a, float) and isinstance(b, float):
        return b - a
    elif isinstance(a, float) and isinstance(b, list):
        return b[-1] - a
    elif isinstance(a, list) and isinstance(b, float):
        return b - a[0]
    elif isinstance(a, list) and isinstance(b, list):
        return b[-1] - a[0]
    else:
        raise Exception("The two entities has to be either float or list")

# function to compute the PPM distance where a < b
def calcPpmDist(a, b):
    if isinstance(a, float) and isinstance(b, float):
        return (b - a) / b * 1E6
    elif isinstance(a, float) and isinstance(b, list):
        return (b[-1] - a) / b[-1] * 1E6
    elif isinstance(a, list) and isinstance(b, float):
        return (b - a[0]) / b * 1E6
    elif isinstance(a, list) and isinstance(b, list):
        return (b[-1] - a[0]) / b[-1] * 1E6
    else:
        raise Exception("The two entities has to be either float or list")

# function to merge two entities
def mergeAB(a, b):
    if isinstance(a, float) and isinstance(b, float):
        return [a, b]
    elif isinstance(a, float) and isinstance(b, list):
        return [a] + b
    elif isinstance(a, list) and isinstance(b, float):
        return a + [b]
    elif isinstance(a, list) and isinstance(b, list):
        return a + b
    else:
        raise Exception("The two entities has to be either float or list")

# function to hclust a sorted list of values to return a list of nested list containing the clusters within a defined distance
def getHclustByDistanceClustersForSortedNumericalList(the_list, thresholdFunc=lambda x: x <= 5, distanceFunc=calcDist):
    while len(the_list) > 1:
        min_dist = None
        min_pos = None
        for i in range(len(the_list)-1):
            dist = distanceFunc(the_list[i], the_list[i+1])
            if min_dist == None or dist < min_dist:
                min_dist = dist
                min_pos = i
        # if the distance is below the threshold, merge the node/cluster
        if thresholdFunc(min_dist):
            x = mergeAB(the_list[min_pos], the_list[min_pos+1])
            del the_list[min_pos:min_pos+2]
            the_list.insert(min_pos, x)
        else:
            break
    return the_list

# function to read original CVS file input
def readScanCsvFile(filename, threshold_ppm=30):
    logging.info(f"Processing CSV file: {filename}")
    with open(filename) as infile:
        # process the data into list of list
        lines = infile.readlines()
        lines = [x.strip().split(",") for x in lines]
        # check the data
        if len(lines[0]) < 4:
            logging.error(
                f"Number of fields in the CSV file {filename} must be 4 or greater.")
            sys.exit(1)
        if lines[0][0] != "m/z":
            logging.error(
                f"The first field in the in the CSV file {filename} has to be m/z: {lines[0][0]}")
            sys.exit(1)
        if lines[0][1] != "Area":
            logging.error(
                f"The second field in the in the CSV file {filename} has to be Area: {lines[0][1]}")
            sys.exit(1)
        if not lines[0][2].startswith("Retention Time = "):
            logging.error(
                f"The third field in the in the CSV file {filename} has to be retention time: {lines[0][2]}")
            sys.exit(1)
        if not lines[0][3].startswith("Parent Ion = "):
            logging.error(
                f"The fourth field in the in the CSV file {filename} has to be parent ion m/z: {lines[0][3]}")
            sys.exit(1)
        scan = {
            "parent_mz": float(lines[0][3].replace("Parent Ion = ", "")),
            "retention_time": float(lines[0][2].replace("Retention Time = ", "")),
            "parent_intensity": 0,
            "fragments": [],
        }
        # threshold for filtering fragments
        threshold = (threshold_ppm/1E6 * scan["parent_mz"]) + scan["parent_mz"]
        for line in lines[1:]:
            if len(line) >= 2:
                fragment = (float(line[0]), float(line[1]))
                if fragment[0] > threshold:
                    logger.info(f"Fragment {fragment[0]:.3f} greater than threshold {threshold:.3f} for scan {scan['parent_mz']:.3f} m/z {scan['retention_time']:.3f} RT")
                else:
                    scan["fragments"].append(fragment)
        return scan

# write out a MS2 file
def writeMs2File(scans, filename):
    with open(filename, "w") as outfile:
        i = 1
        for scan in scans:
            outfile.write(f"S\t{i}\t{i}\t{scan['parent_mz']}\n")
            outfile.write(f"I\tRTime\t{scan['retention_time']}\n")
            outfile.write(f"I\tTIC\t{scan['parent_intensity']}\n")
            for fragment in scan["fragments"]:
                outfile.write(f"{fragment[0]} {fragment[1]}\n")
            i += 1

# read a series of MSP files
def readMspFiles(filenames):
    scans = []
    for filename in filenames:
        logging.info(f"Processing {filename}")
        with open(filename) as infile:
            scan = {
                "parent_mz": 0.0,
                "retention_time": 0.0,
                "parent_intensity": 0.0,
                "fragments": [],
            }
            fragments = False
            for line in infile:
                line = line.strip()
                if fragments:
                    fields = line.split()
                    scan["fragments"].append((float(fields[0]), float(fields[1])))
                elif line.startswith("RETENTIONTIME:"):
                    scan["retention_time"] = float(line.replace("RETENTIONTIME:", "").strip()) / 60.0
                elif line.startswith("PRECURSORMZ:"):
                    scan["parent_mz"] = float(line.replace("PRECURSORMZ:", "").strip())
                elif line.startswith("Num Peaks:"):
                    fragments = True
            scans.append(scan)
    return scans