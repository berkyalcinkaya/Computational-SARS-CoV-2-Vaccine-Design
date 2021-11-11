#/Users/berk/Computational Vaccine Design/MHCIAnalysis.py
#5/18/2021
#In this script, the MHCI binders dataset from the Immune Epitope Dataset will be analyzed. First, the dataset will be cleaned up - duplicates removed and epitopes categorized by protein. Next, MHCI binding predictions will be run on the data and the top epitopes will be selected and subjected to further analysis
#Two datasets had to be utilized in this script. The first - MHCIbinders - contains sequence position data, which will be important to determine from which nonstructural protein the epitope is derived. The second csv file - IEDBMHCI contains only the epitopes that were recognized by two or more studies (analysis done in excel). The epitopes from this list will be referenced with the first list to determine thier start and stop sequences


##IMPORT MODULES
import requests #used to make calls to the IEDB API
import time




##SETUP FUNCTIONS
#the data downloaded contains numerous collumns that are irrelevant for this analysis - this function removes unneccessary data
def clean(list1): #removes unneeded information from the IEDB dataset
    desired_indexes = [0,2,5,6,11] #the indexes that are to be kept
    indexes = [i for i in range(28) if i not in desired_indexes] #indexes to be removed
    for index in sorted(indexes, reverse=True):
        del list1[index]
    return list1




#When SARS-CoV-2 expresses its genome, it first produces one long peptide called open reading frame 1ab (ORF1ab). This peptide is then cleaved into 16 nonstructural proteins by proteases also expressed by the virus. The Immune Epitope Database only shows ORF1ab as an antigen, so the function below assigns a  to all epitopes derived from the orf1ab
#Sequence data: NC_045512.2
def getprotein(list1):
    x = " "
    if int(list1[-2])<=180/3: #NSP1 is known as the leader sequence
        x = "NSP1"
    if int(list1[-2])>180 and int(list1[-2])<=818:
        x = "NSP2"
    if int(list1[-2])>818 and int(list1[-2])<=2763:
        x = "NSP3"
    if int(list1[-2])>2763 and int(list1[-2])<=3263:
        x = "NSP4"
    if int(list1[-2])>3263 and int(list1[-2])<=3569:
        x = "NSP5"
    if int(list1[-2])>3569 and int(list1[-2])<=3859:
        x = "NSP6"
    if int(list1[-2])>3859 and int(list1[-2])<=3942:
        x = "NSP7"
    if int(list1[-2])>3942 and int(list1[-2])<=4140:
        x = "NSP8"
    if int(list1[-2])>4140 and int(list1[-2])<=4253:
        x = "NSP9"
    if int(list1[-2])>4253 and int(list1[-2])<=4392:
        x = "NSP10"
    if int(list1[-2])>4392 and int(list1[-2])<=5324: #for whatever reason, NSP11 does not exist
        x = "NSP12"
    if int(list1[-2])>5324 and int(list1[-2])<=5925:
        x = "NSP13"
    if int(list1[-2])>5925 and int(list1[-2])<=6452:
        x = "NSP14"
    if int(list1[-2])>6452 and int(list1[-2])<=6798:
        x = "NSP15"
    if int(list1[-2])>6798 and int(list1[-2])<=7096:
        x = "NSP16"
    return x




#This function calls the IEDB API to access the Major Histocompatibility Complex I (MHCI) binding affinity prediction algorithm. API documentation can be found here: http://tools.iedb.org/main/tools-api/
HLAI_reference_alleles = "HLA-A*01:01,HLA-A*02:01,HLA-A*02:03,HLA-A*02:06,HLA-A*03:01,HLA-A*11:01,HLA-A*23:01,HLA-A*24:02,HLA-A*26:01,HLA-A*30:01,HLA-A*30:02,HLA-A*31:01,HLA-A*32:01,HLA-A*33:01,HLA-A*68:01,HLA-A*68:02,HLA-B*07:02,HLA-B*08:01,HLA-B*15:01,HLA-B*35:01,HLA-B*40:01,HLA-B*44:02,HLA-B*44:03,HLA-B*51:01,HLA-B*53:01,HLA-B*57:01,HLA-B*58:01"
#print(len(HLAI_reference_alleles.split(",")))
def getMHCi(sequence):
    alleles = []  #list containing each allele that is above the 1.0 percentile binding threshold
    sequence_percentile_ranks = [] #will contain percentile rank values for the alleles under 1 percentile
    for HLAI in HLAI_reference_alleles.split(','): #a call to the API must be made for each allele
        allele = HLAI
        print(allele)
        IEDB_params = {'method':"recommended", 'sequence_text':sequence, 'allele': allele, "length": len(sequence)}
        baseurl = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"
        response = requests.post(baseurl, data=IEDB_params) #make a post to the API
        content = response.text.split("\n")[1].split("\t") #the first line is the header. We only need the second line, hence the [1] index. Response header is 'allele', 'seq_num', 'start', 'end', 'length', 'peptide', 'core', 'icore', 'score', 'percentile_rank'
        print(content)
        percentile_rank = float(content[-1])
        if percentile_rank <= 1.0: #1.0 percentile is the binding percentile rank threshold determined by IEDB to represent a "strong binder". Only the MHC alleles in this top percentile will be considered
            #print(content)
            alleles.append(content[0]) #add the allele to the alleles string
            sequence_percentile_ranks.append(percentile_rank) #
        #else:
            #print("above threshold: " + str(content))

        time.sleep(5) # wait time in order not to overwhelm the API

    average_percentile_rank = sum(sequence_percentile_ranks)/len(sequence_percentile_ranks) #calculate the average value of the percentile ranks below 1 percentile
    return ",".join(alleles), average_percentile_rank #return the alleles and percentile rank as a tuple





##DATA ANALYSIS AND CLEANUP
#open the MHCIbinders file downloaded from IEDB. This file will be henceforth know as the reference file because it contains the sequence position data that will be combined with the epitope data.
with open("MHCIBinders.csv", "r") as MHCIreference:
    MHCIreference_lines = MHCIreference.readlines()

MHCIreference_modified = [] #initialize a list which will contain the "cleaned" data
for line in MHCIreference_lines[1:]:
    line_split = line.split(',')
    MHCIreference_modified.append(clean(line_split)) #append cleaned data
#print("reference set" + str(MHCIreference_modified[0])

#open the data analyzed in excel. This is the binders file which contains the preliminary list of epitopes to be considered for the vaccine. These epitopes were chosen because more than one study found them to be epitopes.
with open("IEDBMHCI.csv", "r") as MHCIbinders:
    MHCIbinders_lines = MHCIbinders.readlines()
    line1 = MHCIbinders_lines[0].strip().split(',') + ["Starting Position"] + ["Ending Position"]
    del MHCIbinders_lines[0] #remove header to make analysis easiers
    print("epitope data set:" + str(line1))
MHCIbinders_lines.sort(key = lambda line: (int(line.split(',')[-2]), int(line.split(',')[-1])), reverse=True) #sort based on number of references
MHCIbinders_modified = [line.strip().split(',') for line in MHCIbinders_lines] #turn each line into its own list, thereby making analysis and comparison simpler and more syntactically clear

#This code takes an epitope from the binders dataset and finds its start and stop position in the reference
for binder in MHCIbinders_modified:
    for reference in MHCIreference_modified:
        if binder[0] == reference[0]: #binder[0] and reference[0] refer to the epitope ID - the parameter used to match epitopes between datasets
            binder.append(reference[2]) #append the start position
            binder.append(reference[3]) #append the stop position

#using the getprotein function to determine to which NSP each orf1ab epitopes belongs
for binder in MHCIbinders_modified:
    antigen = binder[2]
    if antigen == "Replicase polyprotein 1ab": #only the antigens with from ORF1ab need to be further analyzed
        binder[2] = getprotein(binder) #change the antigen name to the correct NSP

##CALL THE MHCI API
for binder in MHCIbinders_modified:
    MHCIresponse = getMHCi(binder[1])
    if len(MHCIresponse[0]) == 0:
        MHCIbinders.remove(binder) #if the epitope has no MHC alleles that bind the MHC molecule with an affinity above the 1.0 percentile threshold, the epitope will be removed from the dataset
    else:
        binder.extend(MHCIresponse) #otherwise, we will record the alleles as well as the average percentile rank in this data

#write the data to a new file
MHCIbinders.sort(key = lambda row: (row[2], row[3]))
with open("MHCIEpitopes.csv", "w") as MHCIeptiopefile:
    MHCIepitopefile.write()
    for i in MHCIbinders:
        MHCIepitopefile.write(i)
        MHCIeptitopefile.write("\n")
print("file written")
