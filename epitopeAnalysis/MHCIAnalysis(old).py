#/Users/berk/Computational Vaccine Design/MHCIAnalysis.py
#5/18/2021
#In this script, the MHCI binders dataset from the Immune Epitope Dataset will be analyzed. First, the dataset will be cleaned up - duplicates removed and epitopes categorized by protein. Next, MHCI binding predictions will be run on the data and the top epitopes will be selected and subjected to further analysis


def clean(list1): #removes unneeded information from the IEDB dataset
    desired_indexes = [2,5,6,11] #the indexes that are to be kept
    indexes = [i for i in range(28) if i not in desired_indexes] #indexes to be removed
    for index in sorted(indexes, reverse=True):
        del list1[index]
    return list1





#open the MHCIbinders file downloaded from IEDB
with open("MHCIBinders.csv", "r") as MHCIbinders:
    MHCIbinders_lines = MHCIbinders.readlines()
    line1 = MHCIbinders_lines[1]
    #print(line1)

MHCIbinders_modified = [] #initialize a list which will contain the "cleaned" data
for line in MHCIbinders_lines[1:]:
    line_split = line.split(',')
    MHCIbinders_modified.append(clean(line_split)) #append cleaned data
MHCIbinders_modified[0].append("frequency")

MHCIbinders_nonduplicates = [MHCIbinders_modified[0], MHCIbinders_modified[1]]
MHCIbinders_nonduplicates[1].append(1)

for epitope1 in MHCIbinders_nonduplicates[1:]:
    for epitope2 in MHCIbinders_modified[1:]:
        if epitope2[0] != epitope1[0]:
            MHCIbinders_nonduplicates.append(epitope2+[1])
            print("append")
        else:
            epitope1[-1] = epitope1[-1] + 1
            print("match")

print(MHCIbinders_nonduplicates)
