import requests
import time


threshold = 20.0


HLAii_reference = "HLA-DRB1*01:01, HLA-DRB1*03:01, HLA-DRB1*04:01, HLA-DRB1*04:05, HLA-DRB1*07:01, HLA-DRB1*08:02, HLA-DRB1*09:01, HLA-DRB1*11:01, HLA-DRB1*12:01, HLA-DRB1*13:02, HLA-DRB1*15:01, HLA-DRB3*01:01, HLA-DRB3*02:02, HLA-DRB4*01:01, HLA-DRB5*01:01, HLA-DQA1*05:01/DQB1*02:01, HLA-DQA1*05:01/DQB1*03:01, HLA-DQA1*03:01/DQB1*03:02, HLA-DQA1*04:01/DQB1*04:02, HLA-DQA1*01:01/DQB1*05:01, HLA-DQA1*01:02/DQB1*06:02, HLA-DPA1*02:01/DPB1*01:01, HLA-DPA1*01:03/DPB1*02:01, HLA-DPA1*01:03/DPB1*04:01, HLA-DPA1*03:01/DPB1*04:02, HLA-DPA1*02:01/DPB1*05:01, HLA-DPA1*02:01/DPB1*14:01"



def getMHCii(sequence):
    IEDB_params = {'method':"recommended", 'sequence_text':sequence, 'allele':HLAii_reference}
    baseurl = "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"

    response = requests.post(baseurl, data=IEDB_params)

    content = response.text
    #content = response.text.split("\n")[1:]
    #content = response.text.split("\n")[1].split("\t")

    return content




def getAlleles(lst1):
    alleles={}

    for allele in lst1[1:-1]:
        #print(allele.split("\t")[5])
        if "Consensus" not in allele.split("\t")[5]:

            if float(allele.split("\t")[-5]) <= threshold:
                    alleles[allele.split("\t")[0]] = float(allele.split("\t")[-5])

        else:

            if float(allele.split("\t")[8]) <= threshold:
                    alleles[allele.split("\t")[0]] = float(allele.split("\t")[8])

    #print(alleles)
    return alleles




def getScore(alleles_dict):
    getAlleles_response = getAlleles(alleles_dict)

    sum_items = len(getAlleles_response.keys())

    return (list(getAlleles_response.keys()), sum_items)





with open("/Users/berk/GrifoniCD4Epitopes.csv","r" ) as results:
    results_lines=results.readlines()

    line_1=results_lines[1]
    print(line_1)
    #print(results_lines[2].split('"')[1])





proteins={}
protein_num=1
for i in results_lines[2:]:
    protein=i.split(",")[0]

    if protein not in proteins:
        proteins[protein]=protein_num
        protein_num = protein_num+1
print(proteins)
print("\n")



line_1_mhcii = getMHCii("CTFEYVSQPFLMDLE").split("\n")[0].split("\t")
print(line_1_mhcii)
print("\n")




sorted_epitopes=[]

for i in results_lines[2:]:

    time.sleep(7)

    if i.split(",")[3] != "HLA class II" and len(i.split('"'))==3:
        sect1=i.split('"')[0]
        sect3=i.split('"')[2]
        #print(sect1)
        #print(sect3)

        if int(sect3.split(',')[1])>=2:
            response = (getMHCii(sect1.split(",")[2])).split("\n")

            sorted_epitopes.append('{}{},{}'.format(sect1, (" ").join(getScore(response)[0]), getScore(response)[1]))

            print('{}{},{}'.format(sect1, (" ").join(getScore(response)[0]), getScore(response)[1]))
            print("\n")
            #sorted_epitopes.append('{},{},{}'.format(i.strip(), (" ").join(getScore(response)[0]), getScore(response)[1]))

            #scored_epitopes.write
            #print(response)
            #print(len(response))
            #print("\n")

    else:

        if int(i.split(',')[4])>=2:
            response = (getMHCii(i.split(",")[2])).split("\n")

            sorted_epitopes.append('{},{},{}'.format(",".join(i.split(",")[0:3]), (" ").join(getScore(response)[0]), getScore(response)[1]))

            print('{},{},{}'.format(",".join(i.split(",")[0:3]), (" ").join(getScore(response)[0]), getScore(response)[1]))
            print("\n")
            #sorted_epitopes.append('{},{},{}'.format(i.strip(), (" ").join(getScore(response)[0]), getScore(response)[1]))
            #print(response)
            #print(len(response))
            #print("\n")

    sorted_epitopes.sort(key = lambda row: (proteins[row.split(",")[0]], float(row.split(",")[-1])))

print(sorted_epitopes)
#sorted_epitopes.sort(key = lambda row: (proteins[row.split(",")[0]], float(row.split(",")[-1])))



with open("/Users/berk/GrifoniCD4EpitopeModified.csv", "w") as scored_epitopes:
    scored_epitopes.write("Protein, Start, Sequence, MHCii alleles, Num of binders")
    scored_epitopes.write("\n")

    for i in sorted_epitopes:
        scored_epitopes.write(i)
        scored_epitopes.write("\n")
