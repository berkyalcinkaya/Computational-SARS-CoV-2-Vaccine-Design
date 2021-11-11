import requests
import time

def getMHCi(sequence,allele):
    IEDB_params = {'method':"recommended", 'sequence_text':sequence, 'allele':allele, "length":len(sequence)}
    baseurl = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"

    response = requests.post(baseurl, data=IEDB_params)

    #content = response.text.split("\n")
    content = response.text.split("\n")[1].split("\t")
    return content


with open("/Users/berk/GrifoniCD8Epitopes.csv","r" ) as results:
    results_lines=results.readlines()
    line_1=results_lines[1]
    print(line_1)
    #print(len(results_lines))


proteins={}
protein_num=1
for i in results_lines[2:]:
    protein=i.split(",")[0]

    if protein not in proteins:
        proteins[protein]=protein_num
        protein_num = protein_num+1


#print(getMHCi(results_lines[2].split(',')[4], ("HLA-"+results_lines[2].split(',')[1])))
#['allele', 'seq_num', 'start', 'end', 'length', 'peptide', 'core', 'icore', 'score', 'percentile_rank']

epitopes_scored=[]
#epitopes_scored.append('{},{},{}'.format(line_1.strip(), "NetMHC score", "percentile rank"))

with open("/Users/berk/GrifoniCD8EpitopeModified.csv", "a") as filtered_epitopes:
    #filtered_epitopes.write('{},{},{}'.format(line_1.strip(), "NetMHC score", "percentile rank"))
    filtered_epitopes.write("\n")


    x=0
    for i in results_lines[299:]:

        MHCiresults = getMHCi(i.split(',')[4], "HLA-"+i.split(',')[1])

        epitopes_scored.append('{},{},{}'.format(i.strip(), MHCiresults[-2], MHCiresults[-1]))

        filtered_epitopes.write('{},{},{}'.format(i.strip(), MHCiresults[-2], MHCiresults[-1]))
        filtered_epitopes.write("\n")

        x=x+1
        print(x)
        print(i)

        time.sleep(10)



#epitopes_scored.sort(key=lambda x: proteins[x[0]] )




#this assinged some arbitrary number to each unique protein in this file, which will be useful to sort the final list by protein
proteins={}
protein_num=1
for i in results_lines[2:]:
    protein=i.split(",")[0]

    if protein not in proteins:
        proteins[protein]=protein_num
        protein_num = protein_num+1
