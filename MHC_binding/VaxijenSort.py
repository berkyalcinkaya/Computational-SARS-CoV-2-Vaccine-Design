def sortCD8 (n, dict1):
    return dict1[n.split(",")[0]]


CD4threshold = 9
CD8threshold = 0.18

CD4_probable_antigens = []

with open("/Users/berk/GrifoniCD4EpitopeModified.csv", "r") as CD4_file:
    CD4_epitopes = CD4_file.readlines()
    line1_CD4 = CD4_epitopes[0]
    print(line1_CD4)

with open("/Users/berk/CD4VaxijenResults.txt", "r") as CD4VaxijenResults:
    CD4_antigenicity_scores = CD4VaxijenResults.read()

    for i in CD4_antigenicity_scores.split(";"):

        if i.strip().endswith("Antigen"):
            if int(CD4_epitopes[int(i.strip().split(",")[0].split("_")[-1])].split(",")[-1])>=CD4threshold:

                CD4_probable_antigens.append((CD4_epitopes[int(i.strip().split(",")[0].split("_")[-1])]).strip() + "," + (i.strip().split(",")[1]).strip())


print(len(CD4_probable_antigens))
CD4_probable_antigens.sort(key = lambda x: x.split(",")[-1], reverse = True)
print(CD4_probable_antigens)

CD4proteins={}
for i in CD4_probable_antigens:
    CD4protein=i.split(",")[0]

    if CD4protein not in CD4proteins:
        CD4proteins[CD4protein]=0

    CD4proteins[CD4protein] = CD4proteins[CD4protein]+1

print(CD4proteins)




CD8_probable_antigens = []

with open("/Users/berk/GrifoniSortedCD8_nonduplicates.csv", "r") as CD8_file:
    CD8_epitopes = CD8_file.readlines()
    line1_CD8 = CD8_epitopes[0]

with open("/Users/berk/CD8VaxijenResults.txt") as CD8VaxijenResults:
    CD8_antigenicity_scores = CD8VaxijenResults.read()

    for i in CD8_antigenicity_scores.split(";"):
        if i.strip().endswith("Antigen"):
            if float((CD8_epitopes[int(i.strip().split(",")[0].split("_")[-1])].split(",")[-1]).strip())<=CD8threshold:
                CD8_probable_antigens.append((CD8_epitopes[int(i.strip().split(",")[0].split("_")[-1])]).strip() + "," + (i.strip().split(",")[1]).strip())

print(len(CD8_probable_antigens))
print(CD8_probable_antigens)


CD8proteins={}
for i in CD8_probable_antigens:
    CD8protein=i.split(",")[0]

    if CD8protein not in CD8proteins:
        CD8proteins[CD8protein]=0

    CD8proteins[CD8protein] = CD8proteins[CD8protein]+1

print(CD8proteins)





CD8_probable_antigens.sort(key = lambda x: (sortCD8(x, CD8proteins), x.split(",")[-1]), reverse = True)


o = 0
for mhciiepitope in CD4_probable_antigens:

    for mhciepitope in CD8_probable_antigens:

        if mhciiepitope.split(",")[0] == mhciepitope.split(",")[0]:

            if mhciepitope.split(",")[4] in mhciiepitope.split(",")[2]:

                    #print(mhciiepitope)
                    #print(mhciepitope)
                    o+=1

print(o)






with open("/Users/berk/GrifoniCD4probableantigens.csv", "w") as CD4vaxijen:
    CD4vaxijen.write(line1_CD4.strip())
    CD4vaxijen.write(", Antigenicity Score")
    CD4vaxijen.write("\n")

    for i in CD4_probable_antigens:
        CD4vaxijen.write(i)
        CD4vaxijen.write("\n")


with open("/Users/berk/GrifoniCD8probableantigens.csv", "w") as CD8vaxijen:
    CD8vaxijen.write(line1_CD8.strip())
    CD8vaxijen.write(", Antigenicity Score")
    CD8vaxijen.write("\n")

    for i in CD8_probable_antigens:
        CD8vaxijen.write(i.strip())
        CD8vaxijen.write("\n")
