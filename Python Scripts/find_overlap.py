

with open("/Users/berk/GrifoniSortedCD8_nonduplicates.csv", "r") as CD8:
    CD8epitopes = CD8.readlines()
    print(CD8epitopes[0])

with open("/Users/berk/GrifoniCD4EpitopeModified.csv", "r") as CD4:
    CD4epitopes = CD4.readlines()
    print(CD4epitopes[0])

x=0
with open("/Users/berk/Overlapping_peptides.csv", "w") as Overlapping_peptides:
    Overlapping_peptides.write("CD4, CD8")
    Overlapping_peptides.write("\n")

    for mhciiepitope in CD4epitopes[1:]:

        for mhciepitope in CD8epitopes:

            if mhciiepitope.split(",")[0] == mhciepitope.split(",")[0]:

                if mhciepitope.split(",")[4] in mhciiepitope.split(",")[2]:

                    if int(mhciiepitope.split(",")[-1]) >= 5:

                        print(mhciiepitope, mhciepitope)

                        Overlapping_peptides.write(mhciiepitope + mhciepitope + "\n")
                        x+=1

print(x)
