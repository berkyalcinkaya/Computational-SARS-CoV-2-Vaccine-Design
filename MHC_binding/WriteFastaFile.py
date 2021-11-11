promiscuity_threshold = 2

with open("/Users/berk/GrifoniCD4EpitopeModified.csv", "r") as CD4_file:
    CD4_epitopes = CD4_file.readlines()[1:]

x=0
with open("/Users/berk/CD4Vaxijen.txt", "w") as CD4fasta:
    for i in CD4_epitopes:
        row = i.split(",")

        if int(row[-1]) >= promiscuity_threshold:
            x = CD4_epitopes.index(i) + 1

            CD4fasta.write(">")
            CD4fasta.write(row[0]+"_")
            CD4fasta.write(str(x))
            CD4fasta.write("\n")
            CD4fasta.write(row[2])
            CD4fasta.write("\n"*2)




with open("/Users/berk/GrifoniSortedCD8_nonduplicates.csv", "r") as CD8_file:
    CD8_epitopes = CD8_file.readlines()[1:]

y=0
with open("/Users/berk/CD8Vaxijen.txt", "w") as CD8fasta:
    for i in CD8_epitopes:
        row = i.split(",")
        x = CD8_epitopes.index(i) + 1

        CD8fasta.write(">")
        CD8fasta.write(row[0]+"_")
        CD8fasta.write(str(x))
        CD8fasta.write("\n")
        CD8fasta.write(row[-2])
        CD8fasta.write("\n"*2)
