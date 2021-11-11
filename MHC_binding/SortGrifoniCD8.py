percentile_threshold = 1.0




def getaverage(list1,pos):
    sum_items=0.0
    total_items=0.0

    for line in list1:
        sum_items = sum_items + float(line.split(",")[pos])
        total_items=total_items+1

    return (sum_items/total_items)





def length(epitope1, epitope2):
    if int(epitope1.split(",")[3]) >= int(epitope2.split(",")[3]):
        r1 = epitope1
        r2 = epitope2.split(",")[1]
    else:
        r1 = epitope2
        r2 = epitope1.split(",")[1]

    if epitope1.split(",")[1] == epitope2.split(",")[1]:
        r2 = ""

    if epitope1.split(",")[-1] <= epitope2.split(",")[-1]:
        r3 = epitope1.split(",")[-1]
    else:
        r3 = epitope2.split(",")[-1]

    return (r1, r2, r3.strip())






with open("/Users/berk/GrifoniCD8EpitopeModified.csv", "r") as filtered_epitopes:
    results_lines=filtered_epitopes.readlines()
    line_1=results_lines[0]
    print(line_1)


proteins={}
protein_num=1
for i in results_lines[1:]:
    protein=i.split(",")[0]

    if protein not in proteins:
        proteins[protein]=protein_num
        protein_num = protein_num+1


epitopes_scored = list(filter(lambda i: float(i.split(",")[-1])<percentile_threshold, results_lines[2:]))


avg_percentilerank = getaverage(epitopes_scored[1:], -1)
print(avg_percentilerank)


epitopes_scored.sort(key = lambda row: (proteins[row.split(",")[0]], float(row.split(",")[-1])))


epitopes_scored_nonduplicates = []
epitopes_scored_nonduplicates_updated = []

for x1 in epitopes_scored:
    i1=x1.split(",")

    for x2 in epitopes_scored:
        i2=x2.split(",")

        if x2 != x1:

            if i1[0] == i2[0] and (i1[4] in i2[4] or i2[4] in i1[4]):

                if x1 not in epitopes_scored_nonduplicates and x2 not in epitopes_scored_nonduplicates:
                    print(x1, x2)
                    epitopes_scored_nonduplicates.append(length(x1, x2)[0])
                    epitopes_scored_nonduplicates_updated.append('{},{},{},{}'.format(length(x1, x2)[0].split(",")[0], (length(x1, x2)[0].split(",")[1] + ' ' + length(x1, x2)[1]), ",".join(length(x1, x2)[0].split(",")[2:5]), length(x1, x2)[2]))


for x1 in epitopes_scored:
    sequence = x1.split(",")[4]
    no_append = False

    for x2 in epitopes_scored_nonduplicates_updated:
        sequence2 = x2.split(",")[4]
        if sequence in sequence2 or sequence2 in sequence:
            no_append = True

    if no_append == False:
        epitopes_scored_nonduplicates_updated.append('{},{}'.format(",".join(x1.split(",")[0:5]), x1.split(",")[-1]))


epitopes_scored_nonduplicates_updated.sort(key = lambda row: (proteins[row.split(",")[0]], float(row.split(",")[-1])))


print(epitopes_scored_nonduplicates_updated)




with open("/Users/berk/GrifoniSortedCD8.csv", "w") as sorted_epitopes:
    sorted_epitopes.write(line_1)
    sorted_epitopes.write("\n")

    for i in epitopes_scored:
        sorted_epitopes.write(i)




with open("/Users/berk/GrifoniSortedCD8_nonduplicates.csv", "w") as nonduplicates_updated:
    nonduplicates_updated.write("protein, alleles, start, length, sequence, percentile")
    nonduplicates_updated.write("\n")

    for i in epitopes_scored_nonduplicates_updated:
        nonduplicates_updated.write(i.strip())
        print(i.strip())
        nonduplicates_updated.write("\n")
