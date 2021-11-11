def sortepi(item,pos1, pos2=None):
    if pos2==None or type(pos2) != "int":
        return float(item.split(",")[pos1])
    else:
        return (float(item.split(",")[pos1]), float(item.split(",")[pos2]))

with open("/Users/berk/CD8Tepi.csv","r" ) as results:

    results_lines=results.readlines()
    #print(len(results_lines))
    line_1=results_lines[0]
    print(line_1)

    with open("/Users/berk/CD8Modified.csv", "w") as epitope_list:
        epitope_list.write("allele,seq_num,start,end,length,peptide,ic50,percentile_rank")
        epitope_list.write("\n")

        results_cutoff=[]
        for allele in results_lines[1:]:
            if float(allele.split(",")[-1])<=1.0 and float(allele.split(",")[-2])<=50.0:
                results_cutoff.append(allele.strip())
        sorted_results = sorted(results_cutoff, key=lambda x: (float(x.split(",")[1]), float(x.split(",")[-1])))


        nucleocapsid = sorted([epitope for epitope in results_cutoff if int(epitope.split(",")[1])==1], key=lambda x: (float(x.split(",")[-2]), float(x.split(",")[-1])))
        membrane = sorted([epitope for epitope in results_cutoff if int(epitope.split(",")[1])==2], key=lambda x: (float(x.split(",")[-2]), float(x.split(",")[-1])))
        envelope =  sorted([epitope for epitope in results_cutoff if int(epitope.split(",")[1])==3], key=lambda x: (float(x.split(",")[-2]), float(x.split(",")[-1])))
        orf8 = sorted([epitope for epitope in results_cutoff if int(epitope.split(",")[1])==4], key=lambda x: (float(x.split(",")[-2]), float(x.split(",")[-1])))
        orf7a =  sorted([epitope for epitope in results_cutoff if int(epitope.split(",")[1])==5], key=lambda x: (float(x.split(",")[-2]), float(x.split(",")[-1])))
        orf6 = sorted([epitope for epitope in results_cutoff if int(epitope.split(",")[1])==6], key=lambda x: (float(x.split(",")[-2]), float(x.split(",")[-1])))
        orf3a = sorted([epitope for epitope in results_cutoff if int(epitope.split(",")[1])==7], key=lambda x: (float(x.split(",")[-2]), float(x.split(",")[-1])))
        nsp3 = sorted([epitope for epitope in results_cutoff if int(epitope.split(",")[1])==8], key=lambda x: (float(x.split(",")[-2]), float(x.split(",")[-1])))
        nsp4 = sorted([epitope for epitope in results_cutoff if int(epitope.split(",")[1])==9], key=lambda x: (float(x.split(",")[-2]), float(x.split(",")[-1])))
        nsp6 = sorted([epitope for epitope in results_cutoff if int(epitope.split(",")[1])==10], key=lambda x: (float(x.split(",")[-2]), float(x.split(",")[-1])))
        nsp12 = sorted([epitope for epitope in results_cutoff if int(epitope.split(",")[1])==11], key=lambda x: (float(x.split(",")[-2]), float(x.split(",")[-1])))

        epitope_list.write("nucleocapsid")
        epitope_list.write("\n")
        for i in nucleocapsid[:5]:
            epitope_list.write(i)
            epitope_list.write("\n")

        epitope_list.write("membrane")
        epitope_list.write("\n")
        for i in membrane[:5]:
            epitope_list.write(i)
            epitope_list.write("\n")

        epitope_list.write("envelope")
        epitope_list.write("\n")
        for i in envelope[:5]:
            epitope_list.write(i)
            epitope_list.write("\n")

        epitope_list.write("orf8")
        epitope_list.write("\n")
        for i in orf8[:5]:
            epitope_list.write(i)
            epitope_list.write("\n")

        epitope_list.write("orf7a")
        epitope_list.write("\n")
        for i in orf7a[:5]:
            epitope_list.write(i)
            epitope_list.write("\n")

        epitope_list.write("orf6")
        epitope_list.write("\n")
        for i in orf6[:5]:
            epitope_list.write(i)
            epitope_list.write("\n")

        epitope_list.write("orf3a")
        epitope_list.write("\n")
        for i in orf3a[:5]:
            epitope_list.write(i)
            epitope_list.write("\n")

        epitope_list.write("nsp3")
        epitope_list.write("\n")
        for i in nsp3[:5]:
            epitope_list.write(i)
            epitope_list.write("\n")

        epitope_list.write("nsp4")
        epitope_list.write("\n")
        for i in nsp4[:5]:
            epitope_list.write(i)
            epitope_list.write("\n")

        epitope_list.write("nsp6")
        epitope_list.write("\n")
        for i in nsp6[:5]:
            epitope_list.write(i)
            epitope_list.write("\n")

        epitope_list.write("nsp12")
        epitope_list.write("\n")
        for i in nsp12[:5]:
            epitope_list.write(i)
            epitope_list.write("\n")
