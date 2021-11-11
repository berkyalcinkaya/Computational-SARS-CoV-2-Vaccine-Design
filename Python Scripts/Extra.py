

while pos_end<=526:
    if break_var==True:
        break

    for i in results_lines[pos_start:pos_end]:
        MHCiresults = getMHCi(i.split(',')[4], "HLA-"+i.split(',')[1])

        if MHCiresults[-1] == '<html><head>':
            print(results_lines.index(i))
            break_var=True
        elif float(MHCiresults[-1])<1.0:
            epitopes_scored.append('{},{},{}'.format(i, MHCiresults[-2], MHCiresults[-1]))

    time.sleep(3)

    if pos_start==2:
        pos_start=pos_start+3
    else:
        pos_start=pos_start+5

    if pos_end == 520:
        pos_end = 526
    else:
        pos_end=pos_end+5





        if len(epitopes_scored)==len(results_lines[2:]):
            print(results_lines.index(i))
            break

        elif MHCiresults == "'', ''":
            print(results_lines.index(i))
            break

        elif float(MHCiresults[-1])<1.0:
            epitopes_scored.append('{},{},{}'.format(i, MHCiresults[-2], MHCiresults[-1]))
            filtered_epitopes.write('{},{},{}'.format(i.strip(), MHCiresults[-2], MHCiresults[-1]))
            filtered_epitopes.write("\n")
