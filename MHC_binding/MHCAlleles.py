with open("/Users/berk/MHCII_Alleles.txt", "r") as MHCii:
    MHCII_alleles = MHCii.read().split(",")

#MHCi_alleles = [HLA-A*68:01, HLA-A*03:01, HLA-A*02:06, HLA-A*02:01, HLA-B*35:01, HLA-B*51:01, HLA-B*35:01, HLA-B*07:02, HLA-A*24:02, HLA-B*57:01, HLA-A*11:01, HLA-A*03:01, HLA-A*68:01, HLA-B*07:02, HLA-B*07:02, HLA-B*08:01, HLA-B*57:01, HLA-A*68:01, HLA-B*35:01, HLA-B*57:01, HLA-A*26:01, HLA-B*15:01, HLA-A*68:01, HLA-A*01:01, HLA-B*57:01, HLA-B*35:01, HLA-A*68:01, HLA-B*08:01, HLA-B*57:01]

MHCII_alleles_unique = []
#MHCi_alleles_unique=[]

for i in MHCII_alleles:
    if i.strip() not in MHCII_alleles_unique:
        MHCII_alleles_unique.append(i.strip())

print(len(MHCII_alleles_unique))
print(MHCII_alleles_unique)

#for i in MHCi_alleles:
#    if i not in MHCi_alleles_unique:
#        MHCi_alleles_unique.append(i)

#print(len(MHCi_alleles_unique))
