import sys

vcf_file=sys.argv[1]

bed=[]

with open(vcf_file, "r") as input:
     vcf=[i.rstrip().split("\t") for i in input if i.rstrip().split("\t")[0].startswith("chr")]

#for i in vcf[:10]:
#    print len(i[-2].split(":")), i[-2]

for i in vcf:
    if len(i[-2].split(":"))==5:
       bed.append((i[0], str(int(i[1])-1), i[1], i[-1].split(":")[2], i[-1].split(":")[0], i[3]+":"+i[4]+":"+i[-1].split(":")[1], i[5]+":"+i[6]))
    else:
       bed.append((i[0], str(int(i[1])-1), i[1], "0", i[-1].split(":")[0], i[3]+":"+i[4],  i[5]+":"+i[6]))

output_file=vcf_file.split(".")[0]+".bed"

with open(output_file,"w") as output:
     for i in bed:
         output.write("\t".join(i))
         output.write("\n")
