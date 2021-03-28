import re

import pandas as pd
import screed


fp = '/Users/phi/data_local/strain_security/VFDB_setB_pro.fas'

pattern = '(VF.*?) \((.*?)\) .*? \[(.*?)\] \[(.*)\]'
# VFG037218(gb|YP_001847229) (basJ) acinetobactin biosynthesis protein BasJ [Acinetobactin (VF0467)] [Acinetobacter baumannii ACICU]

d = {}
with screed.open(fp) as file:
    for i in file:
        r = re.match(pattern, i.name)
        vf, name, description, species = [r.group(j) for j in range(1, 5)]
        d[vf] = (name, description, species)


aln = '/Users/phi/tmp/leclercia/results/aln.m8'

df = pd.read_csv(aln, sep='\t', header=None)
for _, i in df.iterrows():
    print(d[i[1]])



