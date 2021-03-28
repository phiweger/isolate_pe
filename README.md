## README

### TODO

- https://github.com/najoshi/sickle
- antismash
- uv
- prokka
- plasmidfinder
- abricate
- ...


### Run

```bash
# Using data from SRA
nextflow run main.nf --sra PRJEB14100 --results results --minlen 2000 -resume

nextflow run main.nf --genomes input.csv --qc false --results results --minlen 2000 --db VFDB_setB_pro.fas -profile docker -resume
```

