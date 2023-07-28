# fix-sex-chr.smk
This pipeline spawned from an issue with hifiasm/0.16.1 that it cannot resolve/distinguish sex chromosomes with non-trio male diploid samples in dual assembly mode. Hifiasm's developer isn't keen on fixing the algorithm to accommodate and therefore, this pipeline is a workaround to “phasing” the sex chromosome into their rightful XY haplotypes. The output is a modified fasta file and therefore it can be used with e.g. PAV where you'll get correct genotypes and SafFire with strict X and Y chromosomes (ad-hoc way to fix the contig switches).
## Step 1. Prepare directory
```
config
├── config.yaml
└── manifest.tab
runsnake
```
## Step 2. Prepare inputs
### config/config.yaml
```
reference: /path/to/CHM13/T2T/v2.0/T2T-CHM13v2.fasta

# OPTIONAL ARGS
minimap_params: '-x asm20 --secondary=no -s 25000'
manifest: config/manifest.tab
```
### config/manifest.tab
```bash
sample  hap1    hap2
your_sample path/to/hap1_asm    /path/to/hap2_asm
```
## Step 3. Start the analysis!
```bash
ln -s ./runsnake .
./runsnake 30
```