# A step-by-step guide showing how Chrom3D can be used to generate a full 3D genome model
**1.  Download the Hi-C data**
```bash
# Setup the folder to run the analyses
mkdir INC-tutorial
cd INC-tutorial

# Download Hi-C data:
mkdir -p fastq/all
wget -P fastq/all/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR665/000/SRR6657510/SRR6657510_1.fastq.gz
wget -P fastq/all/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR665/000/SRR6657510/SRR6657510_2.fastq.gz
wget -P fastq/all/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR665/001/SRR6657511/SRR6657511_1.fastq.gz
wget -P fastq/all/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR665/001/SRR6657511/SRR6657511_2.fastq.gz
wget -P fastq/all/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR665/002/SRR6657512/SRR6657512_1.fastq.gz
wget -P fastq/all/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR665/002/SRR6657512/SRR6657512_2.fastq.gz
wget -P fastq/all/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR665/003/SRR6657513/SRR6657513_1.fastq.gz
wget -P fastq/all/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR665/003/SRR6657513/SRR6657513_2.fastq.gz
wget -P fastq/all/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR665/004/SRR6657514/SRR6657514_1.fastq.gz
wget -P fastq/all/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR665/004/SRR6657514/SRR6657514_2.fastq.gz
```

**2. Download LAD data**
```bash
mkdir lad
wget -P lad/ https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109924/suppl/GSE109924_lad_D0-rep1.bed.gz
gunzip lad/GSE109924_lad_D0-rep1.bed.gz
```

**3. Setting up HiC-Pro**
```bash
wget https://raw.githubusercontent.com/nservant/HiC-Pro/master/config-hicpro.txt
wget https://raw.githubusercontent.com/nservant/HiC-Pro/master/annotation/chrom_hg19.sizes
wget https://raw.githubusercontent.com/nservant/HiC-Pro/master/annotation/HindIII_resfrag_hg19.bed
wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip
unzip hg19.zip -d hg19
```

**4. Downloading and installing the required processing scripts**
```bash
mkdir processing_scripts
wget https://github.com/Chrom3D/preprocess_scripts/archive/v.1.2.zip
unzip -j -d processing_scripts/ v.1.2.zip preprocess_scripts-v.1.2/NCHG_hic.zip
unzip -j -d processing_scripts/ v.1.2.zip preprocess_scripts-v.1.2/cap_chr_end.py
unzip -j -d processing_scripts/ v.1.2.zip preprocess_scripts-v.1.2/make_diploid_gtrack.py
unzip -j -d processing_scripts/ v.1.2.zip preprocess_scripts-v.1.2/NCHG_fdr_oddratio_calc.py
unzip -j -d processing_scripts/ v.1.2.zip preprocess_scripts-v.1.2/unmappable_blacklist.bed
unzip -j -d processing_scripts/ v.1.2.zip preprocess_scripts-v.1.2/makeGtrack.py
```

**5. Installing NCHG**
```bash
unzip -d processing_scripts/ processing_scripts/NCHG_hic.zip
make -C processing_scripts/NCHG_hic/
```

**6. Adapting the `config-hicpro.txt` file**
- Add the bowtie path to line nr. 39: `BOWTIE2_IDX_PATH =` -> `BOWTIE2_IDX_PATH = [fullpath]/INC-tutorial/hg19]/` where `[fullpath]` is the full path to your current working directory
- Change line nr. 89: `BIN_SIZE = 20000 40000 150000 500000 1000000` -> `BIN_SIZE = 50000 1000000`

**7. Run HiC-Pro (Takes several hours)**
```bash
HiC-Pro --input fastq --output hicpro_results --conf config-hicpro.txt
```

**8. Setup the folder structure for the HiC contacts**
```bash
mkdir -p hic/bedpe/intra
mkdir -p hic/bedpe/inter
mkdir -p hic/matrix
```

**9. Convert intrachromosomal Hi-C to BEDPE an matrix format**
```bash
# BEDPE:
awk 'NR==FNR { map[$4] = $1"\t"$2"\t"$3; next } { print $0,map[$1],map[$2] }' hicpro_results/hic_results/matrix/all/raw/50000/all_50000_abs.bed  hicpro_results/hic_results/matrix/all/raw/50000/all_50000.matrix  | awk '$4==$7' | awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$3>"hic/bedpe/intra/"$4}'

# Matrix format:
for chr in hic/bedpe/intra/*
do
chrname=$(basename $chr)
cut -f 2,5,7 $chr > hic/matrix/$chrname
done
```
**10. Convert interchromosomal Hi-C to BEDPE**
```bash
awk 'NR==FNR { map[$4] = $1"\t"$2"\t"$3; next } { print $0,map[$1],map[$2] }'  hicpro_results/hic_results/matrix/all/raw/1000000/all_1000000_abs.bed  hicpro_results/hic_results/matrix/all/raw/1000000/all_1000000.matrix  | awk '$4<$7' | awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$3>"hic/bedpe/inter/"$4"_"$7}'
```

**11. Running Armatus to call TADs**
```bash
mkdir hic/tads
for chr in hic/matrix/*
do
chrname=$(basename $chr)
armatus-linux-x64 -r 50000 -c $chrname -S -i $chr -g .6 -o hic/tads/$chrname 
done
```

**12. Since Armatus output is end-inclusive, convert to BED by adding 1 to end position**
```bash
ombine TADs per chromosome into a single BED file of TADs:
cat hic/tads/*.txt | awk '{printf("%s\t%i\t%i\n",$1,$2,$3+1)}' | bedtools sort -g chrom_hg19.sizes | bedtools slop -i stdin -g chrom_hg19.sizes -l 0 -r 1 > hic/D0_TADs.bed
```

**13. Convert called TADs into a segmented genome to define Chrom3D beads**
```bash
bedtools complement -i hic/D0_TADs.bed -g chrom_hg19.sizes | cat - hic/D0_TADs.bed | bedtools sort -g chrom_hg19.sizes | awk '$1!="chrY"' | awk '$1!="chrM"' > D0_beads.bed
```

**14. Map intra-chromosomal interactions from Hi-C to the beads defined in the previous step and aggregate the contacts between these beads**
```bash
cat hic/bedpe/intra/chr* | awk '{printf("%s\t%s\t%s\n",$1,$2,$2+1)}' | bedtools intersect -wao -a stdin -b D0_beads.bed | cut -f 4,5,6 > left.tmp
cat hic/bedpe/intra/chr* | awk '{printf("%s\t%s\t%s\t%s\n",$4,$5,$5+1,$7)}' | bedtools intersect -wao -a stdin -b D0_beads.bed | awk '{printf("%s\t%s\t%s\t%s\n",$5,$6,$7,$4)}' > right.tmp

paste left.tmp right.tmp | awk '{a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6] += $7} END{for (i in a) print i"\t"a[i]}' |  awk '$2!=$5' | sort -k 2n,2n > D0_bead_interactions.intra.bedpe
rm left.tmp right.tmp
```

**15. Remove interactions between beads overlapping centromeres**
```bash
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" | gunzip -c | grep acen | bedtools pairtobed -a D0_bead_interactions.intra.bedpe -b stdin -type neither > D0_bead_interactions.intra.nocen.bedpe
```

**16. Identifying statistically significant inter-bead interactions within chromosomes, using the Non-central Hypergeometric distribution (NCHG)**
```bash
processing_scripts/NCHG_hic/NCHG -m 50000 -p D0_bead_interactions.intra.nocen.bedpe > D0_bead_interactions.intra.nocen.NCHG.out

# Correcting for multiple-testing using FDR:
python processing_scripts/NCHG_fdr_oddratio_calc.py D0_bead_interactions.intra.nocen.NCHG.out fdr_bh 2 0.01 > D0_bead_interactions.intra.nocen.NCHG.sig
```
**17. Identifying statistically significant inter-bead interactions between chromosomes**
```bash
# Create a BEDPE file containing interchromosomal interactions where 'blacklisted' regions are removed:
cat hic/bedpe/inter/chr* | bedtools pairtobed -type neither -a stdin -b processing_scripts/unmappable_blacklist.bed | python processing_scripts/cap_chr_end.py chrom_hg19.sizes > D0_bead_interactions.inter.noblist.bedpe

# Run NCHG, like above, but on the interchromosomal interactions: 
processing_scripts/NCHG_hic/NCHG -i -p D0_bead_interactions.inter.noblist.bedpe > D0_bead_interactions.inter.noblist.NCHG.out

python processing_scripts/NCHG_fdr_oddratio_calc.py D0_bead_interactions.inter.noblist.NCHG.out fdr_bh 2 0.01 > D0_bead_interactions.inter.noblist.NCHG.sig

```
**18. Mapping all (intra- and interchromosomal interactions) to the TADs/beads
```bash
awk '{printf("%s\t%s\t%s\n",$1,($2+$3)/2,1+($2+$3)/2)}' D0_bead_interactions.inter.noblist.NCHG.sig | bedtools intersect -wao -a stdin -b D0_beads.bed | cut -f 4,5,6 > left.tmp
awk '{printf("%s\t%s\t%s\t%s\n",$4,($5+$6)/2,1+($5+$6)/2,$7)}' D0_bead_interactions.inter.noblist.NCHG.sig | bedtools intersect -wao -a stdin -b D0_beads.bed | awk '{printf("%s\t%s\t%s\t%s\n",$5,$6,$7,$4)}' > right.tmp

paste left.tmp right.tmp | sort -u -k 2n,2n > D0_bead_interactions.inter.noblist.NCHG.tadwise.sig
rm left.tmp right.tmp

cat D0_bead_interactions.intra.nocen.NCHG.sig D0_bead_interactions.inter.noblist.NCHG.tadwise.sig | awk '$2!=-1 && $5!=-1' > D0_bead_interactions.all.sig
```

**19. Generate the Chrom3D input file in GTrack format, specifying the 3D model setup**
```bash
python processing_scripts/makeGtrack.py D0_bead_interactions.all.sig D0_beads.bed > D0_bead_interactions.gtrack
```

**20. Add LAD information to the GTrack file***
```bash
# Create header
echo -e "##gtrack version: 1.0\n##track type: linked segments\n###seqid\tstart\tend\tid\tradius\tperiphery\tedges" > D0_bead_interactions.lads.gtrack

# Add LAD info:
bedtools intersect -c -a D0_bead_interactions.gtrack -b lad/GSE109924_lad_D0-rep1.bed | awk '{if($7>=1) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t1\t" $6; else  print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t.\t" $6}' >> D0_bead_interactions.lads.gtrack
```

**21. Make the GTrack file define a diploid genome structure (optional)**
```bash
python processing_scripts/make_diploid_gtrack.py D0_bead_interactions.lads.gtrack > D0_bead_interactions.lads.diploid.gtrack
```

**22. Run Chrom3D based on the GTrack file (takes up to 20 hrs)**
```bash
Chrom3D -l 10000 -y 0.15 -r 5.0 -n 3000000 D0_bead_interactions.lads.diploid.gtrack > model_full.cmm 2> model_full.err 
```

**23. Visualizing `model_full.cmm` in ChimeraX**
- If you are running this step-by-step guide on a server, download the `model.cmm` file to your local computer
- The resulting `model.cmm` (and `model_redlad.cmm` from step 19) can be opened in ChimeraX and displays of these turned on and off in the bottom right "Models" panel. To generate tomographic views of models, the command "clip" can be used in the "Command:" field in the bottom panel of ChimeraX. Background color and other graphical adjustments can be performed by clicking the "Graphics" button in the top panel.
- In ChimeraX, the command "shape sphere center 0,0,0 radius 5.0 color #ffc9b5 slab 0.5" can be used in the "Command" field in the bottom panel to display a nucleus structure on top of the model view. To change opacity of the nucleus model, click the colored square called "sphere" in the "Models" panel in the bottom right of the view, and select e.g. 30%. 
- Again, "clip" can be used to clip this to generate tomographic views. The model can also be tilted to allow a better perception of depth in the structures. Figure 4 shows some of the resulting illustrations that can be generated using ChimeraX.


**24. Coloring beads defined by LADs using red color**
```bash
unzip -j -d processing_scripts/ v.1.2.zip preprocess_scripts-v.1.2/color_beads.py

awk '$6==1' D0_bead_interactions.lads.gtrack | cut -f 4 > lads.ids
python3 processing_scripts/color_beads.py model_full.cmm lads.ids 255,0,0 OVERRIDE > model_full_redlad.cmm
```
- Do the steps from step 23., but open and visualize the `model_full_redlad.cmm` file instead
