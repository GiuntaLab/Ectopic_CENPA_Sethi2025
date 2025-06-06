# Data processing
We analyzed spike-in ChIP-seq data for CENP-A protein from three biological replicates per condition in RPE-1 cells treated with siNEG or siEP400. Each replicate included a corresponding input background. The data are available in the NCBI repository (BioProject accession number: PRJNA1182946). Scripts for mapping and spike-in normalization were provided by Tetsuya Hori and adapted by Luca Corda. 

## Index:

```commandline
bowtie2-build --large-index -f RPE1v1.1_dm6.fasta RPE1v1.1_dm6.fasta
```
## Mapping:
Adapter trimming and quality control were performed using Trim Galore, with Cutadapt, and FastQC. MultiQC was used to aggregate quality control reports.  Alignments were performed using Bowtie2 with custom parameters specified within the script, sorting the alignment BAM files by coordinates using samtools.
Running this script will create BAM files and reads count files for each sample.
```commandline
#fastq files must be in the same folder
bash bowtie2_RPE1v1.1_dm6.sh 
```
## Calculating spike-in scale factor (SF)
Now you can navigate to the `Bam/` directory using `cd`. The reads count files for each sample are now available at `../ReadCount/*.txt` and will be used as input for the perl script.\
Running this script will create the scale factors for each sample and condition.
```commandline
#BAM files must be in the same folder and reads count in ../ReadCount/*.txt
perl Step1_ReadCount_combine.pl 
```
## Normalizing mapped reads both by Input and Spike-in genome
To normalize the mapped reads both by input and Spike-in genome, we ran the 3rd script with a proper bin size (e.g. 50bp) and got BigWig files for each sample.

```
#BAM files and ScaleFactor.txt must be in the same folder
perl Step3_BamComp_SpikeIn.pl 50
```
# Data plotting
After completing the data processing steps and generating spike-in-normalized BigWig files (using --operation ratio) representing IP/Input signals, the resulting files were ready for visualization on IGV with Group Autoscale mode. For custom plotting in R, the BigWig files were either merged or converted into a single BedGraph file, which is compatible with the karyoploteR package used to generate genome-wide chromosome plots.

```
#merge multiple bigwig of the same condition into a single bedgraph
bigWigMerge in1.bw in2.bw .. inN.bw out.bedGraph

#convert a single bigwig file to bedgraph format for custom visualization
bigWigToBedGraph in.bigWig out.bedGraph
```

[HumAS-HMMER_for_AnVIL](https://github.com/fedorrik/HumAS-HMMER_for_AnVIL) was used to annotate centromeres in the RPE1v1.1 genome, as described in [Volpe et al.](https://pubmed.ncbi.nlm.nih.gov/38168337/). The centromeric alpha-satellite annotation file was used to intersect the BedGraph signals and retrieve only centromeric or no-centromeric signals.  

#### Fig3a, b and FigS7a

#### Fig3c and FigS6
Genome-wide and chromosome-level plotting of CENP-A IP/Input signals was performed using the custom **Step5_Fig3c_FigS6.R** script in R. 
Input data for the R script were generated using `bedtools` intersect with the centromere annotation file of the RPE1v1.1 genome.

```
#parsing all bedgraph files with bedtools intersect
for i in *bedgraph; do bash bedgraph_intersect.sh $i; done 
```
The resulting files were processed with `awk` to add sample information:

```
#siEP rep1
awk -v OFS="\t" '{print $0, "REP1"}' RPE-gCA1_siEP_CA_SpikeIn_pInp_bin50_sm150_hap1.bw.bedgraph_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siEP"}' | awk -v OFS="\t" '{print $0, "in_centromere"}' | awk -v OFS="\t" '{print $0, "hap1"}' > REP1_siEP_in_AS_HOR_hap1.bed
awk -v OFS="\t" '{print $0, "REP1"}' RPE-gCA1_siEP_CA_SpikeIn_pInp_bin50_sm150_hap1.bw.bedgraph_out_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siEP"}' | awk -v OFS="\t" '{print $0, "out_centromere"}' | awk -v OFS="\t" '{print $0, "hap1"}' > REP1_siEP_out_AS_HOR_hap1.bed
awk -v OFS="\t" '{print $0, "REP1"}' RPE-gCA1_siEP_CA_SpikeIn_pInp_bin50_sm150_hap2.bw.bedgraph_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siEP"}' | awk -v OFS="\t" '{print $0, "in_centromere"}' | awk -v OFS="\t" '{print $0, "hap2"}' > REP1_siEP_in_AS_HOR_hap2.bed
awk -v OFS="\t" '{print $0, "REP1"}' RPE-gCA1_siEP_CA_SpikeIn_pInp_bin50_sm150_hap2.bw.bedgraph_out_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siEP"}' | awk -v OFS="\t" '{print $0, "out_centromere"}' | awk -v OFS="\t" '{print $0, "hap2"}' > REP1_siEP_out_AS_HOR_hap2.bed

#siNG rep1
awk -v OFS="\t" '{print $0, "REP1"}' RPE-gCA1_siNG_CA_SpikeIn_pInp_bin50_sm150_hap1.bw.bedgraph_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siNG"}' | awk -v OFS="\t" '{print $0, "in_centromere"}' | awk -v OFS="\t" '{print $0, "hap1"}' > REP1_siNG_in_AS_HOR_hap1.bed
awk -v OFS="\t" '{print $0, "REP1"}' RPE-gCA1_siNG_CA_SpikeIn_pInp_bin50_sm150_hap1.bw.bedgraph_out_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siNG"}' | awk -v OFS="\t" '{print $0, "out_centromere"}' | awk -v OFS="\t" '{print $0, "hap1"}' > REP1_siNG_out_AS_HOR_hap1.bed
awk -v OFS="\t" '{print $0, "REP1"}' RPE-gCA1_siNG_CA_SpikeIn_pInp_bin50_sm150_hap2.bw.bedgraph_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siNG"}' | awk -v OFS="\t" '{print $0, "in_centromere"}' | awk -v OFS="\t" '{print $0, "hap2"}' > REP1_siNG_in_AS_HOR_hap2.bed
awk -v OFS="\t" '{print $0, "REP1"}' RPE-gCA1_siNG_CA_SpikeIn_pInp_bin50_sm150_hap2.bw.bedgraph_out_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siNG"}' | awk -v OFS="\t" '{print $0, "out_centromere"}' | awk -v OFS="\t" '{print $0, "hap2"}' > REP1_siNG_out_AS_HOR_hap2.bed

#siEP rep2
awk -v OFS="\t" '{print $0, "REP2"}' RPE-gCA2_siEP_CA_SpikeIn_pInp_bin50_sm150_hap1.bw.bedgraph_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siEP"}' | awk -v OFS="\t" '{print $0, "in_centromere"}' | awk -v OFS="\t" '{print $0, "hap1"}' > REP2_siEP_in_AS_HOR_hap1.bed
awk -v OFS="\t" '{print $0, "REP2"}' RPE-gCA2_siEP_CA_SpikeIn_pInp_bin50_sm150_hap1.bw.bedgraph_out_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siEP"}' | awk -v OFS="\t" '{print $0, "out_centromere"}' | awk -v OFS="\t" '{print $0, "hap1"}' > REP2_siEP_out_AS_HOR_hap1.bed
awk -v OFS="\t" '{print $0, "REP2"}' RPE-gCA2_siEP_CA_SpikeIn_pInp_bin50_sm150_hap2.bw.bedgraph_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siEP"}' | awk -v OFS="\t" '{print $0, "in_centromere"}' | awk -v OFS="\t" '{print $0, "hap2"}' > REP2_siEP_in_AS_HOR_hap2.bed
awk -v OFS="\t" '{print $0, "REP2"}' RPE-gCA2_siEP_CA_SpikeIn_pInp_bin50_sm150_hap2.bw.bedgraph_out_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siEP"}' | awk -v OFS="\t" '{print $0, "out_centromere"}' | awk -v OFS="\t" '{print $0, "hap2"}' > REP2_siEP_out_AS_HOR_hap2.bed

#siNG rep2
awk -v OFS="\t" '{print $0, "REP2"}' RPE-gCA2_siNG_CA_SpikeIn_pInp_bin50_sm150_hap1.bw.bedgraph_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siNG"}' | awk -v OFS="\t" '{print $0, "in_centromere"}' | awk -v OFS="\t" '{print $0, "hap1"}' > REP2_siNG_in_AS_HOR_hap1.bed
awk -v OFS="\t" '{print $0, "REP2"}' RPE-gCA2_siNG_CA_SpikeIn_pInp_bin50_sm150_hap1.bw.bedgraph_out_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siNG"}' | awk -v OFS="\t" '{print $0, "out_centromere"}' | awk -v OFS="\t" '{print $0, "hap1"}' > REP2_siNG_out_AS_HOR_hap1.bed
awk -v OFS="\t" '{print $0, "REP2"}' RPE-gCA2_siNG_CA_SpikeIn_pInp_bin50_sm150_hap2.bw.bedgraph_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siNG"}' | awk -v OFS="\t" '{print $0, "in_centromere"}' | awk -v OFS="\t" '{print $0, "hap2"}' > REP2_siNG_in_AS_HOR_hap2.bed
awk -v OFS="\t" '{print $0, "REP2"}' RPE-gCA2_siNG_CA_SpikeIn_pInp_bin50_sm150_hap2.bw.bedgraph_out_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siNG"}' | awk -v OFS="\t" '{print $0, "out_centromere"}' | awk -v OFS="\t" '{print $0, "hap2"}' > REP2_siNG_out_AS_HOR_hap2.bed

#siEP rep3
awk -v OFS="\t" '{print $0, "REP3"}' RPE-gCA3_siEP_CA_SpikeIn_pInp_bin50_sm150_hap1.bw.bedgraph_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siEP"}' | awk -v OFS="\t" '{print $0, "in_centromere"}' | awk -v OFS="\t" '{print $0, "hap1"}' > REP3_siEP_in_AS_HOR_hap1.bed
awk -v OFS="\t" '{print $0, "REP3"}' RPE-gCA3_siEP_CA_SpikeIn_pInp_bin50_sm150_hap1.bw.bedgraph_out_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siEP"}' | awk -v OFS="\t" '{print $0, "out_centromere"}' | awk -v OFS="\t" '{print $0, "hap1"}' > REP3_siEP_out_AS_HOR_hap1.bed
awk -v OFS="\t" '{print $0, "REP3"}' RPE-gCA3_siEP_CA_SpikeIn_pInp_bin50_sm150_hap2.bw.bedgraph_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siEP"}' | awk -v OFS="\t" '{print $0, "in_centromere"}' | awk -v OFS="\t" '{print $0, "hap2"}' > REP3_siEP_in_AS_HOR_hap2.bed
awk -v OFS="\t" '{print $0, "REP3"}' RPE-gCA3_siEP_CA_SpikeIn_pInp_bin50_sm150_hap2.bw.bedgraph_out_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siEP"}' | awk -v OFS="\t" '{print $0, "out_centromere"}' | awk -v OFS="\t" '{print $0, "hap2"}' > REP3_siEP_out_AS_HOR_hap2.bed

#siNG rep3
awk -v OFS="\t" '{print $0, "REP3"}' RPE-gCA3_siNG_CA_SpikeIn_pInp_bin50_sm150_hap1.bw.bedgraph_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siNG"}' | awk -v OFS="\t" '{print $0, "in_centromere"}' | awk -v OFS="\t" '{print $0, "hap1"}' > REP3_siNG_in_AS_HOR_hap1.bed
awk -v OFS="\t" '{print $0, "REP3"}' RPE-gCA3_siNG_CA_SpikeIn_pInp_bin50_sm150_hap1.bw.bedgraph_out_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siNG"}' | awk -v OFS="\t" '{print $0, "out_centromere"}' | awk -v OFS="\t" '{print $0, "hap1"}' > REP3_siNG_out_AS_HOR_hap1.bed
awk -v OFS="\t" '{print $0, "REP3"}' RPE-gCA3_siNG_CA_SpikeIn_pInp_bin50_sm150_hap2.bw.bedgraph_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siNG"}' | awk -v OFS="\t" '{print $0, "in_centromere"}' | awk -v OFS="\t" '{print $0, "hap2"}' > REP3_siNG_in_AS_HOR_hap2.bed
awk -v OFS="\t" '{print $0, "REP3"}' RPE-gCA3_siNG_CA_SpikeIn_pInp_bin50_sm150_hap2.bw.bedgraph_out_AS_HOR.bed | awk -v OFS="\t" '{print $0, "siNG"}' | awk -v OFS="\t" '{print $0, "out_centromere"}' | awk -v OFS="\t" '{print $0, "hap2"}' > REP3_siNG_out_AS_HOR_hap2.bed
```

Before running the R script, make sure to place it in the same directory as the *_AS_HOR_hap1.bed* and *_AS_HOR_hap2.bed* files produced from the previous codes, and manually set the working directory at the beginning of the script using `setwd()`.

```
Rscript Step5_Fig3c_FigS6.R
```

#### Fig3d and FigS7b