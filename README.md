# Ectopic_CENPA_Sethi2025
This repository contains scripts for data processing and visualization related to the study **"Chromatin remodeling activity of EP400 safeguards chromosomal stability by preventing CENP-A mislocalization"** by Sethi et al., 2025 (accepted in principle).\
The analysis was performed on ChIP-seq data with spike-in normalization on RPE-1 cells to assess the localization of CENP-A upon depletion of EP400, which is a subunit of the NuA4 histone acetyltransferase complex. The analysis shows that EP400 depletion promotes CENP-A mislocalization on chromosome arms, confirming previous observation obtained using metaphase chromosome spreads and biochemical assays.

### Main process:
Before mapping, we prepared the concatenated reference genome, RPE1v1.1 (hap1 + hap2) + dm6. Then we ran the 1st script **Step1_bowtie2_RPE1v1.1_dm6.sh** for the alignment and computed the read count for the merged references, RPE1 hap1, RPE1 hap2 and dm6.
After putting the 2nd script **Step2_ReadCount_combine.pl** and the 3rd script **Step3_BamComp_SpikeIn.pl** in the same folder as the BAM files, we calculated the Spike-in scale factor (SF) and created a table named ScaleFactor.txt file by running the 2nd script. To normalize the mapped reads both by input and Spike-in genome, we ran the 3rd script with a proper bin size
(e.g. 50bp) and got BigWig files for each sample. BigWig files are now ready for visualization on IGV with Group Autoscale mode. To quantify signal enrichment, BigWig files were converted to BedGraph, enabling the extraction of the IP/Input ratio in the fourth column for downstream analysis.


### Tools and R packages:
+ FastQC (v0.12.1)
+ MultiQC (v1.24)
+ TrimGalore (v0.6.10)
+ Cutadapt (v2.8)
+ Bowtie2 (v2.4.4)
+ Samtools (v1.12)
+ deepTools (v3.5.1)
+ bedtools (v2.27.1)
+ bigWigToBedGraph (UCSC's utility)
+ bigWigMerge (UCSC's utility)
+ bigwigCompare (v3.5.5)
+ data.table (v1.17.4)
+ ggplot2 (v3.5.2)
+ karyoploteR (v1.28.0)
+ rtracklayer (v1.68.0)
+ dplyr (1.1.4)
+ tidyr (v1.3.1)
+ ggrepel (v0.9.6)
+ gridExtra (v2.3)
+ tidyverse (v2.0.0)

### References:
+ **Chromatin remodeling activity of EP400 safeguards chromosomal stability by preventing CENP-A mislocalization** by Sethi et al., 2025 (accepted in principle)
+ **The complete diploid reference genome of RPE-1 identifies human phased epigenetic landscapes** by Volpe et al., bioRxiv 2023 (currently under revision)
