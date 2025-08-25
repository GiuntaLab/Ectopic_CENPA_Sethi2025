#!/usr/bin/env Rscript
library(karyoploteR)
library(rtracklayer)

#set you working directory
#setwd("/data")
setwd("")

#creating the custom genome object for each haplotype
hap1_genome <- toGRanges(data.frame(chr= c("chr1_hap1", "chr2_hap1", "chr3_hap1", "chr4_hap1", "chr5_hap1", "chr6_hap1", "chr7_hap1", "chr8_hap1", "chr9_hap1", "chr10_hap1", "chr11_hap1", "chr12_hap1", "chr13_hap1", "chr14_hap1", "chr15_hap1", "chr16_hap1", "chr17_hap1", "chr18_hap1", "chr19_hap1", "chr20_hap1", "chr21_hap1", "chr22_hap1", "chrX_hap1"),
                                    start= c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                                    end= c(245130294, 242677997, 198760610, 192622708, 182259134, 172045469, 163647842, 145680301, 137521496, 134988382, 134265560, 133236650, 101016140, 97561019,
                                           95272790, 89874169, 82991830, 76653175, 60904190, 65674320, 38797553, 45250002, 227216509)))

hap2_genome <- toGRanges(data.frame(chr= c("chr1_hap2", "chr2_hap2", "chr3_hap2", "chr4_hap2", "chr5_hap2", "chr6_hap2", "chr7_hap2", "chr8_hap2", "chr9_hap2", "chr10_hap2", "chr11_hap2", "chr12_hap2", "chr13_hap2", "chr14_hap2", "chr15_hap2", "chr16_hap2", "chr17_hap2", "chr18_hap2", "chr19_hap2", "chr20_hap2", "chr21_hap2", "chr22_hap2", "chrX_hap2"),
                                    start= c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                                    end= c(246029683, 242154731, 201810706, 191569848, 181537949, 172484233,
                                           160431239, 145016226, 141000299, 134604704, 135257176, 133212390, 103308910, 95818454,
                                           90244609, 90313584, 82817096, 76658275, 61496696, 66327865, 39805135, 47559405, 154363813)))

#loading the post-processed bedgraph files for the genome-wide visualization
hap1_siEP <- import("RPE1_siEP_bin50_sm150_hap1_merged.bw.over1.bedGraph", format = "bedGraph")
hap2_siEP <- import("RPE1_siEP_bin50_sm150_hap2_merged.bw.over1.bedGraph", format = "bedGraph")
hap1_siNG <- import("RPE1_siNG_bin50_sm150_hap1_merged.bw.over1.bedGraph", format = "bedGraph")
hap2_siNG <- import("RPE1_siNG_bin50_sm150_hap2_merged.bw.over1.bedGraph", format = "bedGraph")

#Fig3a
tiff("hap1_all.tiff", width = 8.27, height = 11.69, units = "in", res = 1200)
pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight <- 0
pp$data2max <- 0
pp$data1max <- 50
pp$data1height <- 200
pp$bottommargin <- 1200
pp$topmargin <- 1200
pp$data1inmargin <- 0
pp$leftmargin <- 0.1
pp$rightmargin <- 0.1
pp$data2outmargin <- 0
pp$data1outmargin <- 0

kp <- plotKaryotype(plot.type=1, genome = hap1_genome, plot.params=pp, cex=0.6, chromosomes="all", cytobands = NULL)
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 2, tick.col="black", cex=0.5, minor.tick.dist = 5000000, minor.tick.col = "black", minor.tick.len = 1, units = "Mb")

kp <- kpPlotDensity(kp, data.panel=1,r0=5, r1=35, data = hap1_siEP, window.size = 3000, col="red2", border = "red2")
kpAxis(kp, data=hap1_siEP, data.panel=1, tick.len = 1e6, r0=5, r1=35, ymax=max(kp$latest.plot$computed.values$max.density), cex=0.4, side = 1)

kp <- kpPlotDensity(kp, data.panel=1,r0=5, r1=35, data = hap1_siNG, window.size = 3000, col="blue2", border = "blue2")
kpAxis(kp, data=hap1_siNG, data.panel=1, tick.len = 1e6, r0=5, r1=35, ymax=max(kp$latest.plot$computed.values$max.density), cex=0.4, side = 1)

dev.off()

#Fig3b
tiff("hap2_all.tiff", width = 8.27, height = 11.69, units = "in", res = 1200)
pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight <- 0
pp$data2max <- 0
pp$data1max <- 50
pp$data1height <- 200
pp$bottommargin <- 1200
pp$topmargin <- 1200
pp$data1inmargin <- 0
pp$leftmargin <- 0.1
pp$rightmargin <- 0.1
pp$data2outmargin <- 0
pp$data1outmargin <- 0

kp <- plotKaryotype(plot.type=1, genome = hap2_genome, plot.params=pp, cex=0.6, chromosomes="all", cytobands = NULL)
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 2, tick.col="black", cex=0.5, minor.tick.dist = 5000000, minor.tick.col = "black", minor.tick.len = 1, units = "Mb")

kp <- kpPlotDensity(kp, data.panel=1,r0=5, r1=35, data = hap2_siEP, window.size = 3000, col="red2", border = "red2")
kpAxis(kp, data=hap2_siEP, data.panel=1, tick.len = 1e6, r0=5, r1=35, ymax=max(kp$latest.plot$computed.values$max.density), cex=0.4, side = 1)

kp <- kpPlotDensity(kp, data.panel=1,r0=5, r1=35, data = hap2_siNG, window.size = 3000, col="blue2", border = "blue2")
kpAxis(kp, data=hap2_siNG, data.panel=1, tick.len = 1e6, r0=5, r1=35, ymax=max(kp$latest.plot$computed.values$max.density), cex=0.4, side = 1)

dev.off()

#loading the post-processed bedgraph files for the chr3 and FHIT visualization
rep1_hap1_siEP <- import("RPE-gCA1_siEP_CA_SpikeIn_pInp_bin50_sm150_hap1.bw.bedgraph.over1.bedGraph", format = "bedGraph")
rep1_hap1_siNG <- import("RPE-gCA1_siNG_CA_SpikeIn_pInp_bin50_sm150_hap1.bw.bedgraph.over1.bedGraph", format = "bedGraph")
rep2_hap1_siEP <- import("RPE-gCA2_siEP_CA_SpikeIn_pInp_bin50_sm150_hap1.bw.bedgraph.over1.bedGraph", format = "bedGraph")
rep2_hap1_siNG <- import("RPE-gCA2_siNG_CA_SpikeIn_pInp_bin50_sm150_hap1.bw.bedgraph.over1.bedGraph", format = "bedGraph")
rep3_hap1_siEP <- import("RPE-gCA3_siEP_CA_SpikeIn_pInp_bin50_sm150_hap1.bw.bedgraph.over1.bedGraph", format = "bedGraph")
rep3_hap1_siNG <- import("RPE-gCA3_siNG_CA_SpikeIn_pInp_bin50_sm150_hap1.bw.bedgraph.over1.bedGraph", format = "bedGraph")
rep1_hap2_siEP <- import("RPE-gCA1_siEP_CA_SpikeIn_pInp_bin50_sm150_hap2.bw.bedgraph.over1.bedGraph", format = "bedGraph")
rep1_hap2_siNG <- import("RPE-gCA1_siNG_CA_SpikeIn_pInp_bin50_sm150_hap2.bw.bedgraph.over1.bedGraph", format = "bedGraph")
rep2_hap2_siEP <- import("RPE-gCA2_siEP_CA_SpikeIn_pInp_bin50_sm150_hap2.bw.bedgraph.over1.bedGraph", format = "bedGraph")
rep2_hap2_siNG <- import("RPE-gCA2_siNG_CA_SpikeIn_pInp_bin50_sm150_hap2.bw.bedgraph.over1.bedGraph", format = "bedGraph")
rep3_hap2_siEP <- import("RPE-gCA3_siEP_CA_SpikeIn_pInp_bin50_sm150_hap2.bw.bedgraph.over1.bedGraph", format = "bedGraph")
rep3_hap2_siNG <- import("RPE-gCA3_siNG_CA_SpikeIn_pInp_bin50_sm150_hap2.bw.bedgraph.over1.bedGraph", format = "bedGraph")

#FigS7a
tiff("hap1_chr3_all.tiff", width = 8.27, height = 11.69, units = "in", res = 1200)

pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight <- 20
pp$data2max <- 0
pp$data1max <- 250
pp$data1height <- 500
pp$bottommargin <- 1200
pp$topmargin <- 1200
pp$data1inmargin <- 0
pp$leftmargin <- 0.3
pp$rightmargin <- 0.3
pp$data2outmargin <- 0
pp$data1outmargin <- 0

kp <- plotKaryotype(plot.type=1, genome = hap1_genome, plot.params=pp, cex=0.6, chromosomes="chr3_hap1", cytobands = NULL)
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 2, tick.col="black", cex=0.5, minor.tick.dist = 5000000, minor.tick.col = "black", minor.tick.len = 1, units = "Mb")

kp <- kpPlotDensity(kp, data.panel=1,r0=10, r1=40, data = rep1_hap1_siNG, window.size = 3000, col="blue2", border = "blue2")
kpAxis(kp, data=rep1_hap1_siNG, data.panel=1, r0=10, r1=40, ymax=max(kp$latest.plot$computed.values$max.density), cex=0.4, side = 1)

kp <- kpPlotDensity(kp, data.panel=1,r0=50, r1=80, data = rep2_hap1_siNG, window.size = 3000, col="blue2", border = "blue2")
kpAxis(kp, data=rep2_hap1_siNG, data.panel=1, r0=50, r1=80, ymax=max(kp$latest.plot$computed.values$max.density), cex=0.4, side = 1)

kp <- kpPlotDensity(kp, data.panel=1,r0=90, r1=120, data = rep3_hap1_siNG, window.size = 3000, col="blue2", border = "blue2")
kpAxis(kp, data=rep3_hap1_siNG, data.panel=1, r0=90, r1=120, ymax=max(kp$latest.plot$computed.values$max.density), cex=0.4, side = 1)

#rep1_hap1_siEP_f <- rep1_hap1_siEP[seqnames(rep1_hap1_siEP) == "chr3_hap1" & start(rep1_hap1_siEP) >= 52060531 & end(rep1_hap1_siEP) <= 53133926]
kp <- kpPlotDensity(kp, data.panel=1,r0=130, r1=160, data = rep1_hap1_siEP, window.size = 3000, col="red2", border = "red2")
kpAxis(kp, data=rep1_hap1_siEP, data.panel=1, r0=130, r1=160, ymax=max(kp$latest.plot$computed.values$max.density), cex=0.4, side = 1)

#rep2_hap1_siEP_f <- rep2_hap1_siEP[seqnames(rep2_hap1_siEP) == "chr3_hap1" & start(rep2_hap1_siEP) >= 52060531 & end(rep2_hap1_siEP) <= 53133926]
kp <- kpPlotDensity(kp, data.panel=1,r0=170, r1=200, data = rep2_hap1_siEP, window.size = 3000, col="red2", border = "red2")
kpAxis(kp, data=rep2_hap1_siEP, data.panel=1, r0=170, r1=200, ymax=max(kp$latest.plot$computed.values$max.density), cex=0.4, side = 1)

#rep3_hap1_siEP_f <- rep3_hap1_siEP[seqnames(rep3_hap1_siEP) == "chr3_hap1" & start(rep3_hap1_siEP) >= 52060531 & end(rep3_hap1_siEP) <= 53133926]
kp <- kpPlotDensity(kp, data.panel=1,r0=210, r1=240, data = rep3_hap1_siEP, window.size = 3000, col="red2", border = "red2")
kpAxis(kp, data=rep3_hap1_siEP, data.panel=1, r0=210, r1=240, ymax=max(kp$latest.plot$computed.values$max.density), cex=0.4, side = 1)
#kp <- kpPlotRegions(kp, data.panel = 1,r0=242, r1=250, data = gene_hap1_fhit, avoid.overlapping = F)

kpRect(kp, chr= "chr3_hap1", x0 = 52060531, x1 = 53133926,  y0=0, y1=245,  r0=5, r1=245, border="black", col=NA, data.panel=1, lty=2)

dev.off()

tiff("hap2_chr3_all.tiff", width = 8.27, height = 11.69, units = "in", res = 1200)

pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight <- 20
pp$data2max <- 0
pp$data1max <- 250
pp$data1height <- 500
pp$bottommargin <- 1200
pp$topmargin <- 1200
pp$data1inmargin <- 0
pp$leftmargin <- 0.3
pp$rightmargin <- 0.3
pp$data2outmargin <- 0
pp$data1outmargin <- 0

kp <- plotKaryotype(plot.type=1, genome = hap2_genome, plot.params=pp, cex=0.6, chromosomes="chr3_hap2", cytobands = NULL)
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 2, tick.col="black", cex=0.5, minor.tick.dist = 5000000, minor.tick.col = "black", minor.tick.len = 1, units = "Mb")

#rep1_hap2_siNG_f <- rep1_hap2_siNG[seqnames(rep1_hap2_siNG) == "chr3_hap2" & start(rep1_hap2_siNG) >= 59684759 & end(rep1_hap2_siNG) <= 61192174]
kp <- kpPlotDensity(kp, data.panel=1,r0=10, r1=40, data = rep1_hap2_siNG, window.size = 3000, col="blue2", border = "blue2")
kpAxis(kp, data=rep1_hap2_siNG, data.panel=1, r0=10, r1=40, ymax=max(kp$latest.plot$computed.values$max.density), cex=0.4, side = 1)

#rep2_hap2_siNG_f <- rep2_hap2_siNG[seqnames(rep2_hap2_siNG) == "chr3_hap2" & start(rep2_hap2_siNG) >= 59684759 & end(rep2_hap2_siNG) <= 61192174]
kp <- kpPlotDensity(kp, data.panel=1,r0=50, r1=80, data = rep2_hap2_siNG, window.size = 3000, col="blue2", border = "blue2")
kpAxis(kp, data=rep2_hap2_siNG, data.panel=1, r0=50, r1=80, ymax=max(kp$latest.plot$computed.values$max.density), cex=0.4, side = 1)

#rep3_hap2_siNG_f <- rep3_hap2_siNG[seqnames(rep3_hap2_siNG) == "chr3_hap2" & start(rep3_hap2_siNG) >= 59684759 & end(rep3_hap2_siNG) <= 61192174]
kp <- kpPlotDensity(kp, data.panel=1,r0=90, r1=120, data = rep3_hap2_siNG, window.size = 3000, col="blue2", border = "blue2")
kpAxis(kp, data=rep3_hap2_siNG, data.panel=1, r0=90, r1=120, ymax=max(kp$latest.plot$computed.values$max.density), cex=0.4, side = 1)

#rep1_hap2_siEP_f <- rep1_hap2_siEP[seqnames(rep1_hap2_siEP) == "chr3_hap2" & start(rep1_hap2_siEP) >= 59684759 & end(rep1_hap2_siEP) <= 61192174]
kp <- kpPlotDensity(kp, data.panel=1,r0=130, r1=160, data = rep1_hap2_siEP, window.size = 3000, col="red2", border = "red2")
kpAxis(kp, data=rep1_hap2_siEP, data.panel=1, r0=130, r1=160, ymax=max(kp$latest.plot$computed.values$max.density), cex=0.4, side = 1)

#rep2_hap2_siEP_f <- rep2_hap2_siEP[seqnames(rep2_hap2_siEP) == "chr3_hap2" & start(rep2_hap2_siEP) >= 59684759 & end(rep2_hap2_siEP) <= 61192174]
kp <- kpPlotDensity(kp, data.panel=1,r0=170, r1=200, data = rep2_hap2_siEP, window.size = 3000, col="red2", border = "red2")
kpAxis(kp, data=rep2_hap2_siEP, data.panel=1, r0=170, r1=200, ymax=max(kp$latest.plot$computed.values$max.density), cex=0.4, side = 1)

#rep3_hap2_siEP_f <- rep3_hap2_siEP[seqnames(rep3_hap2_siEP) == "chr3_hap2" & start(rep3_hap2_siEP) >= 59684759 & end(rep3_hap2_siEP) <= 61192174]
kp <- kpPlotDensity(kp, data.panel=1,r0=210, r1=240, data = rep3_hap2_siEP, window.size = 3000, col="red2", border = "red2")
kpAxis(kp, data=rep3_hap2_siEP, data.panel=1, r0=210, r1=240, ymax=max(kp$latest.plot$computed.values$max.density), cex=0.4, side = 1)

#kp <- kpPlotRegions(kp, data.panel = 1,r0=242, r1=250, data = gene_hap2_fhit, avoid.overlapping = F)
kpRect(kp, chr= "chr3_hap2", x0 = 59684759, x1 = 61192174,  y0=0, y1=245,  r0=5, r1=245, border="black", col=NA, data.panel=1, lty=2)

dev.off()

tiff("hap1_chr3_fhit.tiff", width = 8.27, height = 11.69, units = "in", res = 1200)

pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight <- 20
pp$data2max <- 0
pp$data1max <- 250
pp$data1height <- 500
pp$bottommargin <- 1200
pp$topmargin <- 1200
pp$data1inmargin <- 0
pp$leftmargin <- 0.3
pp$rightmargin <- 0.3
pp$data2outmargin <- 0
pp$data1outmargin <- 0

kp <- plotKaryotype(plot.type=1, genome = hap1_genome, plot.params=pp, cex=0.6, chromosomes="chr3_hap1", cytobands = NULL, zoom = "chr3_hap1:52060531-53133926")
kpAddBaseNumbers(kp, tick.dist = 100000, tick.len = 2, tick.col="black", cex=0.5, minor.tick.dist = 50000, minor.tick.col = "black", minor.tick.len = 1, units = "Mb")

rep1_hap1_siNG_f <- rep1_hap1_siNG[seqnames(rep1_hap1_siNG) == "chr3_hap1" & start(rep1_hap1_siNG) >= 52060531 & end(rep1_hap1_siNG) <= 53133926]
kp <- kpPlotDensity(kp, data.panel=1,r0=10, r1=40, data = rep1_hap1_siNG_f, window.size = 3000, col="white", border = "blue2", ymax = 15)
kpAxis(kp, data=rep1_hap1_siNG_f, data.panel=1, r0=10, r1=40, ymax=15, cex=0.4, side = 1)

rep2_hap1_siNG_f <- rep2_hap1_siNG[seqnames(rep2_hap1_siNG) == "chr3_hap1" & start(rep2_hap1_siNG) >= 52060531 & end(rep2_hap1_siNG) <= 53133926]
kp <- kpPlotDensity(kp, data.panel=1,r0=50, r1=80, data = rep2_hap1_siNG_f, window.size = 3000, col="white", border = "blue2", ymax = 15)
kpAxis(kp, data=rep2_hap1_siNG_f, data.panel=1, r0=50, r1=80, ymax=15, cex=0.4, side = 1)

#rep3_hap1_siNG_f <- rep3_hap1_siNG[seqnames(rep3_hap1_siNG) == "chr3_hap1" & start(rep3_hap1_siNG) >= 52060531 & end(rep3_hap1_siNG) <= 53133926]
kp <- kpPlotDensity(kp, data.panel=1,r0=90, r1=120, data = rep3_hap1_siNG, window.size = 3000, col="white", border = "blue2", ymax = 15)
#kpAxis(kp, data=rep3_hap1_siNG_f, data.panel=1, r0=90, r1=120, ymax=max(kp$latest.plot$computed.values$max.density), cex=0.4, side = 1)

rep1_hap1_siEP_f <- rep1_hap1_siEP[seqnames(rep1_hap1_siEP) == "chr3_hap1" & start(rep1_hap1_siEP) >= 52060531 & end(rep1_hap1_siEP) <= 53133926]
kp <- kpPlotDensity(kp, data.panel=1,r0=130, r1=160, data = rep1_hap1_siEP_f, window.size = 3000, col="white", border = "red2", ymax = 15)
kpAxis(kp, data=rep1_hap1_siEP_f, data.panel=1, r0=130, r1=160, ymax=15, cex=0.4, side = 1)

rep2_hap1_siEP_f <- rep2_hap1_siEP[seqnames(rep2_hap1_siEP) == "chr3_hap1" & start(rep2_hap1_siEP) >= 52060531 & end(rep2_hap1_siEP) <= 53133926]
kp <- kpPlotDensity(kp, data.panel=1,r0=170, r1=200, data = rep2_hap1_siEP_f, window.size = 3000, col="white", border = "red2", ymax = 15)
kpAxis(kp, data=rep2_hap1_siEP_f, data.panel=1, r0=170, r1=200, ymax=15, cex=0.4, side = 1)

rep3_hap1_siEP_f <- rep3_hap1_siEP[seqnames(rep3_hap1_siEP) == "chr3_hap1" & start(rep3_hap1_siEP) >= 52060531 & end(rep3_hap1_siEP) <= 53133926]
kp <- kpPlotDensity(kp, data.panel=1,r0=210, r1=240, data = rep3_hap1_siEP_f, window.size = 3000, col="white", border = "red2", ymax = 15)
kpAxis(kp, data=rep3_hap1_siEP_f, data.panel=1, r0=210, r1=240, ymax=15, cex=0.4, side = 1)

dev.off()

tiff("hap2_chr3_fhit.tiff", width = 8.27, height = 11.69, units = "in", res = 1200)

pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight <- 20
pp$data2max <- 0
pp$data1max <- 250
pp$data1height <- 500
pp$bottommargin <- 1200
pp$topmargin <- 1200
pp$data1inmargin <- 0
pp$leftmargin <- 0.3
pp$rightmargin <- 0.3
pp$data2outmargin <- 0
pp$data1outmargin <- 0

kp <- plotKaryotype(plot.type=1, genome = hap2_genome, plot.params=pp, cex=0.6, chromosomes="chr3_hap2", cytobands = NULL, zoom = "chr3_hap2:59684759-61192174")
kpAddBaseNumbers(kp, tick.dist = 100000, tick.len = 2, tick.col="black", cex=0.5, minor.tick.dist = 50000, minor.tick.col = "black", minor.tick.len = 1, units = "Mb")

rep1_hap2_siNG_f <- rep1_hap2_siNG[seqnames(rep1_hap2_siNG) == "chr3_hap2" & start(rep1_hap2_siNG) >= 59684759 & end(rep1_hap2_siNG) <= 61192174]
kp <- kpPlotDensity(kp, data.panel=1,r0=10, r1=40, data = rep1_hap2_siNG_f, window.size = 3000, col="white", border = "blue2", ymax = 15)
kpAxis(kp, data=rep1_hap2_siNG_f, data.panel=1, r0=10, r1=40, ymax=15, cex=0.4, side = 1)

rep2_hap2_siNG_f <- rep2_hap2_siNG[seqnames(rep2_hap2_siNG) == "chr3_hap2" & start(rep2_hap2_siNG) >= 59684759 & end(rep2_hap2_siNG) <= 61192174]
kp <- kpPlotDensity(kp, data.panel=1,r0=50, r1=80, data = rep2_hap2_siNG_f, window.size = 3000, col="white", border = "blue2", ymax = 15)
kpAxis(kp, data=rep2_hap2_siNG_f, data.panel=1, r0=50, r1=80, ymax=15, cex=0.4, side = 1)

#rep3_hap2_siNG_f <- rep3_hap2_siNG[seqnames(rep3_hap2_siNG) == "chr3_hap2" & start(rep3_hap2_siNG) >= 59684759 & end(rep3_hap2_siNG) <= 61192174]
kp <- kpPlotDensity(kp, data.panel=1,r0=90, r1=120, data = rep3_hap2_siNG, window.size = 3000, col="white", border = "blue2", ymax = 15)
#kpAxis(kp, data=rep3_hap2_siNG_f, data.panel=1, r0=90, r1=120, ymax=max(kp$latest.plot$computed.values$max.density), cex=0.4, side = 1)

rep1_hap2_siEP_f <- rep1_hap2_siEP[seqnames(rep1_hap2_siEP) == "chr3_hap2" & start(rep1_hap2_siEP) >= 59684759 & end(rep1_hap2_siEP) <= 61192174]
kp <- kpPlotDensity(kp, data.panel=1,r0=130, r1=160, data = rep1_hap2_siEP_f, window.size = 3000, col="white", border = "red2", ymax = 15)
kpAxis(kp, data=rep1_hap2_siEP_f, data.panel=1, r0=130, r1=160, ymax=15, cex=0.4, side = 1)

rep2_hap2_siEP_f <- rep2_hap2_siEP[seqnames(rep2_hap2_siEP) == "chr3_hap2" & start(rep2_hap2_siEP) >= 59684759 & end(rep2_hap2_siEP) <= 61192174]
kp <- kpPlotDensity(kp, data.panel=1,r0=170, r1=200, data = rep2_hap2_siEP_f, window.size = 3000, col="white", border = "red2", ymax = 15)
kpAxis(kp, data=rep2_hap2_siEP_f, data.panel=1, r0=170, r1=200, ymax=15, cex=0.4, side = 1)

rep3_hap2_siEP_f <- rep3_hap2_siEP[seqnames(rep3_hap2_siEP) == "chr3_hap2" & start(rep3_hap2_siEP) >= 59684759 & end(rep3_hap2_siEP) <= 61192174]
kp <- kpPlotDensity(kp, data.panel=1,r0=210, r1=240, data = rep3_hap2_siEP_f, window.size = 3000, col="white", border = "red2", ymax = 15)
kpAxis(kp, data=rep3_hap2_siEP_f, data.panel=1, r0=210, r1=240, ymax=15, cex=0.4, side = 1)

dev.off()


