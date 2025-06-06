#!/usr/bin/perl
#Script to create normalized BigWig tracks based on the spike-in scale factors

#Input file with all the scale factors per sample
$SFfile="ScaleFactor.txt";

$bin=$ARGV[0]; #first command-line argument when running the script
$p=10;         #number of threads to use in bamCompare

unless($bin > 0){
	die ("NO bin size!\n");
	}

$sm=$bin*3;   #defines the smoothing window size in BamCompare

#Creating output directories if they don't already exist
$odir="BigWig_SpikeIn/BigWig_bin${bin}";
system("if [ ! -d  BigWig_SpikeIn ];then mkdir BigWig_SpikeIn; fi");
system("if [ ! -d  $odir ];then mkdir $odir; fi");
system("if [ ! -d  BigWig_SpikeIn/Log ];then mkdir BigWig_SpikeIn/Log; fi");

open(IN, "$SFfile");                    #open the file
while(<IN>){
	chomp;
	my @data=split(/\t/, $_);           #loops through each line, removes the newline, and splits by tabs
	my @name=split(/\_/, $data[0]);
	my $out_hap1="$name[0]_$name[1]_$name[2]_SpikeIn_pInp_bin${bin}_sm${sm}_hap1.bw";                                       #extracts names for the outputs from the filenames of the BAM files (Hap1)
	my $out_hap2="$name[0]_$name[1]_$name[2]_SpikeIn_pInp_bin${bin}_sm${sm}_hap2.bw";                                       #extracts names for the outputs from the filenames of the BAM files (Hap2)
	my $sf="$data[4]:$data[5]";                                                                                             #combines the scale factor of IP and Input for bamCompare
	my $opt="--smoothLength $sm --ignoreDuplicates --skipNAs --operation ratio";                                            #parameters for bamCompare
	my $log_hap1="BigWig_SpikeIn/Log/$name[0]_$name[1]_$name[2]_Log_Bam2bw_pInp_Bin${bin}_hap1.txt";                        #extracts names for the log files
	my $log_hap2="BigWig_SpikeIn/Log/$name[0]_$name[1]_$name[2]_Log_Bam2bw_pInp_Bin${bin}_hap2.txt";                        #extracts names for the log files

	print "calculating Compare for Hap1 $data[0]\/$data[1]...\n";

	system("echo $data[0] $data[1] scaleFactor $sf");
	system("bamCompare -b1 $data[0] -b2 $data[1] -bs $bin -p $p --scaleFactors $sf $opt -o $odir/$out_hap1 2>$log_hap1");   #running bamCompare between IP bam and Input bam for Hap1

	print "calculating Compare for Hap2 $data[2]\/$data[3]...\n";

        system("echo $data[2] $data[3] scaleFactor $sf");
        system("bamCompare -b1 $data[2] -b2 $data[3] -bs $bin -p $p --scaleFactors $sf $opt -o $odir/$out_hap2 2>$log_hap2");   #running bamCompare between IP bam and Input bam for Hap1
}

close(IN);