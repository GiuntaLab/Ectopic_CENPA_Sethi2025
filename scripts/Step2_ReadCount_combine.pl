#!/usr/bin/perl
#Script to combine reads count from all the samples

#Getting the current system date and format and use it as a prefix in output filenames
my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
my $Date = sprintf("%04d%02d%02d", $year + 1900, $mon + 1, $mday);

#Defining the BAM filename suffixes for haplotypes
$Fix="RPE1v1.1_hap1_sort.bam";
$Fox="RPE1v1.1_hap2_sort.bam";

#Create the output file for combined read counts and scale factors
open(OUT, ">${Date}_ReadCount_Combine.txt");
print OUT "Sample\tRC_RPE1hap1\tRC_RPE1hap2\tRC_dm6\tSF_dm6\n";

#Read all read count files from the folder
my @f = glob "../ReadCount/*.txt";

#Loop through each read count file and computes the spike-in scale factor for each sample
foreach $f(@f){
	open(IN, "$f");
	while(<IN>){                           #iterates over each line in the file
		next if($_ =~ /Sample/);           #skip header
		chomp;                             #remove the newline character (\n) from the end of a string
		my @data=split(/\t/, $_);          #split line by tab
		my $sf = (1000000/$data[-1]);      #divide 1M by last column = dm6 reads
		print OUT "$_\t$sf\n";             #write original line + new SF_dm6 value

	}
	close(IN);
}

close(OUT);

#Parsed the new combined file to build a hash of scale factors
open(IN, "${Date}_ReadCount_Combine.txt");
while(<IN>){                               #iterates over each line in the file
	next if($_ =~ /Sample/);               #skip header
	chomp;                                 #remove the newline character (\n) from the end of a string
	my @data=split(/\t/, $_);              #split by tab
	my @name=split(/\_/, $data[0]);        #split sample name by underscores
	my $tag="$name[0]_$name[1]_$name[2]";  #rebuild tag
	$sf{$tag}=$data[-1];                   #assign scale factor (last column) to that tag

print  "$tag\t$sf{$tag}\n";
}
close(IN);

#Writing the ScaleFactor.txt file
open(OUT, "> ScaleFactor.txt");

foreach $key(sort keys(%sf)){    #loops through each key in the %sf hash
	my @tag2=split(/\_/, $key);  #splits tag by underscores
	my $tag_inp="$tag2[0]_$tag2[1]_input"; #reconstructs the corresponding input sample name

	#Writes to ScaleFactor.txt a line with BAM file patterns for hap1 and hap2 for both experiment and the input and their corresponding scale factors
	print OUT "${key}_*${Fix}\t${tag_inp}_*${Fix}\t${key}_*${Fox}\t${tag_inp}_*${Fox}\t$sf{$key}\t$sf{$tag_inp}\n";
}

close(OUT);