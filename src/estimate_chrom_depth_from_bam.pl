#!/usr/bin/env perl

use warnings;
use strict;

# debug with stacktrace:
use Carp;
$SIG{__DIE__} = \&Carp::confess;


######################################################
#
# simple library functions:
#
my $LOGFH =\*STDERR;
my $OUTFH =\*STDOUT;

sub errorX($) {
    confess "ERROR: " . $_[0] . "\n";
}

sub logX($) {
    print $LOGFH "INFO: " . $_[0] . "\n";
}

sub checkFile($;$) {
    my ($file,$label) = @_;
    return if(-f $file);
    errorX("Can't find" . (defined($label) ? " $label" : "") . " file: '$file'");
}


#####################################

# optionally point to a samtools binary:
#
my $samtoolsBin=$ENV{'SAMTOOLS'};

if(! defined($samtoolsBin)) {
   $samtoolsBin = "samtools"; 
}

if(scalar(@ARGV)!=1) {
  print STDERR <<ENDE;

Create a simple estimate of mean depth for each chromosome given a WGS BAM file as input.

Usage: $0 bamfile > depth.txt

The file created by this script is designed to be used as input to the 'gatk_to_gvcf' tool.
The script requires samtools to either be in your path, or defined in the environment variable
'SAMTOOLS'. Note this method will not work on BAM files from exome or other targeted
sequencing data.

ENDE
  exit 2;
}

my $bamfile=$ARGV[0];

checkFile($bamfile,"bam");


my $cmd1="$samtoolsBin idxstats $bamfile";
open(my $FP1,"$cmd1 |");

my %chrom;
my @chroms;

while(<$FP1>) {
    my @X =split(/\t/);
    next if($X[0] eq "*");
    $chrom{$X[0]} = [ $X[1] , $X[2] ];
    push @chroms, $X[0];
}

close($FP1);


my $example_chrom = $chroms[0];
my $cmd2="$samtoolsBin view -F 4 -s 0.1 $bamfile $example_chrom";
open(my $FP2,"$cmd2 |");

my $length = 0;
my $count = 0;
while(<$FP2>) {
    my @F = split(/\t/,$_,7);
    next unless($F[5] =~ /[0-9]+M/);
    $length += $1 while($F[5] =~ /([0-9]+)M/g);
    $count++;
    last if($count >= 200000);
}
close($FP2);

if($count <= 100000) {
    die("Unexpected read length approximation results");
}


my $avg_length = ($length/$count);

for my $c (@chroms) {
    next if($chrom{$c}->[0] < $avg_length);
    my $depth = $chrom{$c}->[1]*$avg_length / $chrom{$c}->[0];
    printf "%s\t%.3f\t%s\t%.3f\n",$c,$depth,$count,$avg_length;
}


