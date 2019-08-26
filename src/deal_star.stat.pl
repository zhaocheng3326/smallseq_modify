#!/usr/bin/perl -w
use strict;

my ($total, $umc,$ump,$mmp,$mmc);
while (<>) {
    chomp;
    if ($_=~/Number of input reads\s*\|\s*(.*?)$/) {
        $total=$1;
    }
    if ($_=~/Uniquely mapped reads number\s*\|\s*(.*?)$/) {
        $umc=$1;
    }
    if ($_=~/Uniquely mapped reads %\s*\|\s*(.*?)%$/) {
        $ump=$1;
    }
    if ($_=~/Number of reads mapped to multiple loci\s*\|\s*(.*?)$/) {
        $mmc=$1;
    }
    if ($_=~/% of reads mapped to multiple loci\s*\|\s*(.*?)%$/) {
        $mmp=$1;
    }
}
print "beforemapping\t$total\nUMP_reads\t$umc\nUMP_ratio\t$ump\nMMP_reads\t$mmc\nMMP_ratio\t$mmp";
