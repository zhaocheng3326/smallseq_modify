#!/usr/bin/perl -w 
use strict;

### the stdin is 6 of total region, 6 of select region(intersectbed -a -b -wa -wb ), the argv[1] is the number of windows
my $w=$ARGV[0];
while (<STDIN>) {
    chomp;
    my @a=split "\t",$_;
    my $bg_begin;
    my $bg_end;
    if ($a[5] eq "+") {
        $bg_begin=$a[1];
        $bg_end=$a[2];
    }else{
        $bg_begin=$a[2];
        $bg_end=$a[1];
    }
    my ($b0,$b1);
    $b0=1+int(0.5 + abs( ($a[7]-$bg_begin)/($bg_begin-$bg_end)*$w));
    $b1=1+int(0.5 + abs( ($a[8]-$bg_begin)/($bg_begin-$bg_end)*$w));
    my ($w1,$w2);
    if ($b0 <=$b1) {
        $w1=$b0;
        $w2=$b1;
        if ($w2 >=$w) {$w2=$w};
    }else{
        $w1=$b1;
        $w2=$b0;
        if ($w2 >=$w) {$w2=$w};
    }
    my $key=$a[3]."|".$a[9];
    for (my $n=$w1;$n<=$w2;$n++) {
        print "$key\twin_$n\n";
    }
}

