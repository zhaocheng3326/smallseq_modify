#!/usr/bin/perl -w                                                                                                                                                                                
#die "please input STDIN and one file " if @ARGV !=2;
#
open A,"$ARGV[0]" or die "$ARGV[0] can't be opened";
#
my %stat;

while(<STDIN>) {
    chomp;
    my @array=split "\t",$_;
    $stat{$array[0]}=$array[1];
}

while (<A>) {
    chomp;
    if (not exists $stat{$_}) {
        $stat{$_}=0;
    }
    print "$_\t$stat{$_}\n";
}
