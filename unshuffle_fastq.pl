#!/usr/bin/perl
use Getopt::Long;

GetOptions("f=s"     => \$fastqfile,
           "n=s"     => \$name,
           "o=s"     => \$newPath
    );

unless ($fastqfile || $name || $newPath) {
    die "\nusage: unshuffle_fastq.pl -f FASTQFILE -o OUTPUTPATH -n NAME

";
}

open(FASTQ,$fastqfile) || die "cannot open fastq file: $!";
open(OUT1,">$newPath/$name.1") || die "cannot open output file: $!";
open(OUT2,">$newPath/$name.2") || die "cannot open output file: $!";
while($l1=<FASTQ>) {
#@D0ENAACXX111117:4:1101:10001:111966/1
#TCACCAGCAAACGGTCCGCCAGCAGGCGGGCGACGCCTAAATCATGGGTGACAATCACCACCGCGAGGTTCAGCTCCACCACCAGGCCGCGCAGCAGGTCG
#+
#CCCFFFFFHHHHHJIJJIJJJJIIJJIJIIDCCDDDDDDDDDDDDDDD?BDDDDDDDDDDBDDDDBDD5?@CDDDCDACDBDDDDDDDBBDDDDBDBD?C@
#@D0ENAACXX111117:4:1101:10001:111966/2
#TCGCCTGCTGCGTACCGAATGGGGCGTGGTGCATCAGCATCCACTCGACGGCCTGCGCCGCCAGGTGTCGGCAGGCGGCAATATCGGCGAGGGGGTGGGGG
#+
#C@CFFFFFHHGHHIIGHIGGIIGGGI@FH8CBFGHEEHHGGHBHHGFDFCDBD:BA?B7;@B7<B:@@@@DD>B@9>B#######################
    $l2=<FASTQ>;
    $l3=<FASTQ>;
    $l4=<FASTQ>;
    $l5=<FASTQ>;
    $l6=<FASTQ>;
    $l7=<FASTQ>;
    $l8=<FASTQ>;

    print OUT1 $l1;
    print OUT1 $l2;
    print OUT1 $l3;
    print OUT1 $l4;
    print OUT2 $l5;
    print OUT2 $l6;
    print OUT2 $l7;
    print OUT2 $l8;
}
close FASTQ;
close OUT1;
close OUT2;
