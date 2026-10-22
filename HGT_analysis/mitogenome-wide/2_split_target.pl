use strict;
use warnings;

my %fasta = ();
my $head = "";

open (IN,"<P-tonkinensis_mito.fasta") or die();
open (OUT,">P-tonkinensis_mito.split.fasta") or die();
while (<IN>) {
	$_ =~ s/[\r\n]+//gi;
	if (/>(\S+)/) {$head = $1}
	else {$fasta{$head} .= $_}
}

for (my $i = 0; $i < length $fasta{$head}; $i+=500) {
	print OUT ">",$i+1,"-",$i+1000,"\n";
	print OUT substr($fasta{$head},$i,1000),"\n";
}

if (! -s "P-tonkinensis_mito.split.blastn") {
	system("blastn -query P-tonkinensis_mito.split.fasta -db all_mito -outfmt 6 -perc_identity 70 -evalue 1e-5 -out P-tonkinensis_mito.split.blastn");
}