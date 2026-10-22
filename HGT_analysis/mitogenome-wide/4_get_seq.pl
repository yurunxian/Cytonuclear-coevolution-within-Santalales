use strict;
use warnings;

my $head = "";
my $fasta = ();
open (FAS,"<P-tonkinensis_mito.fasta") or die();
open (IN,"<P-tonkinensis_mito.candidate") or die();
open (OUT,">P-tonkinensis_mito.candidate.fasta") or die();

while (<FAS>) {
	$_ =~ s/[\r\n]+//gi;
	if (/>(\S+)/) {$head = $1}
	else {$fasta .= $_}
}

while (<IN>) {
	chomp;
	my @tmp = split /\t/,$_;
	print OUT ">$tmp[0]\n",substr($fasta,$tmp[2]-1,$tmp[1]),"\n";
}

system("makeblastdb -in P-tonkinensis_mito.candidate.fasta -dbtype nucl -out P-tonkinensis_mito.candidate");
system("blastn -query P-tonkinensis_mito.candidate.fasta -db P-tonkinensis_mito.candidate -outfmt 6 -evalue 1e-5 -perc_identity 100 -out P-tonkinensis_mito.candidate.selfblastn");