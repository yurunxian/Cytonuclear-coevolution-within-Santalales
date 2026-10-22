use strict;
use warnings;

open (IN,"<P-tonkinensis_mito.candidate.final") or die();
open (FAS,"<P-tonkinensis_mito.fasta") or die();
open (MITO,"<all_mito.fasta") or die();

my $head = "";
my %fasta = ();

while (<FAS>) {
	chomp;
	if (/>(\S+)/) {$head = $1}
	else {$fasta{$head} .= $_}
}

while (<MITO>) {
	chomp;
	if (/>(\S+)/) {$head = $1}
	else {$fasta{$head} .= $_}
}

while (<IN>) {
	chomp;
	my @tmp = split /\t/,$_;
	my $name = $tmp[0];
	my @cord = split /, /,$tmp[2];
	my ($start,$end) = split /-/,$cord[0];
	open (OUT,">./P-tonkinensis/$name.fasta") or die();
	open (BLASTO,">>./P-tonkinensis/$name.fasta") or die();
	print OUT ">$name\n",substr($fasta{"P-tonkinensis"},$start-1,$end-$start+1),"\n";
	if (! -s "./P-tonkinensis/$name.blastn") {
		system("blastn -query ./P-tonkinensis/$name.fasta -db all_mito -outfmt 6 -perc_identity 80 -evalue 1e-10 -out ./P-tonkinensis/$name.blastn");
	}
	open (BLAST,"<./P-tonkinensis/$name.blastn") or die();
	open (BLASTF,">./P-tonkinensis/$name.blastn.filter") or die();
	while (<BLAST>) {
		my @blast = split /\t/,$_;
		if ($blast[3] > ($end-$start+1) * 0.8) {
			print BLASTF $_;
			print BLASTO ">$blast[1]_$blast[8]-$blast[9]\n";
			if ($blast[8] < $blast[9]) {
				print BLASTO substr($fasta{$blast[1]},$blast[8]-1,$blast[9]-$blast[8]+1),"\n";
			}
			else {
				print BLASTO revcom(substr($fasta{$blast[1]},$blast[9]-1,$blast[8]-$blast[9]+1)),"\n";
			}
		}
	}
	system("rm ./P-tonkinensis/$name.blastn");
#	if (! -s "./P-tonkinensis/$name.fas") {
#		system("mafft --maxiterate 1000 --localpair --thread 4 ./P-tonkinensis/$name.fasta > ./P-tonkinensis/$name.fas");
#	}

#	if (! -s "./P-tonkinensis/$name.fas") {
#		system("~/software/raxml-ng --all --msa ./P-tonkinensis/$name.fas --model GTR+G4 --tree pars{10} --bs-trees 1000 --threads 20");
#	}
}

sub revcom {
	my $string = shift;
	my $revcom=reverse($string);
	$revcom =~ tr/ACGTacgt/TGCAtgca/;
	return $revcom;
}