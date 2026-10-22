use strict;

open (OG,"<ribosome.txt") or die();
my %total_gene = ();
while (<OG>) {
	chomp;
	if (/Gene name/) {next}
	else {
		my @tmp = split /\t/,$_;
		for (my $i = 1; $i < @tmp; $i++) {
			$tmp[$i] =~ s/\s+//g;
			foreach my $gene (split /,/,$tmp[$i]) {push @{$total_gene{$tmp[0]}},$gene}
		}
	}
}

open (IN,"<all_cds.fa2");
my $head = "";
my %fasta = ();
while (<IN>) {
	chomp;
	if (/>(\S+)/) {$head = $1}
	else {$fasta{$head} .= $_}
}

foreach my $gene (sort keys %total_gene) {

if (! -e "./ribosome/$gene/") {system("mkdir ./ribosome/$gene/")}

if (-s "./ribosome/$gene/$gene.raw.tree") {next}

open (OUT,">./ribosome/$gene/$gene.fa");
foreach my $i (@{$total_gene{$gene}}) {
	print OUT ">$i\n",$fasta{$i},"\n";
}

open(IN,"<./ribosome/$gene/$gene.fa");
open(OUT,">./ribosome/$gene/$gene.safe.fa");
my $head = "";
my %fasta = ();
my @order = ();
while (<IN>) {
	chomp;
	if (/>(\S+)/) {$head = $1; push @order,$head}
	else {$fasta{$head} .= $_}
}
foreach my $i (@order) {
	my $j = 0;
	for ($j = 0; $j < length $fasta{$i}; $j+=3) {
		my $string = substr($fasta{$i},$j,3);
		if ($string =~ /TAA/i or $string =~ /TAG/i or $string =~ /TGA/i) {
			print OUT ">$i\n",substr($fasta{$i},0,$j),"\n";
			last;
		}
	}
	if ($j eq length $fasta{$i}) {
		print OUT ">$i\n",substr($fasta{$i},0,$j),"\n";
	}
}

system("seqkit translate ./ribosome/$gene/$gene.safe.fa > ./ribosome/$gene/$gene.pep.fa");
system("mafft --thread 8 --maxiterate 1000 --quiet --localpair ./ribosome/$gene/$gene.pep.fa > ./ribosome/$gene/$gene.pep.fas");
system("pal2nal.pl ./ribosome/$gene/$gene.pep.fas ./ribosome/$gene/$gene.safe.fa -output fasta -nomismatch -nostderr > ./ribosome/$gene/$gene.cds.fas");
system("trimal -in ./ribosome/$gene/$gene.cds.fas -fasta -gt 0.5> ./ribosome/$gene/$gene.cds.trim.fas");
system("cp ./ribosome/$gene/$gene.cds.trim.fas ./ribosome/$gene/$gene.final.fas");
system("~/software/raxml-ng --all --msa ./ribosome/$gene/$gene.final.fas --model GTR+G4 --threads 8 --redo --tree pars{10} --bs-trees 100");
system("cp ./ribosome/$gene/$gene.final.fas.raxml.support ./ribosome/$gene/$gene.raw.tree");
system("rm ./ribosome/$gene/*raxml*");
}