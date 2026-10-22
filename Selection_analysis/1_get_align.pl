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

foreach my $gene (sort keys %total_gene) {

	system("~/software/raxml-ng --all --msa ./OXPHOS/$gene/$gene.final.fas --model GTR+G4 --threads 8 --redo --tree pars{10} --bs-trees 100");
	system("cp ./OXPHOS/$gene/$gene.final.fas.raxml.support ./OXPHOS/$gene/$gene.raw.tree");
	system("rm ./OXPHOS/$gene/*raxml*");
}