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
	system("~/software/raxml-ng --evaluate --msa ./ribosome/$gene/$gene.cds.trim.fas --model GTR+G4 --tree ./ribosome/$gene/$gene.root.tree --brlen scaled --redo");
	open (IN,"<./ribosome/$gene/$gene.cds.trim.fas.raxml.bestTree") or die();
	while (<IN>) {
		if (/^\(+(\S+?):/) {
			my $root = $1;
			system("nw_reroot ./ribosome/$gene/$gene.cds.trim.fas.raxml.bestTree $root | nw_order -c n - > ./ribosome/$gene/$gene.root.rate.tree");
		}
	}
	system("rm ./ribosome/$gene/*raxml*");
}
