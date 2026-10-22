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
	system("hyphy aBSREL --alignment ribosome/$gene/$gene.cds.trim.fas --tree ribosome/$gene/$gene.root.tree --srv Yes --starting-points 5 --output ribosome/$gene/$gene.ABSREL.json");
}