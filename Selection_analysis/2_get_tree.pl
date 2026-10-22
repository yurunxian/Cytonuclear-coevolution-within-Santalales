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
my $bootstrap_threshold = 80;
foreach my $gene (sort keys %total_gene) {
	system("java -jar Notung-3.0_24-beta.jar -g ./ribosome/$gene/$gene.raw.tree -s Santalales_tree.tree --edgeweights name --rearrange --speciestag prefix --threshold $bootstrap_threshold --treeoutput newick --nolosses --silent");
	system("java -jar Notung-3.0_24-beta.jar -g $gene.raw.tree.rear.0.nwk -s Santalales_tree.tree --edgeweights name --root --speciestag prefix --threshold $bootstrap_threshold --treeoutput newick --nolosses --silent");
	open (TREE,"<$gene.raw.tree.rear.0.root.0.nwk") or die($gene);
	open (TREEO,">tmp.tree") or die();
	while (<TREE>) {
		chomp;
		if (/(\S+)/) {
			my $string = $1;
			$string =~ s/\)n\d+/)/gi;
			$string =~ s/\)r\d+/)/gi;
			$string =~ s/\)\d+\.\d+/)/gi;
			$string =~ s/\:\d+\.\d+//gi;
			$string =~ s/e-\d+//gi;
			print TREEO $string,"\n";
			last;
		}
	}
	close TREE;
	close TREEO;
	system("mv tmp.tree ./ribosome/$gene/$gene.root.tree");
	system("rm ./*nwk");
}