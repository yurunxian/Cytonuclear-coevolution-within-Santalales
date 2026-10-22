use strict;

open(IN,"<OXPHOS.txt");
my @gene = ();
while (<IN>) {
	if (/#/) {next}
	else {
		my @tmp = split /\t/,$_;
		push @gene,$tmp[0];
#		system("hyphy relax --alignment ./OXPHOS/$gene/$gene.final.filtered.fas --tree ./OXPHOS/$gene/$gene.root.rate.TEST.tree -test TEST --srv Yes --starting-points 5 --output ./OXPHOS/$gene/$gene.RELAX.json");
#		system("perl 5_get_omega.pl $gene >> OXPHOS.omega.relax.txt");
#		system("cp ./OXPHOS/$gene/$gene.root.rate.tree ./OXPHOS/$gene/$gene.root.rate.TEST.tree");
#		system("sed -i \"s/)/){TEST}/g\" ./OXPHOS/$gene/$gene.root.rate.TEST.tree");
#		system("hyphy relax --alignment ./OXPHOS/$gene/$gene.final.fas --tree ./OXPHOS/$gene/$gene.root.rate.TEST.tree -test TEST");
	}
}

for (my $i = $ARGV[0]; $i < $ARGV[1]; $i++) {
	my $gene = $gene[$i];
#	print "$gene\n";
#	if (! -s "./ribosome/$gene/$gene.RELAX.VA.json") {system("hyphy relax --alignment ./ribosome/$gene/$gene.cds.trim.fas --tree ./ribosome/$gene/$gene.root.VA.tree -test VA --srv Yes --starting-points 5 --output ./ribosome/$gene/$gene.RELAX.VA.json")}
#	if (! -s "./ribosome/$gene/$gene.RELAX.E.json") {system("hyphy relax --alignment ./ribosome/$gene/$gene.cds.trim.fas --tree ./ribosome/$gene/$gene.root.E.tree -test E --srv Yes --starting-points 5 --output ./ribosome/$gene/$gene.RELAX.E.json")}
#	if (! -s "./ribosome/$gene/$gene.RELAX.L.json") {system("hyphy relax --alignment ./ribosome/$gene/$gene.cds.trim.fas --tree ./ribosome/$gene/$gene.root.L.tree -test L --srv Yes --starting-points 5 --output ./ribosome/$gene/$gene.RELAX.L.json")}
#	if (! -s "./ribosome/$gene/$gene.ABSREL.json") {system("hyphy aBSREL --alignment ribosome/$gene/$gene.cds.trim.fas --tree ribosome/$gene/$gene.root.tree --srv Yes --starting-points 5 --output ribosome/$gene/$gene.ABSREL.json")}
#	system("perl 8_get_test_Tree.OXPHOS.pl $gene");
#	system("perl 2.pl $gene");
#	system("perl bash.ribosome.pl $gene");
	system("perl 5_get_omega.pl $gene");
#	system("perl 9_get_RELAX_data.pl $gene VA");
#	system("perl 9_get_RELAX_data.pl $gene E");
#	system("perl 9_get_RELAX_data.pl $gene L");
#	system("perl 9_get_RELAX_data.pl $gene VATS");
#	system("perl 9_get_RELAX_data.pl $gene VAL");
#	system("perl 9_get_RELAX_data.pl $gene E");
#	system("perl 10_general_descriptive_relax.pl $gene");
}

#my @gene = ("atp1","atp4","atp6","atp8","atp9","ccmC","ccmFC","ccmFN","cob","cox1","cox2","cox3");
#foreach my $gene (@gene) {
#	system("perl 9_get_RELAX_data.pl $gene VA");
#	system("perl 9_get_RELAX_data.pl $gene E");
#	system("perl 9_get_RELAX_data.pl $gene L");
#	system("hyphy busted --alignment ./paml/$gene/$gene.final.fas --tree ./paml/$gene/$gene.root.rate.tree --starting-points 5 --intermediate-fits ./paml/$gene/$gene.file.cache --srv Yes --output ./paml/$gene/$gene.BUSTED.json");
#	system("hyphy busted --alignment ./paml/$gene/$gene.final.fas --tree ./paml/$gene/$gene.root.rate.tree --starting-points 5 --intermediate-fits ./paml/$gene/$gene.file.cache --error-sink Yes --srv Yes --output ./paml/$gene/$gene.BUSTED-E.json");
#	system("hyphy error-filter --json ./paml/$gene/$gene.BUSTED-E.json --output ./paml/$gene/$gene.final.filtered.fas --output-json ./paml/$gene/$gene.filtered-info.json");
#	system("hyphy absrel --alignment ./paml/$gene/$gene.final.filtered.fas --tree ./paml/$gene/$gene.root.rate.tree --srv Yes --starting-points 5 --output ./paml/$gene/$gene.ABSREL.json");
#	system("hyphy relax --alignment ./paml/$gene/$gene.final.filtered.fas --tree ./paml/$gene/$gene.root.rate.TEST.tree -test TEST --srv Yes --starting-points 5 --output ./paml/$gene/$gene.RELAX.json");

#	system("~/software/raxml-ng --evaluate --msa ./paml/$gene/$gene.final.fas --model GTR+G4 --tree ./paml/$gene/Santalales_tree.29.txt --brlen scaled --redo");
#	open (IN,"<./paml/$gene/$gene.final.fas.raxml.bestTree") or die();
#	while (<IN>) {
#		if (/^\(+(\S+?):/) {
#			my $root = $1;
#			system("nw_reroot ./paml/$gene/$gene.final.fas.raxml.bestTree $root | nw_order -c n - > ./paml/$gene/$gene.root.rate.tree");
#		}
#	}
#	system("rm ./paml/$gene/*raxml*");
#	system("hyphy aBSREL --alignment ./paml/$gene/$gene.final.fas --tree ./paml/$gene/$gene.root.rate.tree");
#	system("perl 5_get_omega.mito.pl $gene >> mito.omega.relax.txt");
#	system("cp ./paml/$gene/$gene.root.rate.tree ./paml/$gene/$gene.root.rate.TEST.tree");
#	system("sed -i \"s/)/){TEST}/g\" ./paml/$gene/$gene.root.rate.TEST.tree");
#	system("hyphy relax --alignment ./paml/$gene/$gene.final.fas --tree ./paml/$gene/$gene.root.rate.TEST.tree -test TEST");
#}