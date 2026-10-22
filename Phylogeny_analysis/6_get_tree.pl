use strict;
use warnings;

open (BUSCO,"<embryophyta_odb10.list") or die();

while (<BUSCO>) {
	my $gene = "";
	if (/(\S+)/) {$gene = $1}
	if (! -e "./BUSCO_tree/$gene/") {
		system("mkdir ./BUSCO_tree/$gene/");
		system("cp ./aln/$gene.final.fas ./BUSCO_tree/$gene/$gene.final.cds.fas");
#		system("seqkit translate -x ./BUSCO_tree/$gene/$gene.final.cds.fas > ./BUSCO_tree/$gene/$gene.final.pep.fas");
	}

	if (! -s "./BUSCO_tree/$gene/$gene.final.cds.fas.raxml.support") {
		system("~/software/raxml-ng -all --msa ./BUSCO_tree/$gene/$gene.final.cds.fas --model GTR+G4 --tree pars{10} --bs-trees 100 --threads 8 --redo");
#		system("~/software/raxml-ng -all --msa ./BUSCO_tree/$gene/$gene.final.pep.fas --model LG+G4 --tree pars{10} --bs-trees 100 --threads 8 --redo");
	}

#	if (! -s "./BUSCO_tree/$gene/$gene.final.fas.out") {
#		system("modeltest-ng -i ./BUSCO_tree/$gene/$gene.final.fas --force -h uigf −f ef -T raxml -p 4 -d nt");
#	}

#	if (! -s "./BUSCO_tree/$gene/$gene.final.fas.raxml.support") {
#		open(MODEL,"<./BUSCO_tree/$gene/$gene.final.fas.out") or die();
#		my $raxml_model = "";
#		while (<MODEL>) {
#			chomp;
#			if (/raxml-ng --msa.+ --model (\S+)/) {$raxml_model = $1}
#		}
#		system("~/software/raxml-ng -all --msa ./BUSCO_tree/$gene/$gene.final.fas --model $raxml_model --tree pars{10} --bs-trees 100 --threads 8 --redo");
#	}
}