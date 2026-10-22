use strict;

my $orthogroup = $ARGV[0];

# get final alignment based on the reconciled tree with reduced gene list
if (! -s "./Tree/$orthogroup/$orthogroup.final.fas") {
	open (TREE,"<./Tree/$orthogroup/$orthogroup.50.clean.tree");
	open (OUT,">./Tree/$orthogroup/$orthogroup.final.fas");
	my @tmp = <TREE>;
	my $string = shift @tmp;
	$string =~ s/[(),;]/ /gi;
	my @gene = split /\s+/,$string;
	print scalar(@gene),"\n";
	open (FAS,"<./Tree/$orthogroup/$orthogroup.trim.cds.fas");
	my %fasta = ();
	my $head = "";
	while (<FAS>) {
		chomp;
		if (/>(\S+)/) {$head = $1} else {$fasta{$head} .= $_}
	}
	foreach my $gene (@gene) {if (exists $fasta{$gene}) {print OUT ">$gene\n",$fasta{$gene},"\n"}}
}

# optimize branch lengths in the reconciled tree
if (! -s "./Tree/$orthogroup/$orthogroup.branch_opt.tree") {
	system("~/software/raxml-ng --evaluate --msa ./Tree/$orthogroup/$orthogroup.final.fas --model GTR+G4 --tree ./Tree/$orthogroup/$orthogroup.50.clean.tree --brlen scaled --redo");
	system("mv ./Tree/$orthogroup/$orthogroup.final.fas.raxml.bestTree ./Tree/$orthogroup/$orthogroup.branch_opt.tree");
	system("rm ./Tree/$orthogroup/*raxml*");
}

# root the final tree 
if (! -s "./Tree/$orthogroup/$orthogroup.branch_opt.root.tree") {
	system("java -jar Notung-3.0_24-beta.jar -g ./Tree/$orthogroup/$orthogroup.branch_opt.tree -s Santalales_tree.txt --edgeweights name --root --speciestag prefix --treeoutput newick --nolosses --silent");
	open (TREE,"<$orthogroup.branch_opt.tree.root.0.nwk") or die($orthogroup);
	open (TREEO,">tmp.tree") or die();
	while (<TREE>) {
		chomp;
		if (/(\S+)/) {
			my $string = $1;
			$string =~ s/\)n\d+/)/gi;
			$string =~ s/\)r\d+/)/gi;
			$string =~ s/\)\d+\.\d+/)/gi;
	#		$string =~ s/\:\d+\.\d+//gi;
	#		$string =~ s/e-\d+//gi;
			print TREEO $string,"\n";
			last;
		}
	}
	close TREE;
	close TREEO;
	system("mv tmp.tree ./Tree/$orthogroup/$orthogroup.branch_opt.root.tree");
	system("rm ./*nwk");
}