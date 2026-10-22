use strict;

my $dir = "/data/yurx/work/Amphorogynaceae/dNdS/ribosome";

foreach my $i (glob "$dir/*") {
	if ($i =~ /$dir\/(\S+?)$/) {
		my $gene = $1;
		if (-s "$i/$gene.free.codeml") {
			print $gene,"\n";
			open (IN,"<$i/$gene.free.codeml") or die();
			open (OUT,">tmp.dNdS.tree") or die();
			my @codeml = <IN>;
			my $output = "";
			for (my $j = 0; $j < @codeml; $j++) {
				if ($codeml[$j] eq "w ratios as labels for TreeView:\n") {$output = $codeml[$j+1]}
			}
			$output =~ s/ #/:/gi;
			$output =~ s/([a-zA-Z]+)/$1_$gene/gi;
			print OUT $output;
			system("perl get_branch_len.pl -i Santalales_tree.tree -g tmp.dNdS.tree -n $gene -b Santalales_tree.txt -o ribosome.dNdS.txt");
		}
	}
}