
use IO::All;

my ($start,$end) = ($ARGV[0],$ARGV[1]);
my @txt = io("Orthogroups.list")->chomp->getlines;

for (my $i = $start; $i < $end+1; $i++) {
	my $orthogroup = $txt[$i];
	if ($orthogroup =~ /#/) {next}
	print $orthogroup,"\n";
#	system("perl 1_align_trim_tree.pl $orthogroup");
#	system("perl 2_remove_long_branch.pl -i $orthogroup -b 50");
#	system("perl 3_duplication.pl -i $orthogroup -s 2 -print_lineage");
#	system("perl 4_reconciliation_optimization.pl $orthogroup");
#	system("perl 5_get_branch_len.pl -i Santalales_tree.txt -g ./Reconciled_tree/$orthogroup.branch_opt.root.tree -n $orthogroup -b Santalales_tree.node -o OG.Branch_length.txt");
}