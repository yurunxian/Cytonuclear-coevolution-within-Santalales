#!/usr/bin/env perl
use Bio::TreeIO;
use Modern::Perl;
use IO::All;
use Data::Dumper;
use Array::Utils qw(:all);
use List::Util qw(min max);
use List::MoreUtils qw(:all);
use Getopt::Long;

#define the parameters
my $bootstrap_threshold = 50;
my $lineage_number = 13;
my $force_VA = "";
my $orthogroup = "";
my $shared = 2;
my $print_lineage = "";
my @AV = ("Kor","Val","Vov","Ach","Api","Dae","Den","Pri","Pto");
my @order = ("Amb","Ari","Ara","Ery","Mal","Rph","Bya","Hac","Mco","Tma","Scu","Tax","Cha","The","Pyr","Scl","Ecu","Ene","San","Osy","Dae","Den","Pri","Pto","Kor","Val","Vov","Ach","Api");
my %clade = ("Amb"=>"1","Ari"=>"1","Ara"=>"1","Ery"=>"2","Mal"=>"3","Rph"=>"4","Bya"=>"4","Hac"=>"5","Mco"=>"5","Tma"=>"5","Scu"=>"5","Tax"=>"5","Cha"=>"6","The"=>"7","Pyr"=>"7","Scl"=>"7",
	"Ecu"=>"8","Ene"=>"8","San"=>"8","Osy"=>"8","Dae"=>"9","Den"=>"9","Pri"=>"9","Pto"=>"9","Kor"=>"10","Val"=>"11","Vov"=>"11","Ach"=>"12","Api"=>"13");

GetOptions ("force_VA" => \$force_VA,
			"i=s" => \$orthogroup,
			"b=i" => \$bootstrap_threshold,
			"s=i" => \$shared,
			"print_lineage" => \$print_lineage);
if ($force_VA eq 1) {$force_VA = 1} else {$force_VA = 0}
if ($print_lineage) {open (GENE,">>./Duplication_analysis/lineage.node14.$shared.txt")}
if ($orthogroup !~ /OG/) {die("Please specify orthogroup!")}
open (OUT,">>./Duplication_analysis/duplication.node14.$shared.txt");


#if (! -s "/data/yurx/work/Amphorogynaceae/Tree_0414/Tree/$orthogroup/$orthogroup.$bootstrap_threshold.clean.tree") {
#	system("java -jar Notung-3.0_24-beta.jar -g /data/yurx/work/Amphorogynaceae/Tree_0414/Tree/$orthogroup/$orthogroup.clean.tree -s Santalales_tree.txt --edgeweights name --rearrange --speciestag prefix --threshold $bootstrap_threshold --treeoutput newick --nolosses --silent");
#	system("java -jar Notung-3.0_24-beta.jar -g $orthogroup.clean.tree.rear.0.nwk -s Santalales_tree.txt --edgeweights name --root --speciestag prefix --threshold $bootstrap_threshold --treeoutput newick --nolosses --silent");
#	open (TREE,"<$orthogroup.clean.tree.rear.0.root.0.nwk") or die($orthogroup);
#	open (TREEO,">tmp.tree") or die();
#	while (<TREE>) {
#		chomp;
#		if (/(\S+)/) {
#			my $string = $1;
#			$string =~ s/\)n\d+/)/gi;
#			$string =~ s/\)r\d+/)/gi;
#			$string =~ s/\)\d+\.\d+/)/gi;
#			$string =~ s/\:\d+\.\d+//gi;
#			$string =~ s/e-\d+//gi;
#			print TREEO $string,"\n";
#			last;
#		}
#	}
#	close TREE;
#	close TREEO;
#	system("mv tmp.tree /data/yurx/work/Amphorogynaceae/Tree_0414/Tree/$orthogroup/$orthogroup.$bootstrap_threshold.clean.tree");
#	system("rm ./*nwk");
#}

#read the newick tree
if (! -e "./Reconciled_tree/$orthogroup.branch_opt.root.tree") {die($orthogroup,"\tno file\n")}
my $treeio = Bio::TreeIO->new(-format => 'newick',
			      -file => "./Reconciled_tree/$orthogroup.branch_opt.root.tree",
			      -internal_node_id => 'bootstrap');
my $tree = $treeio->next_tree;

#reroot the tree at a random leaf tip of the OUTTEST species
my @leaf_all = ();
my %leaf_id = ();
my $root_clade = 0;
foreach my $node ($tree->get_leaf_nodes) {push @leaf_all,$node->id; $leaf_id{$node->id} = $node}
foreach my $clade (@order){
	if (grep {$_ =~ /^$clade/} @leaf_all) {
		my @outgroup = grep {$_ =~ /^$clade/} @leaf_all;
		@outgroup = sort {$a cmp $b} @outgroup;
		$tree->reroot($leaf_id{$outgroup[0]});
#		say "Reroot the tree at ",$outgroup[0];
		$root_clade = $clade{get_sp_name($outgroup[0])};
		last;
	}
}

#get nodes with a >= "$bootstrap_threshold" bootstrap
#my @confident_nodes = grep {(($_->bootstrap)and($_->bootstrap > $bootstrap_threshold))} $tree->get_nodes;
my @confident_nodes = $tree->get_nodes;

#main function
my %duplicated_gene = ();
my @duplication_result = ();
for (my $i = 0; $i < $lineage_number-1; $i++) {$duplication_result[$i] = 0}
foreach my $node (@confident_nodes){
	#find confident nodes with at least four Descendents which include leaves from at least two clades
	my @children = $node->get_all_Descendents;
	my ($clade_number_first,$min_clade_first) = get_clade_number($node);
	if (($min_clade_first <= $root_clade)or($clade_number_first < 2)) {next}
	my @confident_children = intersect(@confident_nodes,@children);
	my @new_confident_children = ();
	for (my $i = 0; $i < scalar @confident_children; $i++) {
		my ($clade_number,$min_clade) = get_clade_number($confident_children[$i]);
		#avoid paraphylic outgroup caused by the reroot processing
		if (($min_clade > $root_clade)and($clade_number >= 2)) {push @new_confident_children,$confident_children[$i]}
	}
	if (scalar @new_confident_children eq 0) {next}

	#find confident child nodes whose leaves share genes from at least two sp with any other leaves
	#find the node whose leaves contain the most ancient lineage
	my @duplication_result_tmp = ();
	for (my $i = 0; $i < $lineage_number-1; $i++) {$duplication_result_tmp[$i] = 0}
	my %duplicate_child_clade = ();
	foreach my $child (@new_confident_children) {
		my @child_Descendents = $child->get_all_Descendents;
		my @other_Descendents = array_minus(@children,@child_Descendents);
		my ($duplicate_clade,$earliest_clade,$shared_sp_number) = duplication(\@child_Descendents,\@other_Descendents);
		if ($shared_sp_number >= $shared and $duplicate_clade >= $force_VA) {
			$duplication_result_tmp[$earliest_clade-1]++;
			push @{$duplicate_child_clade{$earliest_clade-1}},$child;
		}
	}

	#reduce the redundant node
	my $earliest_duplicate_clade = 0;
	for (my $i = 0; $i < $lineage_number-1; $i++) {
		if ($duplication_result_tmp[$i] > 0) {$earliest_duplicate_clade = $i;last;}
	}
	#find internal node
	my @internal_node = ();
	foreach my $child (@{$duplicate_child_clade{$earliest_duplicate_clade}}) {
		my $path = get_path_between_nodes($node,$child);
		my @confident_internal_node = intersect(@$path,@confident_children);
		if (@confident_internal_node eq 0) {
			my ($node_early,$child_early,$both_B_V) = get_lineage_number($node,$child,$earliest_duplicate_clade,0);
#			if ($node_early eq $child_early or $node_early+$child_early eq $lineage_number * 2 - 1) {
			if ($node_early eq $child_early or $node_early+$child_early >= 22) {
#				if ($both_B_V > 0) {
					$duplication_result[$earliest_duplicate_clade] += 1;
					if ($print_lineage) {get_lineage_number($node,$child,$earliest_duplicate_clade,1)}
#				}
			}
			last;
		}
	}
}

#output
print $orthogroup,"\n";
print OUT $orthogroup;
for (my $i = 11; $i < $lineage_number-1; $i++) {print OUT "\t",$duplication_result[$i]}
print OUT "\n";

sub get_sp_name {
	my $gene = shift;
	if ($gene =~ /^(\S+?)\_/) {return $1}
}

sub get_clade_number {
	my $node = shift;
	my @children = $node->get_all_Descendents;
	my %children_leaf = ();
	my $leaf_number = 0;
	for (@children) {
		if ($_->id) {
			$leaf_number++;
			$children_leaf{$clade{get_sp_name($_->id)}}++;
		}
	}
	my @keys = keys %children_leaf;
	my ($clade_number,$min_clade) = (0,0,);
	$clade_number = scalar @keys;
	if ($clade_number > 0) {
		my @tmp = sort {$a <=> $b} keys %children_leaf;
		$min_clade = $tmp[0];
	}
	return ($clade_number,$min_clade);
}

sub duplication {
	my ($tmp1,$tmp2) = @_;
	my @nearest = @$tmp1;
	my @other = @$tmp2;
	my (%nearest_leaf,%other_leaf);
	for (@nearest) {if ($_->id) {$nearest_leaf{get_sp_name($_->id)}++}}
	for (@other) {if ($_->id) {$other_leaf{get_sp_name($_->id)}++}}
	my ($duplicate_clade,$earliest_clade,$shared_sp_number) = (0,0,0);
	for(@AV) {if (exists $nearest_leaf{$_}) {$duplicate_clade += 1; last;}}
	for(@AV) {if (exists $other_leaf{$_}) {$duplicate_clade += 1; last;}}
	for(@order) {if ((exists $nearest_leaf{$_})and(exists $other_leaf{$_})) {$shared_sp_number++}}
	for(@order) {if (exists $nearest_leaf{$_}) {$earliest_clade = $clade{$_}; last;}}
	return ($duplicate_clade,$earliest_clade,$shared_sp_number);
}

sub get_tips {
	my $node = shift;
	my @tips = $node->get_all_Descendents;
	my @output = ();
	for (@tips) {if ($_->id) {push @output,$_->id}}
	return \@output;
}

sub get_path_between_nodes {
	# by default, node1 is deeper than node2
	my ($node1,$node2) = @_;
	my @path = ();
	my @tmp = ();
	while (!grep {$_ eq $node2} @tmp) {
		@tmp = $node1->each_Descendent;
		foreach my $i (@tmp) {
			if (grep {$_ eq $node2} $i->get_all_Descendents) {
				push @path,$i;
				$node1 = $i;
			}
		}
	}	
	return \@path;
}

sub get_lineage_number {
	my ($node,$child,$earliest_duplicate_clade,$print) = @_;
	my @child_result = ();
	my @node_result = ();
	my $child_tip = get_tips($child);
	my @child_tip = @$child_tip;
	my $node_tip = get_tips($node);
	my @node_tip = array_minus(@$node_tip,@child_tip);
	my %child_tip = ();
	my %node_tip = ();
	my @child_tip_gene = ();
	my @node_tip_gene = ();
	for (my $i = 0; $i < scalar @node_tip; $i++) {push @node_tip_gene,$node_tip[$i]; $node_tip[$i] = get_sp_name($node_tip[$i])}
	for (my $i = 0; $i < scalar @child_tip; $i++) {push @child_tip_gene,$child_tip[$i]; $child_tip[$i] = get_sp_name($child_tip[$i])}
	@child_tip = uniq(@child_tip);
	@node_tip = uniq(@node_tip);
	foreach my $i (@node_tip) {$node_tip{$clade{$i}}++}
	foreach my $i (@child_tip) {$child_tip{$clade{$i}}++}
	for (my $i = 1; $i < $lineage_number+1; $i++) {
		if (exists $node_tip{$i}) {push @node_result, $node_tip{$i}}
		else{push @node_result, "0"}
		if (exists $child_tip{$i}) {push @child_result, $child_tip{$i}}
		else{push @child_result, "0"}
	}
	if ($print eq 1) {
		if ($earliest_duplicate_clade+3 eq 14) {
			print GENE $orthogroup,"\t",$earliest_duplicate_clade+3,"\t",scalar(@node_tip),"\t",join("-",@node_result),"\t",join(", ",sort @node_tip_gene),"\t",scalar(@child_tip),"\t",join("-",@child_result),"\t",join(", ",sort @child_tip_gene),"\n";
		}		
	}
	else{
		my $node_early = 0;
		my $child_early = 0;
		my $both_B_V = 0;
		for (my $i = 0; $i < scalar @node_result; $i++) {if ($node_result[$i] > 0) {$node_early = $i+1; last}}
		for (my $i = 0; $i < scalar @child_result; $i++) {if ($child_result[$i] > 0) {$child_early = $i+1; last}}
		if ($child_result[-1] > 0 and $child_result[-2] > 0) {$both_B_V++}
		return $node_early,$child_early,$both_B_V;
	}
#	return join("-",@node_result),join("-",@child_result);
}