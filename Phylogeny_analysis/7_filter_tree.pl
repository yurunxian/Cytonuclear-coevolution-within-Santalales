use strict;
use warnings;
use Bio::TreeIO;
use Modern::Perl;
use Data::Dumper;
use Array::Utils qw(:all);
use List::Util qw(min max);
use List::MoreUtils qw(:all);
use Getopt::Long;

my $bootstrap_threshold = 70;
my @order = ("Amb","Ara","Ari","Ery","Mal","Bya","Rph","Hac","Mco","Tma","Scu","Tax","Cha","The","Pyr","Scl",
	"Ecu","Ene","Sya","Osy","Ach","Api","Kor","Val","Vov","Dae","Den","Pri","Pto");
my %clade = ("Amb"=>"0","Ara"=>"0","Ari"=>"0","Ery"=>"1","Mal"=>"1","Mco"=>"2","Tma"=>"2","Scu"=>"2","Tax"=>"2",
	"Cha"=>"3","The"=>"3","Pyr"=>"3","Scl"=>"3","Ecu"=>"4","Ene"=>"4","Sya"=>"4","Osy"=>"4",
	"Ach"=>"5","Api"=>"5","Kor"=>"5","Val"=>"5","Vov"=>"5","Dae"=>"6","Den"=>"6","Pri"=>"6","Pto"=>"6");

open (BUSCO,"<embryophyta_odb10.list") or die();

print "Gene\tSantales\tLoranthaceae\tSantalaceae\tViscaceae\tAmphorogynaceae\n";

while (<BUSCO>) {
	my $gene = "";
	if (/(\S+)/) {$gene = $1}
	print "$gene\t";

	my $treeio = Bio::TreeIO->new(-format => 'newick',
				      -file => "./BUSCO_tree/$gene/$gene.final.cds.fas.raxml.support",
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
			$root_clade = $clade{$outgroup[0]};
			last;
		}
	}

	#Santalales
	my $Santales = "Fail";
	my @node = ();
	my $count = 0;
	foreach my $node ($tree->get_leaf_nodes) {if ($clade{$node->id} > 0) {push @node, $node}}
	my $lca = $tree->get_lca(-nodes => \@node);
	foreach my $node ($lca->get_all_Descendents) {if ($node->id and $clade{$node->id} eq 0) {$count++}}
	if ($count eq 0 and $lca->bootstrap >= $bootstrap_threshold) {$Santales = "Pass"}

	#Loranthaceae
	my $Loranthaceae = "Fail";
	@node = ();
	$count = 0;
	foreach my $node ($tree->get_leaf_nodes) {if ($clade{$node->id} eq 2) {push @node, $node}}
	if (scalar @node < 2) {$Loranthaceae = "Pass"}
	else {
		$lca = $tree->get_lca(-nodes => \@node);
		foreach my $node ($lca->get_all_Descendents) {if ($node->id and $clade{$node->id} ne 2) {$count++}}
		if ($count eq 0 and $lca->bootstrap >= $bootstrap_threshold) {$Loranthaceae = "Pass"}
	}


	#Santalaceae
	my $Santalaceae = "Fail";
	@node = ();
	$count = 0;
	foreach my $node ($tree->get_leaf_nodes) {if ($clade{$node->id} eq 4) {push @node, $node}}
	if (scalar @node < 2) {$Santalaceae = "Pass"}
	else {
		$lca = $tree->get_lca(-nodes => \@node);
		foreach my $node ($lca->get_all_Descendents) {if ($node->id and $clade{$node->id} ne 4) {$count++}}
		if ($count eq 0 and $lca->bootstrap >= $bootstrap_threshold) {$Santalaceae = "Pass"}
	}

	#Viscaceae
	my $Viscaceae = "Fail";
	@node = ();
	$count = 0;
	foreach my $node ($tree->get_leaf_nodes) {if ($clade{$node->id} eq 5) {push @node, $node}}
	if (scalar @node < 2) {$Viscaceae = "Pass"}
	else {
		$lca = $tree->get_lca(-nodes => \@node);
		foreach my $node ($lca->get_all_Descendents) {if ($node->id and $clade{$node->id} ne 5) {$count++}}
		if ($count eq 0 and $lca->bootstrap >= $bootstrap_threshold) {$Viscaceae = "Pass"}
	}

	#Amphorogynaceae
	my $Amphorogynaceae = "Fail";
	@node = ();
	$count = 0;
	foreach my $node ($tree->get_leaf_nodes) {if ($clade{$node->id} eq 6) {push @node, $node}}
	if (scalar @node < 2) {$Amphorogynaceae = "Pass"}
	else {
		$lca = $tree->get_lca(-nodes => \@node);
		foreach my $node ($lca->get_all_Descendents) {if ($node->id and $clade{$node->id} ne 6) {$count++}}
		if ($count eq 0 and $lca->bootstrap >= $bootstrap_threshold) {$Amphorogynaceae = "Pass"}
	}

	print "$Santales\t$Loranthaceae\t$Santalaceae\t$Viscaceae\t$Amphorogynaceae\n";
}