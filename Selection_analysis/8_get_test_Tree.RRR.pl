use Bio::TreeIO;
use Modern::Perl;
use IO::All;
use Data::Dumper;
use Array::Utils qw(:all);
use List::Util qw(min max);
use List::MoreUtils qw(:all);
use Getopt::Long;

my $gene = $ARGV[0];

my @VA = ("Vov", "Val", "Pri", "Api", "Ach", "Kor", "Pto", "Den", "Dae");
my @E = ("Ene","Ecu");
my @L = ("Scu","Tax");
my @VA_node = ();
my @E_node = ();
my @L_node = ();

############### VA ###############
my $treeio_VA = Bio::TreeIO->new(-format => 'newick',
			      -file => "./RRR/$gene/$gene.root.rate.tree",
			      -internal_node_id => 'bootstrap');
my $tree_VA = $treeio_VA->next_tree;
foreach my $node ($tree_VA->get_nodes) {
	my @gene = ();
	foreach my $i ($node->get_all_Descendents) {if ($i->id) {push @gene,get_sp_name($i->id)}}
	@gene = uniq(@gene);
	if (scalar @gene > 1) {
		my @intersect = intersect(@gene,@VA);
		if (scalar @intersect eq scalar @gene) {
			#print join(", ",@intersect),"\n";
			push @VA_node,$node;
		}
	}
}
foreach my $node (@VA_node) {
	$node->bootstrap("{VA}");
	foreach my $i ($node->get_all_Descendents) {
		$i->bootstrap("{VA}");
		if ($i->id) {
			my $string = $i->id;
			if ($string !~ /{VA}/) {
				$string .= "{VA}";
				$i->id($string);
			}
		}
	}
}
my $out_VA = Bio::TreeIO->new(-format => 'newick',
                           -file   => ">./RRR/$gene/$gene.root.rate.VA.tree",
                           -internal_node_id => 'bootstrap');
$out_VA->write_tree($tree_VA);

############### E ###############
my $treeio_E = Bio::TreeIO->new(-format => 'newick',
			      -file => "./RRR/$gene/$gene.root.rate.tree",
			      -internal_node_id => 'bootstrap');
my $tree_E = $treeio_E->next_tree;
foreach my $node ($tree_E->get_nodes) {
	my @gene = ();
	foreach my $i ($node->get_all_Descendents) {if ($i->id) {push @gene,get_sp_name($i->id)}}
	@gene = uniq(@gene);
	if (scalar @gene > 1) {
		my @intersect = intersect(@gene,@E);
		if (scalar @intersect eq scalar @gene) {
			#print join(", ",@intersect),"\n";
			push @E_node,$node;
		}
	}
}
foreach my $node (@E_node) {
	$node->bootstrap("{E}");
	foreach my $i ($node->get_all_Descendents) {
		$i->bootstrap("{E}");
		if ($i->id) {
			my $string = $i->id;
			if ($string !~ /{E}/) {
				$string .= "{E}";
				$i->id($string);
			}
		}
	}
}
my $out_E = Bio::TreeIO->new(-format => 'newick',
                           -file   => ">./RRR/$gene/$gene.root.rate.E.tree",
                           -internal_node_id => 'bootstrap');
$out_E->write_tree($tree_E);

############### L ###############
my $treeio_L = Bio::TreeIO->new(-format => 'newick',
			      -file => "./RRR/$gene/$gene.root.rate.tree",
			      -internal_node_id => 'bootstrap');
my $tree_L = $treeio_L->next_tree;
foreach my $node ($tree_L->get_nodes) {
	my @gene = ();
	foreach my $i ($node->get_all_Descendents) {if ($i->id) {push @gene,get_sp_name($i->id)}}
	@gene = uniq(@gene);
	if (scalar @gene > 1) {
		my @intersect = intersect(@gene,@L);
		if (scalar @intersect eq scalar @gene) {
			#print join(", ",@intersect),"\n";
			push @L_node,$node;
		}
	}
}
foreach my $node (@L_node) {
	$node->bootstrap("{L}");
	foreach my $i ($node->get_all_Descendents) {
		$i->bootstrap("{L}");
		if ($i->id) {
			my $string = $i->id;
			if ($string !~ /{L}/) {
				$string .= "{L}";
				$i->id($string);
			}
		}
	}
}
my $out_L = Bio::TreeIO->new(-format => 'newick',
                           -file   => ">./RRR/$gene/$gene.root.rate.L.tree",
                           -internal_node_id => 'bootstrap');
$out_L->write_tree($tree_L);


sub get_sp_name {
	my $gene = shift;
	if ($gene =~ /^(\S+?)\_/) {return $1}
}