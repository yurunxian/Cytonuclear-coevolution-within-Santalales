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
my $outdir = "/data/yurx/work/Amphorogynaceae/Tree_0414/Tree";
my $bootstrap_threshold = 50;
my $orthogroup = "";
my @order = ("Amb","Ari","Ara","Ery","Mal","Rph","Bya","Hac","Mco","Tma","Scu","Tax","Cha","The","Scl","Pyr","Ene","Ecu","Osy","San",
	"Dae","Den","Pri","Pto","Kor","Ach","Api","Val","Vov");
my %clade = ("Amb"=>"1","Ari"=>"1","Ara"=>"1","Ery"=>"2","Mal"=>"2","Rph"=>"2","Bya"=>"2","Hac"=>"3","Mco"=>"3","Tma"=>"3","Scu"=>"3","Tax"=>"3",
	"Cha"=>"4","The"=>"4","Pyr"=>"4","Scl"=>"4","Ecu"=>"5","Ene"=>"5","San"=>"5","Osy"=>"5","Dae"=>"6","Den"=>"6","Pri"=>"6","Pto"=>"6",
	"Kor"=>"7","Val"=>"7","Vov"=>"7","Ach"=>"7","Api"=>"7");
my %type = ("Amb"=>"genome","Ari"=>"genome","Ara"=>"genome","Ery"=>"genome","Mal"=>"genome","Rph"=>"transcriptome","Bya"=>"transcriptome",
	"Hac"=>"transcriptome","Mco"=>"transcriptome","Tma"=>"transcriptome","Scu"=>"genome","Tax"=>"genome","Cha"=>"transcriptome","The"=>"transcriptome",
	"Pyr"=>"transcriptome","Scl"=>"transcriptome","Ecu"=>"transcriptome","Ene"=>"transcriptome","San"=>"genome","Osy"=>"genome","Dae"=>"transcriptome",
	"Den"=>"transcriptome","Pri"=>"transcriptome","Pto"=>"transcriptome","Kor"=>"transcriptome","Val"=>"transcriptome","Vov"=>"transcriptome",
	"Ach"=>"transcriptome","Api"=>"transcriptome");

GetOptions ("i=s" => \$orthogroup,
			"b=i" => \$bootstrap_threshold);
if ($orthogroup !~ /OG/) {die("Please specify orthogroup!")}

#read the newick tree
if (! -e "$outdir/$orthogroup/$orthogroup.raw.tree") {die($orthogroup,"\tno file\n")}
my $treeio = Bio::TreeIO->new(-format => 'newick',
			      -file => "$outdir/$orthogroup/$orthogroup.raw.tree",
			      -internal_node_id => 'bootstrap');
my $tree = $treeio->next_tree;
open (FAS,"<$outdir/$orthogroup/$orthogroup.trim.cds.fas") or die();
open (OUT,">$outdir/$orthogroup/$orthogroup.trim.clean.cds.fas") or die();
open (TXT,">$outdir/$orthogroup/$orthogroup.trim.clean.cds.txt") or die();
my $head = "";
my @fasta_order = ();
my %fasta = ();
while (<FAS>) {chomp; if (/>(\S+)/) {$head = $1; push @fasta_order,$head} else {$fasta{$head} .= $_}}
my %leaf_id = ();
my %leaf_len = ();
my %removed_gene = ();
foreach my $node ($tree->get_leaf_nodes) {$leaf_id{$node->id} = $node; $leaf_len{$node->id} = $node->branch_length}
my @confident_nodes = grep {(($_->bootstrap)and($_->bootstrap > $bootstrap_threshold))} $tree->get_nodes;

#find redundancy according to the topology
foreach my $node (@confident_nodes) {
	my $all_tips = get_tips($node);
	my %species = ();
	my %len = ();
	foreach my $gene (@$all_tips) {my $sp = get_sp_name($gene); $species{$sp}++}
	if (scalar keys %species ne 1) {next}
	foreach my $gene (@$all_tips) {
		my $string = $fasta{$gene};
		$string =~ s/-//g;
		$len{$gene} = length ($fasta{$gene}) - length ($string);
	}
	my @removed_gene = sort {$len{$a} <=> $len{$b}} keys %len;
	shift @removed_gene;
#	print join(", ",@removed_gene),"\n";
	foreach my $gene (@removed_gene) {$removed_gene{$gene} = "Species-specific duplicated gene or redundant transcripts"}
}

#find contamination
foreach my $node (@confident_nodes) {
	my $all_tips = get_tips($node);
	my %sp_clade = ();
	my %data_type = ();
	my @gene = ();
	foreach my $gene (@$all_tips) {if (!exists $removed_gene{$gene}) {push @gene, $gene; my $sp = get_sp_name($gene); $sp_clade{$clade{$sp}}++; $data_type{$type{$sp}}++}}
	if (! exists $data_type{"genome"} and scalar keys %sp_clade >= 4) {
		foreach my $gene (@gene) {if (!exists $removed_gene{$gene}) {$removed_gene{$gene} = "Very possible contamination";print "$orthogroup\t$gene\tContamination!\n"}}
	}
#	print "\n";
}


#find tips with exceptionally large branch length
my @len = ();
foreach my $node (sort keys %leaf_id) {if (!exists $removed_gene{$node}) {push @len,$leaf_len{$node}}}
@len = sort {$a <=> $b} @len;
#my $len_threshold = $len[int(scalar(@len) * 0.75)] * 2.5 - $len[int(scalar(@len) * 0.25)] * 1.5;
my $len_threshold = $len[int(scalar(@len) * 0.75)] * 3;
#print $len[int(scalar(@len) * 0.75)],"\t",$len_threshold,"\n";
foreach my $node (sort keys %leaf_id) {
	if (!exists $removed_gene{$node} and $leaf_len{$node} >= $len_threshold and $clade{get_sp_name($node)} > 1) {$removed_gene{$node} = "Exceptionally large branch length"}
}


#output
foreach my $gene (@fasta_order) {if (!exists $removed_gene{$gene}) {print OUT ">$gene\n$fasta{$gene}\n"} else {print TXT $gene,"\t",$removed_gene{$gene},"\n"}}
close FAS;
close OUT;

#tree reconstruction
if (! -s "$outdir/$orthogroup/$orthogroup.clean.tree") {
	system("~/software/raxml-ng -all --msa $outdir/$orthogroup/$orthogroup.trim.clean.cds.fas --model GTR+G --tree pars{10} --bs-trees 100 --threads 8 --redo");
	system("mv $outdir/$orthogroup/$orthogroup.trim.clean.cds.fas.raxml.support $outdir/$orthogroup/$orthogroup.clean.tree");
	system("rm $outdir/$orthogroup/*raxml*");
}

sub get_sp_name {
	my $gene = shift;
	if ($gene =~ /^(\S+?)\_/) {return $1}
}

sub get_tips {
	my $node = shift;
	my @tips = $node->get_all_Descendents;
	my @output = ();
	for (@tips) {if ($_->id) {push @output,$_->id}}
	return \@output;
}
