use strict;

my $gene = $ARGV[0];

open (ABSREL,"<./paml/$gene/$gene.ABSREL.json") or die();
open (RELAX,"<./paml/$gene/$gene.RELAX.json") or die();
open (CODEML,"<./paml/$gene/$gene.codeml") or die();
open (TREE,">./paml/$gene/$gene.codeml.rate.tree") or die();
open (TREEO,">./paml/$gene/$gene.codeml.omega.tree") or die();

my $head = "";
my $absrel_tree = "";
my $rate_tree = "";
my $omega_tree = "";
my %relax_k = ();
my %p_value = ();
my %MG94_omega = ();
my %adaptive_dS = ();
my %adaptive_dN = ();
my %node_gene = ();
my %species_node = ();
my %branch_length = ();
my %codeml_omega = ();

while (<ABSREL>) {
	$_ =~ s/[\r\n]+//g;
	if (/"(\S+)":{/) {$head = $1}
	elsif (/"Baseline MG94xREV omega ratio":(\S+),/) {$MG94_omega{$head} = $1}
	elsif (/"Corrected P-value":(\S+),/) {$p_value{$head} = $1}
	elsif (/-synonymous subs\/site\)":(\S+),/) {$adaptive_dN{$head} = $1}
	elsif (/\(synonymous subs\/site\)":(\S+),/) {$adaptive_dS{$head} = $1}
	if (/"0":"(\S+)"/) {$absrel_tree = $1}
}

while (<RELAX>) {
	$_ =~ s/[\r\n]+//g;
	if (/\"(\S+)\":{/) {$head = $1}
	elsif (/k \(general descriptive\)":(\S+)/) {$relax_k{$head} = $1}
}

$absrel_tree =~ s/:\d+\.*\d+//gi;
$absrel_tree =~ s/:\de-\d+//gi;

foreach my $i (sort keys %p_value) {
	if ($i =~ /Node/) {
		my @tmp = split /$i[),]/,$absrel_tree;
		my $left = 0;
		my $right = 0;
		my $length = 0;
		until ($left eq $right and $left > 0 and $right > 0) {
			my $letter = substr($tmp[0],length($tmp[0])-$length,1);
			if ($letter eq "(") {$left++}
			if ($letter eq ")") {$right++}
			$length++;
		}
		my $string = substr($tmp[0],-$length,$length);
		$string =~ s/[(),]/ /g;
		$string =~ s/Node\d+/ /g;
		$string =~ s/^\s+//g;
		$string =~ s/\s+$//g;
		my @gene = split /\s+/,$string;
		$node_gene{join(", ",sort @gene)} = $i;
#		print "$i\t",$MG94_omega{$i},"\t",$adaptive_dN{$i}/$adaptive_dS{$i},"\t",$p_value{$i},"\n";
	}
}

my @codeml = <CODEML>;
for (my $i = 0; $i < @codeml; $i++) {
	if ($codeml[$i] =~ /^tree length =/) {$rate_tree = $codeml[$i+4]}
	if ($codeml[$i] =~ /^w ratios as /) {
		my $string = $codeml[$i+1];
		$string =~ s/ #/:/g;
		$omega_tree = $string;
	}
}

$rate_tree =~ s/([a-zA-Z]+)/$1_$gene/g;
$omega_tree =~ s/([a-zA-Z]+)/$1_$gene/g;

print TREE $rate_tree;
print TREEO $omega_tree;

system("perl get_branch_len.pl -g ./paml/$gene/$gene.codeml.rate.tree -o ./paml/$gene/$gene.codeml.rate.txt -i Santalales_tree.tree -n $gene -b Santalales_tree.txt");
system("perl get_branch_len.pl -g ./paml/$gene/$gene.codeml.omega.tree -o ./paml/$gene/$gene.codeml.omega.txt -i Santalales_tree.tree -n $gene -b Santalales_tree.txt");

open (RATE,"<./paml/$gene/$gene.codeml.rate.txt");
while (<RATE>) {
	chomp;
	my @tmp = split /\t/,$_;
	if (/Branch_length/) {next}
	else {
		$tmp[3] =~ s/_$gene//g;
		my @gene = split /, /,$tmp[3];
		my $string = join(", ",sort @gene);
		$species_node{$string} = $tmp[1];
		$branch_length{$string} = $tmp[2];
	}
}

open (OMEGA,"<./paml/$gene/$gene.codeml.omega.txt");
while (<OMEGA>) {
	chomp;
	my @tmp = split /\t/,$_;
	if (/Branch_length/) {next}
	else {
		$tmp[3] =~ s/_$gene//g;
		my @gene = split /, /,$tmp[3];
		my $string = join(", ",sort @gene);
		$codeml_omega{$string} = $tmp[2];
	}
}

foreach my $i (sort keys %node_gene) {
	if (exists $species_node{$i} and $branch_length{$i} > 0.001) {
		if (!$codeml_omega{$i}) {$codeml_omega{$i} = "<0.001"}
		if ($codeml_omega{$i} > 2) {$codeml_omega{$i} = ">2"}
#		if (!$relax_k{$node_gene{$i}}) {$relax_k{$node_gene{$i}} = "N.A."}
		print $gene,"\t",$species_node{$i},"\t",$branch_length{$i},"\t",$codeml_omega{$i},"\t";
		print $relax_k{$node_gene{$i}},"\t";
		print $p_value{$node_gene{$i}},"\t",$i,"\n";
	}
	
#	print $MG94_omega{$node_gene{$i}},"\t",$adaptive_dN{$node_gene{$i}},"\t",$adaptive_dS{$node_gene{$i}},"\t",$p_value{$node_gene{$i}},$i,"\n";
#	print $node_gene{$i},"\t",$i,"\n";
}