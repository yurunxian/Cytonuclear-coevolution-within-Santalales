use strict;

my $gene = $ARGV[0];

open (ABSREL,"<./OXPHOS/$gene/$gene.ABSREL.json") or die($gene);
open (RELAX,"<./OXPHOS/$gene/$gene.RELAX.VA.json") or die($gene);
open (CODEML,"<./OXPHOS/$gene/$gene.codeml") or die();
open (TREE,">./OXPHOS/$gene/$gene.root.rate.tree") or die();
open (TREEO,">./OXPHOS/$gene/$gene.codeml.omega.tree") or die();

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
	else {
		$node_gene{$i} = $i;
	}
}

my @codeml = <CODEML>;
for (my $i = 0; $i < @codeml; $i++) {
	$codeml[$i] =~ s/\./_/gi;
	if ($codeml[$i] =~ /^tree length =/) {$rate_tree = $codeml[$i+4]}
	if ($codeml[$i] =~ /^w ratios as /) {
		my $string = $codeml[$i+1];
		$string =~ s/ #/:/g;
		$omega_tree = $string;
	}
}

print TREE $rate_tree;
print TREEO $omega_tree;

#system("perl get_branch_len.pl -g ./OXPHOS/$gene/$gene.root.rate.tree -o ./OXPHOS/$gene/$gene.root.rate.txt -i Santalales_tree.tree -n $gene -b Santalales_tree.txt");
#system("perl get_branch_len.pl -g ./OXPHOS/$gene/$gene.codeml.omega.tree -o ./OXPHOS/$gene/$gene.codeml.omega.txt -i Santalales_tree.tree -n $gene -b Santalales_tree.txt");

open (RATE,"<./OXPHOS/$gene/$gene.root.rate.txt") or die();
while (<RATE>) {
	chomp;
	my @tmp = split /\t/,$_;
	if (/Branch_length/) {next}
	else {
		my @gene = ();
		my $string = "";
		if ($tmp[3] =~ /,/) {@gene = split /, /,$tmp[3]; $string = join(", ",sort @gene)}
		else {$string = $tmp[3]}
		$string =~ s/\./_/g;
		$species_node{$string} = $tmp[1];
		$branch_length{$string} = $tmp[2];
	}
}

open (OMEGA,"<./OXPHOS/$gene/$gene.codeml.omega.txt") or die();
while (<OMEGA>) {
	chomp;
	my @tmp = split /\t/,$_;
	if (/Branch_length/) {next}
	else {
		my @gene = ();
		my $string = "";
		if ($tmp[3] =~ /,/) {@gene = split /, /,$tmp[3]; $string = join(", ",sort @gene)}
		else {$string = $tmp[3]}
		$string =~ s/\./_/g;
		$codeml_omega{$string} = $tmp[2];
	}
}

foreach my $i (sort keys %node_gene) {
	if (exists $species_node{$i} and $branch_length{$i} > 0.001) {
		if (!$codeml_omega{$i}) {$codeml_omega{$i} = "<0.001"}
		if ($codeml_omega{$i} > 2) {$codeml_omega{$i} = ">2"}
#		if (!$relax_k{$node_gene{$i}}) {$relax_k{$node_gene{$i}} = "N.A."}
		my $node = 0;
		if ($species_node{$i} =~ /(\d+)/) {$node = $1}
		if ($node > 4 and $node < 22) {$node += 1}
		elsif ($node >= 22 and $node < 25) {$node += 2}
		elsif ($node >= 25 and $node < 46) {$node += 3}
		elsif ($node >= 46) {$node += 3}
		print $gene,"\t",$node,"\t",$branch_length{$i},"\t",$codeml_omega{$i},"\t";
#		print $relax_k{$node_gene{$i}},"\t";
		print $p_value{$node_gene{$i}},"\t",$i,"\n";
	}
#	print $MG94_omega{$node_gene{$i}},"\t",$adaptive_dN{$node_gene{$i}},"\t",$adaptive_dS{$node_gene{$i}},"\t",$p_value{$node_gene{$i}},$i,"\n";
#	print $node_gene{$i},"\t",$i,"\n";
}