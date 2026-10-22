use strict;
use warnings;
use List::MoreUtils qw(uniq);
use Array::Utils qw(:all);

open (OG,"<Orthogroups.4sp.tsv") or die();
open (OXPHOS,"<OXPHOS.txt") or die();
open (SPLIC,"<PPR.txt") or die();
open (CP,"<ribosome_cp.txt") or die();

my @genome = ("Amb","Ari","Ara","Ery","Mal","Scu","Tax","San","Osy");
my @species = ();
my @OG = ();
my %OG = ();
my @PPR = ();
my @CP_gene = ();
my %MT = ();
my %PT = ();
my @OXPHOS = ();

while (<OG>) {
	chomp;
	my @tmp = split /\t/,$_;
	if (/Orthogroup/) {
#	if (/Homogroup/) {
		for (my $i = 1; $i < @tmp; $i++) {push @species,$tmp[$i]}
	}
	else {
		push @OG,$tmp[0];
		for (my $i = 1; $i < @tmp; $i++) {
			foreach my $gene (split /, /,$tmp[$i]) {
				if ($gene =~ /(\S+_c\d+_g\d+)_/) {push @{$OG{$tmp[0]}{$species[$i-1]}},$1}
				else {push @{$OG{$tmp[0]}{$species[$i-1]}},$gene}
			}
		}
		foreach my $sp (sort keys %{$OG{$tmp[0]}}) {
			@{$OG{$tmp[0]}{$sp}} = uniq(@{$OG{$tmp[0]}{$sp}});
		}
	}
}

while (<SPLIC>) {
	$_ =~ s/[\r\n]+//g;
	$_ =~ s/, /\t/g;
	my @tmp = split /\t/,$_;
	foreach my $gene (@tmp) {
		if ($gene =~ /(\S+)/) {push @PPR,$1}
	}
}

while (<CP>) {
	$_ =~ s/[\r\n]+//g;
	$_ =~ s/, /\t/g;
	my @tmp = split /\t/,$_;
	foreach my $gene (@tmp) {
		if ($gene =~ /(\S+)/) {push @CP_gene,$1}
	}
}
@PPR = uniq(@PPR);
@CP_gene = uniq(@CP_gene);

while (<OXPHOS>) {
	$_ =~ s/[\r\n]+//g;
	$_ =~ s/, /\t/g;
	my @tmp = split /\t/,$_;
	foreach my $gene (@tmp) {
		if ($gene =~ /(\S+)/) {push @OXPHOS,$1}
	}
}

foreach my $sp (@species) {
	if ($sp eq "Ara") {
		open (SUBA,"<./Target/SUBA5.txt") or die();
		while (<SUBA>) {
			chomp;
			my @tmp = split /\t/,$_;
			$tmp[0] =~ s/\.\d+//g;
			$tmp[0] = "Ara_".$tmp[0];
			if ($tmp[2] ne "" or $tmp[3] ne "") {
				my $string = $tmp[2].$tmp[3];
				if ($tmp[1] =~ /plastid/ and $string =~ /plastid/) {push @{$PT{$sp}},$tmp[0]}
				if ($tmp[1] =~ /mito/ and $string =~ /mito/) {push @{$MT{$sp}},$tmp[0]}
			}
		}
	}
	else {
		open (CP,"<./Target/$sp.cp.txt") or die();
		open (MT,"<./Target/$sp.mt.txt") or die();
		while (<CP>) {if (/(\S+_c\d+_g\d+)_/) {push @{$PT{$sp}},$1} elsif (/(\S+)/) {push @{$PT{$sp}},$1}}
		while (<MT>) {if (/(\S+_c\d+_g\d+)_/) {push @{$MT{$sp}},$1} elsif (/(\S+)/) {push @{$MT{$sp}},$1}}
	}
}




foreach my $OG (@OG) {
	my ($PT,$MT,$PPR,$OXPHOS,$cp) = (0,0,0,0,0);
	my ($PT_genome,$MT_genome,$PPR_genome,$OXPHOS_genome,$cp_genome) = (0,0,0,0,0);
	print $OG;
	foreach my $sp (@species) {
		if (exists $OG{$OG}{$sp}) {
			my @PPR2 = intersect(@{$OG{$OG}{$sp}},@PPR);
			if (scalar @PPR2 > 0) {$PPR++}
			my @cp = intersect(@{$OG{$OG}{$sp}},@CP_gene);
			if (scalar @cp > 0) {$cp++}
			my @PT = intersect(@{$OG{$OG}{$sp}},@{$PT{$sp}});
			if (scalar @PT > 0) {$PT++}
			my @MT = intersect(@{$OG{$OG}{$sp}},@{$MT{$sp}});
			if (scalar @MT > 0) {$MT++}
			my @OXPHOS = intersect(@{$OG{$OG}{$sp}},@OXPHOS);
			if (scalar @OXPHOS > 0) {$OXPHOS++}
#			print "\t",join(", ",@{$OG{$OG}{$sp}});
		}
#		else {print "\t"}
	}
	foreach my $sp (@genome) {
		if (exists $OG{$OG}{$sp}) {
			my @PPR2 = intersect(@{$OG{$OG}{$sp}},@PPR);
			if (scalar @PPR2 > 0) {$PPR_genome++}
			my @cp = intersect(@{$OG{$OG}{$sp}},@CP_gene);
			if (scalar @cp > 0) {$cp_genome++}
			my @PT = intersect(@{$OG{$OG}{$sp}},@{$PT{$sp}});
			if (scalar @PT > 0) {$PT_genome++}
			my @MT = intersect(@{$OG{$OG}{$sp}},@{$MT{$sp}});
			if (scalar @MT > 0) {$MT_genome++}
			my @OXPHOS = intersect(@{$OG{$OG}{$sp}},@OXPHOS);
			if (scalar @OXPHOS > 0) {$OXPHOS_genome++}
#			print "\t",join(", ",@{$OG{$OG}{$sp}});
		}
#		else {print "\t"}
	}

	my $sp_num = scalar keys %{$OG{$OG}};
	if ($PT /$sp_num > 0.7 or $PT_genome >= 5) {print "\tPT"} else {print "\t"}
	if ($MT /$sp_num > 0.7 or $MT_genome >= 5) {print "\tMT"} else {print "\t"}
	if ($PPR > 0) {print "\tPPR"} else {print "\t"}
	if ($cp > 0) {print "\tcp"} else {print "\t"}
	if ($OXPHOS > 0) {print "\tOXPHOS"} else {print "\t"}
	print "\n";
}