use strict;

my $gene = $ARGV[0];
my $type = $ARGV[1];

open(IN,"<./ribosome/$gene/$gene.RELAX.$type.json");

my $p_value = "";
my $k = "";
my $model = "";
my $branch = "";
my %omega = ();
my %ratio = ();
my %LRT = ();
my %AICc = ();

while (<IN>) {
	if (/"p-value":(\S+),/) {$p_value = $1}
	if (/"relaxation or intensification parameter":(\S+)/) {$k = $1}
	if (/"RELAX alternative":{/) {$model = "alternative"}
	if (/"RELAX null":{/) {$model = "null"}

	if (/"Reference":{/) {$branch = "ref"}
	if (/"Test":{/) {$branch = "test"}
	if (/"Synonymous site-to-site rates":{/) {$branch = "Synonymous"}

	if (/"AIC-c":(\S+),/ and ! exists $AICc{$model}) {$AICc{$model} = $1}
	if (/"Log Likelihood":(\S+),/ and ! exists $LRT{$model}) {$LRT{$model} = $1}

	if (/"omega":(\S+),/) {push @{$omega{$model}{$branch}},$1}
	if (/"proportion":(\S+)/) {push @{$ratio{$model}{$branch}},$1}
}

my $mean_omega_null = $omega{"null"}{"ref"}[0]*$ratio{"null"}{"ref"}[0] + $omega{"null"}{"ref"}[1]*$ratio{"null"}{"ref"}[1] + $omega{"null"}{"ref"}[2]*$ratio{"null"}{"ref"}[2];
my $mean_omega_test = $omega{"alternative"}{"test"}[0]*$ratio{"alternative"}{"test"}[0] + $omega{"alternative"}{"test"}[1]*$ratio{"alternative"}{"test"}[1] + $omega{"alternative"}{"test"}[2]*$ratio{"alternative"}{"test"}[2];
my $mean_omega_ref = $omega{"alternative"}{"ref"}[0]*$ratio{"alternative"}{"ref"}[0] + $omega{"alternative"}{"ref"}[1]*$ratio{"alternative"}{"ref"}[1] + $omega{"alternative"}{"ref"}[2]*$ratio{"alternative"}{"ref"}[2];

print "$type\t$gene\t$k\t$p_value\t",$LRT{"null"},"\t",$AICc{"null"},"\t$mean_omega_null\t",$LRT{"alternative"},"\t",$AICc{"alternative"},"\t$mean_omega_ref\t$mean_omega_test\n";