use strict;

my $gene = $ARGV[0];

my %clade = ("Amb"=>"1","Ari"=>"1","Ara"=>"1","Rph"=>"2","Bya"=>"2","Ery"=>"2","Mal"=>"3","Hac"=>"4","Mco"=>"4","Tma"=>"4","Scu"=>"4","Tax"=>"4","Cha"=>"5","The"=>"6","Pyr"=>"6","Scl"=>"6",
	"Ecu"=>"7","Ene"=>"7","San"=>"7","Osy"=>"7","Dae"=>"8","Den"=>"8","Pri"=>"8","Pto"=>"8","Kor"=>"9","Val"=>"9","Vov"=>"9","Ach"=>"9","Api"=>"9");

my @VA = ();
my @TS = ();
my @auto = ();
my @L = ();
my @other = ();

open (FAS,"<./RRR/$gene/$gene.final.fas") or die();
open (VA,">./RRR/$gene/VA.txt");
open (TS,">./RRR/$gene/TS.txt");
open (AUTO,">./RRR/$gene/auto.txt");
open (OTHER,">./RRR/$gene/other.txt");
open (LO,">./RRR/$gene/L.txt");

while (<FAS>) {
	if (/>(\S+)/) {
		my $string = $1;
		my $sp = get_sp_name($string);
		if ($clade{$sp} <= 2) {push @auto,$string}
		elsif ($clade{$sp} >= 8) {push @VA,$string}
		elsif ($sp eq "Scu" or $sp eq "Tax") {push @TS,$string}
		else {push @other,$string}
		if ($clade{$sp} eq 4) {push @L,$string}
	}
}

print VA join("\n",@VA),"\n";
print TS join("\n",@TS),"\n";
print AUTO join("\n",@auto),"\n";
print OTHER join("\n",@other),"\n";
print LO join("\n",@L),"\n";

system("cat ./RRR/$gene/VA.txt ./RRR/$gene/TS.txt > ./RRR/$gene/VATS.txt");
system("cat ./RRR/$gene/VA.txt ./RRR/$gene/L.txt > ./RRR/$gene/VAL.txt");

#system("sed -r \'s/_p/.p/gi\' ./RRR/$gene/$gene.final.filtered.fas > ./RRR/$gene/$gene.final.filtered.fas2");

system("hyphy label-tree --tree ./RRR/$gene/$gene.root.rate.tree --list ./RRR/$gene/VA.txt --label VA --output ./RRR/$gene/$gene.root.rate.VA.tree");
system("hyphy label-tree --tree ./RRR/$gene/$gene.root.rate.tree --list ./RRR/$gene/TS.txt --label TS --output ./RRR/$gene/$gene.root.rate.TS.tree");
system("hyphy label-tree --tree ./RRR/$gene/$gene.root.rate.tree --list ./RRR/$gene/L.txt --label L --output ./RRR/$gene/$gene.root.rate.L.tree");
system("hyphy label-tree --tree ./RRR/$gene/$gene.root.rate.tree --list ./RRR/$gene/VATS.txt --label VATS --output ./RRR/$gene/$gene.root.rate.VATS.tree");
system("hyphy label-tree --tree ./RRR/$gene/$gene.root.rate.tree --list ./RRR/$gene/VAL.txt --label VAL --output ./RRR/$gene/$gene.root.rate.VAL.tree");

system("nohup hyphy relax --alignment ./RRR/$gene/$gene.final.fas --tree ./RRR/$gene/$gene.root.rate.VA.tree --test VA --srv Yes --output ./RRR/$gene/$gene.RELAX.VA.json &");
system("nohup hyphy relax --alignment ./RRR/$gene/$gene.final.fas --tree ./RRR/$gene/$gene.root.rate.VATS.tree --test VATS --srv Yes --output ./RRR/$gene/$gene.RELAX.VATS.json &");
system("nohup hyphy relax --alignment ./RRR/$gene/$gene.final.fas --tree ./RRR/$gene/$gene.root.rate.VAL.tree --test VAL --srv Yes --output ./RRR/$gene/$gene.RELAX.VAL.json &");
system("nohup hyphy relax --alignment ./RRR/$gene/$gene.final.fas --tree ./RRR/$gene/$gene.root.rate.L.tree --test L --srv Yes --output ./RRR/$gene/$gene.RELAX.L.json &");
system("hyphy relax --alignment ./RRR/$gene/$gene.final.fas --tree ./RRR/$gene/$gene.root.rate.TS.tree --test TS --srv Yes --output ./RRR/$gene/$gene.RELAX.TS.json");
#system("hyphy label-tree --tree ./RRR/$gene/$gene.root.rate.tree --list ./RRR/$gene/VA.txt --label VA --output ./RRR/$gene/tmp.1.tree");
#system("hyphy label-tree --tree ./RRR/$gene/tmp.1.tree --list ./RRR/$gene/TS.txt --label TS --output ./RRR/$gene/tmp.2.tree");
#system("hyphy label-tree --tree ./RRR/$gene/tmp.2.tree --list ./RRR/$gene/auto.txt --label auto --output ./RRR/$gene/tmp.3.tree");
#system("hyphy label-tree --tree ./RRR/$gene/tmp.3.tree --list ./RRR/$gene/other.txt --label other --output ./RRR/$gene/tmp.4.tree");
#system("nohup hyphy relax --alignment ./RRR/$gene/$gene.final.fas --tree ./RRR/$gene/tmp.4.tree --reference-group auto --srv Yes --output MSH1.RELAX.auto.other.json --mode 'Group mode' &");
#system("nohup hyphy relax --alignment ./RRR/$gene/$gene.final.fas --tree ./RRR/$gene/tmp.3.tree --reference-group auto --srv Yes --output MSH1.RELAX.auto.json --mode 'Group mode' &");

sub get_sp_name {
	my $gene = shift;
	if ($gene =~ /^(\S+?)\_/) {return $1}
}