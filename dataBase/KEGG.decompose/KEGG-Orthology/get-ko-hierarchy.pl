#perl get-ko-hierarchy.pl kegg-T00030.txt ko00001.keg >dme.keg.txt
$f=shift;
open F,"<$f";
while(<F>){
	s/\s+$//;
	$a=(split/\s+/)[0];
	$b=substr($a,3);
	push @keep, $b;
	$dme{$b}=$_
#print "$a $b\n";
}
close F;

# parse ko00001.keg tree
while(<>){
	next if(/^D/);
	s/\s+$//;
	if(/^A/){
		$lay1=$_;
		next
	}
	if(/^B/){
		@a=(split/\s+/);
		next if(scalar(@a)==1);
		$lay2=$_;
	}
	if(/^C/){
		@a=(split/\s+/);
		shift @a;
		$id=shift @a;
		$line=join(' ',@a);
		$h->{$id}->{'lay2'}=$lay2;
		$h->{$id}->{'lay1'}=$lay1;
	}
}
@id=keys %{$h};
print '#id ',scalar(@id),"\n";
#for(@id){
#	print "$_, $h->{$_}->{'lay2'}, $h->{$_}->{'lay1'}\n";
#}

for(@keep){
	$b=$_;
	#print "dme$_: $h->{$_}->{'lay2'}: $h->{$_}->{'lay1'}\n";
	print "$dme{$b}: $h->{$_}->{'lay2'}: $h->{$_}->{'lay1'}\n";
}
