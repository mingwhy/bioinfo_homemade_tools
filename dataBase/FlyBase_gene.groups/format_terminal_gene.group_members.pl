while(<>){ #terminal_gene.group_members.txt
	s/\s+$//;
	@a=(split/\t/);
	#print $a[2],$a[3],"\n";
	$gene=pop @a;
	@genes=(split/\s+/,$gene);
	for(@genes){
		print "$a[0]\t$a[1]\t$a[2]\t$_\n";
	}
}

