$f=shift; #terminal_gene.group_members.txt
open F,"<$f";
while(<F>){
	$gg=(split/\s+/)[0];
	$h{$gg}=1;
}
close F;

while(<>){ #gene_to_all_its_groups.txt
	($gene,$group)=(split/\s+/)[0,1];
	if($h{$group}==1){
		$h1{$gene}++
	}
}
@genes=sort { $h1{$a} <=> $h1{$b} } keys  %h1;
for(@genes){
	print "$_\t$h1{$_}\n";
}

	

