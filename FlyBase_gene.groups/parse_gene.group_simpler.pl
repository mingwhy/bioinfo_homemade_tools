#$perl parse_gene.group.pl gene_group_data_fb_2021_02.tsv >terminal_gene.group_members.txt
# perl parse_gene.group.pl gene_group_data_fb_2021_02.tsv|less

## i) the member genes are only associated with the terminal subgroups
## ii) the immediate parent of any subgroup is identified in the â<U+0080><U+0098>Parent_FB_group_id' and 'Parent_FB_group_symb
while(<>){
    next if(/^#/);
    s/\s+$||^\s+//g;
    next if($_ eq '');
    @all=(split/\t/);
    #print $all[2],"\n";
    $symbol=$all[1];
    $name=$all[2];
    $anno{$all[0]}=$symbol."\t".$name;
    if($symbol eq '' || $name eq ''){
        #die("check here $all[0] $symboal $name @all\n")
    }
    $a=$_;
    @a=($a=~m/(FBgg\d+)/g);
    @b=($a=~m/(FBgn\d+)/g);   
    if(@b>1){
        die('more than one gene per row!\n')
    }

    #print "@a @b\n";
    if(@a==1){
        #a top gene group layer
    }elsif(@a==2 && @b==1){
        # yes a terminal group
        # has a child->parent relationship
        $link->{$a[0]}->{$a[1]}=1; #child -> parent
        $member->{$a[0]}->{$b[0]}=1; #terminal group gene member
        print "$a[0]\t$symbol\t$name\t$b[0]\n";
    }elsif(@b==1){
        # yes a terminal group
        # no parent term, but contain genes
        $link->{$a[0]}=1;
        $member->{$a[0]}->{$b[0]}=1;
        print "$a[0]\t$symbol\t$name\t$b[0]\n";
    }else{
        #not a terminal group
        $link->{$a[0]}->{$a[1]}=1; #not a terminal group, only parent-parent relationship
    }
}


