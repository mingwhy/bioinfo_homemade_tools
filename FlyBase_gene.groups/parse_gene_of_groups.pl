# perl parse_gene_of_groups.pl gene_group_data_fb_2021_02.tsv >gene_to_all_its_groups.txt

## i) the member genes are only associated with the terminal subgroups
## ii) the immediate parent of any subgroup is identified in the â<U+0080><U+0098>Parent_FB_group_id' and 'Parent_FB_group_symb

while(<>){
    next if(/^#/);
    s/\s+$||^\s+//g;
    next if($_ eq '');
    @all=(split/\t/);
    $child=$all[0];
    $symbol=$all[1];
    $name=$all[2];
     $anno{$all[0]}=$symbol.":".$name;
    
     if($symbol eq '' || $name eq ''){
        #die("check here $all[0] $symboal $name @all\n")
    }
    
    $a=$_;
    @a=($a=~m/(FBgg\d+)/g);
    @b=($a=~m/(FBgn\d+)/g);   
    if(@b>1){
        die('more than one gene per row!\n')
    }
    $h{$a[0]}=$anno;
    #print "@a @b\n";
    if(@a==1){
        #a top gene group layer or group directly attached to root
        $link->{$a[0]}=2; #root
    }elsif(@a==2 && @b==1){
        # yes a terminal group
        # has a child->parent relationship
        $link->{$a[0]}->{$a[1]}=1; #child -> parent
        $gene->{$b[0]}->{$a[0]}=1; #gene-> group
    }elsif(@b==1){
        # yes a terminal group
        # no parent term, but contain genes
        $link->{$a[0]}=1;
        $gene->{$b[0]}->{$a[0]}=1;
    }else{
        #not a terminal group
        $link->{$a[0]}->{$a[1]}=1; #not a terminal group, only parent-parent relationship
    }
}

@genes= keys %{$gene};
for(@genes){
    $name=$_;
    @groups=keys %{$gene->{$name}}; #gene to terminal groups
    #print "check $name\t@groups\n";
    for(@groups){
        $gg=$_;
        print "$name\t$gg\t$anno{$gg}\n";
        if($link->{$gg}==1 || $link->{$gg}==2){
            next
        }else{
            $here=$gg;
            @keys=keys %{$link->{$here}};
            while(@keys>0){
                $here=$keys[0];
                @keys=keys %{$link->{$here}};
                print "$name\t$here\t$anno{$here}\n";
                if($link->{$here}==1 || $link->{$here}==2){last}
            }
        }
    }
}


