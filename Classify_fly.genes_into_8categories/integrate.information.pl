open F,"<ribosomal_proteins.txt";
while(<F>){
    s/\s+$//g;
    $h{$_}='ribosomal genes';
    print "$_\tribosomal genes\n";
}
close F;

open F,"transcription_factors.txt";
while(<F>){
    s/\s+$//g;
    if($h{$_} eq ''){
        $h{$_}="transcription factors";
        print "$_\ttranscription factors\n";
    }
}
close F;

open F,"<translation_factors.txt";
while(<F>){
    s/\s+$//;
    $h{$_}='translation factors';
}
close F;
open F,"<tRNA_genes.txt";
while(<F>){
    s/\s+$//;
    $h{$_}='tRNA genes';
}
close F;

open F,"<GO_RNA_binding_proteins.txt";
while(<F>){
    s/\s+$//;
    $gene=(split/\s+/)[0];
    if($h{$gene} eq ''){
        $h{$gene}='RNA-binding proteins';
        print "$gene\tRNA-binding proteins\n";
    }
}
close F;

open F,"<GO_non-coding_RNA.genes.txt";
while(<F>){
    s/\s+$//;
    $gene=(split/\s+/)[0];
    if($h{$gene} eq ''){
        $h{$gene}='non-coding RNA genes';
        print "$gene\tnon-coding RNA genes\n";
    }
}
close F;

open F,"<CAM_cell.adhesion.molecules.txt";
while(<F>){
    s/\s+$//;
    if($h{$_} eq ''){
        $h{$_}='cell adhesion molecules';
        print "$_\tcell adhesion molecules\n";
    }
}
close F;

open F,"<transmembrane_receptors.txt";
while(<F>){
    s/\s+$//;
    if($h{$_} eq ''){
        $h{$_}='receptor and ligands';
        print "$_\treceptor and ligands\n";
    }
}
close F;
open F,"<receptor_ligands.txt";
while(<F>){
    s/\s+$//;
    if($h{$_} eq ''){
        $h{$_}='receptor and ligands';
        print "$_\treceptor and ligands\n";
    }
}
close F;


open F,"<ion_channels.txt";
while(<F>){
    s/\s+$//;
    if($h{$_} eq ''){
        $h{$_}='ion channels';
        print "$_\tion channels\n"
    }
}
close F;

open F,"<GO_synaptic.genes.txt";
while(<F>){
    s/\s+$//;
    $gene=(split/\s+/)[0];
    if($h{$gene} eq ''){
        $h{$gene}='synaptic genes';
        print "$gene\tsynaptic genes\n";
    }
}
close F;



