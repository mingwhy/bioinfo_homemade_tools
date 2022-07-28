#</td><td>
while(<>){
    next unless(/FBgn/);
    s/\s+$//;
    $_=~s/\<span.+?\>//g;
    $_=~s/\<\/span\>//g;
    $_=~s/\<a href.+?\>//g;
    $_=~s/\<\/a\>//g;
    $_=~s/\<u\>||\<\/u\>//g;
    @a=($_=~m/\<td\>(.+?)\<\/td\>/g);
    ##@a=($_=~m/\>(.+?)\</g);
    #print scalar(@a),"\n"; #all 16 elements
    #print "@a\n";
    $gene=$a[2];
    $symbol=$a[3];
    $go=pop @a;
    $XC=$a[11];
    print "$gene\t$symbol\t$XC\t$go\n";
}

