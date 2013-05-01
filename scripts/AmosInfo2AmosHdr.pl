$ocontig = "";
while (<>) {
    ($contig, $clen, $read, $roffset, @rest) = split;
    if ($ocontig and $ocontig ne $contig) { 
	print "##$ocontig ", scalar(@reads), " $oclen bases, 0 checksum.\n"; 
	print(join("\n", @reads), "\n"); 
	@reads = ();
    } 
    $ocontig = $contig; 
    $oclen = $clen;
    push @reads, "#$read($roffset) ".join(" ", @rest);
}
print "##$ocontig ", scalar(@reads), " $oclen bases, 0 checksum.\n";
print(join("\n", @reads), "\n");
