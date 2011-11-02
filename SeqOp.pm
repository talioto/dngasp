package SeqOp;
use strict;
use File::Basename;
use Bio::Seq;
use Bio::DB::Fasta;
############################################################
# 
# Author: T. Alioto
# Date: 060717
############################################################

sub chr_subseq{

    my ( $file, $start, $end, $str, $len) = @_;
    my $end_diff = 0;
    my $start_diff =0;
    my $length;
    #$len = 1;
    if (defined $len and $len > 0) {
	$length =  fasta_length($file);
	if ($end > $length) {
	    $end_diff = $end - $length;
	    $end = $length;
	}
	if ($start < 1) {
	    $start_diff = 1 - $start;
	    $start = 1;
	}
    }
    my $command = "/home/ug/talioto/bin/chr_subseq $file $start $end |";
    open( SEQ, $command ) || die("Error running command $command");
    my $seq = <SEQ>;

    if ($seq) {
	chomp $seq;
	close( SEQ );
	if ($str eq '-' || $str eq '-1' ) {
	    $seq =~tr/ACGTacgt/TGCAtgca/;
	    return "N" x $end_diff . reverse($seq) . "N" x $start_diff;
	} else {
	    return "N" x $start_diff .($seq). "N" x $end_diff;
	}
    } else {
	print STDERR "no seq found\n";
	return "";
    }
}
sub seqobj_subseq{

    my ( $sobj, $start, $end, $str) = @_;

    my $seq = $sobj->subseq($start,$end);

    if ($seq) {
	if ($str eq '-' || $str eq '-1' ) {
	    $seq =~tr/ACGTacgt/TGCAtgca/;
	    return reverse($seq);
	} else {
	    return ($seq);
	}
    } else {
	print STDERR "no seq found\n";
	return "";
    }
}
sub fasta_fetch{

    my ( $multifasta, $seqid, $indexdir, $outdir) = @_;
    my ($name,$path,$suffix) = fileparse($multifasta,qr{\.fa \.fasta \.fsa \.mfa});
    if (! -e "$indexdir/$name.index"){`fastaindex -f $multifasta -i $indexdir/$name.index`;}
    my $fasta_string = `fastafetch -i $indexdir/$name.index -f $multifasta -q $seqid`;
    if ($outdir){
	open OUT, ">$outdir/$seqid.fa";
	print OUT $fasta_string;
	close OUT;
    }else{
	return $fasta_string;
    }
}

sub fasta_length{
    my $fasta_file = shift;
    open FALEN, "fastalength -f $fasta_file |";
    my $lenin = <FALEN>;
    close FALEN;
    $lenin =~ m/(\d+)\s/;
    my $len = $1;
    return $len
}
sub get_seq_BioDBFasta{
    my $d = shift;
    #print STDERR "$d\n";
    my $c = shift;
    my $s = shift;
    my $e = shift;
    my $strand = shift;
    my $tmpend;
    my $rc = 0;
    if (($strand eq "-")||($strand eq "-1")){
	$tmpend = $s;
	$s = $e;
	$e = $tmpend;
	$rc = 1 if $s==$e;
    }
    my $seq =$d->seq($c,$s,$e);
    if ($rc){$seq =~ tr/ACGTacgt/TGCAtgca/;}
    return $seq;
}
sub get_genome{
    my $d = shift;
    my $c = shift;
    my $s = shift;
    my $e = shift;
    my $strand = shift;
    my $tmpend;
    if (($strand eq "-")||($strand < 0)){
	$tmpend = $s;
	$s = $e;
	$e = $tmpend;
    }
    my $seq = '';
    #print STDERR "get-genome -d $d '$c:$s-$e' \n";
    open GG, "get-genome -d $d '$c:$s-$e' |";
    while (<GG>){
	next if m/^>/;
	chomp;
	$seq .= $_;
    }
    close GG;
    return $seq;
}

1;
