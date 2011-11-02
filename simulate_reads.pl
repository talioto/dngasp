#!/usr/bin/env perl
use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use Math::Random qw(:all);
use Bio::DB::Fasta;
use File::Basename;
use lib "~/myperlmods/";
use SeqOp;
use Time::HiRes;  

my $fqfile = '';
my $N = 1;
my $p1file ='';	#'/home/devel/talioto/proj-dev-tal/GAG/assembly_evaluation/sample_qual/lynx621D2AAXX_p1.quals.txt';
my $p2file ='';	#'/home/devel/talioto/proj-dev-tal/GAG/assembly_evaluation/sample_qual/lynx621D2AAXX_p2.quals.txt';
#my $genome_mfa = '/scratch/devel/talioto/assembly/dnGASP/sim/sim_genome_final.fa';
#test sequence /project/devel/talioto/GAG/assembly_evaluation/simulate_reads/test.fa

my $genome_mfa = '/home/devel/talioto/data/sim_genome_final.fa';
my $xcov = 50;
my $insertsize_avg = 0;
my $mp_insertsize_avg = 5000;
my $mp_insert_sd = 50;
my $type = 'pe';		# other option is mp;
my $mp=0;
my $pe=1;
my $rl=100; #dnGASP 114
my $mprl=100; #dnGASP 36
my $minpereadlength = 50;
my $minreadlength = $mprl;
my $sdpct = 0.08;
my $printGC = 0;
my $bound = 0;
my $pehist = '/project/devel/talioto/GAG/assembly_evaluation/simulate_reads/insert_hist.txt';
my $minpelen = 150;
my $filetag = $$;
my $mp_contamination_rate = 0.05; #0.05;
my $qdir = '/project/devel/talioto/GAG/assembly_evaluation/sample_qual/lynx';
my $BP=0;
my $genome_size = 0;
my $loc_in_id = 0;
my $label_contam = 1;
my $bias = 1;
my $error = 1;
my $btrim = 1;
my $printObservedGC = 0;
my $pairfile = 0;
GetOptions(
	   'bp' =>\$BP,
	   'qdir:s'=>\$qdir,
	   'qp1:s'=>\$p1file,
	   'qp2:s'=>\$p2file,
	   'g:s'=>\$genome_mfa,
	   'x:s'=>\$xcov,
	   'i:s'=>\$insertsize_avg,
	   'mi:s'=>\$mp_insertsize_avg,
	   'sd:s'=>\$mp_insert_sd,
	   'mp'=>\$mp,
	   'pe'=>\$pe,
	   'rl|readlength:s'=>\$rl,
	   'mprl|mpreadlength:s'=>\$mprl,
	   'mpcr=i'=> \$mp_contamination_rate,
	   'gc'=>\$printGC,
	   'gcbias=i'=>\$bias,
	   'pehist:s'=>\$pehist,
	   'bound' =>\$bound,
	   't|tag:s'=>\$filetag,
	   'l|loc!' =>\$loc_in_id,
	   'k!' => \$label_contam,
	   'e=s'=> \$error,
	   'btrim!' => \$btrim,
	   'p'=>\$pairfile
	  );

####### PROB OF SEQ given GC% #########
my %gcprob=(
	    '0.00'=>0.567676367893438,
	    '0.02'=>0.531251483324315,
	    '0.04'=>0.678531907523418,
	    '0.06'=>0.75754322576544,
	    '0.08'=>0.777638961590955,
	    '0.10'=>0.867007003050834,
	    '0.12'=>0.878788514920345,
	    '0.14'=>0.982633893251847,
	    '0.16'=>1,
	    '0.18'=>0.986966679736818,
	    '0.20'=>0.952755336891541,
	    '0.22'=>0.921717398411219,
	    '0.24'=>0.90140050089748,
	    '0.26'=>0.868171349693093,
	    '0.28'=>0.85352242071229,
	    '0.30'=>0.83507105155154,
	    '0.32'=>0.811850888007334,
	    '0.34'=>0.80436291301415,
	    '0.36'=>0.794178198567866,
	    '0.38'=>0.77321879224528,
	    '0.40'=>0.72685326523463,
	    '0.42'=>0.674024803953745,
	    '0.44'=>0.621663146247214,
	    '0.46'=>0.564046144192179,
	    '0.48'=>0.510152506933564,
	    '0.50'=>0.454927475026177,
	    '0.52'=>0.403678136416014,
	    '0.54'=>0.365044704216372,
	    '0.56'=>0.326333223936909,
	    '0.58'=>0.300428213413583,
	    '0.60'=>0.272803155255271,
	    '0.62'=>0.245935147740802,
	    '0.64'=>0.221283960074104,
	    '0.66'=>0.191294651000942,
	    '0.68'=>0.152176693692474,
	    '0.70'=>0.119656870662224,
	    '0.72'=>0.0801352939439834,
	    '0.74'=>0.0486969167123914,
	    '0.76'=>0.0254712345407072,
	    '0.78'=>0.0125323699824664,
	    '0.80'=>0.00483856703647052,
	    '0.82'=>0.00194324175700659,
	    '0.84'=>0.000923268767291863,
	    '0.86'=>0.000506424480288979,
	    '0.88'=>0.000330958361204984,
	    '0.90'=>0.000378621335717881,
	    '0.92'=>0.000771047209259158,
	    '0.94'=>0.00234666430871269,
	    '0.96'=>0.00392899320529253,
	    '0.98'=>0.00392899320529253,
	    '1.00'=>0.00392899320529253,
	   );
if ($rl < $minpereadlength) {
  $rl = $minpereadlength; print STDERR "PE read length set to minimum of $minpereadlength\n";
}
if ($mp_contamination_rate < 0){$mp_contamination_rate = 0;}elsif($mp_contamination_rate > 1){$mp_contamination_rate = 1;} #mpcr>1 not implemented
my %observed_gc;
foreach my $gcbin (sort {$a<=>$b} keys %gcprob) {
  $observed_gc{$gcbin}=0;
}

my $uid = sprintf("%s_%s", $filetag,generate_random_password(3));
opendir (QDIR,$qdir) or die "Couldn't open $qdir: $!\n";
my @p1files = grep /pair1/, readdir QDIR;
close QDIR;
my ($p1file,$p2file) = getQualFile($qdir,\@p1files);

$type = 'mp' if $mp;
my $SCRATCH = $ENV{'TMPDIR'};
#`mkdir -p $SCRATCH/$uid`;
my ($base,$path,$ext) = fileparse($genome_mfa,qw(\.fa \.fasta \.mfa));
if (!-e "$SCRATCH/$base$ext") {
  `cp $genome_mfa $SCRATCH`;
}
my $genome_mfa_orig = $genome_mfa;
$genome_mfa = "$SCRATCH/$base$ext";

print STDERR "$genome_mfa\n";
if (!-e "$genome_mfa.index") {
  print STDERR "indexGenome.pl -f $genome_mfa\n";
  `indexGenome.pl -f $genome_mfa`;
}
die "$genome_mfa.index does not exist!\n" if !-e "$genome_mfa.index";
my $db = Bio::DB::Fasta->new($genome_mfa);

print STDERR "Reading PE distribution from $pehist...\n";
open(HIST, "<$pehist") or die "$!\n";
my @pelengths;
my $csum = 0;
while (<HIST>) {
  my @f = split " ",$_;
  my $len = $f[0];
  next if $len < 50;
  next if $len > 2000;
  my $count = $f[1];
  $csum += $count;
  my $p = $f[2];
  push (@pelengths,{'len'=>$len,'csum'=>$csum});
}
close HIST;
my $hist_insertsize_median = 400;

my $rval = int($csum/2);
foreach my $length (@pelengths) {
  if ($rval < $length->{'csum'}) {
    $hist_insertsize_median = $length->{'len'};
    last;
  }
}

print STDERR "median insert in $pehist: $hist_insertsize_median\n";
$insertsize_avg = $hist_insertsize_median if !$insertsize_avg;
my $run = "$type$insertsize_avg"."_$uid";
if ($type eq 'mp') {
  $run = "$type$mp_insertsize_avg"."_$uid";
}
my %seqs;
if ($BP) {
  print STDERR "Loading genome into memory...\n";
  my $seqin = Bio::SeqIO->new(-file=>$genome_mfa , '-format' => 'Fasta' );
  while ( my $seqobj = $seqin->next_seq()) {
    print STDERR "\t",$seqobj->id,"\n";
    $seqs{$seqobj->id}=$seqobj;
  }
}
print STDERR "Getting chr lengths\n";
# get chr lengths
my %chrlen;
my @seqids = sort $db->ids;
die if ! @seqids;
my $slcsum = 0;
my @seqlengths;
my $minpelen = (($type eq 'pe')?($hist_insertsize_median - int($hist_insertsize_median/4)):($mp_insertsize_avg - int($mp_insertsize_avg/4)));
my $min_contam_pelen = ($hist_insertsize_median - int($hist_insertsize_median/4));
print STDERR "Minimum pe/mp length: $minpelen\n";
open(SAM,">$run.reads.sam") or die "$!\n";
print SAM "\@HD\tVN:1.0\n\@PG\tID:SIMREADS\n";
foreach my $id (@seqids) {
  $chrlen{$id}=$db->length($id);
  print SAM '@SQ'."\tSN:$id\tLN:$chrlen{$id}\n";
  $genome_size += $chrlen{$id};
  print STDERR "$id\t$chrlen{$id}\n";
  next if $chrlen{$id} < (($type eq 'pe')?$hist_insertsize_median:$mp_insertsize_avg);
  my $offset = $slcsum;
  my $sampling_length = $chrlen{$id} - 2*($minpelen);
  $slcsum += $sampling_length;
  push (@seqlengths,{'sid'=>$id,'csum'=>$slcsum,'offset'=>$offset,'real_length'=>$chrlen{$id}});
}


my %normal = ('A'=> ['C','G','T'],
	      'C'=> ['A','G','T'],
	      'G'=> ['C','A','T'],
	      'T'=> ['C','G','A'],
	      'N'=> ['C','G','A','T']
	     );
my %blist = ('A'=> ['C','G','T','C','G','T','C','G','T','C','G','T','C','G','T','C','G','T','N'],
	     'C'=> ['A','G','T','A','G','T','A','G','T','A','G','T','A','G','T','A','G','T','N'],
	     'G'=> ['C','A','T','C','A','T','C','A','T','C','A','T','C','A','T','C','A','T','N'],
	     'T'=> ['C','G','A','C','G','A','C','G','A','C','G','A','C','G','A','C','G','A','N'],
	     'N'=> ['C','G','A','T','C','G','A','T','C','G','A','T','C','G','A','T','N']
	    );

#index probability according to ascii character
my %prob_from_ascii;
for (my $q = 2; $q<=40; $q++) {
  $prob_from_ascii{chr($q+64)}=exp(($q*log(10))/-10);
  #print STDERR chr($q+64) ."=>". exp(($q*log(10))/-10)."\n";
}
my %prob_from_qual;
for (my $q = 2; $q<=40; $q++) {
  $prob_from_qual{$q}=exp(($q*log(10))/-10);
  #print STDERR chr($q+64) ."=>". exp(($q*log(10))/-10)."\n";
}

#foreach my $k (keys %blist){print "$k\n";}
die "must provide -g genome.fa" if !-e $genome_mfa;

my ($base,$path,$ext) = fileparse($genome_mfa,qw(.fa .fasta));
if ($printGC) {
  open(GC,">$run.gc.txt") or die "$!\n";
}

print STDERR "opening $p1file and $p2file...\n";
open(my $p1fq, "<$p1file") or die "$!\n";
open(my $p2fq, "<$p2file") or die "$!\n";
my $p1fqout;
my $p2fqout;
my $pairfileout;
if ($pairfile){
  open($pairfileout, ">$run"."_pair.txt") or die "$!\n";
}else{
  open($p1fqout, ">$run"."_1.fastq") or die "$!\n";
  open($p2fqout, ">$run"."_2.fastq") or die "$!\n";
}
if ($type eq 'pe') {
  $N = int(($xcov*$genome_size)/($rl*2))+1;
} else {
  $N = int(($xcov*$genome_size)/($mprl*2))+1;
}
print STDERR "For $xcov"."x coverage you will need $N read pairs\n";
random_set_seed_from_phrase($genome_mfa);
my @peis_result = ();
my @mpis_result = ();
my @contam_peis_result = ();
print STDERR "Generating random pe insert sizes with average $insertsize_avg\n";
my $offset = $insertsize_avg - $hist_insertsize_median;
for (my $i = 0; $i < int($N*(1.1)); $i++) {
  #print "$i\n" if !($i%100000);
  my $rval = int(rand($csum));
  foreach my $length (@pelengths) {
    if ($rval < $length->{'csum'}) {
      push (@peis_result,($length->{'len'}+$offset));
      last;
    }
  }
}

if ($type eq 'mp') {
  # mp
  print STDERR "Generating random $type insert sizes with average $mp_insertsize_avg and SD ".int($sdpct*$mp_insertsize_avg)."\n";
    
  @mpis_result = (random_normal($N+int($N*(1.1)),$mp_insertsize_avg,($sdpct*$mp_insertsize_avg)));

  ### for contaminant PEs in MP prep
  my $ncontam = int($N * $mp_contamination_rate * 1.25);
  for (my $i = 0; $i < $ncontam; $i++) {
    #print "$i\n" if !($i%100000);
    my $rval = int(rand($csum));
    foreach my $length (@pelengths) {
      if ($rval < $length->{'csum'}) {
	push (@contam_peis_result,($length->{'len'}+$offset));
	last;
      }
    }
  }
}
print STDERR "Generating reads...\n";
my $contam_count = 0;
my $j = $N;
my $begintime = 0;
my $getseqtime = 0;
my $makefatime = 0;
my $makesamtime = 0;
my $num_passed_gc_filter = 0;
my $i = 0;
while ($i < $N) {
  my $t0 = Time::HiRes::time();
  my $t1;
  my $t2;
  my $t3;
  my $t4;
  my $t5;
  if ($i && (!($i%10000) )) {
    print STDERR "$i\t";
    printf STDERR "BeginTime: %.3f\tGetSeqTime: %.3f\tMakeFastQTime: %.3f\tMakeSAMTime: %.3f\n",$begintime,$getseqtime,$makefatime,$makesamtime;
  }
  my $pedist = $peis_result[$i];
  my $mpdist = 0;
  if ($type eq 'mp') {
    $mpdist = int($mpis_result[$i]);
  }

  my $s1 = 0;
  my $s2 = 0;
  my $bigfrag = 0;
  my $id = "$run:";		#.($i+1);
  my $loc1 = '';
  my $loc2 = '';
  my $sid = "";
  if ($type eq 'pe') {
    if ($bound) {
      my $upperbound = int($insertsize_avg + ($insertsize_avg/2));
      $pedist = $upperbound if $pedist > $upperbound;
      my $lowerbound = $insertsize_avg - int($insertsize_avg/1.5);
      $pedist = $lowerbound if $pedist < $lowerbound;
	  
    }
    #print STDERR "PEDIST: $pedist\n";
    my $fwd_strand = int(rand(2));
    my $first = int(rand(2));
    my $genomepos = int(rand($slcsum))+1;
    ###now determine which sequence and coordinate this is
    my $break = 1;
    my $forstart = 1;
    $sid = "";
    foreach my $chrom (@seqlengths) {
      if ($genomepos <= $chrom->{'csum'}) {
	$break = $genomepos -  $chrom->{'offset'} + $minpelen;
	$sid = $chrom->{'sid'};
	last;
      }
    }
    my $slen = $chrlen{$sid};
    my $revend = $slen; 
    if ($first) {
      $forstart = $break; 
      $revend = $break + $pedist - 1;
      $revend = $slen if $revend > $slen;
    } else {
      $revend = $break; 
      $forstart = $break - $pedist + 1;
      $forstart = 1 if $forstart < 1;
    }
    $t1 = Time::HiRes::time();
    my $fragseq = "";
    if ($fwd_strand) {
      if ($BP) {
	$s1  = SeqOp::seqobj_subseq($seqs{$sid},$forstart,$forstart + $rl -1,'+');
	$s2  = SeqOp::seqobj_subseq($seqs{$sid},$revend - $rl + 1,$revend,'-');
	$fragseq = SeqOp::seqobj_subseq($seqs{$sid},$forstart,$revend,'+');
      } else {
	$s1  = SeqOp::get_seq_BioDBFasta($db, $sid,$forstart,$forstart + $rl -1,'+');
	$s2  = SeqOp::get_seq_BioDBFasta($db, $sid,$revend - $rl + 1,$revend,'-');
	$fragseq = SeqOp::get_seq_BioDBFasta($db, $sid,$forstart,$revend,'+');
      }
      if (!$s1) {
	die "$db, $sid,$forstart,$forstart + $rl -1,+\n";
      }
      if (!$s2) {
	die "$db, $sid,$revend - $rl + 1,$revend,-\n";
      }
      my $ncount = 0;
      while ($s1=~/[^ACGTacgt]/g) {
	$ncount++;
      }
      if ($ncount > ($rl/2)) {
	#$i++;
	next;
      }
      $ncount = 0;
      while ($s2=~/[^ACGTacgt]/g) {
	$ncount++;
      }
      if ($ncount > ($rl/2)) {
	#$i++;
	next;
      }
      ### GC Filter ###
      if ($bias) {
	my $gc = gc_pct($fragseq);
	my $prob = (1 - ($bias * (1 - $gcprob{$gc})));
	#if (rand() <= $gcprob{$gc}) {
	if (rand() <= $prob ) {
	  #my $subseq = substr($s1,0,50);
	  #my $lsub = length($subseq);
	  #my $gcp = gc_pct($subseq);	
	  $num_passed_gc_filter++;
	  #$id.=$num_passed_gc_filter;
	} else {
	  next;
	}
      }
      if ($printObservedGC){
	$observed_gc{gc_pct(substr($s1,0,50))}++;
	$observed_gc{gc_pct(substr($s2,0,50))}++;
      }
      $loc1 = "$sid:$forstart:".($forstart + $rl -1).":+";
      $loc2 = "$sid:".($revend - $rl + 1).":".($revend).":-";
    } else {
      if ($BP) {
	$s2  = SeqOp::seqobj_subseq($seqs{$sid},$forstart,$forstart + $rl -1,'+');
	$s1 = SeqOp::seqobj_subseq($seqs{$sid},$revend - $rl + 1,$revend,'-');
	$fragseq = SeqOp::seqobj_subseq($seqs{$sid},$forstart,$revend,'+');
      } else {
	$s2  = SeqOp::get_seq_BioDBFasta($db, $sid,$forstart,$forstart + $rl -1,'+');
	$s1  = SeqOp::get_seq_BioDBFasta($db, $sid,$revend - $rl + 1,$revend,'-');
	$fragseq = SeqOp::get_seq_BioDBFasta($db, $sid,$forstart,$revend,'+');
      }
      if (!$s2) {
	die "$db, $sid,$forstart,$forstart + $rl -1,+\n";
      }
      if (!$s1) {
	die "$db, $sid,$revend - $rl + 1,$revend,-\n";
      }
      my $ncount = 0;
      while ($s1=~/[^ACGTacgt]/g) {
	$ncount++;
      }
      if ($ncount > ($rl/2)) {
	#$i++;
	next;
      }
      $ncount = 0;
      while ($s2=~/[^ACGTacgt]/g) {
	$ncount++;
      }
      if ($ncount > ($rl/2)) {
	#$i++;
	next;
      }
      ### GC Filter ###
      if ($bias) {
	my $gc = gc_pct($fragseq);
	my $prob = (1 - ($bias * (1 - $gcprob{$gc})));
	#if (rand() <= $gcprob{$gc}) {
	if (rand() <= $prob ) {	
	  $num_passed_gc_filter++;
	  #$id.=$num_passed_gc_filter;
	} else {
	  next;
	}
      }
      if ($printObservedGC){
	$observed_gc{gc_pct(substr($s1,0,50))}++;
	$observed_gc{gc_pct(substr($s2,0,50))}++;
      }
      $loc2 = "$sid:$forstart:".($forstart + $rl -1).":+";
      $loc1 = "$sid:".($revend - $rl + 1).":".($revend).":-";
    }
    $t2 = Time::HiRes::time();
  } elsif ($type eq 'mp') {	#Type = MP
    if ($bound) {
      my $upperbound = int($mp_insertsize_avg + ($mp_insertsize_avg/4));
      $mpdist = $upperbound if $mpdist > $upperbound;
      my $lowerbound = ($mp_insertsize_avg - int($mp_insertsize_avg/2));
      $mpdist = $lowerbound if $mpdist < $lowerbound;
    }
    my $fwd_strand = int(rand(2));
    my $contam_strand_fwd = int(rand(2));
    my $add_contam = 0;
    if ($mp_contamination_rate) {
      $add_contam =(int(rand(1/$mp_contamination_rate)) == 0);
    }
    my $first = int(rand(2));
    my $break = 1;
    my $genomepos = int(rand($slcsum))+1;
    ###now determine which sequence and coordinate this is
    my $fragstart = 1;
    $sid = "";
    foreach my $chrom (@seqlengths) {
      if ($genomepos <= $chrom->{'csum'}) {
	$break = $genomepos -  $chrom->{'offset'} + $minpelen;
	$sid = $chrom->{'sid'};
	last;
      }
    }
    my $slen = $chrlen{$sid};
    my $fragend = $slen; 
    if ($first) {
      $fragstart = $break; 
      $fragend = $break + $mpdist - 1;
      $fragend = $slen if ($fragend > $slen);
    } else {
      $fragend = $break; 
      $fragstart = $break - $mpdist + 1;
      $fragstart = 1 if $fragstart < 1;
    }
    $mpdist = $fragend - $fragstart + 1;

    ###############################################
    $bigfrag = '';
    $t1 = Time::HiRes::time();
    if ($BP) {
      $bigfrag = SeqOp::seqobj_subseq($seqs{$sid},$fragstart,$fragend,'+');
    } else {
      $bigfrag = SeqOp::get_seq_BioDBFasta($db, $sid,$fragstart,$fragend,'+');
      die "Couldn't get mp frag with arguments: db: $db  sid: $sid fragstart: $fragstart fragend: $fragend\n" if !$bigfrag;
    }
    
    $t2 = Time::HiRes::time();
    my $doublefrag = $bigfrag x 2;
    my $leftsize = int(rand($pedist));
    if ($fwd_strand) {
      my $coord = ($mpdist-$leftsize-1);
      $s1 = substr($doublefrag,($mpdist-$leftsize-1),$mprl);
      my $rseq = substr($doublefrag,($mpdist+$pedist-$leftsize-$mprl+1),$mprl);
      $rseq =~tr/ACGTacgt/TGCAtgca/;
      $s2 = reverse($rseq);
      my $lds = length($doublefrag);
      die "fragstart: $fragstart fragend: $fragend mpd:$mpdist ped:$pedist ls:$leftsize coord:$coord lds:$lds mprl:$mprl \n\t$s1\n\t$s2\n" if (!$s1 || !$s2);
      my $ncount = 0;
      while ($s1=~/[^ACGTacgt]/g) {
	$ncount++;
      }
      if ($ncount > ($rl/2)) {
	$i++;next;
      }
      $ncount = 0;
      while ($s2=~/[^ACGTacgt]/g) {
	$ncount++;
      }
      if ($ncount > ($rl/2)) {
	$i++;next;
      }
      ### GC Filter ###
      if ($bias) {
	my $pefragseq = substr($doublefrag,($mpdist-$leftsize-1),$pedist);
	my $gc = gc_pct($pefragseq);
	my $prob = (1 - ($bias * (1 - $gcprob{$gc})));
	#if (rand() <= $gcprob{$gc}) {
	if (rand() <= $prob ) {
	  $num_passed_gc_filter++;
	  #$id.=$num_passed_gc_filter;
	} else {
	  next;
	}
      }
      if ($printObservedGC){
	$observed_gc{gc_pct(substr($s1,0,50))}++;
	$observed_gc{gc_pct(substr($s2,0,50))}++;
      }
      if ($leftsize < $mprl) {
	$loc2 = "$sid:".($fragstart + $pedist - $leftsize - $mprl).":".($fragstart + $pedist - $leftsize).":-";
	$loc1 = "$sid:".($fragend -$leftsize).":$fragend:+|$sid:$fragstart:".($fragstart + $mprl -1 - $leftsize).":+";
      } elsif (($pedist - $leftsize) < $mprl) {
	$loc2 = "$sid:$fragstart:".($fragstart + ($pedist - $leftsize)).":-|$sid:".($fragend - ($mprl - ($pedist - $leftsize)) +1).":$fragend:-";
	$loc1 = "$sid:".($fragend -$leftsize).":".($fragend -$leftsize + $mprl - 1).":+";
      } else {
	$loc2 = "$sid:".($fragstart + $pedist - $leftsize - $mprl + 1).":".($fragstart + $pedist - $leftsize).":-";
	$loc1 = "$sid:".($fragend -$leftsize).":".($fragend -$leftsize + $mprl -1).":+";
      }
    } else {
      $s2 = substr($doublefrag,($mpdist-$leftsize-1),$mprl);
      my $rseq = substr($doublefrag,($mpdist+$pedist-$leftsize-$mprl+1),$mprl);
      $rseq =~tr/ACGTacgt/TGCAtgca/;
      $s1 = reverse($rseq);
      my $coord = $mpdist-$leftsize-1;
      my $lds = length($doublefrag);
      die "fragstart: $fragstart fragend: $fragend mpd:$mpdist ped:$pedist ls:$leftsize coord:$coord lds:$lds mprl:$mprl \n\t$s1\n\t$s2\n" if (!$s1 || !$s2);
      my $ncount = 0;
      while ($s1=~/[^ACGTacgt]/g) {
	$ncount++;
      }
      if ($ncount > ($rl/2)) {
	$i++;next;
      }
      $ncount = 0;
      while ($s2=~/[^ACGTacgt]/g) {
	$ncount++;
      }
      if ($ncount > ($rl/2)) {
	$i++;next;
      }
      if ($bias) {
	my $pefragseq = substr($doublefrag,($mpdist-$leftsize-1),$pedist);
	my $gc = gc_pct($pefragseq);
	my $prob = (1 - ($bias * (1 - $gcprob{$gc})));
	#if (rand() <= $gcprob{$gc}) {
	if (rand() <= $prob ) {
	  $num_passed_gc_filter++;
	  #$id.=$num_passed_gc_filter;
	} else {
	  next;
	}
      }
      if($printObservedGC){
	$observed_gc{gc_pct(substr($s1,0,50))}++;
	$observed_gc{gc_pct(substr($s2,0,50))}++;
      }
      if ($leftsize < $mprl) {
	$loc1 = "$sid:".($fragstart + $pedist - $leftsize - $mprl).":".($fragstart + $pedist - $leftsize).":-";
	$loc2 = "$sid:".($fragend -$leftsize).":$fragend:+|$sid:$fragstart:".($fragstart + $mprl -1 - $leftsize).":+";
      } elsif (($pedist - $leftsize) < $mprl) {
	$loc1 = "$sid:$fragstart:".($fragstart + ($pedist - $leftsize)).":-|$sid:".($fragend - ($mprl - ($pedist - $leftsize)) +1).":$fragend:-";
	$loc2 = "$sid:".($fragend -$leftsize).":".($fragend -$leftsize + $mprl -1).":+";
      } else {
	$loc1 = "$sid:".($fragstart + $pedist - $leftsize - $mprl + 1).":".($fragstart + $pedist - $leftsize).":-";
	$loc2 = "$sid:".($fragend -$leftsize).":".($fragend -$leftsize + $mprl -1).":+";
      }
    }
    if ($add_contam) {
      my $skip = 0;
      my $contamid = "$run:";	#.($j++);
      if ($label_contam) {
	$contamid = "contam_".$contamid;
      }
      my $fstart = $pedist - $leftsize + 1;
      my $fend = $mpdist - $leftsize;
      my $break = $min_contam_pelen + int(rand($mpdist - $pedist - 2*($min_contam_pelen)))+($fstart);
      my $contam_pedist = $contam_peis_result[$contam_count++];
      my $l = $break;
      my $r = $break + $contam_pedist - 1;
      my $first = int(rand(2));
      if ($first) {
	$l = $break; 
	$r = $break + $contam_pedist - 1;
	$r = $fend if ($r > $fend);
      } else {
	$r = $break; 
	$l = $break - $contam_pedist + 1;
	$l = $fstart if $l < $fstart;
      }
      $contam_pedist = $r - $l + 1;
      my $cs1 = '';
      my $cs2 = '';
      my $cloc1 = '';
      my $cloc2 = '';
      if ($bias) {
	my $pefragseq = substr($bigfrag,$l-1,$contam_pedist);
	my $gc = gc_pct($pefragseq);
	my $prob = (1 - ($bias * (1 - $gcprob{$gc})));
	#if (rand() <= $gcprob{$gc}) {
	if (rand() <= $prob ) {
	  $j++;
	  $contamid.=$j;
	} else {
	  $skip = 1;
	}
      }
      if (!$skip) {
	if ($contam_strand_fwd) {
	  $cs1 = substr($bigfrag,$l-1,$mprl);
	  my $rseq = substr($bigfrag,$r-$mprl,$mprl);
	  $rseq =~tr/ACGTacgt/TGCAtgca/;
	  $cs2 = reverse($rseq);
	  my $ncount = 0;
	  while ($cs1=~/[^ACGTacgt]/g) {
	    $ncount++;
	  }
	  if ($ncount > ($rl/2)) {
	    $i++;$skip=1;
	  }
	  $ncount = 0;
	  while ($cs2=~/[^ACGTacgt]/g) {
	    $ncount++;
	  }
	  if ($ncount > ($rl/2)) {
	    $i++;$skip=1;
	  }
	  if (!$skip) {
	    $cloc1 = "$sid:".($fragstart + $l -1).":".($fragstart + $l + $mprl -2).":+";
	    $cloc2 = "$sid:".($fragstart +$r - $mprl).":".($fragstart + $r -1).":-";
	  }
	} else {
	  $cs2 = substr($bigfrag,$l-1,$mprl);
	  my $rseq = substr($bigfrag,$r-$mprl,$mprl);
	  $rseq =~tr/ACGTacgt/TGCAtgca/;
	  $cs1 = reverse($rseq);
	  my $ncount = 0;
	  while ($cs1=~/[^ACGTacgt]/g) {
	    $ncount++;
	  }
	  $skip = 1 if $ncount > ($rl/2);
	  $ncount = 0;
	  while ($cs2=~/[^ACGTacgt]/g) {
	    $ncount++;
	  }
	  $skip = 1 if $ncount > ($rl/2);
	  if (!$skip) {
	    $cloc2 = "$sid:".($fragstart + $l -1).":".($fragstart + $l + $mprl -2).":+";
	    $cloc1 = "$sid:".($fragstart +$r - $mprl).":".($fragstart + $r -1).":-";
	  }
	}
	if (length($cs1)<$mprl) {
	  print STDERR "fraglen: ".($fragend-$fragstart+1)." or ".length($bigfrag)." fstart:$fstart fend:$fend l:$l r:$r \n";
	  die "cs1 $cs1 is less than $mprl nt long\n";
	}
	if (length($cs2)<$mprl) {
	  print STDERR "fraglen: ".($fragend-$fragstart+1)." or ".length($bigfrag)." fstart:$fstart fend:$fend l:$l r:$r \n";
	  die "cs2 $cs2 is less than $mprl nt long\n";
	}
	if (!$skip) {
	  if($printObservedGC){
	    $observed_gc{gc_pct(substr($cs1,0,50))}++;
	    $observed_gc{gc_pct(substr($cs2,0,50))}++;
	  }
	  if (eof($p1fq)) {
	    close $p1fq;
	    close $p2fq;
	    ($p1file,$p2file) = getQualFile($qdir,\@p1files);
	    open($p1fq, "<$p1file") or die "$!\n";
	    open($p2fq, "<$p2file") or die "$!\n";
	
	  }
	  my $p1qline = <$p1fq>;
	  my $p2qline = <$p2fq>;
	  chomp($p1qline);
	  chomp($p2qline);
	  my $fq1 = make_fastq_from_template("$contamid#0/1",$cs1,$p1qline,$p1fqout);
	  my $fq2 = make_fastq_from_template("$contamid#0/2",$cs2,$p2qline,$p2fqout);
	  #print SAM using $cloc1 and $cloc2
	  next if !($fq1->[1] && $fq2->[1]);
	  if($pairfile){
	    print $pairfileout "$contamid#0\t",$fq1->[1],"\t",$fq1->[2],"\t",$fq2->[1],"\t",$fq2->[2],"\n";
	  }else{
	    print $p1fqout "\@",$fq1->[0],"\n",$fq1->[1],"\n+\n",$fq1->[2],"\n";
	    print $p2fqout "\@",$fq2->[0],"\n",$fq2->[1],"\n+\n",$fq2->[2],"\n";
	  }

	  my $lseq1=length($fq1->[1]);
	  my $lseq2=length($fq2->[1]);
	  ###MATE 1
	  if ($lseq1 > 0) {
	    my ($seqid,$start,$end,$strand) = split ":",$cloc1;
	    my $segment_length = $end-$start+1;
	    if ($segment_length > $lseq1) {
	      if ($strand eq '+') {
		$end=$start+$lseq1-1;
	      } else {
		$start=$end-$lseq1+1;
	      }
	      $segment_length = $end-$start+1;
	    }
	    my $cigar = $segment_length."M";
	    my $PRIM = 1;
	    die "segment length is $segment_length = $end-$start+1\n$cloc1\n" if $segment_length < 1;
	    print SAM formatSAM({'primary'=>$PRIM,'length'=>$segment_length,'pair'=>"$contamid#0", 'mate'=>1,'seq'=>substr($fq1->[1],0,$segment_length), 'qual'=>substr($fq1->[2],$offset,$segment_length),'ref'=>$sid, 'split'=>0,'start'=>$start,'end'=>$end,'strand'=>$strand, 'cigar'=>$cigar});
	    $lseq1 = $lseq1-$segment_length;
	  }

	  if ($lseq2 > 0) {
	    my ($seqid,$start,$end,$strand) = split ":",$cloc2;
	    my $segment_length = $end-$start+1;
	    if ($segment_length > $lseq2) {
	      if ($strand eq '+') {
		$end=$start+$lseq2-1;
	      } else {
		$start=$end-$lseq2+1;
	      }
	      $segment_length = $end-$start+1;
	    }
	    die "segment length is $segment_length = $end-$start+1\n$cloc2\n" if $segment_length < 1;
	    my $cigar = $segment_length."M";
	    my $PRIM = 1;
	    print SAM formatSAM({'primary'=>$PRIM,'length'=>$segment_length,'pair'=>"$contamid#0", 'mate'=>2,'seq'=>substr($fq2->[1],0,$segment_length), 'qual'=>substr($fq2->[2],$offset,$segment_length),'ref'=>$sid, 'split'=>0,'start'=>$start,'end'=>$end,'strand'=>$strand, 'cigar'=>$cigar});
	    $lseq2 = $lseq2-$segment_length;
	  }
  
	  $i++;
	}
      }
    }
  }
  
  my $id = "$run:".($i+1);
  die "Failed to get sequences\n$bigfrag\n" if !($s1 && $s2);
  if ($printGC) {
    my $gc = 0;
    my $at = 0;
    $gc++ while $s1=~/[GC]/g;
    $at++ while $s1=~/[AT]/g;
    if ($gc + $at) {
      my $gccontent = sprintf( "%1.3f\n",$gc/($gc+$at));
      print GC "$gccontent";
    }
    $gc = 0;
    $at = 0;
    $gc++ while $s2=~/[GC]/g;
    $at++ while $s2=~/[AT]/g;
    if ($gc + $at) {
      my $gccontent = sprintf( "%1.3f\n",$gc/($gc+$at));
      print GC "$gccontent";
    }
  }
  
  
  if (eof($p1fq)) {
    close $p1fq;
    close $p2fq;
    ($p1file,$p2file) = getQualFile($qdir,\@p1files);
    print STDERR "Opening $p1file and $p2file\n";
    open($p1fq, "<$p1file") or die "$!\n";
    open($p2fq, "<$p2file") or die "$!\n";
	
  }
  my $p1qline = <$p1fq>;
  my $p2qline = <$p2fq>;
  chomp($p1qline);
  chomp($p2qline);
  $t3 = Time::HiRes::time();
  my $fq1 = make_fastq_from_template("$id#0/1",$s1,$p1qline,$p1fqout);
  my $fq2 = make_fastq_from_template("$id#0/2",$s2,$p2qline,$p2fqout);
  next if !($fq1->[1] && $fq2->[1]);
  if($pairfile){
    print $pairfileout "$id#0\t",$fq1->[1],"\t",$fq1->[2],"\t",$fq2->[1],"\t",$fq2->[2],"\n";
  }else{
    print $p1fqout "\@",$fq1->[0],"\n",$fq1->[1],"\n+\n",$fq1->[2],"\n";
    print $p2fqout "\@",$fq2->[0],"\n",$fq2->[1],"\n+\n",$fq2->[2],"\n";
  }
  $t4 = Time::HiRes::time();
  # $fq1->[0]  => newid
  # $fq1->[1]   => newseq
  # $fq1->[2]  => qualities

  my @locs1= split /\|/,$loc1;
  my $lseq1=length($fq1->[1]);
  my $lseq2=length($fq2->[1]);
  my @locs2= split /\|/,$loc2;

  ###MATE 1
  my $offset = 0;
  my $remaining_lseq1 = $lseq1;
  foreach my $l (@locs1) {
    if ($remaining_lseq1 > 0) {
      #print STDERR "$l\n";
      my ($seqid,$start,$end,$strand) = split ":",$l;
      my $segment_length = $end-$start+1;
      if ($segment_length > $remaining_lseq1) {
	if ($strand eq '+') {
	  $end=$start+$remaining_lseq1-1;
	} else {
	  $start=$end-$remaining_lseq1+1;
	}
	$segment_length = $end-$start+1;
      }
      my $cigar = $segment_length."M";
      # Do hardmasking. none if $lseq1 == $segment_length;
      if ($lseq1 != $segment_length) {
	if (!$offset) {
	  $cigar = $cigar . ($lseq1-$segment_length).'H';
	} else {
	  $cigar = $offset.'H' . $cigar;
	}
      }
      my $PRIM = 1;
      print SAM formatSAM({'primary'=>$PRIM,'length'=>$segment_length,'pair'=>"$id#0", 'mate'=>1,'seq'=>substr($fq1->[1],$offset,$segment_length), 'qual'=>substr($fq1->[2],$offset,$segment_length),'ref'=>$sid, 'split'=>0,'start'=>$start,'end'=>$end,'strand'=>$strand, 'cigar'=>$cigar});
      
      $remaining_lseq1 = $remaining_lseq1-$segment_length;
      $offset = $offset + $segment_length;
    }
  }

  ###MATE 2
  $offset = 0;
  my $remaining_lseq2 = $lseq2;
  foreach my $l (@locs2) {
    if ($remaining_lseq2 > 0) {
      my ($seqid,$start,$end,$strand) = split ":",$l;
      my $segment_length = $end-$start+1;
      if ($segment_length > $remaining_lseq2) {
	if ($strand eq '+') {
	  $end=$start+$remaining_lseq2-1;
	} else {
	  $start=$end-$remaining_lseq2+1;
	}
	$segment_length = $end-$start+1;
      }
      my $cigar = $segment_length."M";
      if ($lseq2 != $segment_length) {
	if (!$offset) {
	  $cigar = $cigar . ($lseq2-$segment_length).'H';
	} else {
	  $cigar = $offset.'H' . $cigar;
	}
      }
      my $PRIM = 1;

      print SAM formatSAM({'primary'=>$PRIM,'length'=>$segment_length,'pair'=>"$id#0", 'mate'=>2,'seq'=>substr($fq2->[1],$offset,$segment_length), 'qual'=>substr($fq2->[2],$offset,$segment_length),'ref'=>$sid, 'split'=>0,'start'=>$start,'end'=>$end,'strand'=>$strand, 'cigar'=>$cigar});
      $remaining_lseq2 = $remaining_lseq2-$segment_length;
      $offset = $offset + $segment_length;
    }
  }

  $t5 = Time::HiRes::time();
  
  $begintime += $t1-$t0;
  $getseqtime += $t2-$t1;
  $makefatime += $t4-$t3;
  $makesamtime += $t5-$t4;
  $i++;
}
printf STDERR "BeginTime: %.3f\tGetSeqTime: %.3f\tMakeFastQTime: %.3f\tMakeSAMTime: %.3f\n",$begintime,$getseqtime,$makefatime,$makesamtime;

close GC if $printGC;
close SAM;
close $p1fq;
close $p2fq;
if ($pairfile){
  close $pairfileout;
}else{
  close $p1fqout;
  close $p2fqout;
}


`samtools view -S -b $run.reads.sam | samtools fixmate - $run.reads.fixmate.bam`;
`samtools sort $run.reads.fixmate.bam $run.reads`;
`samtools index $run.reads.bam`;
unlink "$run.reads.sam";
unlink "$run.reads.fixmate.bam";
#`rm -r /scratch_tmp/$uid`;
my $fq1 = $run."_1.fastq";
my $fq2 = $run."_2.fastq";
`bzip2 $fq1`;
`bzip2 $fq2`;
if ($printObservedGC){
  open(GCOBS,">$run.observed_gc.txt") or die "$!\n";
  print GCOBS "GC\tobserved\n";
  foreach my $k (sort {$a<=>$b} keys %observed_gc) {
    print GCOBS "$k\t",($observed_gc{$k}/$N),"\n";
  }
  close GCOBS;
}
#run samtools to make bam and then fixmate and then sort by position and index
print STDERR "samtools tview $run.reads.bam $genome_mfa_orig\n";
sub make_fastq_from_template
  {
    my $id = shift;
    my $seqstring = shift;
    my $qstring_temp = shift;
    my $qstring = substr($qstring_temp,0,length($seqstring));
    $qstring=~s/B+$//;
    my $qlen = length($qstring);
    if ($qlen < $minreadlength) {
      return [$id,0,0];
    }
    my $fhout = shift;
    my $seqlen = length($seqstring);
    $seqlen = $qlen if $qlen < $seqlen;
    chomp $seqstring;
    #uc($seqstring);
    $seqstring=~s/([ACGTacgt])([^ACGTacgt])([ACGTacgt])/$1.($normal{'N'}->[int(rand 3)]).$3/gi;
    my @nts = split '',$seqstring;
    my @ascii_quals = split '',$qstring;
    #my @numeric_quals = map(qual_from_ascii($_), @ascii_quals);
    my $newseq = '';
    my $phase = int(rand 2);
    for (my $n = 0; $n < $seqlen; $n++) {
      my $nt = uc($nts[$n]);
      my $QC = 'B';
      my $p = 1;
      if ($nt !~/[ACGT]/) {
	#$Q=2;
	substr($qstring,$n,1,'B');
	$newseq .= $nt;
      } else {
	$QC = $ascii_quals[$n];
	if ($QC eq 'B') {
	  $p = $prob_from_qual{2+int(rand(8))};
	} else {
	  $p = $prob_from_ascii{$QC}; #exp(($Q*log(10))/-10);
	}
	if (!make_error($p,$error)) {
	  $newseq .= $nt;
	} else {
	  #print STDERR "Making an error at $n ($nt)\n";
	  if ($QC ne 'B') {
	    $newseq .= $normal{$nt}->[int(rand 3)];
	  } else {		### B illumina QC clipping or for Ns
	    #print STDERR "Q == 2\n";
	    if (($n == ($seqlen - 1)) || ($ascii_quals[$n + 1] ne 'B')) {
	      $phase = 1;	#  $newseq .= 'N';
	    }			#else {
				#print STDERR "$nt $blist{$nt}\n";
	    if ($phase) {
	      my $prev = uc($nts[$n-1]);
	      my @possibilities = (@{$blist{$nt}},$prev,$prev,$prev,$prev,$prev,$prev,$prev,$prev,$prev,$prev,$prev,$prev,$prev,$prev,$prev,$prev);
	      $newseq .= $possibilities[int(rand (scalar @possibilities))];
	    } else {  
	      my $next = uc($nts[$n+1]);
	      my @possibilities = (@{$blist{$nt}},$next,$next,$next,$next,$next,$next,$next,$next,$next,$next,$next,$next,$next,$next,$next,$next);
	      $newseq .= $possibilities[int(rand (scalar @possibilities))];
	    }
	    #}
	  } 
	} 
      }
    }
    #print $fhout "\@$id\n$newseq\n+\n$qstring\n";
    return [$id,$newseq,solexa2phred($qstring)];

  }

sub qual_from_ascii
  {
    my $ascii = shift;
    return ord($ascii) - 64;
  }
sub make_error
  {
    my $prob = shift;
    my $error_factor = shift;
    if ($error_factor) {
      $prob = $prob * $error_factor;
      my $random_number = int(rand((1/$prob)));
      return $random_number?0:1;
    } else {
      return 0;
    }
  }

sub generate_random_password 
  {
    my $passwordsize = shift;
    my @alphanumeric = ('a'..'z', 'A'..'Z', 0..9);
    my $randpassword = join '', 
      map $alphanumeric[rand @alphanumeric], 0..$passwordsize;

    return $randpassword;
  }

sub getQualFile{

  my $qd = shift;
  my $p1filesref = shift;
  my $qf1 = $p1filesref->[int(rand(scalar(@$p1filesref)))];
  my $qf2 = $qf1;
  $qf2 =~ s/pair1/pair2/;
  return ("$qdir/$qf1","$qdir/$qf2");
}

sub formatSAM{
  
  # No. Name Description
  # 1 QNAME Query NAME of the read or the read pair
  # 2 FLAG Bitwise FLAG (pairing, strand, mate strand, etc.)
  # 3 RNAME Reference sequence NAME
  # 4 POS 1-Based leftmost POSition of clipped alignment
  # 5 MAPQ MAPping Quality (Phred-scaled)
  # 6 CIGAR Extended CIGAR string (operations: MIDNSHP)
  # 7 MRNM Mate Reference NaMe (‘=’ if same as RNAME)
  # 8 MPOS 1-Based leftmost Mate POSition
  # 9 ISIZE Inferred Insert SIZE
  # 10 SEQ Query SEQuence on the same strand as the reference
  # 11 QUAL Query QUALity (ASCII-33=Phred base quality
  # 12 OPT variable OPTional fields in the format TAG:VTYPE:VALUE

  # Flag	Chr	Description
  # 0x0001	p	the read is paired in sequencing
  # 0x0002	P	the read is mapped in a proper pair
  # 0x0004	u	the query sequence itself is unmapped
  # 0x0008	U	the mate is unmapped
  # 0x0010	r	strand of the query (1 for reverse)
  # 0x0020	R	strand of the mate
  # 0x0040	1	the read is the first read in a pair
  # 0x0080	2	the read is the second read in a pair
  # 0x0100	s	the alignment is not primary (secondary)
  # 0x0200	f	the read fails platform/vendor quality checks
  # 0x0400	d	the read is either a PCR or an optical duplicate
  my %flag = (
	      'p'=>hex('0x0001'), #paired
	      'P'=>hex('0x0002'), #mapped in Proper pair
	      'u'=>hex('0x0004'), #unmapped
	      'U'=>hex('0x0008'), #mate unmapped
	      'r'=>hex('0x0010'), #this read is on minus strand and seq is revcomp
	      'R'=>hex('0x0020'), #mate is on minus strand and seq is revcomp
	      '1'=>hex('0x0040'), #read 1
	      '2'=>hex('0x0080'), #read 2
	      's'=>hex('0x0100'), #secondary
	      'f'=>hex('0x0200'), #fail
	      'd'=>hex('0x0400')  #dup
	     );
  my $no_store_seq = 0;
  my $bitwise_flag = 0;
  my $num_maps = scalar @_;
  my $record;
  my $mate;
  $record = shift;
  my $rnext = '*';
  my $pnext = 0;
  my $seq = '*';
  my $qual = '*';
  if ($record->{primary}) {
    $seq = $record->{seq};
    $qual = $record->{qual};

    #print STDERR "$seq\n";
  }
  my $opt = '';			#'XG:Z:'.$record->{'gem-match'};
  if ($num_maps>1) {		#mated
    $mate = shift;
    $bitwise_flag += $flag{'p'};
    $bitwise_flag += $flag{'P'};
    $bitwise_flag += $flag{'s'} if !$record->{primary};
    $bitwise_flag += $flag{'r'} if $record->{strand} eq '-';
    $bitwise_flag += $flag{'R'} if $mate->{strand} eq '-';
    if ($record->{mate} == 1) {
      $bitwise_flag += $flag{'1'};
    } else {
      $bitwise_flag += $flag{'2'};
    } 
    $rnext = '=';
    $pnext = $mate->{start};
    my $isize = 0;
    if ($mate->{start} > $record->{end}) {
      $isize =  $mate->{end} - $record->{start} + 1;
    } else {
      $isize =  $mate->{start} - $record->{end} - 1;
    }
    if ($no_store_seq) {
      $seq = '*';
      $qual = '*';
    } else {
      if ($record->{strand} eq '-') {
	my $rseq = reverse($seq);
	$rseq =~ tr/ACGT/TGCA/;
	$seq = $rseq;
	my $rqual = reverse($qual);
	$qual = $rqual;
      }
    }
    return (join("\t",($record->{pair},$bitwise_flag,$record->{ref},$record->{start},'255',$record->{cigar},$rnext,$pnext,$isize,$seq,$qual))."\n");
  
  } else {			#unmated
    #$bitwise_flag += $flag{'p'};
    $bitwise_flag += $flag{'p'};
    $bitwise_flag += $flag{'P'};
    
    $bitwise_flag += $flag{'s'} if !$record->{primary};
    $bitwise_flag += $flag{'r'} if $record->{strand} eq '-';
    if ($record->{mate} == 1) {
      $bitwise_flag += $flag{'1'};
    } else {
      $bitwise_flag += $flag{'2'};
    } 
    my $isize = 0;
    
    if ($no_store_seq) {
      $seq = '*';
      $qual = '*';
    } else {
      if ($record->{strand} eq '-') {
	my $rseq = reverse($seq);
	$rseq =~ tr/ACGT/TGCA/;
	$seq = $rseq;
	my $rqual = reverse($qual);
	$qual = $rqual;
      }
    }
    return (join("\t",("$record->{pair}",$bitwise_flag,$record->{ref},$record->{start},'255',$record->{cigar},$rnext,$pnext,$isize,$seq,$qual))."\n");
  
  }

}

sub solexa2phred{
  my $qin = shift;
  my @ch = split '',$qin;
  my $qout = '';
  foreach my $c (@ch) {
    $qout.=chr(ord($c)-31);
  }                             #
  return $qout;
}
sub gc_pct{
  my $seq = shift;
  #print STDERR "$seq\n";
  my $gc = 0;
  my $at = 0;
  $gc++ while $seq=~/[GC]/ig;
  $at++ while $seq=~/[AT]/ig;
  if ($gc + $at) {
    my $gccontent = sprintf( "%1.2f",(2*(sprintf( "%1.0f",(100*($gc/($gc+$at)))/2)))/100);
    if ($gccontent eq '0.57') {
      die "$gc / $gc + $at\n".(100*($gc/($gc+$at)))." ".(sprintf( "%1.0f",(100*($gc/($gc+$at)))/2))." ".((2*(sprintf( "%1.0f",(100*($gc/($gc+$at)))/2)))/100)."\n";
    }
    return $gccontent;
  } else {
    return "1.00";
  }
}
