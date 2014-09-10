#!/usr/bin/perl

##
##
##  Create CAGEscan cluster from paired-end data in BED12, BAM or SAM format
##  See usage() below for a detaillled descrition of the script
##
##  Remark :
##    * mandatory use of BedTools intersectBed version 2.9.0
##    * in part optional use of BedTools mergeBed version 2.9.0
##    *                         (home-brewed) pairedBamToBed12 (re-enginering of BedTools bamToBed)
##                              samtools version: 0.1.7 (r510)
##    * default paths to intersectBed, bamToBed and mergeBed and samtools are hardcoded
##    * read/write access to /tmp/ optional, but recommanded
##
## Nicolas Bertin <nbertin@gsc.riken.jp>
##       Created : Wed Feb  2 16:53:59 JST 2011
## Last modified : Wed Feb 27 13:05:19 JST 2013
##
##

use strict;
use Getopt::Long;
use File::Basename;


my $samtools_bin = "/usr/bin/samtools";
my $bamtobed_bin = "/usr/bin/pairedBamToBed12";
my $intersectbed_bin = "/usr/bin/intersectBed";
my $mergebed_bin = "/usr/bin/mergeBed";
my $tmp_bed6_tss_input_file = '/tmp/' . join(".", basename($0, ".pl"), 'bed6_tss_formatted_input', $$, 'tmp');
my $tmp_bed6_cluster_file = '/tmp/' . join(".", basename($0,".pl"), 'bed6_induced_cluster', $$, 'tmp');
my @valid_input_format = qw(bed  sam  bam  bed.gz  sam.gz  bam.gz);

my $help;
my $input_format;
my $input_file;
my $predefined_cluster_file;
my $presorted_bed6_tss_input_file;
my $min_map_quality;
my $denovo_clustering_distance = 0;
my $stored_denovo_cluster_file;
my $stored_uprocessed_bam_file;
my $usage;
my $verbose = 0;

my $opt = GetOptions(
		     "help"              => \$help,            
		     "usage"             => \$usage,
		     "verbose=i"         => \$verbose,
		     "cluster_file=s"               => \$predefined_cluster_file,
		     "input=s"                      => \$input_file,
		     "format=s"                     => \$input_format,
		     "min_map_quality|q=s"          => \$min_map_quality,
		     "denovo_cluster_distance=s"    => \$denovo_clustering_distance,
		     "stored_denovo_cluster_file|s=s"  => \$stored_denovo_cluster_file,
		     "stored_uprocessed_bam_file|x=s"  => \$stored_uprocessed_bam_file,
		     "bin_samtools=s"     => \$samtools_bin,
		     "bin_bamtobed=s"     => \$bamtobed_bin,
		     "bin_mergebed=s"     => \$mergebed_bin,
		     "bin_intersectbed=s" => \$intersectbed_bin,
		     "tmp_bed6_tss_input_file=s"  => \$tmp_bed6_tss_input_file ,
		     "tmp_bed6_cluster_file=s" => \$tmp_bed6_cluster_file ,
                     "presorted_bed6_tss_input_file" => \$presorted_bed6_tss_input_file,
		    );

sub usage(){
  die <<HERE;
usage: $0 -i <input bed12, bam or sam file to be clustered>
       type “perldoc $0” for details
HERE
}

=pod

=head1 NAME

CAGEscan-Clustering - Create CAGEscan cluster from BED12, BAM or SAM pairs

=head1 SYNOPSIS

CAGEscan-Clustering.pl -i <input bed12, bam or sam file to be clustered>

=head1 DESCRIPTION

Create CAGEscan cluster from paired-end data in BED12, BAM or SAM
format.  Input data can be sent to STDIN, in which case the format
(C<-f>) of the input must be specified, or read from a file C<-i>, in
which case the format of the data can be specified using the C<-f>
option, in the absence of which it is guessed from the suffix of the
input file.  Valid suffixes / format are: 'bed', 'sam' 'bam' or gziped
version of those formats).

Importantly if the data is sent in BAM or SAM format, B<it must be
sorted by name>.

An optional BED6 formatted cluster file (C<-c>) can be used to guide
the construction of CAGEscan clusters. Note that only paired-end data
whose TSS overlap with at least one cluster will be grouped.

In the absence of such an optional BED6 formatted cluster file,
"guiding" clusters are created from the paired-end data themselves
using a simple "single-linkage clustering" approach using the distance
C<-d> (ie clustering together all the data based on their TSS being
distant of C<-d> basepair). The result of the TSS single-linkage-based
clustering can be stored (BED6) if a file path (C<-s>) is provided or
is otherwise discarded.

Output data will be BED12 formatted CAGEscan cluster.

=over

=item

the name of each cluster will correspond to the name of the guiding
cluster (BED6 formatted cluster name), or, if clusters are created
from the paired-end data, will correspond to:
C<Lx_chr_strand_start_end> with the chr,strand,start and end of the
TSS single-linkage derived clusters.

=item

The score of each cluster will correspond to the number of input
paired-end tags composing the resulting CAGEscan cluster.

=item

The hallmark of the TSS/guiding cluster on the CAGEscan cluster will
be encoded by a 3'UTR region with :

=over

=item

The I<cdsStart> will correspond to either the start of CAGEscan cluster
or end of the guiding cluster (when on the plus strand), or to the
start of the guiding cluster or end of the CAGEScan cluster (when on
the minus strand).

=item

Similarly, I<cdsEnd> will correspond to either the start of the
guiding cluster or end of the CAGEscan cluster (when on the plus
strand), or to the start of the CAGEscan cluster or end of the guiding
cluster (when on the minus strand).

=back

=back

=head2 Note about using guiding cluster

(instead of de-novo clustering from the paired-end TSS)

Whenever the boundaries of the guiding cluster extend beyond those of
the intersected paired-end input data, the reported CAGEscan cluster
boundaries will be those of the paired-end data, not that of the
guiding cluster.  As it can be problematic to integrate multiple run
with various input paired-end data, keep in mind that the cluster name
is that of the guide and shall be used.

=head2 Note about dealing with very large input

The rate limiting step of this script is the preparation of the
paired-end input data that need to be transformed into a BED6 with
I<start> and I<end> corresponding to the 0-based pe_tss position, name
corresponding to the semicolon delimited concatenated full BED12 line
and (more importantly and time-intensive) sorted by TSS position.

This step can be `precomputed` and skipped providing the Boolean flag
C<-p> which will informs the script that the provided input file/stdin
already correspond to the precomputed data content, format and
sorting. Also since CAGEscan cluster are by definition, assemblies of
reads located on the same chromosome and strand, I recommend for large
files to split up the input paired-end reads by chromosome and strand
and precompute the sorted BED6 internal temporary file (with I<start>
and I<end> corresponding to the 0-based pe_tss position, I<name>
corresponding to the semicolon delimited concatenated full BED12 line
and sorted by TSS position).

=head1 OPTIONS

=over 8

=item B<-h>

Show this message

=item B<-v>

Verbose level (1-9)

=back

=head2 BASIC OPTIONS

=over 8

=item B<-i,--input>

Input paired-end data file. Can be a BAM, SAM, BED or gzip compressed
BED12 file.  It is recommended to use sensible suffix for the detection
of the format of the file, although this can also be specified using
the C<-f> option.  Note that this parameter is not mandatory, as input
data can also be read from STDIN.

=item B<-q,--min_map_quality>

The minumum "mapping quality value" for including the input data into
a CAGEscan cluster.  If the input is SAM or BAM this value corresponds
to the sum of each pair member MapQ.  If the input is BED12 , this
value is extracted from the score column.  This value is inclusive
(aka "bigger or equal than").  [Default : undefined, aka no filtering]

=item B<-f,--format>

Format of the input (paired-end) data. Valid format are: C<bed>,
C<sam>, C<bam>, C<bed.gz>, C<sam.gz>, or C<bam.gz>.  C<-f> is
mandatory if the input comes from STDIN but optionnal if the input
(C<-i>) is a file ending with one of the valid formats.

=item B<-c,--cluster_file>

Path to a BED6 formatted cluster file directing the construction of
CAGEscan clusters.  Note that only paired-end data whose TSS overlap
with at least one cluster will be grouped.  If none provided, the input
data itself will be used.

=item B<-d,--denovo_cluster_distance>

The distance for paired-end data TSS single linkage clustering.  If no
BED6 formatted cluster file directing the construction of CAGEscan
clusters were provided, the paired-end data themselves using a simple
"single-linkage clustering" approach using this distance to group
neighboring TSS.

=item B<-s,--stored_denovo_cluster_file>

The path of the BED6 file where (TSS) cluster file directing the
construction of CAGEscan clusters is to be stored (otherwise stored in
a temporary file).

=item B<-x,--stored_uprocessed_bam_file>

The path where non properly paired or pairs not fulfilling the
min_map_quality criteria of BAM or SAM input can be stored.

=item B<-p,--presorted_bed6_tss_input_file>

Boolean flag skipping the transformation of the input paired-end data
into a BED6 with I<start> and I<end> correspond to the 0-based pe_tss
position name correspond to the semicolon delimited concatenated full
BED12 line sort by TSS position.  B<Importantly>, it of course means
that the input file (C<-i>) B<must> conform to the above mentionned
specifications

=back

=head2 OTHER OPTIONS (path to binaries and temp files)

=over 8

=item B<--tmp_bed6_tss_input_file>

Path to the location where a temporary BED12 formatted file holding
the input data might be created.  This is in particular the case when
"guiding" clusters are created from the paired-end data themselves.
[default '/tmp/CAGEscan-Clustering.bed12_formatted_input.PID.tmp']

=item B<--tmp_bed6_cluster_file>

Path to the location where a temporary BED6 formatted file holding the
"guiding" clusters can be created.  [default
'/tmp/CAGEscan-Clustering.bed6_formatted_induced_cluster.PID.tmp']

=item B<--bin_samtools>

Path to C<samtools> binary ("samtools view -Sb ..." used only when
input is in SAM format).

=item B<--bin_bamtobed>

Path to C<pairedBamToBed12> (re-enginering of BEDTools C<bamToBed>)
binary which transform properly paired BAM alignments into a single
BED entry (used only when input is in BAM or SAM format).

=item B<--bin_mergebed>

Path to BedTools::mergeBed binary (used only when "guiding" clusters
must be generated from single-linkage clustering of the the paired-end
data themselves).

=item B<--bin_intersectbed>

Path to BedTools C<intersectBed> binary (always used).

=back

=cut

## print usage and exit
if (defined $usage){ usage(); exit;}
if (defined $help) { usage(); exit;}

## check the input format
unless ($input_format){
  if (defined $input_file){
    my ($name,$path,$suffix) = fileparse($input_file,@valid_input_format);
    $input_format = $suffix;
  }
  else{
    usage();
    die $0,"*****ERROR: you must provide the format of the stdin-ed data via the option '--format' or a file via the option '--input file', the format of which will be guessed from the its suffix\n";
  }
}

my $cmd;
if (defined $input_file){
  if    ($input_format =~ /.gz$/i ){ $cmd .= "zcat $input_file | "; }
  else                             { $cmd .= "cat  $input_file | "; }
}

if   ($input_format =~ /bed/i ){
  $cmd .= " awk -v minQ=$min_map_quality '{if (\$5 >= minQ) print \$0}' | " if ($min_map_quality);
}
elsif ($input_format =~ /sam|bam/i ){
  if ($presorted_bed6_tss_input_file){
     usage();
     die $0,"*****ERROR: incompatible input format '$input_format', given you provided a presorted bed6 tss input\n";
  }
  if ($input_format =~ /sam/i )   { $cmd .= "$samtools_bin view -Sb - | "; }
  $cmd .= " $bamtobed_bin ";
  if ($min_map_quality)           { $cmd .= " -qual $min_map_quality "; }
  if ($verbose == 0)              { $cmd .= " -quiet "; }
  if ($stored_uprocessed_bam_file){ $cmd .= " -x $stored_uprocessed_bam_file "; }
  $cmd .= " | ";
}
else {
  usage();
  die $0,"*****ERROR: unknown input format '$input_format', valid values are '", join("','", @valid_input_format), "'\n";
}

unless ($presorted_bed6_tss_input_file){
  ## transform the input bed12 formatted pe data into a bed6 with
  ##   start and end correspond to the 0based pe_tss position
  ##   name correspond to the semicolon delimited concatenated full bed12 line
  ##   sort by tss position
  $cmd .= " sed -e 's/\t/;/g'                                          \\
               | awk 'BEGIN{FS=\";\"}{OFS=\"\t\"}                      \\
                   {                                                   \\
                   if (\$6 == \"+\") print \$1,\$2,\$2+1,\$0,\".\",\$6;   \\
                   else           print \$1,\$3-1,\$3,\$0,\".\",\$6;   \\
                   }'                                                  \\
               |  sort -k1,1d -k6,6d -k2,2n                            \\
               |  ";
}

if (defined $predefined_cluster_file){
  ## predefined_cluster_file is provided no need to built it from the input
  ## intersect the transformed input bed12 formatted pe data with the bed6
  ## cluster, reporting both the "pe bed12" and "cluster bed6" content
  $cmd .= " $intersectbed_bin -s -wo -a stdin -b $predefined_cluster_file | sort -k10,10d";
  if ($verbose > 0){ print STDERR "merging predefined cluster in $predefined_cluster_file with $input_format input data ...\n"; }
  if ($verbose > 2){ print STDERR $cmd, "\n"; }
}
else{
  ## predefined_cluster_file is not provided; must built it from the input
  ## cluster names follow the format Lx_<chr>_<strand>_<start>_<end>
  my $bed6_cluster_file = (defined $stored_denovo_cluster_file)?$stored_denovo_cluster_file:$tmp_bed6_cluster_file;
  my $bed6_tss_input_file = $tmp_bed6_tss_input_file;
  $cmd .= " tee $bed6_tss_input_file                                         \\
               | $mergebed_bin -s -i stdin -d $denovo_clustering_distance    \\
               | awk '{OFS=\"\t\"}{print \$1,\$2,\$3,\"Lx_\"\$1\"_\"\$4\"_\"\$2\"_\"\$3,\"1\",\$4}'> $bed6_cluster_file";
  if ($verbose > 0){ print STDERR "computing single linkage cluster -d $denovo_clustering_distance from $input_format input data (this may take some time)...\n"; }
  if ($verbose > 2){ print STDERR $cmd, "\n"; }
  system($cmd);
  ## intersect the transformed input bed12 formatted pe data with the de novo
  ## bed cluster, reporting both the "pe bed12" and "cluster bed6" content
  $cmd = "$intersectbed_bin -s -wo -a $bed6_tss_input_file -b $bed6_cluster_file | sort -k10,10d";
  if ($verbose > 0){ print STDERR "now, merging de-novo single linkage cluster with $input_format input data...\n"; }
  if ($verbose > 2){ print STDERR $cmd, "\n"; }
}


my @curr_cluster;
my $curr_cluster_key;
my %bsize;
my $pe_count=0;
## read each "pe bed12"+"cluster bed6" reported intersectbed output
## use the "cluster bed6 as the key to group pe data content
open(CLUST, "$cmd |");
while(<CLUST>){
  chomp;
  my @data = split /\t/;
  my @cluster = @data[6,7,8,9,11];
  my $cluster_key = join(";", @cluster);
  my $cluster_start = $cluster[1];
  my $cluster_end = $cluster[2];

  my @bed12 = split /;/, $data[3];
  my $pe_start = $bed12[1];
  my $pe_end = $bed12[2];
  my $pe_block_count = $bed12[9];
  my @pe_block_size = split /,/, $bed12[10];
  my @pe_block_start = split /,/, $bed12[11];

  unless (scalar @curr_cluster){
    @curr_cluster = @cluster;
    $curr_cluster_key = $cluster_key;
    %bsize = ();
    $pe_count = 0;
  }

  if ($cluster_key eq $curr_cluster_key){
    ## this pe is a part of the current cluster
    ## fill out the sparse array bsize{'absolute_pe_bloc_start_pos'}='(max_seen)pe_block_size'
    for my $i (0..$pe_block_count-1){
      if ($verbose > 3){ print STDERR join(" ", 'pe_block in cluster', $i, $pe_block_start[$i], $pe_block_size[$i]), "\n";  }
      if ( $pe_block_size[$i] > $bsize{$pe_block_start[$i]+$pe_start} ) {
	$bsize{$pe_block_start[$i]+$pe_start} = $pe_block_size[$i];
	if ($verbose > 3){ print STDERR "\t", join(" ", 'selected start and len', $i, $pe_block_start[$i], $pe_block_size[$i]), "\n"; }
      }
    }
    $pe_count++;
  }
  else{
    ## this pe is a part of the next cluster
    ## gather current cagescan cluster data to be exported
    ## construct the list of cagecluster block sizes and (absolute) start positions
    if ($verbose > 3){   print STDERR join(" ", 'pe_block in next cluster process previous'), "\n"; }
    my @istart = sort {$a<=>$b} keys %bsize;
    my @ostart;
    my @osize;
    push @ostart , shift @istart;
    push @osize , $bsize{$ostart[$#ostart]};
    if ($verbose > 3){ print STDERR join("\t", 'first block', @ostart, @osize), "\n"; }

    while (defined (my $istart = shift @istart)){
      if ($verbose > 3){ print STDERR join("\t", 'next block', $istart, $bsize{$istart}), "\n";  }
      if ($istart > $ostart[$#ostart] +  $osize[$#osize]){
	if ($verbose > 3){ print STDERR "\t", join("\t", 'this input block defines a new dest block', $ostart[$#ostart], $istart, $bsize{$istart}), "\n"; }
	push @ostart, $istart;
	push @osize , $bsize{$istart};
	if ($verbose > 3){ print STDERR "\t", join("\t",'last dest block', $ostart[$#ostart] , $osize[$#osize]), "\n" }
      }
      else{
	if ($verbose > 3){ print STDERR "\t", join("\t", 'this input block DO NOT defines a new dest block', $ostart[$#ostart], $istart, $bsize{$istart}), "\n"; }
	if ($istart + $bsize{$istart} > $ostart[$#ostart] + $osize[$#ostart]){
	  if ($verbose > 3){ print STDERR "\t", join("\t", 'its end point ', $istart + $bsize{$istart}, 'is further modify end point', $ostart[$#ostart] + $osize[$#ostart] ), "\n"; }
	  $osize[$#osize] =  $istart + $bsize{$istart} - $ostart[$#ostart];
	  if ($verbose > 3){ print STDERR "\t", join("\t",'last dest block', $ostart[$#ostart] , $osize[$#osize],  $ostart[$#ostart] + $osize[$#osize]), "\n" }
	}
	else{
	  if ($verbose > 3){ print STDERR "\t", join("\t", 'its end point ', $bsize{$istart}, 'is NOT further the previous  end point',  $osize[$#ostart] ), "\n"; }
	  if ($verbose > 3){ print STDERR "\t", join("\t",'last dest block', $ostart[$#ostart] , $osize[$#osize]), "\n" }
	}
      }
    }

    ## gather the data to be stdout-ed
    my $ochr = $curr_cluster[0];
    my $ostart = $ostart[0];
    my $oend = $ostart[$#ostart] + $osize[$#osize];
    my $oname = $curr_cluster[3];
    my $oscore = $pe_count;
    my $ostrand = $curr_cluster[4];
    ## bed12 cdsStart|cdsStart is used to mark the 3'extend of the cluster
    my ($oostart, $ooend) = ($curr_cluster[4] eq '+')?($ostart, $curr_cluster[2]):($curr_cluster[1], $oend);
    my $ocolor = '255,0,0';
    my $obcount =  scalar @osize;
    ## correct the cagecluster block starts to be relative to the 1st cagecluster block
    foreach (@ostart){ $_ = $_ - $ostart }

        ## make sure that cdsEnd is at most the end of the CAGEscan cluster (proper cdsStart has been dealt with
        ## when building the cluster)
        $ooend = $oend if ($ooend > $oend);
        $oostart = $ostart if ($oostart < $ostart);
    ## stdout the current cagescan cluster data
    print STDOUT join("\t", $ochr,  $ostart, $oend, $oname, $oscore, $ostrand, $oostart, $ooend, $ocolor, $obcount, join(",", @osize), join(",", @ostart)), "\n";

    ## re-initialize the new current cluster
    @curr_cluster = @cluster;
    $curr_cluster_key = join(";", @curr_cluster);
    %bsize = ();
    $pe_count = 0;
    ## this pe is a part of the (newly initialize) current cluster
    ## fill out the sparse array bsize{'absolute_pe_bloc_start_pos'}='pe_block_size'
    for my $i (0..$pe_block_count-1){
      if ( $pe_block_size[$i] > $bsize{$pe_block_start[$i]+$pe_start} ) {
	$bsize{$pe_block_start[$i]+$pe_start} = $pe_block_size[$i];
	#if ($verbose > 3){ print STDERR join("\t", 'bsize', $bsize{$pe_block_start[$i]+$pe_start}, $pe_block_size[$i]), "\n"; }
      }
    }
    $pe_count++;
  }
}

## last cluster
## construct the list of cagecluster block sizes and (absolute) start positions
if ($verbose > 3){ print STDERR join(" ", 'process last cluster'), "\n";}
my @istart = sort {$a<=>$b} keys %bsize;
my @ostart;
my @osize;
if ($verbose > 3){ for my $start (@istart){  print STDERR "start ", $start," of size ", $bsize{$start}, "\n"; }}
push @ostart , shift @istart;
push @osize , $bsize{$ostart[$#ostart]};

if ($verbose > 3){ print STDERR join("\t", 'first block', @ostart, @osize), "\n"; }
while (defined (my $istart = shift @istart)){
  if ($verbose > 3){ print STDERR join("\t", 'next block', $istart, $bsize{$istart}), "\n"; }
  if ($istart > $ostart[$#ostart] +  $osize[$#osize]){
    if ($verbose > 3){ print STDERR "\t", join("\t", 'this input block defines a new dest block', $ostart[$#ostart], $istart, $bsize{$istart}), "\n"; }
    push @ostart, $istart;
    push @osize , $bsize{$istart};
    if ($verbose > 3){ print STDERR "\t", join("\t",'last dest block', $ostart[$#ostart] , $osize[$#osize]), "\n" }
  }
  else{
    if ($verbose > 3){ print STDERR "\t", join("\t", 'this input block DO NOT defines a new dest block', $ostart[$#ostart], $istart, $bsize{$istart}), "\n"; }
    if ($istart + $bsize{$istart} > $ostart[$#ostart] + $osize[$#ostart]){
      if ($verbose > 3){ print STDERR "\t", join("\t", 'its end point ', $istart + $bsize{$istart}, 'is further modify end point', $ostart[$#ostart] + $osize[$#ostart] ), "\n"; }
      $osize[$#osize] =  $istart + $bsize{$istart} - $ostart[$#ostart];
      if ($verbose > 3){ print STDERR "\t", join("\t",'last dest block', $ostart[$#ostart] , $osize[$#osize],  $ostart[$#ostart] + $osize[$#osize]), "\n" }
    }
    else{
      if ($verbose > 3){ print STDERR "\t", join("\t", 'its end point ', $bsize{$istart}, 'is NOT further the previous  end point',  $osize[$#ostart] ), "\n"; }
      if ($verbose > 3){ print STDERR "\t", join("\t",'last dest block', $ostart[$#ostart] , $osize[$#osize]), "\n" }
    }
  }
}
## gather the data to be stdout-ed
my $ochr = $curr_cluster[0];
my $ostart = $ostart[0];
my $oend = $ostart[$#ostart] + $osize[$#osize];
my $oname = $curr_cluster[3];
my $oscore = $pe_count;
my $ostrand = $curr_cluster[4];
## bed12 cdsStart|cdsStart is used to mark the 3'extend of the cluster
my ($oostart, $ooend) = ($curr_cluster[4] eq '+')?($ostart, $curr_cluster[2]):($curr_cluster[1], $oend);
my $ocolor = '255,0,0';
my $obcount =  scalar @osize;
## correct the cagecluster block starts to be relative to the 1st cagecluster block
foreach (@ostart){ $_ = $_ - $ostart }

## make sure that cdsEnd is at most the end of the CAGEscan cluster (proper cdsStart has been dealt with
## when building the cluster)
$ooend = $oend if ($ooend > $oend);
$oostart = $ostart if ($oostart < $ostart);
print STDOUT join("\t", $ochr,  $ostart, $oend, $oname, $oscore, $ostrand, $oostart, $ooend, $ocolor, $obcount, join(",", @osize), join(",", @ostart)), "\n";

## clean up the tmp file if created
if ($verbose > 1){
  print STDERR "$0 WARNING: temp file $tmp_bed6_tss_input_file has not been deleted\n" if (-e $tmp_bed6_tss_input_file);
  print STDERR "$0 WARNING: temp file $tmp_bed6_cluster_file has not been deleted\n" if (-e $tmp_bed6_cluster_file);
}
else {
  system("rm $tmp_bed6_tss_input_file") if (-e $tmp_bed6_tss_input_file);
  system("rm $tmp_bed6_cluster_file") if (-e $tmp_bed6_cluster_file);
}
