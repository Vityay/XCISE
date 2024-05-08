#!/usr/bin/perl -w
use strict;

# Input parameters
# -b bam (can be multiple)
# -s snv
# -r chromosome
# -o output name
# -t number of tries
# -u min number of UMIs per SNV
# -m allowed MAF per SNV
# -p discordant penalty

warn "Reading parameters ...\n"; 
my ( $sample, $chromosome, $snv_file);
my $tries = 100;
my $min_maf = 0.25;
my $min_umis = 10;
my $discordant_penalty = 5;
my %bams = ();
my $ele = 0;
while ( $ele <= $#ARGV ) {
    if ( $ARGV[$ele] eq '-b' ) { #BAM files
        while ( $ele < $#ARGV and $ARGV[$ele+1] !~ m/^\-/) {
            die 'BAM file '.$ARGV[$ele+1].' does not exist' unless -e $ARGV[$ele+1];
            $bams{$ARGV[$ele+1]}=1;
            $ele++;
        }
        $ele++;
    }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-o' ) { #Output prefix
        $sample = $ARGV[$ele+1];
        $ele+=2;
    }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-s' ) { #SNV file
        $snv_file = $ARGV[$ele+1];
        $ele+=2;
    }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-r' ) { #Chromosome
        $chromosome = $ARGV[$ele+1];
        $ele+=2;
    }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-t' ) { #Number of tries
        $tries = $ARGV[$ele+1];
        $ele+=2;
    }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-u' ) { #Min number of UMIs per SNV
        $min_umis = $ARGV[$ele+1];
        $ele+=2;
    }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-m' ) { #Allowed MAF per SNV
        $min_maf = $ARGV[$ele+1];
        $ele+=2;
    }
    elsif ( $ele < $#ARGV and $ARGV[$ele] eq '-p' ) { #Discrodant penalty
        $discordant_penalty = $ARGV[$ele+1];
        $ele+=2;
    }
    else {
        die "Unexpected/incomplete parameter:$ARGV[$ele]";
    }
}
my $usage = 'Usage: perl '.$0.' -o <output_prefix> -s <vcf_file> -b <bam_file1> [bam_file2] [-r <chromosome>] [-t <tries>] [-u <min_umis>] [-m <min_maf>] [-p <penalty>]'."\n";
die "No Sample name\n".$usage unless $sample;
warn "Sample: ", $sample, "\n";
die "No VCF(SNV) file\n".$usage unless $snv_file;
warn "VCF file: ", $snv_file, "\n";
die "No BAM file\n".$usage unless keys %bams;
warn "BAM files: ", scalar(keys %bams),' (',join(',', keys %bams), ")\n";
$chromosome = 'X' unless $chromosome;
warn "Chromosome: ", $chromosome, "\n";
$tries = 100 if $tries <= 0;
warn "Number of tries: ", $tries, "\n";
$min_umis = 10 if $min_umis <= 0;
warn "Min number of UMIs per SNV: ", $min_umis, "\n";
$min_maf = 0.25 if $min_maf <= 0 or $min_maf >= 1;
warn "Min MAF per SNV: ", $min_maf, "\n";
$discordant_penalty = 5 if $discordant_penalty < 1;
warn "Discordant penalty: ", $discordant_penalty, "\n";

warn "Reading SNVs ...\n";
my %snp_info = ();
open F, $snv_file;
open F, 'gunzip -c '.$snv_file.' |' if $snv_file =~ m/\.gz$/;
while ( <F> ) {
    next if m/^\#/;
    chomp;
    my @arr = split /\t/;
    next unless $arr[0] eq $chromosome;
    $snp_info{$arr[1]}{'rs'} = $arr[2];
    $snp_info{$arr[1]}{'ref'} = $arr[3];
    $snp_info{$arr[1]}{'alt'} = $arr[4];
}
close F;
warn scalar( keys %snp_info), " SNVs were loaded from VCF file $snv_file\n";
die "Zero SNVs, run terminated" unless keys %snp_info;

warn "Reading BAM files ...\n";
my %umis = ();
my %cb2allele = ();
my $inf_reads = 0;
my $inf_alleles = 0;
my ( $hard_umis, $soft_umis ) = ( 0, 0 );
foreach my $file ( sort keys %bams ) {
    warn "    Reading $file ...\n";
    open F, 'samtools view -F 256 '.$file.' '.$chromosome.' |';
    while ( my $line = <F> ) {
        next unless $line =~ m/\tvG\:B\:i\,(\d+)/;
        my ( $vGs ) = $line =~ m/\tvG\:B\:i\,([\d+\,]+)/;
        my @vG = split /\,/, $vGs;
        my ( $vAs ) = $line =~ m/\tvA\:B\:c\,([\d+\,]+)/;
        my @vA = split /\,/, $vAs;
        my $cb;
        if ( $line =~ m/\tCB\:Z\:(\S+)/ or $line =~ m/\tRG\:Z\:(\S+)/ ) {
            $cb = $1;
            next if $cb eq '-';
        }
        else {
            die "Found no CB, no RG tags in BAM file in read $line";
        }
        my $umi;
        if ( $line =~ m/\tUB\:Z\:(\S+)/ ) {
            $umi = $1;
            next if $umi eq '-';
            $hard_umis++;
        }
        else {
            my ( $read, $flag, $chr, $pos, $mapq, $cigar, $chr2, $pos2, $tlen ) = split /\t/, $line;
            $umi = join('_', $cb, $chr, $pos, $flag, $cigar, $tlen );
            $soft_umis++;
            #warn $umi, "\n";
        }
        $inf_reads++;
        foreach my $ele ( 0 .. $#vG ) {
            my $pos = $vG[$ele]+1;
            next unless exists($snp_info{$pos});
            my $allele = $vA[$ele];
            next unless $allele == 1 or $allele == 2;
            $inf_alleles++;
            $umis{$pos}{$allele}{$umi."\t".$cb}++;
            $cb2allele{$cb}{$pos}{$umi}=$allele;
        }
    }
    close F;
    warn "    Reads with allelic info: $inf_reads Alleles: $inf_alleles UMIs:$hard_umis\/$soft_umis (hard/soft)\n";
}
warn scalar keys %umis, " SNV positions were covered in ", scalar keys %cb2allele, " cells were loaded from BAM files\n";
die "Zero SNVs, run terminated" unless keys %umis;
die "Zero cells, run terminated" unless keys %cb2allele;

warn "Checking AFs after WASP genotyping ...\n";
my @phased_pos = ();
my @phased_dir = ();
my $low_umis = 0;
my $low_maf = 0;
foreach my $pos ( sort {$a<=>$b} keys %umis ) {
    if ( scalar(keys %{$umis{$pos}{1}}) + scalar(keys %{$umis{$pos}{2}} ) < $min_umis ) {
        $low_umis++;
        next;
    }
    my $af = scalar(keys %{$umis{$pos}{1}}) / ( scalar(keys %{$umis{$pos}{1}}) + scalar(keys %{$umis{$pos}{2}} ) );
    if ( $af < $min_maf or $af > ( 1 - $min_maf ) ) {
        $low_maf++;
        next;
    }
    push @phased_pos, $pos;
    push @phased_dir, int(rand(3))-1;
}
die "No SNVs to phase\n" unless @phased_pos;
warn "Excluded $low_umis SNVs having less than $min_umis reads/UMIs\n" if $low_umis;
warn "Excluded $low_maf SNVs having less than $min_maf minor allele frequency\n" if $low_maf;

warn "Starting XCI calling\n";
my $global_best = q{};
foreach my $try ( 1 .. 100 ) {

    warn "Try #$try, randomizing order and XCI status ...\n";
    foreach my $ele ( 0 .. $#phased_pos ) {
        my $ele1 = int(rand(@phased_pos));
        my $ele2 = int(rand(@phased_pos));
        @phased_pos[$ele1,$ele2] = @phased_pos[$ele2,$ele1] if $ele1 != $ele2;
        $phased_dir[$ele] = int(rand(3))-1;
    }
    
    warn "    Calculating initial score ...\n";
    my ( $t, $d, $c, $hap2, $hap0, $hap1 ) = (0,0,0,0,0,0);
    my %cb2phase = ();
    foreach my $ele ( 0 .. $#phased_pos ) {
        foreach my $allele ( keys %{$umis{$phased_pos[$ele]}} ) {
            foreach my $bc ( keys %{$umis{$phased_pos[$ele]}{$allele}} ) {
                my ( $umi, $cb ) = split /\t/, $bc;
                $cb2phase{$cb}{$allele}++ if $phased_dir[$ele] == 1;
                $cb2phase{$cb}{3-$allele}++ if $phased_dir[$ele] == -1;
            }
        }
    }
    foreach my $cb ( keys %cb2phase ) {
        my $min = exists($cb2phase{$cb}{1}) ? $cb2phase{$cb}{1} : 0; 
        my $max = exists($cb2phase{$cb}{2}) ? $cb2phase{$cb}{2} : 0; 
        $t += $min + $max;
        ( $min, $max ) = ( $max, $min ) if $min > $max;
        next unless $max;
        $d += $min;
        $c += $max - 1;
    }
    my $score = $c - $discordant_penalty * $d;
    warn "    SNVs: ",scalar(@phased_pos)," Initial score: $score, discordance rate : ",100*$d/($d+$c),"\n";

    my $best_score = $score;
    my ( $pass, $imps, $last_imp ) = ( 0, 1, 0 );
    while ( $imps ) {
        $imps = 0;
        $pass++;
        foreach my $ele ( 0 .. $#phased_dir ) {
            foreach my $new ( -1, 0, 1 ) {
                next if $new == $phased_dir[$ele];
                my $old = $phased_dir[$ele];

                my ( $total, $disc, $conc ) = ( 0, 0, 0 );
                foreach my $cb ( keys %cb2phase ) {
                    my $min = exists($cb2phase{$cb}{1}) ? $cb2phase{$cb}{1} : 0; 
                    my $max = exists($cb2phase{$cb}{2}) ? $cb2phase{$cb}{2} : 0; 
                    $total += $min + $max;
                    ( $min, $max ) = ( $max, $min ) if $min > $max;
                    next unless $max;
                    $disc += $min;
                    $conc += $max - 1;
                }
                my $score_before = $conc - $discordant_penalty * $disc;
#                print 'Before change:', $score_before, "\n";
                die 'score_before ne best_score' if $score_before != $best_score;

                foreach my $allele ( keys %{$umis{$phased_pos[$ele]}} ) {
                    foreach my $bc ( keys %{$umis{$phased_pos[$ele]}{$allele}} ) {
                        my ( $umi, $cb ) = split /\t/, $bc;
                        $cb2phase{$cb}{$allele}-- if $old == 1;
                        $cb2phase{$cb}{$allele}++ if $new == 1;
                        $cb2phase{$cb}{3-$allele}-- if $old == -1;
                        $cb2phase{$cb}{3-$allele}++ if $new == -1;
                    }
                }

                ( $total, $disc, $conc ) = ( 0, 0, 0 );
                foreach my $cb ( keys %cb2phase ) {
                    my $min = exists($cb2phase{$cb}{1}) ? $cb2phase{$cb}{1} : 0; 
                    my $max = exists($cb2phase{$cb}{2}) ? $cb2phase{$cb}{2} : 0; 
                    $total += $min + $max;
                    ( $min, $max ) = ( $max, $min ) if $min > $max;
                    next unless $max;
                    $disc += $min;
                    $conc += $max - 1;
                }
                my $score_after = $conc - $discordant_penalty * $disc;
#                print 'After change:', $score_after, "\n";

                if ( $score_after > $best_score or ( $score_after == $best_score and $new == 0 ) ) {
                    $imps++;
                    $last_imp = $ele;
                    #print "$sample\@$chromosome Try: $try Pass $pass Pos: $ele Score: $best_score -> $score_after imps: $imps\n";
                    $phased_dir[$ele] = $new;
                    $best_score = $score_after;
                }
                else {
                    #Restore hash
                    foreach my $allele ( keys %{$umis{$phased_pos[$ele]}} ) {
                        foreach my $bc ( keys %{$umis{$phased_pos[$ele]}{$allele}} ) {
                            my ( $umi, $cb ) = split /\t/, $bc;
                            $cb2phase{$cb}{$allele}++ if $old == 1;
                            $cb2phase{$cb}{$allele}-- if $new == 1;
                            $cb2phase{$cb}{3-$allele}++ if $old == -1;
                            $cb2phase{$cb}{3-$allele}-- if $new == -1;
                        }
                    }
                    my ( $total_, $disc_, $conc_ ) = ( 0, 0, 0 );
                    foreach my $cb ( keys %cb2phase ) {
                        my $min = exists($cb2phase{$cb}{1}) ? $cb2phase{$cb}{1} : 0; 
                        my $max = exists($cb2phase{$cb}{2}) ? $cb2phase{$cb}{2} : 0; 
                        $total_ += $min + $max;
                        ( $min, $max ) = ( $max, $min ) if $min > $max;
                        next unless $max;
                        $disc_ += $min;
                        $conc_ += $max - 1;
                    }
                    my $score_restored = $conc_ - $discordant_penalty * $disc_;
#                    print 'After restore:', $score_restored, "\n";
                    die 'score_restored ne score_before' if $score_before != $score_restored;
                } # if best score
            } # foreach -1,0,1 
            last if $imps == 0 and $last_imp - $ele == 1; #No imps so far and the ext element is where we made improvement in the last pass
        } # foreach ele
        warn "    Pass: $pass Improvements: $imps Score: $best_score  Global best: $global_best\n";
    } # while imps
    warn "    Finished at Pass: $pass, SNV: $last_imp/",scalar(@phased_dir), " Final score: $best_score Global best: $global_best\n"; 
    next if $global_best and $global_best >= $best_score;
    $global_best = $best_score;

    my ( $total0, $total1, $total2, $total3, $total4 );
    do {
        my ( $ref1, $alt1, $noninf ) = ( 0, 0, 0 );
        foreach my $ele ( 0 .. $#phased_pos ) {
            if ( $phased_dir[$ele] == 1 ) {
                $ref1++;
            }
            elsif ( $phased_dir[$ele] == -1 ) {
                $alt1++;
            }
            else {
                $noninf++;
            }
        }
        my ( $total, $disc, $conc ) = ( 0, 0, 0 );
        foreach my $cb ( keys %cb2phase ) {
            my $min = exists($cb2phase{$cb}{1}) ? $cb2phase{$cb}{1} : 0; 
            my $max = exists($cb2phase{$cb}{2}) ? $cb2phase{$cb}{2} : 0; 
            $total += $min + $max;
            ( $min, $max ) = ( $max, $min ) if $min > $max;
            next unless $max;
            $disc += $min;
            $conc += $max - 1;
        }

        ( $total0, $total1, $total2, $total3, $total4 ) = ( 0, 0, 0, 0, 0 );
        open FLOG, '>', $sample.'_chr'.$chromosome.'_XCISE_summary.txt';
        print FLOG 'Best score  : ', $best_score, "\n";
        print FLOG 'Total UMIs  : ', $total, "\n";
        print FLOG 'Concordant  : ', $conc, "\n";
        print FLOG 'Discordant  : ', $disc, "\n";
        print FLOG 'Discordancy : ', $disc/($conc+$disc), "\n";
        print FLOG 'Total SNVs  : ', scalar(@phased_dir), "\n";
        print FLOG 'SNV Ref/Alt : ', $ref1, "\n";
        print FLOG 'SNV Alt/Ref : ', $alt1, "\n";
        print FLOG 'Non-inf SNVs: ', $noninf, "\n";
        print FLOG 'XCI_inf SNVs: ', $ref1+$alt1, "\n";

        warn '    Total XCI-informative UMIs   : ', $total, "\n";
        warn '    SNVs XCI-informative/non-inf : ', $ref1+$alt1,'/',$noninf, "\n";
        warn '    Final discordancy rate       : ', 100*$disc/($conc+$disc), "\n";

        warn "    Outputting VCF...\n";
        open FVCF, '>', $sample.'_chr'.$chromosome.'_XCISE.vcf';
        print FVCF "##fileformat=VCFv4.2\n";
        print FVCF join( "\t", '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'INFO' ), "\n";
        foreach my $ele ( sort {$phased_pos[$a]<=>$phased_pos[$b]} ( 0 .. $#phased_dir ) ) {
            my $pos = $phased_pos[$ele];
            print FVCF join( "\t", $chromosome, $pos, $snp_info{$pos}{'rs'}, $snp_info{$pos}{'ref'}, $snp_info{$pos}{'alt'}, scalar(keys %{$umis{$pos}{1}})+scalar(keys %{$umis{$pos}{2}}), 'PASS', 'DIR='.$phased_dir[$ele].';HAP1='.scalar(keys %{$umis{$pos}{1}}).';HAP2='.scalar(keys %{$umis{$pos}{2}}) ),"\n";  
        }
        close FVCF;

        warn "    Outputting barcode to phase data...\n";
        open F, '>', $sample.'_chr'.$chromosome.'_XCISE_barcode_to_phase.txt';
        foreach my $cb ( sort keys %cb2allele ) {
            my ( $hap1, $hap2 ) = ( 0, 0 );
            foreach my $ele ( 0 .. $#phased_dir ) {
                my $pos = $phased_pos[$ele];
                next unless exists( $cb2allele{$cb}{$pos} );
                foreach my $umi ( keys %{$cb2allele{$cb}{$pos}} ) {
                    $hap1++ if $cb2allele{$cb}{$pos}{$umi} == 1 and $phased_dir[$ele] == 1;
                    $hap1++ if $cb2allele{$cb}{$pos}{$umi} == 2 and $phased_dir[$ele] == -1;
                    $hap2++ if $cb2allele{$cb}{$pos}{$umi} == 1 and $phased_dir[$ele] == -1;
                    $hap2++ if $cb2allele{$cb}{$pos}{$umi} == 2 and $phased_dir[$ele] == 1;
                }
            }
            my $hap = '?';
            if ( $hap1 == 0 and $hap2 == 0 ) {
                $total0++;
                $hap = 'UNKN';
            }
            elsif ( $hap1 >= 2 and $hap1/($hap1+$hap2) >= 0.9 ) {
                $total1++;
                $hap = 'HAP1';
            } 
            elsif ( $hap2 >= 2 and $hap2/($hap1+$hap2) >= 0.9 ) {
                $total2++;
                $hap = 'HAP2';
            } 
            elsif ( $hap1 > 0 and $hap2>0 ) {
                $total3++;
                $hap = 'BOTH';
            } 
            else {
                $total4++;
                $hap = 'LOWC';
            }
            print F join( "\t", $cb, $hap1, $hap2, $hap ), "\n";
        }
        close F;
        my $grand_total = $total0+$total1+$total2+$total3+$total4;

        warn "    Outputing summary...\n";
        printf FLOG "Haplotype 1 : %5d( %3.2f %% )\n", $total1, 100*$total1/$grand_total;
        printf FLOG "Haplotype 2 : %5d( %3.2f %% )\n", $total2, 100*$total2/$grand_total;
        printf FLOG "Unknown hap : %5d( %3.2f %% )\n", $total0, 100*$total0/$grand_total;
        printf FLOG "Both H seen : %5d( %3.2f %% )\n", $total3, 100*$total3/$grand_total;
        printf FLOG "LowCoverage : %5d( %3.2f %% )\n", $total4, 100*$total4/$grand_total;
        close FLOG;
        warn "    X1/X2/Both/LowC/Unknown cells: ", join( ' / ', $total1, $total2, $total3, $total4, $total0 ), "\n";
        warn "    Switching X1 and X2, so that there are more X2 cells than X1 cells ...\n" if $total1 > $total2;
        foreach my $ele ( 0 .. $#phased_dir ) { $phased_dir[$ele] = -$phased_dir[$ele] }
    } until ( $total1 <= $total2 );
}


=head
my $sample = shift or die 'No sample';
my $ctg = shift;
$ctg = 'X' unless $ctg;
my $min_af = 0.01; #0 Will skip SNVs that are rare in GNOMAD
my $min_maf = 0.25; #0.25 Will skip SNVs that do not loom heterozygous
my $max_discordance = 0.1; #0.1 Will skip candidate SNVs for first SNV that show high rate of biallelic calls in individual cells 
my $min_switch_disc_ratio_diff = 0.01; #0.005 Switching of alleles should make difference
my $min_umi = 10; # SNP should be covered by enough UMIs to be considered
my $sc_mode = '10xreg'; # can be 10x, 10xreg or plate
my $discordant_penalty = 5;

my %known = ();
warn "Preloading known cells\n";
open F,'01_bc2celltype.txt';
<F>;
while ( <F> ) {
    chomp;
    my ( $sample2, $region2, $celltype2, $bc2 ) = split /\t/;
    $known{$region2.'_'.$bc2} = $celltype2 if $sample2 eq $sample;
}
close F;
warn "Total known cells:", scalar keys %known, "\n";

warn "Counting UMIs per SNV\n";
my %umis = ();
my %cb2allele = ();
my $i = 0;
my $j = 0;
foreach my $file ( <$sample/$sample*Aligned.sortedByCoord.out.bam> ) {
    warn "$file\n";
    my $lib = $file;
    $lib =~ s/Aligned\.sortedByCoord\.out\.bam$//;
    $lib =~ s/_WASP$//;
    $lib =~ s/\_R1\_all$//;
    $lib =~ s/^.*\///g;
    $lib =~ s/^plate//;
    
    warn "Parsing lib $lib\n";
    open F, 'samtools view -F 256 '.$file.' '.$ctg.' |';
    while ( <F> ) {
        if ( $sc_mode eq '10x' ) {
            next unless m/\tvG\:B\:i\,(\d+)/;
            my $pos = $1+1;
            next unless m/\tvA\:B\:c\,(\d)/;
            my $allele = $1;
            next unless m/\tCB\:Z\:(\w+)/;
            $j++;
            my $cb = $1;
            next unless m/\tUB\:Z\:(\w+)/;
            my $umi = $1;
            $umis{$pos}{$allele}{$umi."\t".$cb}{$lib}++;
            $cb2allele{$cb}{$pos}{$umi}=$allele;
            $i++;
        }
        elsif ( $sc_mode eq '10xreg' ) {
            my ( $none, $reg ) = $lib =~ m/^(IB|CTR|FTD)\d+\_([^\_]+)/;
            next unless m/\tvG\:B\:i\,(\d+)/;
            my $pos = $1+1;
            next unless m/\tvA\:B\:c\,(\d)/;
            my $allele = $1;
            next unless m/\tCB\:Z\:(\w+)/;
            $j++;
            my $cb = $1;
            next unless exists( $known{$reg.'_'.$cb} );
            next unless m/\tUB\:Z\:(\w+)/;
            my $umi = $1;
            $umis{$pos}{$allele}{$umi."\t".$reg.'_'.$cb}{$lib}++;
            $cb2allele{$reg.'_'.$cb}{$pos}{$umi}=$allele;
            $i++;
        }
        elsif ( $sc_mode eq 'plate' ) {
            next unless m/\tvG\:B\:i\,(\d+)/;
            my $pos = $1+1;
            next unless m/\tvA\:B\:c\,(\d)/;
            my $allele = $1;
            $j++;
            my ( $cb ) = m/RG\:Z\:(\S+)/;
            $cb = $lib unless $cb;
            $cb =~ s/_R1_all_WASP//;
            $cb =~ s/^plate//;
            die 'No read group in SAM string: '.$_ unless $cb;
            my @arr = split /\t/;
            my $umi = join( ':', @arr[3,1,5] ); 
            $umis{$pos}{$allele}{$umi."\t".$cb}{$lib}++;
            $cb2allele{$cb}{$pos}{$umi}=$allele;
            $i++;
        }
        else {
            die "Unsuppoted single-cell mode: $sc_mode";
        }
    }
    close F;
    warn "Lines: $j Alleles: $i\n";
}
warn "Sorting cell barcodes\n";
my @cbs = sort keys %cb2allele;

warn "Counting UMI per allele for SNVs\n";
my @phased_pos = ();
my @phased_dir = ();
my %poss = ();
my %discordance_per_pos = ();
my %cb2phase = ();
open F, '>', '05b_chr'.$ctg.'_SNVs_concordance_'.$sample.'.txt';
foreach my $pos ( sort {scalar(keys %{$umis{$b}{1}})+scalar(keys %{$umis{$b}{2}}) <=> scalar(keys %{$umis{$a}{1}})+scalar(keys %{$umis{$a}{2}})} keys %umis ) {
    next unless exists( $afs{$pos} );
    $poss{$pos}{1} = scalar keys %{$umis{$pos}{1}};
    $poss{$pos}{2} = scalar keys %{$umis{$pos}{2}};
    next unless $poss{$pos}{1}+$poss{$pos}{2}; #Skip SNVs that do not have support
    next if $poss{$pos}{1} / ($poss{$pos}{1}+$poss{$pos}{2}) < $min_maf or $poss{$pos}{1} / ($poss{$pos}{1}+$poss{$pos}{2}) > ( 1-$min_maf );
    @phased_pos = ( $pos );
    @phased_dir = ( 1 );
    %cb2phase = ();
    foreach my $ele ( 0 .. $#phased_pos ) {
        foreach my $allele ( keys %{$umis{$phased_pos[$ele]}} ) {
            foreach my $bc ( keys %{$umis{$phased_pos[$ele]}{$allele}} ) {
                my ( $umi, $cb ) = split /\t/, $bc;
                $cb2phase{$cb}{$allele}++ if $phased_dir[$ele] == 1;
                $cb2phase{$cb}{3-$allele}++ if $phased_dir[$ele] == -1;
            }
        }
    }
    my ( $t0, $d0, $c0 ) = ( 0, 0, 0 );
    foreach my $cb ( keys %cb2phase ) {
        my $min = exists($cb2phase{$cb}{1}) ? $cb2phase{$cb}{1} : 0; 
        my $max = exists($cb2phase{$cb}{2}) ? $cb2phase{$cb}{2} : 0; 
        $t0 += $min + $max;
        ( $min, $max ) = ( $max, $min ) if $min > $max;
        next unless $max;
        $d0 += $min;
        $c0 += $max - 1;
    }
    $discordance_per_pos{$pos}{disc} = $d0;
    $discordance_per_pos{$pos}{conc} = $c0;
    my $discordance = ($d0+$c0) ? $d0/($d0+$c0) : '?';
}
close F;


my $global_best;
if ( -s '05b_chr'.$ctg.'_haplotype_inf_SNVs_'.$sample.'.log' ) {
     open F, '05b_chr'.$ctg.'_haplotype_inf_SNVs_'.$sample.'.log';
     my $line = <F>;
     $global_best = $1 if $line =~ m/^Best\s+score\s+\:\s+(\d+)/;
     close F;
     warn "Existsing global best:$global_best\n" if $global_best;
}
foreach my $try ( 1 .. 100 ) {

warn "Initiating phasing\n";
@phased_pos = ();
@phased_dir = ();
foreach my $pos ( sort {$a<=>$b} keys %umis ) {
    next unless exists( $afs{$pos} );
    next unless exists($poss{$pos}) and ( $poss{$pos}{1}+$poss{$pos}{2} ) >= $min_umi;
    next if $poss{$pos}{1} / ($poss{$pos}{1}+$poss{$pos}{2}) < $min_maf or $poss{$pos}{1} / ($poss{$pos}{1}+$poss{$pos}{2}) > ( 1-$min_maf );
    my $discordance_rate = ( $discordance_per_pos{$pos}{disc} + $discordance_per_pos{$pos}{conc}) ? $discordance_per_pos{$pos}{disc} / ( $discordance_per_pos{$pos}{disc} + $discordance_per_pos{$pos}{conc}) : '?';
    next unless $discordance_rate eq '?' or $discordance_rate <= $max_discordance;
    push @phased_pos, $pos;
    push @phased_dir, int(rand(3))-1;
}
die "No SNVs to phase\n" unless @phased_pos;
warn "Randomizing order\n";
foreach ( 0 .. scalar(@phased_pos) ) {
    my $ele1 = int(rand(@phased_pos));
    my $ele2 = int(rand(@phased_pos));
    @phased_pos[$ele1,$ele2] = @phased_pos[$ele2,$ele1] if $ele1 != $ele2;
}

warn "Calculating score\n";
my ( $t, $d, $c, $hap2, $hap0, $hap1 ) = (0,0,0,0,0,0);
%cb2phase = ();
foreach my $ele ( 0 .. $#phased_pos ) {
    if ( $phased_dir[$ele] == 1 ) {
        $hap1++;
    }
    elsif ( $phased_dir[$ele] == -1 ) {
        $hap2++;
    }
    else {
        $hap0++;
    }
    foreach my $allele ( keys %{$umis{$phased_pos[$ele]}} ) {
        foreach my $bc ( keys %{$umis{$phased_pos[$ele]}{$allele}} ) {
            my ( $umi, $cb ) = split /\t/, $bc;
            $cb2phase{$cb}{$allele}++ if $phased_dir[$ele] == 1;
            $cb2phase{$cb}{3-$allele}++ if $phased_dir[$ele] == -1;
        }
    }
}
foreach my $cb ( keys %cb2phase ) {
    my $min = exists($cb2phase{$cb}{1}) ? $cb2phase{$cb}{1} : 0; 
    my $max = exists($cb2phase{$cb}{2}) ? $cb2phase{$cb}{2} : 0; 
    $t += $min + $max;
    ( $min, $max ) = ( $max, $min ) if $min > $max;
    next unless $max;
    $d += $min;
    $c += $max - 1;
}
my $score = $c - $discordant_penalty * $d;
warn "Total: $t  Concordant: $c  Discordant: $d  Score: $score\n";

my $best_score = $score;
my ( $pass, $imps ) = ( 0, 1 );
while ( $imps ) {
    $imps = 0;
    $pass++;
#    open FLAST, '>', '05b_chr'.$ctg.'_last_round_scores_SNVs_'.$sample.'.log';
    foreach my $ele ( 0 .. $#phased_dir ) {
        foreach my $new ( -1, 0, 1 ) {
            next if $new == $phased_dir[$ele];
            my $old = $phased_dir[$ele];

            my ( $total, $disc, $conc ) = ( 0, 0, 0 );
            foreach my $cb ( keys %cb2phase ) {
                my $min = exists($cb2phase{$cb}{1}) ? $cb2phase{$cb}{1} : 0; 
                my $max = exists($cb2phase{$cb}{2}) ? $cb2phase{$cb}{2} : 0; 
                $total += $min + $max;
                ( $min, $max ) = ( $max, $min ) if $min > $max;
                next unless $max;
                $disc += $min;
                $conc += $max - 1;
            }
            my $score_before = $conc - $discordant_penalty * $disc;
#            print 'Before change:', $score_before, "\n";
            die 'score_before ne best_score' if $score_before != $best_score;

            foreach my $allele ( keys %{$umis{$phased_pos[$ele]}} ) {
                foreach my $bc ( keys %{$umis{$phased_pos[$ele]}{$allele}} ) {
                    my ( $umi, $cb ) = split /\t/, $bc;
                    $cb2phase{$cb}{$allele}-- if $old == 1;
                    $cb2phase{$cb}{$allele}++ if $new == 1;
                    $cb2phase{$cb}{3-$allele}-- if $old == -1;
                    $cb2phase{$cb}{3-$allele}++ if $new == -1;
                }
            }

            ( $total, $disc, $conc ) = ( 0, 0, 0 );
            foreach my $cb ( keys %cb2phase ) {
                my $min = exists($cb2phase{$cb}{1}) ? $cb2phase{$cb}{1} : 0; 
                my $max = exists($cb2phase{$cb}{2}) ? $cb2phase{$cb}{2} : 0; 
                $total += $min + $max;
                ( $min, $max ) = ( $max, $min ) if $min > $max;
                next unless $max;
                $disc += $min;
                $conc += $max - 1;
            }
            my $score_after = $conc - $discordant_penalty * $disc;
#            print 'After change:', $score_after, "\n";

            if ( $score_after > $best_score or ( $score_after == $best_score and $new == 0 ) ) {
                print "$sample\@$ctg Try: $try Pass $pass Pos: $ele Score: $best_score -> $score_after imps: ", ++$imps, "\n";
                $phased_dir[$ele] = $new;
                $best_score = $score_after;
            }
            else {
                #Restore hash
                foreach my $allele ( keys %{$umis{$phased_pos[$ele]}} ) {
                    foreach my $bc ( keys %{$umis{$phased_pos[$ele]}{$allele}} ) {
                        my ( $umi, $cb ) = split /\t/, $bc;
                        $cb2phase{$cb}{$allele}++ if $old == 1;
                        $cb2phase{$cb}{$allele}-- if $new == 1;
                        $cb2phase{$cb}{3-$allele}++ if $old == -1;
                        $cb2phase{$cb}{3-$allele}-- if $new == -1;
                    }
                }
                my ( $total_, $disc_, $conc_ ) = ( 0, 0, 0 );
                foreach my $cb ( keys %cb2phase ) {
                    my $min = exists($cb2phase{$cb}{1}) ? $cb2phase{$cb}{1} : 0; 
                    my $max = exists($cb2phase{$cb}{2}) ? $cb2phase{$cb}{2} : 0; 
                    $total_ += $min + $max;
                    ( $min, $max ) = ( $max, $min ) if $min > $max;
                    next unless $max;
                    $disc_ += $min;
                    $conc_ += $max - 1;
                }
                my $score_restored = $conc_ - $discordant_penalty * $disc_;
#                print 'After restore:', $score_restored, "\n";
                die 'score_restored ne score_before' if $score_before != $score_restored;
            }
        }
    }
#    close FLAST;
    warn "Chr: $ctg Pass: $pass  Improvements: $imps Global best: $global_best\n";
}

next if defined($global_best) and $global_best >= $best_score;
$global_best = $best_score;

my ( $fwd, $rev, $mix ) = ( 0, 0, 0 );
foreach my $ele ( 0 .. $#phased_pos ) {
    if ( $phased_dir[$ele] == 1 ) {
        $fwd++;
    }
    elsif ( $phased_dir[$ele] == -1 ) {
        $rev++;
    }
    else {
        $mix++;
    }
}
my ( $total, $disc, $conc ) = ( 0, 0, 0 );
foreach my $cb ( keys %cb2phase ) {
    my $min = exists($cb2phase{$cb}{1}) ? $cb2phase{$cb}{1} : 0; 
    my $max = exists($cb2phase{$cb}{2}) ? $cb2phase{$cb}{2} : 0; 
    $total += $min + $max;
    ( $min, $max ) = ( $max, $min ) if $min > $max;
    next unless $max;
    $disc += $min;
    $conc += $max - 1;
}
open FLOG, '>', '05b_chr'.$ctg.'_haplotype_inf_SNVs_'.$sample.'.log';
print FLOG 'Best score  : ', $best_score, "\n";
print FLOG 'Total UMIs  : ', $total, "\n";
print FLOG 'Concordant  : ', $conc, "\n";
print FLOG 'Discordant  : ', $disc, "\n";
print FLOG 'Discordancy : ', $disc/($conc+$disc), "\n";
print FLOG 'Total SNVs  : ', scalar(@phased_dir), "\n";
print FLOG 'Forward SNVs: ', $fwd, "\n";
print FLOG 'Reverse SNVs: ', $rev, "\n";
print FLOG 'Mixed SNVs  : ', $mix, "\n";
print FLOG 'Useful SNVs : ', $fwd+$rev, "\n";

warn "Outputting VCF...\n";
open FVCF, '>', '05b_chr'.$ctg.'_haplotype_inf_SNVs_'.$sample.'.vcf';
print FVCF "##fileformat=VCFv4.2\n";
foreach my $ele ( 0 .. $#phased_dir ) {
    my $pos = $phased_pos[$ele];
    print FVCF join( "\t", $ctg, $pos, $snp_info{$pos}{'rs'}, $snp_info{$pos}{'ref'}, $snp_info{$pos}{'alt'}, $poss{$pos}{1}+$poss{$pos}{2}, 'PASS', 'DIR='.$phased_dir[$ele].';HAP1='.$poss{$pos}{1}.';HAP2='.$poss{$pos}{2} ),"\n";  
}
close FVCF;

warn "Outputting barcode to phase\n";
my ( $total0, $total1, $total2, $total3, $total4 ) = ( 0, 0, 0, 0, 0 );
open F, '>', '05b_chr'.$ctg.'_barcode_to_phase_'.$sample.'.txt';
#foreach my $cb ( sort keys %cb2allele ) {
foreach my $cb ( @cbs ) {
    my ( $hap1, $hap2 ) = ( 0, 0 );
    foreach my $ele ( 0 .. $#phased_dir ) {
        my $pos = $phased_pos[$ele];
        next unless exists( $cb2allele{$cb}{$pos} );
        foreach my $umi ( keys %{$cb2allele{$cb}{$pos}} ) {
            $hap1++ if $cb2allele{$cb}{$pos}{$umi} == 1 and $phased_dir[$ele] == 1;
            $hap1++ if $cb2allele{$cb}{$pos}{$umi} == 2 and $phased_dir[$ele] == -1;
            $hap2++ if $cb2allele{$cb}{$pos}{$umi} == 1 and $phased_dir[$ele] == -1;
            $hap2++ if $cb2allele{$cb}{$pos}{$umi} == 2 and $phased_dir[$ele] == 1;
        }
    }
    my $hap = '?';
    if ( $hap1 == 0 and $hap2 == 0 ) {
        $total0++;
        $hap = 'UNKN';
    }
    elsif ( $hap1 >= 2 and $hap1/($hap1+$hap2) >= 0.9 ) {
        $total1++;
        $hap = 'HAP1';
    } 
    elsif ( $hap2 >= 2 and $hap2/($hap1+$hap2) >= 0.9 ) {
        $total2++;
        $hap = 'HAP2';
    } 
    elsif ( $hap1 > 0 and $hap2>0 ) {
        $total3++;
        $hap = 'BOTH';
    } 
    else {
        $total4++;
        $hap = 'LOWC';
    }
    print F join( "\t", $cb, $hap1, $hap2, $hap ), "\n";
}
close F;

warn "Finishing log\n";
print FLOG 'Haplotype 1: ', $total1, ' (', 0.01 * int(0.5+10000*$total1/($total0+$total1+$total2+$total3+$total4)), "\% )\n";
print FLOG 'Haplotype 2: ', $total2, ' (', 0.01 * int(0.5+10000*$total2/($total0+$total1+$total2+$total3+$total4)), "\% )\n";
print FLOG 'Unknown hap: ', $total0, ' (', 0.01 * int(0.5+10000*$total0/($total0+$total1+$total2+$total3+$total4)), "\% )\n";
print FLOG 'Both H seen: ', $total3, ' (', 0.01 * int(0.5+10000*$total3/($total0+$total1+$total2+$total3+$total4)), "\% )\n";
print FLOG 'LowCoverage: ', $total4, ' (', 0.01 * int(0.5+10000*$total4/($total0+$total1+$total2+$total3+$total4)), "\% )\n";
close FLOG;
warn "Done\n";
} # 1..100

