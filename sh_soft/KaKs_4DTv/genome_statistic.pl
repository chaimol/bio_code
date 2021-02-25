#!/usr/bin/perl
use strict;
use File::Basename;
my $usage = "Usage:\n\tperl genome_statistic.pl genomeFasta\n\n";
if (@ARGV==0) {die $usage}

foreach ( @ARGV ) {
	my $name = basename( $_ );
	print "\n==>>$name<<==\n\n";
	my ($genome_size, $contig_size, $identifier, %scaffold_length, %scaffold_seq, $N_num, $G_num, $C_num, $T_num, $A_num);
	open GENOME, '<',"$_" or die "Can't open file $_!($!)";
	while ( <GENOME> ) {
		chomp;
		if ( m/^>(\S+)/ ) { $identifier = $1 }
		else {
			my $length = length $_;
			$genome_size += $length;
			$scaffold_length{$identifier} += $length;
			$scaffold_seq{$identifier} .= $_;
			$N_num += ($_ =~ tr/Nn/Nn/);
			$G_num += ($_ =~ tr/Gg/Gg/);
			$C_num += ($_ =~ tr/Cc/Cc/);
			$T_num += ($_ =~ tr/Tt/Tt/);
			$A_num += ($_ =~ tr/Aa/Aa/);
		}
	}
	close GENOME;
	
	my @contig_length;
	foreach (keys %scaffold_seq) {
		my $scaffold_seq = $scaffold_seq{$_};
		#my @contigs = split /N{10,}/, $scaffold_seq;
		my @contigs = split /[Nn]+/, $scaffold_seq;
		foreach (@contigs) {
			my $length = length $_;
			$contig_size += $length;
			push @contig_length, $length
		}
	}
	@contig_length = sort { $b <=> $a } @contig_length;
        my $longest_contig = $contig_length[0];
        my $small_contig = $contig_length[$#contig_length];
	my @scaffold_length = values %scaffold_length;
	@scaffold_length = sort { $b <=> $a } @scaffold_length;
	
	my $longest_fragment = $scaffold_length[0];
	my $shortest_fragment = $scaffold_length[$#scaffold_length];
	my $sequence_number = @scaffold_length;
	my $contig_number = @contig_length;
	print "the genome scaffolds number is $sequence_number\nthe genome contigs number is $contig_number\nthe longest contig is $longest_contig\nthe smallest contig is$small_contig\nthe longest scaffold length is $longest_fragment\nthe shortest scaffold length is $shortest_fragment\n";

	my $rate_of_N = $N_num / $genome_size;
	my $rate_of_GC =  ( $G_num + $C_num ) / ($G_num + $C_num + $T_num + $A_num);
	print "the genome scaffolds size is $genome_size\nthe genome contig size is $contig_size\n";
	print "the rate of N is $rate_of_N\n";
	print "the rate of GC is $rate_of_GC\n";
        my $genome_size10= $genome_size * 0.1; my $contig_size10 = $contig_size * 0.1;
	my $genome_size20 = $genome_size * 0.2; my $contig_size20 = $contig_size * 0.2;
        my $genome_size30 = $genome_size * 0.3; my $contig_size30 = $contig_size * 0.3;
	my $genome_size40 = $genome_size * 0.4; my $contig_size40 = $contig_size * 0.4;
	my $genome_size50= $genome_size * 0.5; my $contig_size50 = $contig_size * 0.5;
        my $genome_size60 = $genome_size * 0.6; my $contig_size60 = $contig_size * 0.6;
        my $genome_size70 = $genome_size * 0.7; my $contig_size70 = $contig_size * 0.7;
        my $genome_size80 = $genome_size * 0.8; my $contig_size80 = $contig_size * 0.8;
	my $genome_size90 = $genome_size * 0.9; my $contig_size90 = $contig_size * 0.9;
	my $aclength = 0;
        my $sca=0;
        my $aclength=0;
	my %check;
	foreach (@scaffold_length) {
		$aclength += $_;
		$sca++;	
		if ($aclength >= $genome_size10 and $aclength <= $genome_size20) {
			if(!$check{$genome_size10}){print "the scaffold N10 is $_\t$sca\n";$check{$genome_size10}=1;}
		}elsif ($aclength >= $genome_size20 and $aclength <= $genome_size30) {
                        if(!$check{$genome_size20}){print "the scaffold N20 is $_\t$sca\n";
			$check{$genome_size20}=1;}
                }elsif ($aclength >= $genome_size30 and $aclength <= $genome_size40) {
                        if(!$check{$genome_size30}){print "the scaffold N30 is $_\t$sca\n";
			$check{$genome_size30}=1;}	
                }elsif ($aclength >= $genome_size40 and $aclength <= $genome_size50) {
                        if(!$check{$genome_size40}){print "the scaffold N40 is $_\t$sca\n";
			$check{$genome_size40}=1;}
                }elsif ($aclength >= $genome_size50 and $aclength <= $genome_size60) {
                        if(!$check{$genome_size50}){print "the scaffold N50 is $_\t$sca\n";
			$check{$genome_size50}=1;}
                }elsif ($aclength >= $genome_size60 and $aclength <= $genome_size70) {
                       if(!$check{$genome_size60}){ print "the scaffold N60 is $_\t$sca\n";
			$check{$genome_size60}=1;}
                }elsif ($aclength >= $genome_size70 and $aclength <= $genome_size80) {
                       if(!$check{$genome_size70}){ print "the scaffold N70 is $_\t$sca\n";
			$check{$genome_size70}=1;}
                }elsif ($aclength >= $genome_size80 and $aclength <= $genome_size90) {
                       if(!$check{$genome_size80}){ print "the scaffold N80 is $_\t$sca\n";
			$check{$genome_size80}=1;}
                }elsif ($aclength >= $genome_size90) {
                       if(!$check{$genome_size90}){ print "the scaffold N90 is $_\t$sca\n";
			$check{$genome_size90}=1;}
                }
	}
	$sca = 0;
	$aclength = 0;
	foreach (@contig_length) {
		$aclength += $_;
		$sca++;
	if ($aclength >= $contig_size10 and $aclength <= $contig_size20) {
                        if(!$check{$contig_size10}){print "the contig N10 is $_\t$sca\n";$check{$contig_size10}=1;}
                }elsif ($aclength >= $contig_size20 and $aclength <= $contig_size30) {
                        if(!$check{$contig_size20}){print "the contig N20 is $_\t$sca\n";
                        $check{$contig_size20}=1;}
                }elsif ($aclength >= $contig_size30 and $aclength <= $contig_size40) {
                        if(!$check{$contig_size30}){print "the contig N30 is $_\t$sca\n";
                        $check{$contig_size30}=1;}
                }elsif ($aclength >= $contig_size40 and $aclength <= $contig_size50) {
                        if(!$check{$contig_size40}){print "the contig N40 is $_\t$sca\n";
                        $check{$contig_size40}=1;}
                }elsif ($aclength >= $contig_size50 and $aclength <= $contig_size60) {
                        if(!$check{$contig_size50}){print "the contig N50 is $_\t$sca\n";
                        $check{$contig_size50}=1;}
                }elsif ($aclength >= $contig_size60 and $aclength <= $contig_size70) {
                       if(!$check{$contig_size60}){ print "the contig N60 is $_\t$sca\n";
                        $check{$contig_size60}=1;}
                }elsif ($aclength >= $contig_size70 and $aclength <= $contig_size80) {
                       if(!$check{$contig_size70}){ print "the contig N70 is $_\t$sca\n";
                        $check{$contig_size70}=1;}
                }elsif ($aclength >= $contig_size80 and $aclength <= $contig_size90) {
                       if(!$check{$contig_size80}){ print "the contig N80 is $_\t$sca\n";
                        $check{$contig_size80}=1;}
                }elsif ($aclength >= $contig_size90) {
                       if(!$check{$contig_size90}){ print "the contig N90 is $_\t$sca\n";
                        $check{$contig_size90}=1;}
                }

	}
	my ($length1000Num,$length2000Num,$length3000Num,$total_length3000,$total_length2000_3000,$total_length1000_2000,$total_length2000,$total_length1000,$length2000_3000Num,$length1000_2000Num);
	foreach (@scaffold_length) {
		if ($_ >= 3000) {
			$length3000Num ++;
			$total_length3000 += $_
		}elsif ($_ >= 2000) {
			$length2000_3000Num ++;
			$total_length2000_3000 += $_
		}elsif ($_ >= 1000) {
			$length1000_2000Num ++;
			$total_length1000_2000 += $_
		}
		$length2000Num = $length3000Num + $length2000_3000Num;
		$length1000Num = $length2000Num + $length1000_2000Num;
		$total_length2000 = $total_length3000 + $total_length2000_3000;
		$total_length1000 = $total_length2000 + $total_length1000_2000;
	}
	print "the number of sequences >= 1kb is $length1000Num\ttotal length is $total_length1000\n";
	print "the number of sequences >= 2kb is $length2000Num\ttotal length is $total_length2000\n";
	print "the number of sequences >= 3kb is $length3000Num\ttotal length is $total_length3000\n";
	print "\n";
}
