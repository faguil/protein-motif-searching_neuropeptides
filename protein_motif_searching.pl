#!/usr/bin/env perl
#
## Code by: Felipe Aguilera
## April 2017; Modified March 2018
## Script that search for user-defined motif protein sequences in fasta files (protein sequences)
#
## Assummes pfscan and psa2msa software are located system-wide (e.g., /usr/local/bin) 
## Also, assumes all additional files (i.e., prosite.dat, Prosite.pm) are in the same folder

#####################################
##
## The MIT License
##
## Copyright (c) 2017 Felipe Aguilera
##
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in
## all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
## THE SOFTWARE.
##
######################################
#

use 5.005_03; # uses perl >= v5.5
use IO::File;
use Carp qw(confess cluck);
use vars qw(@ISA $VERSION $errpos $errstr);
use strict;


################################################################################
# integrated subs taken from Prosite.pm module:


# scan a sequence with a perl pattern
sub scanPattern {
    my ($pattern, $sequence, $behavior, $max_x, $opt_miniprofiles) = @_;
    $behavior ||= 0;
    my $allowOverlap = !($behavior & 1);
    my $allowInclude = $behavior & 2 ;
    $max_x ||= 0;
    my @hits;
    my $pos = 0;
    my @comb = $pattern =~ /\((.*?)\)/g;
    if ($pattern) {
        my $prevstop = -1;
        my @tok;
        my $nter_anchor = $pattern=~/^\^/ ? 1 : 0;
        while (@tok = (substr($sequence, $pos) =~ /^(.*?)($pattern)/)) {
            my $prematch = shift @tok;
            my $subseq = shift @tok;
            my $length = length $subseq;
            my $number_x = 0;
            if (@tok == @comb and @tok) {
                $subseq = "";
                for (my $i = 0; $i < @tok; $i++) {
                    if ($comb[$i] =~ /\./) {
                        $tok[$i] =~ tr/A-Z/a-z/;
                    }
                    elsif ($comb[$i] =~ /^\[\^/) {
                    }
                    elsif (my $x_count = $tok[$i] =~ tr/Xx/Xx/) {
                        $number_x += $x_count;
                    }
                    my $tok = $comb[$i] =~ /\./ ? lc($tok[$i]) : $tok[$i];
                    $subseq .= $tok;
                    if (my @numbers = $comb[$i] =~ /(\d+)/g) {
                        my $biggest = pop @numbers;
                        $subseq .= "." x ($biggest - length $tok);
                    }
                }
            }
            elsif (@tok != @comb) {
                cluck "Internal error with regular expression $pattern\n";
            }
            my $shift = ($length || 1) - 1;
            my $stop = $pos + $length + length $prematch;
            $pos = $stop;
            $pos -= $shift if $allowOverlap;
            last if $pos > length $sequence;
            if ($length) {
                if ($allowInclude or $stop > $prevstop) {
                    if ($number_x <= $max_x or $max_x<0) {
                        push @hits, [$subseq, $stop - $length + 1,
                            $stop, undef,
                            ($opt_miniprofiles && !@hits ? $sequence : undef) ];
                        $prevstop = $stop;
                    }
                }
            }
            else {
                $pos++;
            }
            last if $nter_anchor and $pos;
        }
    }
    return \@hits;
}

sub scanProfiles {
    my ($file, $level_min) = @_;
    $level_min = -99 unless defined $level_min;
    my $results = {};
    local $/ = "\n";
    my $must_open_file = !UNIVERSAL::isa($file, "GLOB");
    my $pfscan_h;
    if ($must_open_file) {
        $pfscan_h = new IO::File $file or confess "Could not open $file : $!";
    }
    else {
        $pfscan_h = $file;
    }
    my $lastac = "";
    my $last_entry;
    while (defined (local $_=<$pfscan_h>)) {
        if (my ($id1, $level, $levelna, $nscore, $rawscore, $from, $to,
            $pffrom, $pfto, $repregion, $repnumber, $ac, $de) = m/
        (
        >\S*
        )
        (?<!L=\d) #These three lines are to fix a bug in the output of pfscan if >99 matches are found: the id and the level are then pasted together with an intervening "1"
        (?<!L=[-\d]\d)
        (?<!L=[-\d]\d\d)
        \s*
        (?:L=(\S+)|(NA))? 
        \s*
        (?:(\S+)\s+)? 
        (\S+) #raw score
        \s+
        pos\.
        \s*
        (\d+) #seq from
        \s*-\s*
        (\d+) #seq to
        \s*
        (?:
        \[\s*
        (\d+) #profile from
        ,\s*
        (-\d+) #profile to
        \]
        \s*
        )? 
        (REGION\d+\s)? 
        (\d+REP\d+\s)? 
        (\S+) 
        (?:\s*(.+))? 
            /x) {
            # fix bug in pfscan which report "******" if nscore>999.999
            $nscore = "999.999" if $nscore eq "*******";
            $level = $level_min if !defined($level) and defined $levelna;
            my $entry = ["", $from, $to, $ac, $pffrom, $pfto, $rawscore, $nscore, $level, undef, $de, []];
            if ($repnumber) { push @{$results->{$lastac}->[-1]->[11]}, $entry }
            else { push @{$results->{$ac}}, $entry }
            $lastac = $ac;
            $last_entry = $entry;
        }
        else {
            next unless $last_entry;
            $last_entry->[0] .= $_;
        }
    }
    if ($must_open_file) {
        close $pfscan_h or confess "Error $? closing $file";
    }
    return $results;
}


sub prositeToRegexp_old {
    local $_ = shift;
    my $notGreedy = shift;
    my $ungreed = $notGreedy ? "?" : "";
    s/\.$//;
    s/-//g;
    1 while s/\{([^\}^\^]*)([^\^^\}][^\}]*\})/\{^$1$2/;
    s/\}/]/g;
    s/\{/[/g;
    s/\(/{/g;
    s/\)/}$ungreed/g;
    s/x/./ig;
    s/B/[ND]/g;
    s/Z/[QE]/g;
    s/\[([^\[\]]*)([<>])([^\[\]]*)\]/(?:[$1$3]|$2)/g;
    s/</^/g;
    s/>/\$/g;
    s/ (\[[^\]]*\]|[\w.]) ( \{ \d+(,\d+)? \} )/)($1$2)(/xg;
    return "($_)";
}

# same, using tokenizing parser
sub prositeToRegexp {
    my $string = shift;
    my $notGreedy = shift;
    my $ungreed = $notGreedy ? "?" : "";
    my $preventX = shift;
    $errstr = undef;
    $errpos = undef;
    my $pushback = "";
    my $ntok = 0;
    my $regexp = "";
    my $get = sub {
        $ntok++;
        return ($pushback =~ s/(.)// ? $1 : ($string =~ s/(.)// ? $1 : undef));
    };
    while (defined (my $tok = &$get)) {
        my $state;
        my $not;
        if ($tok eq "-") {
        }
        elsif ($tok eq "[") {
            while(defined (my $tok = &$get)) {
                last if $tok eq "]";
                $state .= $tok;
            }
        }
        elsif ($tok eq "{") {
            $not = 1;
            while(defined(my $tok = &$get)) {
                last if $tok eq "}";
                $state .= $tok;
            }
        }
        elsif ($tok =~ /[A-Za-z]/) {
            $state = $tok;
        }
        elsif ($tok eq "<") {
            $regexp .= "^";
        }
        elsif ($tok eq ">") {
            $regexp .= '$';
        }
        else {
            $errstr = "Parsing error";
            $errpos = $ntok-1;
            return undef;
        }
        if (defined $state) {
        # read range, e.g. "x(2,5)"
            my $range;
            my $range_char;
            if (defined(my $tok = &$get)) {
                if ($tok eq "(") {
                    while (defined(my $tok = &$get)) {
                        last if $tok eq ")";
                        $range .= $tok;
                    }
                }
                elsif ($tok eq "*") {# support e.g. "<{C}*>"
                    $range_char = $tok;
                }
                else {
                    $pushback .= $tok;
                    $ntok--;
                }
            }

            if ($state =~ /x/i) {
                $state = ".";
            }
            else {
            # handle B/Z unsure amino acids both in pattern and sequence
                if ($not) {
                    $state =~ s/B/NDB/g;
                    $state =~ s/Z/QEZ/g;
                }
                else {
                    $state =~ s/B/NDB/g or $state =~ s/([ND])/$1B/g;
                    $state =~ s/Z/QEZ/g or $state =~ s/([QE])/$1Z/g;
                    $state .= "X" unless $preventX;
                }
            }
            my $mod = $1 if $state =~ s/([<>])//g;
            # "$" is not valid in a character range, we have to convert this
            # to (?:[GH]|$)
            $regexp .= "(";
            $regexp .= "(?:" if $mod;
            $regexp .= "[" if length($state) > 1 or $not;
            $regexp .= "^" if $not;
            $regexp .= $state;
            $regexp .= "]" if length($state) > 1 or $not;
            $regexp .= "|" . ($mod eq "<" ? "^" : '$') . ")" if $mod;
            $regexp .= "{$range}$ungreed" if defined $range;
            $regexp .= "$range_char" if defined $range_char;
            $regexp .= ")";
        }
    }
    return $regexp;
}


################################################################################
# protein_motif_searching.pl core:


BEGIN {
   $VERSION = '1.1';
}

eval { require IPC::Open2 };
my $NO_DIRECT_PIPE=$? || $^O eq "MSWin32" ? 1 : 0;
my $PFSCAN  = 'pfscan';
my $PSA2MSA = 'psa2msa';
my $errcode = 0;
my $MOTIF_AC_REGEXP = '\w+\d+';#'PS\d{5}';
$|= 1;
use Getopt::Long;
use Data::Dumper;
use IO::File;
my @formats = qw(msa gff sequence);
my $formats_string = join " | ", @formats;

sub usage {
    my $progname =  $0;
    $progname    =~ s/.*[\\\/]//;
    print <<EOF;
$progname [options] sequence-file(s)
protein_motif_searching version $VERSION options:
-h : this help screen

Input/Output:
  -p <string> : specify a protein motif pattern based on prosite patterns
  -o <string> : specify output format : gff

Other options (Use under your own risk, without warranty of any kind):
  -e <string> : specify the ID or AC (based on PROSITE database) of an entry in sequence-file
  --reverse   : randomize the sequence database by taking the reverse sequence of each individual entry

Note:
  * The sequence-file must be in FASTA format
  * There may be several -p arguments

EOF
    exit 1;
}

my $opt_noprofiles;
my $opt_onlyprofiles;
my $opt_skip;
my $opt_skiponly;
my $opt_help;
my $opt_format;
my $opt_max_x;
my $opt_nongreedy;
my $opt_nooverlaps;
my $opt_miniprofiles;
my $opt_includes;
my $opt_level = 0;
my $opt_rep_pp_4allprofiles;
my $opt_pfsearch;
my $opt_cutoff;
my $opt_raw;
my $opt_minhits;
my $opt_maxhits;
my $opt_filterheader;
my $opt_reverse;
my $opt_shuffle;
my $opt_no_postprocessing;


my @prosite_files;
my @motifAC_or_userpattern;
my @entries;
my @followpp;
my @external_gff_files;

my $SLASH  = $^O eq "MSWin32" ? "\\" : "\/";
my $TMPDIR = ".";
for my $dir ( $ENV{TMPDIR}, $ENV{SP_TEMP}, $ENV{TMP}, $ENV{TEMP},
    "/tmp", "c:\\temp", "c:\\tmp" ) {
    if ( defined($dir) and -d $dir ) {
        $TMPDIR = $dir;
        last;
    }
}
my $TMP_COUNTER = 1;
my $TMP_BASE    = int( rand( 1000000 ) ) + int( substr( abs( $$ ), -6 ) );

my $scan_profiles;
my $scan_pattern;

my $last_profile_tmp_filename;

Getopt::Long::Configure ("bundling", "no_ignorecase");
GetOptions (
    "r"   => \$opt_noprofiles,
    "m"   => \$opt_onlyprofiles,
    "s"   => \$opt_skip,
    "h"   => \$opt_help,
    "v"   => \$opt_nooverlaps,
    "i"   => \$opt_includes,
    "g"   => \$opt_nongreedy,
    "b:s" => \$opt_miniprofiles,
    "x=i" => \$opt_max_x,
    "l=i" => \$opt_level,
    "o=s" => \$opt_format,
    "d=s" => \@prosite_files,
    "p=s" => \@motifAC_or_userpattern,
    "e=s" => \@entries,
    "f=s" => \@followpp,
    "a"   => \$opt_rep_pp_4allprofiles,
    "w=s" => \$opt_pfsearch,
    "C=f" => \$opt_cutoff,
    "R"   => \$opt_raw,
    "pfscan=s"       => \$PFSCAN,
    "psa2msa=s"      => \$PSA2MSA,
    "minhits=i"      => \$opt_minhits,
    "maxhits=i"      => \$opt_maxhits,
    "filterheader=s" => \$opt_filterheader,
    "reverse"        => \$opt_reverse,
    "shuffle=i"      => \$opt_shuffle,
    "skipflag-only"  => \$opt_skiponly,
    "nopp"           => \$opt_no_postprocessing,
    "gff=s"          => \@external_gff_files
) or &usage;

&usage if $opt_help;
&usage if !@ARGV && -t STDIN && !@external_gff_files;
if ( $opt_pfsearch ) {
    die "OPTION CONFLICT: can't use option".
    " -reverse together with -e option"
    &usage if !@prosite_files &&
        !grep {/^$MOTIF_AC_REGEXP$/} @motifAC_or_userpattern;
    $opt_raw = $1 if defined($opt_cutoff) and $opt_cutoff=~ /^(\d+)$/mg;
}
my $use_pfsearchV3 = index( $opt_pfsearch, 'pfsearchV3') != -1 ? 1 : 0;


my $scan_behavior;
$scan_behavior |= 1 if $opt_nooverlaps;
$scan_behavior |= 2 if $opt_includes;

$opt_format =~ tr/A-Z/a-z/;
die "ERROR:Output format must be one of $formats_string\n"
    unless grep {$_ eq $opt_format} @formats;
my $opt_psa_or_msa = $opt_format eq "msa";

# user patterns specified with -p option
my @userpat;
my $specifiedPrositeMotifByAc={};
map {
    ( /^$MOTIF_AC_REGEXP$/ ?
        $specifiedPrositeMotifByAc->{$_} = 1 : push @userpat,$_ )
} @motifAC_or_userpattern;
keys(%$specifiedPrositeMotifByAc) or $specifiedPrositeMotifByAc = undef;

# find default prosite.dat file
if ( !@prosite_files && !@userpat ) {
    if (defined $ENV{PROSITE} and -e "$ENV{PROSITE}/prosite.dat") {
        @prosite_files = "$ENV{PROSITE}/prosite.dat";
    }
    elsif ( defined $ENV{PROSITE} and -e $ENV{PROSITE} ) {
        @prosite_files = $ENV{PROSITE};
    }
    elsif ( defined $ENV{SPROT} and -e "$ENV{SPROT}/prosite/prosite.dat" ) {
        @prosite_files = "$ENV{SPROT}/prosite/prosite.dat";
    }
    elsif ( -e "prosite.dat" ) {
        @prosite_files = "prosite.dat";
    }
    else {
        die "prosite.dat file not found, please use the -d option";
    }
}

my %SkipFlag;
my %KnownFalsePos;
my @MotifInfo;

my $hideMotifByPSAC          = {};
my $files_miniprofiles       = {};
my $postProcessingByPSAC     = {};
my $motifRank4PostProcessing = {};

# dispatch table for post processing of (matching) target
my $postProcessDispatchTable = {
    'COMPETES_SEQ_WITH' => sub {
        my ( $hit_set_a,$hit_set_b ) = @_;
        return unless $hit_set_a and $hit_set_b;
        ( my $ac_a = $hit_set_a->[0]->[3] || '' ) =~ s/\|\w+$//;
        ( my $ac_b = $hit_set_b->[0]->[3] || '' ) =~ s/\|\w+$//;

        my $max_a = 0.0; map { $max_a=$_->[7] if $_->[7]>$max_a } @$hit_set_a;
        my $max_b = 0.0; map { $max_b=$_->[7] if $_->[7]>$max_b } @$hit_set_b;
        if ( $max_a>$max_b ) {
            map {
                $_->[8]--;
                $_->[9] = "OUTCOMPETED_SEQ_BY_$ac_a"
            } @$hit_set_b
        }
        else {
            map {
                $_->[8]--;
                $_->[9] = "OUTCOMPETED_SEQ_BY_$ac_b"
            } @$hit_set_a
        }
    },

    'FTREP' => sub {
        my ( $target, $effector, $paramstr ) = @_;
        return unless $target;

        my ( $max_score, $ac_cutoff, $size, $gf, $filter ) =
            split( '\|', $paramstr );

        my $motif_param = {
            'motif_id'   => 'fake',
            'motif_acc'  => 'fake',
            'motif_len'  => $size,
            'maxscore'   => $max_score,
            'filter'     => undef,
            'normscore1' => undef,
        };

        my $matchlist = [ map {
            my $estart = $_->[1] - ($_->[4]||1)  +1;
            my $estop  = $_->[2] - ($_->[5]||-1) -1;
            my $size  = 1 + $estop - $estart;
            {
                'match_start'  => $_->[1],
                'match_end'    => $_->[2],
                'match_nscore' => $_->[7],
                'match_rscore' => $_->[6],
                'match_pos'    => $size/2.0 + $estart,
                'match_size'   => $size
            }
        } @$target ];

        my $groups = group_matches( $matchlist, $gf * $size );
        return unless @$groups;

        my ( $ft, $fa ) = FT( $matchlist, $groups, $motif_param, 1 );

        my $ac = AC( $matchlist, $groups, $fa );

        select_match( $matchlist, $groups, $ac, $ac_cutoff, $filter );

        for ( my $i = 0 ; $i < @$target ; $i++ ) {
            my $match      = $target->[ $i ];
            my $ftrepmatch = $matchlist->[ $i ];
            if ( $match->[ 8 ] >= 0 ) {
            }
            elsif ( $ftrepmatch->{ selected } ) {
                $match->[ 8 ] = 0;
                $match->[ 9 ] = 'PROMOTED_FTREP'
            }
            else {
                $match->[ 8 ] = -2;
                $match->[ 9 ] = 'DEMOTED_FTREP'
            }
        }
    }

};

my $allowBidirectionalPP = {
    'PROMOTED_BY' => 1,
};


main();
exit $errcode;


##############################################################################
# Methods

sub tmpnam {
    my $tmp;
    do {
        $tmp = $TMPDIR.$SLASH."ps${TMP_BASE}-".$TMP_COUNTER++.".tmp";
    } while (-e $tmp);
    return $tmp;
}

# -------------------------- output methods --------------------------

# format a field with a certain width
sub pf {
    return $_[0] . ( " " x ( $_[1] - length $_[0] ) )
}

my $HIT_COUNT = 0;
sub dispHits {# display hits (for 1 sequence, 1 motif)
    my ( $header, $sq, $hits, $seqid, $de, $aclist, $psac, $psid, $psde ) = @_;
    return if !$hits || !$psac || $hideMotifByPSAC->{ $psac };

    $sq ||= ''; $seqid ||= ''; $de ||= ''; $psid ||= ''; $psde ||= '';
    ( $de = $hits->[0]->[10] || '' ) =~ s/\.\s*$// if !$de && $opt_pfsearch;
    $de =~ s/[\n\r\t]/ /g;

    my $visible_hits;
    @$visible_hits =
        grep {
            !defined( $_->[7] ) || $_->[8] >= $opt_level
        } @$hits;
    @$visible_hits or return;

    unroll_hits( $visible_hits );
    if ( defined $opt_maxhits ) {
        splice( @$visible_hits, $opt_maxhits - $HIT_COUNT );
        exit 0 if $HIT_COUNT >= $opt_maxhits ;
    }
    my $hit_count  = @$visible_hits;
    $HIT_COUNT    += $hit_count;
    return unless $hit_count;
    return if defined ( $opt_minhits ) and $opt_minhits > $hit_count;

    $visible_hits->[0]->[4] = undef unless defined( $visible_hits->[0]->[7] );

    print $header if defined $header;

    if ( $opt_format eq "fasta" or $opt_psa_or_msa ) {
        for my $hit ( @$visible_hits ) {
            my ( $subseq, $from, $to, $_psac, $pffrom, $pfto, $rawscore,
                 $nscore, $leveln, $levelt, $seqde ) = @$hit;
            my $print_level = defined( $levelt ) ? "L=$levelt " :
                defined($leveln) ? "L=$leveln " : "";
            if ( $opt_pfsearch ) {
                $opt_psa_or_msa = 1;
                print ">$seqid/$from-$to : $de : $psac $print_level\n";
            }
            else {
                print ">$seqid/$from-$to : $psac $psid $print_level\n";
            }
            if ( $subseq and $opt_psa_or_msa ) {
                $subseq =~ s/\n?$/\n/;
                print $subseq;
            }
            else {
                $subseq = substr( $sq, $from-1, $to-$from+1 );
                while ( $subseq =~ /(.{1,60})/g ) {
                    print "$1\n";
                }
                print "\n" if $subseq eq "";
            }
        }
    }
    elsif ( $opt_format eq "gff" ) {
        for my $hit ( @$visible_hits ) {
            my ( $subseq, $from, $to, $_psac, $pffrom, $pfto, $rawscore,
                 $nscore, $leveln, $levelt, $seqde ) = @$hit;
            print join( "\t", $seqid, "protein_motif_searching|v$VERSION", $psac,
                                $from, $to, $nscore || ".", ".", "." );
            my @attr;
            if ( $psid ) {
                $psid =~ s/.*\|//;
                push @attr, "Name \"$psid\"";
            }
            push @attr, "AccessionNumbers " .
                join( " ", map {"\"$_\""} @$aclist )
                if defined $aclist and @$aclist;
            push @attr, "Level $leveln" if defined $leveln;
            push @attr, "LevelTag \"$levelt\"" if defined $levelt;
            push @attr, "RawScore $rawscore" if defined $rawscore;
            push @attr, "FeatureFrom $pffrom"
                if defined $pffrom and $pffrom =~ /\d+/;
            push @attr, "FeatureTo $pfto" if defined $pfto;
            $subseq =~ s/\s//g;
            push @attr, "Sequence \"$subseq\"" if defined $subseq;
            push @attr, "SequenceDescription \"$de\"" if $de;
            push @attr, "SkipFlag 1" if $SkipFlag{$psac};
            push @attr, "KnownFalsePos $KnownFalsePos{$psac}"
            if exists $KnownFalsePos{ $psac };
            print "\t", join " ; ", @attr if @attr;
            print "\n";
        }
    }
    elsif ( $opt_format eq "sequence" ) {
        print ">$seqid : $de : $psac\n", map { "$_\n" }
            $sq =~ /(.{1,60})/g if @$visible_hits;
    }
    else {# default output format
        # fasta like header
        if ( $opt_pfsearch ) {
            print ">$seqid : $de : $psac\n";
        }
        else {
            print ">$seqid : $psac $psid $psde\n";
        }
        for my $hit ( @$visible_hits ) {
            my ( $subseq, $from, $to, $_psac, $pffrom, $pfto, $rawscore,
                 $nscore, $leveln, $levelt, $seqde ) = @$hit;
            my $print_level = defined( $levelt ) ? " L=$levelt" :
                defined( $leveln ) ? " L=$leveln" : "";
            my $fromto = "$from - $to";
            print " " x ( 13-length $fromto ), $fromto;
            if ( $subseq ) {
                $subseq =~ s/\n?$/\n/;# add \n to scanPattern output
                $subseq =~ s/^(?<!\A)(.*)/ $1/mg;
                $subseq =~
                    s/^(.*)/$1 . (" " x (60-length $1)) . $print_level/e
                        if $print_level;
                print "  ", $subseq;
            }
            else {
                $subseq = substr( $sq, $from - 1, $to-$from + 1 );
                my $notfirst;
                while ( $subseq =~ /(.{1,60})/g ) {
                    print " " x 13 if $notfirst++;
                    print " $1\n";
                }
                print "\n" if $subseq eq "";
            }
        }
    }
}

# -------------------------- scan methods --------------------------

# scan entries in specified sequence file 'collection' against motif
sub scanSeqFile {
    my $seqfile = shift or return;

    my $entry = "";
    my $opt_fasta;
    my @motifs_4_normal_scan;
    my @motifs_4_pfsearch;
    my $psinfo_bypsac = {};
    map {
        $_->[0] and $psinfo_bypsac->{$_->[0]} = $_;
        $opt_pfsearch && $_->[2] && $_->[2] eq 'MATRIX' ?
            push @motifs_4_pfsearch, $_ : push @motifs_4_normal_scan, $_
    } @MotifInfo;
    my $all_hits_bypsac_byseqid = {};
    if (grep {$_->[7]} @motifs_4_pfsearch) {
        foreach my $prf_info (@motifs_4_pfsearch) {
            next if $prf_info->[5];
            my $psac = $prf_info->[0] or next;
            my $profile_tmp_file = $prf_info->[7] or next;
            my $psid = $prf_info->[1] || '';
            my $pfsearchhits =
                do_profile_scan($profile_tmp_file, $seqfile, undef, $prf_info);
            foreach my $seq_id (keys %$pfsearchhits) {
                foreach my $hit (@{$pfsearchhits->{$seq_id}}) {
                    $hit->[3] = $psid;
                    push @{$all_hits_bypsac_byseqid->{$seq_id}->{$psac}},
                        $hit;
                }
            }
        }

         foreach my $seq_id ( keys %$all_hits_bypsac_byseqid ) {
            my $pfhits_bypsac = $all_hits_bypsac_byseqid->{ $seq_id } or next;
            postProcess( $pfhits_bypsac ) unless $opt_no_postprocessing;
            foreach my $psac ( keys %$pfhits_bypsac ) {
                my $psid = $psinfo_bypsac->{$psac}->[1] || $psac;
                 my $psde = $psinfo_bypsac->{$psac}->[3] || '';
                dispHits( undef, '',
                    $all_hits_bypsac_byseqid->{$seq_id}->{$psac},
                    $seq_id, '', undef, $psac, $psid, $psde );
            }
        }
    }
    # 'NORMAL' SCAN (every sequence against profile collection)
    if ( @motifs_4_normal_scan ) {# scan entries in sequence file collection
    	my $seqfile_h = new IO::File $seqfile
            or die "Cannot open $seqfile: $!";
        while ( <$seqfile_h> ) {# read sequence file
            my $line = $_;
            if ( /^>(.*)/ ) {# entry is in fasta
                scanFromFastaEntry( $entry );
                $opt_fasta = 1;
                $entry     = '';
            }
            $entry .= $line;
            if ( /^\/\// ) {
                $opt_fasta = 0;
                my @id = $entry =~ /^\s*ID\s+(\w+)/mg;
                my $ac_lines;
                $ac_lines .= $1 while $entry =~ /^\s*AC\s+(.*)/mg;
                my @ac;
                while ( $ac_lines =~ /(\w+)/g ) { push @ac, $1; }
                if ( @id ) {
                    if ( not (@entries) or
                        grep { my $ent = $_;
                                grep{ $_ eq $ent } @id, @ac } @entries) {
                        my $id = $id[ 0 ];
                        my $de =
                            $entry =~ /^DE   (?:Rec|Sub)Name: Full=(.+);/m ?
                                $1 : '';
                        unless ( $de ) {
                            my @de = $entry =~ /^\s*DE\s+(.+)/mg;
                            my $add_space = 0;
                            for ( @de ) {
                                $de .= " " if $add_space;
                                $de .= $_;
                                $add_space = !/-$/;
                            }
                        }

                        if ($entry =~ /^\s*SQ\s+SEQUENCE\b.*\n((.+\n)+)/m) {
                            my $sq = $1;
                            $sq =~ tr/A-Z//cd;
                            scanSeq($id, \@ac, $de || $id, $sq);
                        }
                        else {
                            warn "No sequence found in entry $id";
                        }
                    }
                }
                elsif ($entry =~ /^\s*id /m) {
                }
                elsif ($entry =~ /(.*\S.*)/) {
                    warn "Bad sequence found in file, first line: $1\n";
                    $errcode = 1;
                }
                $entry = "";
            }
        }
        close $seqfile_h;
        if ($entry =~ /^>/) {# process last fasta entry
            scanFromFastaEntry($entry);
        }
        elsif ($entry =~ /(.*\S.*)/) {
            warn "Bad sequence found in file, first line : $1\n";
            $errcode = 1
        }
    }
}

# scan from 1 fasta entry (from a sequence entry collection)
sub scanFromFastaEntry {
    my $entry = shift or return;
    return unless $entry =~ s/^>((\S*).*)\n//;
    my ($fasta_header, $primary_id) = ($1, $2);
    return if defined($opt_filterheader)
        and $fasta_header !~ /$opt_filterheader/o;
    if (not (@entries) or grep {$_ eq $primary_id} @entries) {
        $entry =~ tr/A-Z//cd;
        return scanSeq($primary_id, [], $fasta_header, $entry);
    }
}

# scans one sequence against prosite motifs 
sub scanSeq {
    my ($id, $aclist, $de, $sq) = @_;
    if ($opt_reverse) {
        $sq = reverse $sq;
    }
    if ($opt_shuffle) {
        srand 0;
        my @seq = grep {$_ ne "\n"} split(//,$sq);
        $sq = "";
        for (my $start_win = 0; $start_win < @seq;
                $start_win += $opt_shuffle) {
            my $stop_win = $start_win + $opt_shuffle - 1;
            $stop_win = @seq - 1 if $stop_win >= @seq;
            my @residues = @seq[$start_win..$stop_win];
            while (@residues) {
                $sq .= splice(@residues, int(rand(scalar @residues)) ,1);
            }
        }
    }
    my $all_pfhits = {};
    if ($scan_profiles && !$opt_pfsearch) {
        my $file_source = ($last_profile_tmp_filename ?
            [$last_profile_tmp_filename] : \@prosite_files);
        foreach my $profile_file (@$file_source) {#
            my $pfhits = do_profile_scan($profile_file, undef, $sq);
            foreach my $ps_acid (keys %$pfhits) {
                (my $psac = $ps_acid) =~ s/\|.+//;
                $all_pfhits->{$psac} = $pfhits->{$ps_acid};
            }
        }
    }
    my $hits_by_motif_ac = {};
    foreach (@MotifInfo) {
    	my ($psac, $psid, $type, $psde, $pat, $skip) = @$_;
        next if $skip;
        my $hits;
        if ($type eq "MATRIX") {
            next if $opt_pfsearch;
            $hits = $hits_by_motif_ac->{$psac} = $all_pfhits->{$psac} or next;
        }
        else {
            warn("Empty pattern for $psac\n"), next unless $pat;
            $hits = scanPatternWrapper([$pat, $sq, $scan_behavior,
                        $opt_max_x, $opt_miniprofiles], $psid);
            $hits_by_motif_ac->{$psac} = $hits if @$hits;
        }
        dispHits(undef, $sq, $hits, $id, $de, $aclist, $psac, $psid, $psde)
            if $opt_no_postprocessing;
    }
    unless ($opt_no_postprocessing) {
        postProcess($hits_by_motif_ac);
        foreach (@MotifInfo) {
            my ($psac, $psid, $type, $psde, $pat, $skip) = @$_;
            my $hits = $hits_by_motif_ac->{$psac} or next;
            dispHits(undef, $sq, $hits, $id, $de, $aclist, $psac, $psid, $psde);
        }
    }
    return $hits_by_motif_ac;
}


sub unroll_hits {
    my ($hits) = @_;
    return unless(grep {$_->[11] && @{$_->[11]}} @$hits);
    for (my $i = 0; $i < @$hits; $i++) {
        my ($subseq, $from, $to, $pfid, $pffrom, $pfto, $rawscore,
            $nscore, $leveln, $levelt, $seqde, $subhits) = @{$hits->[$i]};
        next unless($subhits && @$subhits);
        map {$_->[8] = $leveln; $_->[9] = $levelt; $_->[10] = $seqde}
            @$subhits;
        splice @$hits, $i--, 1, @$subhits;
    }
}

sub scanPatternWrapper {
    my ($args, $id) = @_;
    my $out = scanPattern(@$args);
    if ($opt_format eq "psa") {
        $_->[0] =~ s/\.//g for @$out;
    }
    if ($id) {
        $_->[3] = $id for @$out;
    }
    return $out;
}

sub prositeToRegexpWrapper {
    my $out = prositeToRegexp(@_);
    unless (defined $out) {
        print STDERR "protein_motif_searching.pl: Syntax error in pattern".
        " at position $Prosite::errpos\n";
        print STDERR "$_[0]\n";
        print STDERR " " x $Prosite::errpos, "^--- $Prosite::errstr\n";
        exit 1;
    }
    return $out;
}

# -------------------------- post processing --------------------------

sub pp_scan {
    my ($hits, $all_psac_in_seq, $intra_not_inter) = @_;
    return unless($hits && $all_psac_in_seq);
    my $seen = {};
    foreach my $target_psac (@$all_psac_in_seq) {
        my $hit_target_set = $hits->{$target_psac} or next;
        @$hit_target_set or next;
        foreach my $potential_pp (grep {($intra_not_inter ? $_->{effector} eq
            $target_psac : $_->{effector} ne $target_psac) }
            @{$postProcessingByPSAC->{$target_psac}}) {
            my $effector_psac = $potential_pp->{effector} or next;
            my $pp_type = $potential_pp->{type} or next;
            next if $seen->{$pp_type}->{$target_psac}->{$effector_psac};
             $seen->{$pp_type}->{$effector_psac}->{$target_psac} = 1
                unless $allowBidirectionalPP->{$pp_type};

            my $hit_effector_set = $hits->{$effector_psac};
            next unless $hit_effector_set or
                $potential_pp->{allow_no_effector_matches};
            next if $hit_effector_set && !@$hit_effector_set;

            my $pp_value = $potential_pp->{value};

            if ($postProcessDispatchTable->{$pp_type}) {
                &{$postProcessDispatchTable->{$pp_type}}
                    ($hit_target_set,$hit_effector_set,$pp_value);

            } else {
                print STDERR "unknown post-processing key [$pp_type]: ignored";
            }
        }
    }
}

# post process results (inside same entry)
sub postProcess {
    my $hits = shift or return;
    my @all_psac_in_seq = sort {
        ($motifRank4PostProcessing->{$b}||0) <=>
            ($motifRank4PostProcessing->{$a}||0)
    } keys %$hits;
    pp_scan($hits, \@all_psac_in_seq, 1);
    pp_scan($hits, \@all_psac_in_seq, 0);
}


#############################################################################
# CORE

my $isPPLinkedTo = {};
my $isInSpecifiedPPGroup = {};

sub addPPLinkedMotifs2FetchStruct {
    foreach my $psfile (@prosite_files) {
        open PSFILE, $psfile or die "Cannot open $psfile : $!";
        my $ac = '';
        while (<PSFILE>) {# read prosite file
            $ac = $1 if m/^AC (\w+);/;
            $ac = '' if m/^\/\//;
            map {
                s/\(.+\)$//;
                $isPPLinkedTo->{$ac}->{$_} = 1, $isPPLinkedTo->{$_}->{$ac} = 1
                if ($_ && $ac ne $_)
            } split(';\s*',$1)
                if $ac && /^PP \/?\w+:\s*([^\r\n]+)/;
        }
    }
    close PSFILE;
    # collect all profile 'linked by pp' to specified psac
    my @followpp_stack = @followpp;
    while (my $psac = pop @followpp_stack) {# for each specified psac
        $specifiedPrositeMotifByAc->{$psac} = 1;
        # add motif ac to 'motif to fetch/use' struct
        foreach my $linked_motifac (keys(%{$isPPLinkedTo->{$psac}})) {
            # add linked profile ac to profile to be hidden at display:
            $hideMotifByPSAC->{$linked_motifac} = 1
                unless grep {$linked_motifac eq $_} @followpp,
                @motifAC_or_userpattern;
            # add linked profile ac to list of profile to be 'pp-followed'
            # (unless already in specifiedPrositeMotifByAc struct):
            push @followpp_stack, $linked_motifac
                unless $specifiedPrositeMotifByAc->{$linked_motifac};
        }
    }
}


sub main {

    ######################### Motif data parsing/processing #####################
    addPPLinkedMotifs2FetchStruct() if !$opt_no_postprocessing && @followpp;
    my $rank = 0;
    for my $psfile ( @prosite_files ) {
        open PSFILE, $psfile or die "Cannot open $psfile : $!";
        my $ps_entry = ""; my $ac = ""; my $id = ""; my $type = "";
        my $pos = 0;

        PROSITE: while ( <PSFILE> ) {# read prosite file
            $ps_entry .= $_;
            $ac=$1 if !$ac && m/^AC   (\w+);/;
            $id=$1, $type=$2
                if !$id && !$type && m/^ID   (\w+);\s+(\w+)\./;

            if (/^\/\//) {
                my $use_motif = 1;
                $use_motif = 0 if (!$id || !$ac);# skip 'bad' entries
                if ($use_motif && $opt_noprofiles && $type eq "MATRIX") {
                    $use_motif = 0;
                    die "ERROR: profile(s) specified".
                        " with -p or -f option, but -r".
                        " (do not scan profile) option used"
                        if $specifiedPrositeMotifByAc;
                }
                if ($use_motif && $opt_onlyprofiles && $type ne "MATRIX") {
                    $use_motif = 0;
                    die "ERROR: pattern(s) specified with -p or -f option,".
                        " but -m (only scan profile) option used"
                        if $specifiedPrositeMotifByAc;
                }
                $use_motif = 0 if $use_motif && $specifiedPrositeMotifByAc &&
                                    !$specifiedPrositeMotifByAc->{$ac};
                delete ($specifiedPrositeMotifByAc->{$ac})
                    if $specifiedPrositeMotifByAc;
                $motifRank4PostProcessing->{$ac} =++ $rank
                    if $use_motif && !$opt_no_postprocessing;
                 processMotif( $ps_entry ) if $use_motif;
                $ps_entry = $ac = $id = $type = "";
                $pos = tell PSFILE;
            }
        }
        close PSFILE;
    }
    close PROFILE_TMP if $last_profile_tmp_filename;
    if ( $specifiedPrositeMotifByAc ) {
        my @notfound;
        foreach my $ac_id_not_found ( keys( %$specifiedPrositeMotifByAc ) ) {
            push @notfound, $ac_id_not_found;
        }
        die "Prosite entry [@notfound] not found in".
            " specified prosite file(s)\n"
                if @notfound;
    }
    # process -p option user defined patterns
    my $user_ctr = 1;
    for ( @userpat ) {# transform userpat into correct array representation
        my $i = "0" x ( 3 - length( $user_ctr ) ) . $user_ctr++;
        push @MotifInfo, [ "USER$i", undef, "USER", undef,
                        prositeToRegexpWrapper( $_, $opt_nongreedy, 1 ), 0 ];
    }

    ###################### Perform scan (and show results) #####################

    if ( @external_gff_files ) {
        process_external_gff();
    }
    elsif ( $opt_format ne "matchlist" ) {# SCAN (non matchlist output format)
        unshift( @ARGV, '-' ) unless @ARGV;
        while ( my $seqfile = shift @ARGV ) {
            scanSeqFile($seqfile);
        }
    }
    else {# SCAN (matchlist output format) (fugly!?)
        showMatchList();
    }
}
