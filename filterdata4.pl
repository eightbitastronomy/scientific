#!/usr/bin/perl
#
###############################################################################################################################
# filterdata5.pl
# author: eightbitastronomy (eightbitastronomy@protonmail.com)
# License: none. 
#          This source code may be used freely and without constraint.
# Last updated: 25 oct 2019 for commenting.
#
# Don't judge >-( This is a strictly utilitarian script with absolutely nothing clever in it :-\
#
# script to filter data by
# a) isolating data points with a given value within a RELATIVE margin.
# b) screening for less than or greater than a threshold
# c) finding a minimum or maximum value
# data must be one record per line, fields space-delineated
# sorting & subsorting may be applied to cases (a) and (b) above
# Usage:
#   -p  position of field to be matched, 0 to n-1
#   -v  value to be matched
#   -r  fractional (relative) error in matching. Default value is 1e-04.
#   -d  discard n lines from beginning of input
#   -g  filter by greater than (-v) value
#   -l  filter by less than (-v) value
#   -M  filter for maximum value
#   -m  filter for minimum value
#   -1  primary sort key field index, 0 to n-1. 1. If not specified, will not sort
#   -2  secondary sort key field index, 0 to n-1. If not specified, sub-sorting will not be used.
# 31 oct 2017
##############################################################################################################################


use Getopt::Std;


my $encoding = ":encoding(ASCII)";
my $precision;
my $master=[];
my $first;
my $second;
my $func;
my $discard;        #number of lines at beginning of input to be discarded. Will be defined.
my $cmp;            #comparison function
my $process;        #filter function


my $usage = "Usage:\n  ' [OPTS] file1 file2 ...' or use a pipe with OPTS\n   -p  position of field to be matched, 0 to n-1\n   -v  value to be matched\n   -r  fractional (relative) error in matching. Default value is 1e-04.\n   -1  primary sort key field index, 0 to n-1. Default value is the first (0th) field.\n   -d  discard n lines at beginning of input\n   -g  filter by greater than (-v) value\n   -l  filter by less than (-v) value\n   -M  filter for maximum value\n   -m  filter for minimum value\n   -1  primary sort key field index, 0 to n-1. If not specified, will not sort\n   -2  secondary sort key field index, 0 to n-1. If not specified, sub-sorting will not be used.";



sub nonextremum
{
    # Arguments: FILE, $position, $matchvalue, \@masterlist
    my $ff = *{$_[0]};
    my $line;
    for $i (0..$discard-1)
    {
	$line=<$ff>;
    }
    while ($line = <$ff>)
    {
	my @buf = split(' ',$line);
	if ( &$cmp(@buf[$_[1]],$_[2]) )
	{
	    push(@{$_[3]},\@buf);   
	}
    } 

}


sub extremum
{
    # Arguments: FILE, $position, [useless], \@masterlist
    my $ff = *{$_[0]};
    my $line;
    my @best;
    for $i (0..$discard-1)
    {
	$line=<$ff>;
    }
    if ( $line = <$ff> )
    {
	@best = split(' ',$line);
    }
    while ($line = <$ff>)
    {
	my @buf = split(' ',$line);
	if ( &$cmp(@buf[$_[1]],@best[$_[1]]) )
	{
	    @best = @buf;
	}
    } 
    push(@{$_[3]},\@best);   
}


sub near
{
    #@ARGV[0] is being passed in as second parameter
    if ( (abs($_[0]-$_[1])/abs($_[1]))<$precision )
    {
	return 1;
    }
    else
    {
	return 0;
    }
}


sub compareDouble
{
    #$a and $b will be pre-defined by sort
    if ($a->[$first] < $b->[$first])
    {
	return -1;
    }
    if ($a->[$first] == $b->[$first] )
    {
	if ($a->[$second] < $b->[$second])
	{
	    return -1;
	}
	if ($a->[$second] > $b->[$second])
	{
	    return 1;
	    }
	return 0;
    }
    return 1;
}


sub compareSingle
{
    #$a and $b will be pre-defined by sort
    if ($a->[$first] < $b->[$first])
    {
	return -1;
    }
    if ($a->[$first] == $b->[$first] )
    {
	return 0;
    }
    return 1;
}


sub greaterThan
{
    return ( $_[0] > $_[1] );
}


sub lessThan
{
    return ( $_[0] < $_[1] );
}



# ----- process cmd-line args -----

my %args = ( 'p' => undef, 
	     'v' => undef,
	     'r' => 1e-04,
	     'h' => undef,
	     'l' => undef,
	     'g' => undef, 
	     'd' => undef,
	     'M' => undef,
	     'm' => undef,
             '1' => undef,    
             '2' => undef      );
getopts('p:v:r:d:12lgMmh',\%args);

defined $args{'h'} and die "$usage\n";
defined $args{'p'} or die "Positional key required for filtering multi-field data.\n$usage\n";

# below: if neither greater-than nor less-than, then filter by 'near', then if neither min nor max, 'near' is still the default. Process is nonextremum unless changed to nonextremum.
$process = \&nonextremum;
if ( defined $args{'g'} )
{
    defined $args{'l'} and die "Cannot use -l and -g simultaneously.\n$usage\n";
    ( (defined $args{'M'}) or (defined $args{'m'}) ) and die "Cannot use -l or -g with a min/max filter.\n$usage\n";
    $cmp = \&greaterThan;
}
else 
{
    if ( defined $args{'l'} )
    {
	( (defined $args{'M'}) or (defined $args{'m'}) ) and die "Cannot use -l or -g with a min/max filter.\n$usage\n";
	$cmp = \&lessThan;
    }
    else
    {
	$cmp = \&near;
    }
}
if ( defined $args{'M'} )
{
    defined $args{'m'} and die "Cannot use -M and -m simultaneously.\n$usage\n";
    ( (defined $args{'l'}) or (defined $args{'g'}) ) and die "Cannot use -m or -M with a less-than/greater-than filter.\n$usage\n";
    $cmp = \&greaterThan;
    $process = \&extremum;
}
else 
{
    if ( defined $args{'m'} )
    {
	( (defined $args{'l'}) or (defined $args{'g'}) ) and die "Cannot use -m or -M with a less-than/greater-than filter.\n$usage\n";
	$cmp = \&lessThan;
	$process = \&extremum;
    }
    else
    {
	defined $args{'v'} or die "Key for comparison value required for all but -m or -M searches\n$usage\n";
    }
}

# below, these vars are not pretty, but they reduce the number of hash calls. This, for a large data file (a million lines), adds up to time saved.
$precision = $args{'r'};
if (defined $args{'2'})
{
    if ( (defined $args{'M'}) or (defined $args{'m'}) )
    {
	print("Sorting will not be applied.\n");
    }
    else
    {
	if ( defined $args{'1'} )
	{
	    $first = $args{'1'};
	    $second = $args{'2'};
	    $func = \&compareDouble;
	}
	else
	{
	    die "Cannot sub-sort (-2) without primary sort (-1).\n$usage\n";
	}
    }
}
else
{
    if ( defined $args{'1'} )
    {
	if ( (defined $args{'M'}) or (defined $args{'m'}) )
	{
	    print("Sorting will not be applied.\n");
	}
	else
	{
	    $first = $args{'1'};
	    $func = \&compareSingle;
	}
    }
}
if ( defined $args{'d'} )
{
    $discard=$args{'d'};
}
else
{
    $discard=0;
}



# ----- filter data -----

if ( -p STDIN )
{
    ( @ARGV==0 ) or die "Ambiguous usage or data sources: pipe detected with excess command-line arguments.\n";
    &$process( STDIN, $args{'p'}, $args{'v'}, $master );
}
else
{
    @ARGV>0 or die "No input to process.\n";
    while ( my $name = shift @ARGV )
    {
	
	open(f,"< $encoding",$name) or die "Can't open $name: $!";
	&$process( f, $args{'p'}, $args{'v'}, $master );
	close(f);
    }
}



# ----- sort data and output -----


if ( defined $func )
{
    @{$master} = sort $func @{$master};
}
foreach my $out (@{$master})
{
    #my $vval = abs($out->[1] / $out->[2]);
    #print "@{$out}[0] @{$out}[2] $vval @{$out}[3] @{$out}[4]\n";
    print "@{$out}\n";
}



# ----- clean up -----

undef($master);

