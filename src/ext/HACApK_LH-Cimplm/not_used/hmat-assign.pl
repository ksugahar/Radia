#/bin/env perl
use strict;
use IO::Handle;
use Getopt::Long qw(:config posix_default no_ignore_case gnu_compat);

sub help {
    print <<'EOS';
    
    hmat-assign.pl: Read <inputfile> or input files whose names are prefixed by <prefix> and generate an eps file that illustrates how leaf matrices in an H-matrix are assigned to MPI processes and worker threads.
        
    Usage:
        
        perl hmat-assign.pl [options not including --mproc] <inputfile>
        perl hmat-assign.pl [options including --mproc]
        
    Options:
        
        --keep, -k: keep all plot files and a gnuplot script after execution.
	--help, -h: show this message.
	--mproc=<prefix>, -m <prefix>:
           Read all input files named <prefix><nnnn> such that all of <prefix>0000, <prefix>0001, ..., and <prefix><nnnn> exist. The postfix number in each input file name is considered as a process number and leaf matrix assignment in all the processes are illustrated in a single figure with process boundaries. 

        <prefix> and <inputfile> cannot be specified at the same time.

    Input file:
    
        Each input file is a text file, in which each line denotes the thread number to which the matrix is assigned, location, size, and type (full or Rk) of a leaf matrix as follows:

          <thr>, <x0>, <y0>, <x1>, <y1>, <mattype>

	where

          <thr>: thread number (>=0),
          <x0>, <y0>: index of the most upper left element,
          <x1>, <y1>: index of the most lower right element,
	  <mattype>: 1 -> Rk-matrix, 2 -> Full-matrix.

EOS
    exit (1);
}

sub y_or_n_p {
    my $message = shift;
    while (1) {
	print "$message (y/n)?: ";
	my $input = <STDIN>;
	chomp ($input);
	if ($input eq 'y') {
	    return 1;
	} elsif ($input eq 'n') {
	    return 0;
	} else {
	    print "Please input y or n.\n";
	}
    }
}

### Get Command-line options
# Keep plot files and gnuplot script? (-k option)
my $Keep_files = 0;
# prefix
my $Prefix;
# debug
my $Debug;
# Help
my $help = 0;

GetOptions (
    'keep|k!' => \$Keep_files,
    'mproc|m=s' => \$Prefix,
    'debug|d!' => \$Debug,
    'help|h' => \$help,
    ) or help();
help() if $help;

if ( $Keep_files ) {
    print "Keep-file option is enabled.\n";
}
 
my $Inputfile0 = shift (@ARGV);
if ($Inputfile0 and $Prefix) {
    print "<inputfile> and <prefix> cannot be specified at the same time.\n";
    help();
}
unless ($Inputfile0 or $Prefix) {
    print "Neither <inputfile> nor <prefix> is specified.\n";
    help();
}

# $Inputfile[$iproc] = Input file for the $iproc-th process
my @Inputfile;
if ( $Inputfile0 ) {
    @Inputfile = ($Inputfile0);
} else {
    foreach my $iproc (0..9999) {
	my $inp = sprintf ("%s%04d", $Prefix, $iproc);
	if ( -e $inp ) {
	    push (@Inputfile, $inp);
	} else {
	    last;
	}
    }
    unless (@Inputfile) {
	print "File $Prefix"."0000 not found.\n";
	help ();
    }
}

# Directory to save all generated files
my $Dirname = ($Inputfile0)?"$Inputfile0-dir":"$Prefix-dir";

# The max thread id number.
my $Thread_max = -1;
# $Plot_file[thread][mattype] = plot file for each <thread, mattype>
my @Plot_file = ();
# $Plot_file[thread][mattype] = file handler for $Plot_file[thread][mattype]
my @Plot_file_fp = ();
# Returns file handler of plot file for matrices of $mattype in $ithr-th thread
# (create and open the file if needed) and update $Thread_max
sub get_plot_file_fp {
    my $ithr = shift;
    my $mattype = shift;
    unless ( $mattype == 1 || $mattype == 2 ) {
	die "Undefined mattype in get_plot_file_fp: $mattype";
    }
    unless ( $Plot_file_fp[$ithr][$mattype] ) {
	my $pltfile = "thread-$ithr-$mattype.plt";
	open ($Plot_file_fp[$ithr][$mattype], "> $Dirname/$pltfile")
	    || die "Failed to open $Dirname/$pltfile";
	$Plot_file[$ithr][$mattype] = $pltfile;
	if ($ithr > $Thread_max) { $Thread_max = $ithr; }
    }
    return $Plot_file_fp[$ithr][$mattype];
}

# Plot file for process boundaries
my $Plot_file_bd = "boundaries.plt";
# File handler for $Plot_file_bd
my $Plot_file_bd_fp;
sub get_plot_file_bd_fp {
    unless ( $Plot_file_bd_fp ) {
	open ($Plot_file_bd_fp, "> $Dirname/$Plot_file_bd")
	    || die "Failed to open $Dirname/$Plot_file_bd";
    }
    return $Plot_file_bd_fp;
}

# eps file to be created
my $Fname_eps = ($Inputfile0)?"$Inputfile0.eps":"$Prefix.eps";

# Gnuplot script file
my $Fname_gp = "plot.gnuplot";

# Range
my $Xmax = 0;
my $Ymax = 0;


### Main
# Make directory
if ( -e $Dirname ) {
    unless (y_or_n_p "$Dirname already exists. Overwrite") {
	print "Abort.";
	exit 99;
    }
} else {
    mkdir $Dirname or die "Failed to generate $Dirname";
}

# Read inputfiles and write to plot files.
foreach my $inp (@Inputfile) {
    my $inp_fp;
    my %x0_vline; my %x1_vline; my %vline_check;
    my %y0_hline; my %y1_hline; my %hline_check;
    open ($inp_fp, "< $inp") or die "Failed to open $inp.";
    while (my $line = <$inp_fp>) {
	if ($line =~ /^\s*(\S+)\s*,\s*(\S+)\s*,\s*(\S+)\s*,\s*(\S+)\s*,\s*(\S+)\s*,\s*(\S+)/) {
	    my ($thr, $x0, $y0, $x1, $y1, $mattype)
		= ($1, $2, $3, $4, $5, $6);
	    my $fp = get_plot_file_fp ($thr, $mattype);
	    my ($x1p, $y1p) = ($x1+1, $y1+1);
	    my ($xc, $yc) = (($x0+$x1p)/2, ($y0+$y1p)/2);
	    my ($xd, $yd) = ($x1p-$xc, $y1p-$yc);
	    $fp->print("$xc $yc $xd $yd\n");
#           # Update the maximum index numbers
	    if ( $x1 > $Xmax ) { $Xmax = $x1; }
	    if ( $y1 > $Ymax ) { $Ymax = $y1; }
#           # entry vertical/horizontal lines
	    $vline_check{$x0} = 1;
	    push (@{$x0_vline{$x0}},   [$y0, $y1+1]);
	    $vline_check{$x1+1} = 1;
	    push (@{$x1_vline{$x1+1}}, [$y0, $y1+1]);
	    $hline_check{$y0} = 1;
	    push (@{$y0_hline{$y0}},   [$x0, $x1+1]);
	    $hline_check{$y1+1} = 1;
	    push (@{$y1_hline{$y1+1}}, [$x0, $x1+1]);
	} else {
	    print "Format error:\n$line\n";
	}
    }
    close ($inp_fp);
    # Draw boundaries
    foreach my $x (sort {$a<=>$b} (keys %vline_check)) {
	my @ranges = ranges_xor ($x0_vline{$x}, $x1_vline{$x});
	if ($Debug) {
	    print "x=$x\n";
	    print "x0: "; print_ranges ($x0_vline{$x});
	    print "x1: "; print_ranges ($x1_vline{$x});
	    print "-> "; print_ranges (\@ranges);
	}
	foreach my $rng (@ranges) {
	    my $fp = get_plot_file_bd_fp ();
	    $fp->print("$x ".($rng->[0])."\n");
	    $fp->print("$x ".($rng->[1])."\n");
	    $fp->print("\n");
	}
    }
    foreach my $y (sort {$a<=>$b} (keys %hline_check)) {
	my @ranges = ranges_xor ($y0_hline{$y}, $y1_hline{$y});
	if ($Debug) {
	    print "y=$y\n";
	    print "y0: "; print_ranges ($y0_hline{$y});
	    print "y1: "; print_ranges ($y1_hline{$y});
	    print "-> "; print_ranges (\@ranges);
	}
	foreach my $rng (@ranges) {
	    my $fp = get_plot_file_bd_fp ();
	    $fp->print(($rng->[0])." $y \n");
	    $fp->print(($rng->[1])." $y \n");
	    $fp->print("\n");
	}
    }
}

sub print_ranges {
    my $ranges = shift;
    foreach my $rng (@{$ranges}) {
	print "[".($rng->[0]).",".($rng->[1])."] "
    }
    print "\n";
}

# Take two references to range lists and returns a list of ranges,
# which are included in just one input ranges.
# Each range is represented as an array that contains two elements [$r0, $r1] ($r0<$r1).
# Ranges in each list should be disjoint.
sub ranges_xor {
    my @switch_points;
    foreach my $i (1..2) {
	my $x_rlist = shift;
	foreach my $rng (@{$x_rlist}) {
	    push (@switch_points, $rng->[0], $rng->[1]);
	}
    }
    @switch_points = sort {$a<=>$b} @switch_points;
    my @ranges = ();
    my $pen_down = 0;
    my $pen_start;
    foreach my $sw (@switch_points) {
	if ($pen_down) {
	    $pen_down = 0;
	    if ($sw > $pen_start) {
		push (@ranges, [$pen_start, $sw]);
	    }
	} else {
	    $pen_down = 1;
	    $pen_start = $sw;
	}
    }
    return @ranges;
}

# Close the plot files
foreach my $ithr (0..$Thread_max) {
    foreach my $mattype (1, 2) {
	if ($Plot_file[$ithr][$mattype])
	{ close ($Plot_file_fp[$ithr][$mattype]); }
    }
}
if ($Plot_file_bd_fp) {
    close ($Plot_file_bd_fp); 
}

# Run gnuplot to generate an eps file and make a gnuplot script file
# with which the same eps file can be generated again.
open (GP, "| gnuplot" ) or die "Failed to exeute gnuplot";
open (FP_GP, "> $Dirname/$Fname_gp") or die "Failed to open $Dirname/$Fname_gp to write";
sub print_gp {
    my $str = shift;
    print GP $str;
    print FP_GP $str;
    print $str;
}
#print_gp "set key left top font 'Helvetica Neue,10' width 2"."\n";
print_gp "set key outside"."\n";
print_gp "unset zeroaxis"."\n";
print_gp "set border 15 back lw 0"."\n";
print_gp "set size ratio ".($Ymax/$Xmax)."\n";

print_gp "set xrange [0:".($Xmax+1)."]"."\n";
#print_gp "set noxtics"."\n";
print_gp "set xtics (0, $Xmax)"."\n";
#print_gp 'set xlabel "time [s]"'."\n";

print_gp "set yrange [".($Ymax+1).":0]"."\n";
#print_gp "set noytics"."\n";
print_gp "set ytics (0, $Ymax)\n";
#print_gp 'set ylabel "worker ID"'."\n";

print_gp 'set terminal postscript eps enhanced color'."\n";
print_gp "set output \"$Fname_eps\""."\n";

my $replot = 0;
my @Thread_color;
{
    my $col = 0;
    my %col_excluded;
    $col_excluded{4} = 1;
    $col_excluded{6} = 1;
    foreach my $ithr (0..$Thread_max) {
	while ($col_excluded{++$col}) {}
	$Thread_color[$ithr] = $col;
    }
}
# Rk matrices
foreach my $ithr (0..$Thread_max) {
    my $file = "$Dirname/$Plot_file[$ithr][1]";
    if (-f $file) {
	print_gp (($replot==0)?'plot ':", \\\n");
	print_gp "\"$file\" w boxxyerrorbars lc ".($Thread_color[$ithr])." lt 1 linewidth 1 fill solid 0.1 title \"Thread $ithr Rk\"";
	$replot = 1;
    }
}
# Full matrices
foreach my $ithr (0..$Thread_max) {
    my $file = "$Dirname/$Plot_file[$ithr][2]";
    if (-f $file) {
	print_gp (($replot==0)?'plot ':", \\\n");
	print_gp "\"$file\" w boxxyerrorbars lc ".($Thread_color[$ithr])." lt 1 linewidth 1 fill solid 0.5 title \"Thread $ithr full\"";
	$replot = 1;
    }
}
# Process boundaries
{
    my $file = "$Dirname/$Plot_file_bd";
    if (-f $file) {
	print_gp (($replot==0)?'plot ':", \\\n");
	print_gp "\"$file\" w lines lc 4 lt 1 linewidth 3 title \"Proc boundaries\"";
    }
}
print_gp "\n";

# Close gnuplot process and script
close (GP);
close (FP_GP);

unless ($Keep_files) {
# remove plot files
    my @rmv_list;
    foreach my $ithr (0..$Thread_max) {
	foreach my $mattype (1,2) {
	    my $file = "$Dirname/$Plot_file[$ithr][$mattype]";
	    if (-f $file) {
		push (@rmv_list, $file);
	    }
	}
    }
    {
	my $file = "$Dirname/$Plot_file_bd";
	if (-f $file) { push (@rmv_list, $file); }
    }
    foreach my $rmv (@rmv_list) {
	unlink "$rmv"
	    or warn "Failed to remove $rmv";
    }
# remove gnuplot script
    unlink "$Dirname/$Fname_gp"
	or warn "Failed to remove $Dirname/$Fname_gp";
# remove directory
    rmdir $Dirname
	or warn "Failed to remove $Dirname";
}

# Show summary
print "* * * Summary * * *\n";
if ( -e $Fname_eps ) {
    print "An eps file $Fname_eps generated.\n";
}
if ( $Keep_files && -d $Dirname ) {
    print "Plot files and a gnuplot script file are stored in $Dirname/.\n";
    print "You can generate $Fname_eps again by:\n";
    print "\$ gnuplot $Dirname/$Fname_gp\n";
}
