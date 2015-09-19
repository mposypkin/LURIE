#!/usr/bin/perl

# the number of an atom to start with in each layer
my $afrom = 0;

# the number an atom to stop in each layer
my $ato = 10;

while(<>) {
  if(/\[(.*)\]/) {
     my @a = split(/,/, $1);
     print join(':', @a), "\n";    
     my $s = @a;
     print "array of size ", $s, "\n";
     my $nl = $s / 3;
     my $j;
     for($j = $afrom; $j <= $ato; $j ++) { 
       my $i;
       for($i = 0; $i < $nl; $i ++) {
         my $y;
         if($i == 0) {
           $y = 0;
         } else {
           $y = $a[3 * $i]
         }
         my $d = $a[3 * $i + 1];
         my $s = $a[3 * $i + 2];
         my $x = $d + $j * $s;
         print $x, " ", $y, " ";
       }
       print "\n";
     }
  }
}