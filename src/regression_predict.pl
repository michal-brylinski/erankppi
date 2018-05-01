#!/project/michal/apps/perl/bin/perl

use strict;
use warnings ;
use File::Slurp;



 my @matrix = read_file($ARGV[0]);
 #my ($y,$x1,$x2,$x3,$x4);
 foreach my $a1 (@matrix)
   {

   $a1 =~ s/^\s+|\s+$//g;
   my @split1 = split (/\s+/,$a1);
   my $y_actual=$split1[0];
   my $x1=substr($split1[2],2,8)*1;
   my $x2=substr($split1[3],2,8)*1;
   my $x3=substr($split1[4],2,8)*1;
   my $x4=substr($split1[5],2,8)*1;
   #$reg->include( $y, [ 1.0, $x1, $x2, $x3, $x4 ] );
   #my $y_predicted = ( 957.5815  - $x1*40.4366 - $x2*282.7437 + $x3*33.9166 + $x4*94.5994 ); #regular 
   my $y_predicted = (  959 + $x1*171.9103 - $x2*891.8082 + $x3*12.6889 - $x4*2.1780 ); #high mcc

   print "$y_predicted\n";
   }

