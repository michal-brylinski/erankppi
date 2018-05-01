#!/project/michal/apps/perl/bin/perl -w
 
 use strict;
 use warnings;
 use File::Basename;
 use File::Path;
 use Cwd;
 use Getopt::Long;
 use Pod::Usage;
 use Time::Local;
 use File::Slurp;
 use File::Copy;
 use File::Temp qw/ tempfile tempdir /;
 use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
 use File::Glob;
 use Algorithm::NeedlemanWunsch;
 use List::Util qw(sum);
## $ARGV[0] : zdock out file 
## $ARGV[1] : any dimer file for ialign
## $ARGV[2] : receptor chain id
## $ARGV[3] (x): cut off for ialign (default was 4.5 , we change it to 5A) 
## $ARGV[4] (y): min number of interface residues (default 20 , we change it to 1 to reduce number of bad models)
## $ARGV[5] : efindsiteppi file 
## $ARGV[6] : HOMO - 0 , HETERO - 1


print "------------------------------------------------------------\n";
 print "                        eRankPPI_input\n";
 print "                        version 1.0\n";
 print "       report bugs and issues to smahes2\@tigers.lsu.edu\n";
 print "                                 michal\@brylinski.org\n";
 print "------------------------------------------------------------\n\n";

 die "PPI_IALIGN_BIN is not set\n" if !( $ENV{'PPI_IALIGN_BIN'} );
 die "PPI_CREATLIG is not set\n" if !( $ENV{'PPI_CREATLIG'} );
 

 my $env_bin_path = $ENV{'PPI_IALIGN_BIN'};
 my $env_creatlig = $ENV{'PPI_CREATLIG'};

 #export PPI_IALIGN_BIN="/home/surabhi/apps/ialign/bin"
 #export PPI_CREATLIG="/home/surabhi/apps/zdock3.0.2_linux_x64"

 if ($#ARGV < 7)
 {
  print "        -z <zdock output file>\n";
  print "        -d <any dimer pdb file>\n";
  print "        -c <receptor chain id>\n";
  print "        -p <efindsitePPI prediction file for receptor>\n\n";
  print "additional options:\n\n";
  print "        -x <interface distance cutoff for ialign, default 5.0 Angstrom>\n";
  print "        -y <min number of interface residues, default 1>\n";
  print "        -t < HOMODIMER - 0, HETERODIMER -1, default 0>\n";
  print "        -o < output file name, default erankppi_input.txt>\n";
  die "\n";
 }
 my $i; my $zfile=''; my $dfile=''; my $chain_id='';my $efppifile='';my $thresh_ix=5;my $thresh_iy=1;my $homo_hetero=0;  my $outfile_name = "erankppi_input.txt";
 for ($i = 0; $i <= $#ARGV; $i++)
 {
     if ($ARGV[$i] eq '-z') {$zfile = $ARGV[$i+1] ;}
  elsif ($ARGV[$i] eq '-d') {$dfile = $ARGV[$i+1] ;}
  elsif ($ARGV[$i] eq '-c') {$chain_id = $ARGV[$i+1] ;}
  elsif ($ARGV[$i] eq '-p') {$efppifile = $ARGV[$i+1] ;}
  elsif ($ARGV[$i] eq '-x') {$thresh_ix = $ARGV[$i+1] ;}
  elsif ($ARGV[$i] eq '-y') {$thresh_iy = $ARGV[$i+1] ;}
  elsif ($ARGV[$i] eq '-t') {$homo_hetero = $ARGV[$i+1] ;}
  elsif ($ARGV[$i] eq '-o') {$outfile_name = $ARGV[$i+1] ;}
 }


############ Global Parameters ####################################################################
  local $| = 1;

 my $CONT_CUTOFF ;    #### 4.5 Angstrom distance cutoff for contact calculations
 my $MIN_INT_RES ;     #### minimum number of protein residues of a interface
 my $MIN_NUM_RES  = 25;     #### minimum number of protein residues, ignore short peptide in a PDB file 
 my $VERSION = '1.0b7';     ### ialign version number of this package
 my $pdb_file1; my $pdb_file2;
 my $pchains1 = ''; my $pchains2 = '';
 
 my $pdb_path; my $out_file; my $set_file1; my $set_file2;
 my $tra_flag; my $aln_flag = 2; my $sel_flag;
 my $measure  = 'is';
 
######## READ _DATA_ ################################################################################


 my %matrix;
 my @row1= qw (ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL NH CO);
 while (my $wdat1 = <DATA>)
 {
   chomp $wdat1;
   my $check = substr($wdat1,0,2);
   if ($check eq 'XX' && $wdat1 ne "XX PDPs ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL NH CO")
    {
     my @split_data = split (/\s+/,$wdat1);
     my $nwr1 = $split_data[1];
     for ( my $xg = 0; $xg < 22; $xg++ ) 
      {
       my $xy=$xg+2;
       $matrix{$nwr1.$row1[$xg]} = ($split_data[$xy]*1);
       }
     }
  }

## sequence alignment ################################################################################

 use vars qw( %nwmat4 @nwseq3 @nwseq4 $nwseq3o $nwseq4o );
 my %nwmat1 = (); 
 my @nwmat3 = qw(A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X);
 
 my @nwseq1 = ();my @nwseq2 = ();my $score;
 my $nwseq1o = '';my $nwseq2o = '';
 while ( my $wdat1 = <DATA> )
 {
  chomp $wdat1;my $check = substr($wdat1,0,2);
  if ( length($wdat1) == 70 && $check ne 'XX' and $wdat1 ne '   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X' )
  {
   my $nwr1 = substr($wdat1, 0, 1);   
   for ( my $xg = 0; $xg < 23; $xg++ )
   {
    $nwmat1{$nwr1.$nwmat3[$xg]} = substr($wdat1, 1 + $xg * 3, 3) * 1;
   }
  }
 }
 my $matcher1 = Algorithm::NeedlemanWunsch->new(\&blosum);
 $matcher1->gap_open_penalty(-10);
 $matcher1->gap_extend_penalty(-2);


##Read Zdock.outfile################################################################################

my $outfile = $zfile;
 # open the zdock output file
 open (ZDOUT, "$outfile") || die "Unable to open file $outfile!\n";
 my @zdout_lines = <ZDOUT>;
 chomp(@zdout_lines);
 close(ZDOUT);

##Parse Zdock output file###########################################################################

# parse the header of the zdock output file
(my $n, my $spacing, my $switch_num)=split(" ", $zdout_lines[0]);

 my $line_num = 1;
 my $rec_rand1 = 0.0, my $rec_rand2 = 0.0, my $rec_rand3 = 0.0;
 if ($switch_num ne "")
  {
    ($rec_rand1, $rec_rand2, $rec_rand3)= split(" ", $zdout_lines[$line_num++]);
  }
 (my $lig_rand1, my $lig_rand2, my $lig_rand3)=split(" ", $zdout_lines[$line_num++]);
 (my $rec, my $r1, my $r2, my $r3) = split (" ", $zdout_lines[$line_num++]);
 (my $lig, my $l1, my $l2, my $l3) = split (" ", $zdout_lines[$line_num++]);

 if ($switch_num eq "1")
  {
    my $temp_name = $rec;
    $rec = $lig;
    $lig = $temp_name;
  }
 
 my $dir1 = getcwd(); my $dir2 = tempdir(CLEANUP => 1);

 my $total_lines = $#zdout_lines; my $pred_num = 1;
  
 my %jacc_score;my %interface_model;my %zdock_score;my %pdp_score;
 my @pdp_values; my @zdock_values;

 $line_num++; my $flag=0; my @bad_models;

 my $count=0;my @isscore;my @target = ();my $filename;

 for (my $i= 5; $i <= $total_lines ; $i++)
  {
     $count++;
    (my $angl_x, my $angl_y, my $angl_z, my $tran_x, my $tran_y, my $tran_z, my $score) = split ( " ", $zdout_lines[$i] );
     my ($fh1, $newligfile) = tempfile( DIR => $dir2,SUFFIX => ".$count", UNLINK => 1);
     my $create_cmd = "$env_creatlig/create_lig $newligfile $lig $lig_rand1 $lig_rand2 $lig_rand3 $r1 $r2 $r3 $l1 $l2 $l3 $angl_x $angl_y $angl_z $tran_x $tran_y $tran_z $n $spacing\n";
     if ($switch_num ne "")
     {
        $create_cmd = "$env_creatlig/create_lig $newligfile $lig $switch_num $rec_rand1 $rec_rand2 $rec_rand3 $lig_rand1 $lig_rand2 $lig_rand3 $r1 $r2 $r3 $l1 $l2 $l3 $angl_x $angl_y $angl_z $tran_x $tran_y $tran_z $n $spacing\n";
     }
     my ($fh0, $dimer_target) = tempfile( DIR => $dir2,SUFFIX => ".pdb", UNLINK => 1);
     my ($fh3, $dimer_target1) = tempfile( DIR => $dir2,SUFFIX => ".pdb", UNLINK => 1);
     system($create_cmd);
    
    open(INFILE1, "$rec") or die "ERROR:canont open file $!";
    open(INFILE2, "$newligfile") or die "ERROR:canont open file $!";
    open(OUTFILE1, ">$dimer_target") or die "ERROR:cannot create file $!"; 
    open(OUTFILE2, ">$dimer_target1") or die "ERROR:cannot create file $!"; 
    foreach(<INFILE1>)
     {
      my $x = substr($_,0,55);
      my $y = sprintf("$x\%5.2f \%5.2f \n",1,0);
      print OUTFILE1  "$y"; print OUTFILE2  "$y";
     } 
    print OUTFILE1 "TER\n";print OUTFILE2  "TER\n";
    foreach(<INFILE2>)
     {
      my $x = substr($_,0,55);
      my $y = sprintf("$x\%5.2f \%5.2f \n",1,0);
      print OUTFILE1  "$y";print OUTFILE2  "$y";
     } 
    print OUTFILE1 "TER\n";print OUTFILE2  "TER\n";
    close INFILE1;
    close INFILE2;
    close OUTFILE1;close OUTFILE2;
    $zdock_score{$count}= log $score;
    $pdp_score{$count}= pdp($dimer_target);
    
   push(@pdp_values,$pdp_score{$count});
   push(@zdock_values,$zdock_score{$count});

    my ($fh4, $ialign_out) = tempfile( DIR => $dir2);#, UNLINK => 1);
    $filename = basename($dimer_target);
    $filename=~ s/.pdb//; 

    $CONT_CUTOFF  = $thresh_ix*1;    #### 4.5 Angstrom distance cutoff for contact calculations
    $MIN_INT_RES  = $thresh_iy*1;     #### minimum number of protein residues of a interface
    ialign ($dfile,$dimer_target,$ialign_out,$dir2);
    my @ialign_out1 = read_file($ialign_out);
    my @grep1 = grep(/IS-score/,@ialign_out1);#print "@ialign_out1\n";
    my @target_interface=();
 #print "@grep1\n"; 
    if (exists $grep1[0])
         { 
           my $int_name1 = $dir2."/".$filename; 
           my @int_name = glob "$int_name1??_int.pdb";
           my @cont_list = glob "$int_name1??_con.lst";
           ##### START : interface file generation for svm and nb score ########
           my @interface_file = read_file($int_name[0]);
           my $i=0; my $j=0;  
           for (0..$#interface_file)
            { 
             if ($interface_file[$_] =~ (/TER/) ) 
              {
              $j=0;
              next;
              }
             $i=$i+1;
             $j=$j+1;
             my @split_l1 = split(/ +/,$interface_file[$_]);
             if ($split_l1[4] eq $chain_id) 
             { 
              push(@target_interface,$split_l1[5]);
             }
            }
           $interface_model{$count}=[@target_interface];
          }
          ##### END :interface file generation for svm and nb score ########

     ##### START : JACC Calculation ###################################

     $CONT_CUTOFF  = 10; #$ARGV[3]*1;    #### 4.5 Angstrom distance cutoff for contact calculations
     $MIN_INT_RES  = 2; #$ARGV[4]*1;     #### minimum number of protein residues of a interface
     my($fh5, $ialign_outj1) = tempfile( DIR => $dir2);#, UNLINK => 1);
     ialign ($dfile,$dimer_target1,$ialign_outj1,$dir2);
     my @ialign_out2 = read_file($ialign_outj1);#print "@ialign_out2\n";
     my @grep2 = grep(/IS-score/,@ialign_out2);
     my @pdb1 = read_file($rec);my @pdb2 = read_file($newligfile);
    # print "@grep2\n"; 
     if (exists $grep2[0] && $grep1[0])
         { 
          my $filename1 = basename($dimer_target1);
          $filename1 =~ s/.pdb//; 
          my $int_name1_j = $dir2."/".$filename1; 
          my @int_name_j = glob "$int_name1_j??_int.pdb";
          my @cont_list_j = glob "$int_name1_j??_con.lst";
          #my @pdb1 = read_file($rec);my @pdb2 = read_file($newligfile);
          my @rec_resi_no=();
          foreach (@pdb1)
           {
            my $atom_name = substr($_, 12,  4);$atom_name =~ s/^\s+//;$atom_name =~ s/\s+//;
            if ( length($_) > 53 && substr($_, 0, 6) eq "ATOM  " && $atom_name eq "CA")
            {
            my $resi_no = substr($_,22,4); $resi_no =~ s/^\s+//;
            push (@rec_resi_no,$resi_no);
            }
           }
           my $length_receptor = $#rec_resi_no +1;
           my $seq1 = pdb2seq(@pdb1);my $seq2 = pdb2seq(@pdb2);
           my ($seq1a,$seq2a) = seq_align($seq1,$seq2);

           my $map_out=map_alignments($seq1a,$seq2a,$length_receptor);
           my %alignlig2rec1 =%{$map_out};
  
           my ($mcc1,$jac1)=find_jacc($int_name_j[0],$cont_list_j[0],\%alignlig2rec1,\@rec_resi_no);
           $jacc_score{$count}= sprintf "%0.5f",$jac1; #print "$count $jac1\n";
          }
         system "rm $newligfile\n";   
         system "rm $ialign_out\n";
 
  }

########SVM & NB Calculation########################################################

####Parse efppi.sites.dat#### 
 my @efppi1=  read_file($efppifile);
 my @efppi2 = grep (/^RESIDUE/,@efppi1);
 my %ppi_svm;my %ppi_nbc;
 foreach (@efppi2)
  {
   my @split_21 = split (/\s+/,$_);
   my @array21=();
   my $key2 = substr ($_,22,4)*1;
   $ppi_svm{$key2}= substr ($_,37,7);
   $ppi_nbc{$key2}= substr ($_,46,8); 
   }

 my %svm ; my %nb; 
 my $sum_svm=0; my $sum_nb=0;
 my @svm_values =();my @nb_values =();
 foreach my $modelno_w1(sort  {$a<=>$b} keys %interface_model)
  {
   foreach my $resno (@{$interface_model{$modelno_w1}})
    {
     if (exists $ppi_svm{$resno})
      {
        $sum_svm = $sum_svm + $ppi_svm{$resno};
        $sum_nb = $sum_nb + $ppi_nbc{$resno};
       } 
     }
    $svm{$modelno_w1}=$sum_svm;
    $nb{$modelno_w1}=$sum_nb;
    push (@svm_values,$sum_svm);
    push (@nb_values,$sum_nb);
    $sum_svm=0;$sum_nb=0;$count=0;
   }

  my $sum_svm_all= sum (@svm_values); my $sum_nb_all= sum (@nb_values);
  my $n1 = scalar keys %svm;
  my $mean_svm = $sum_svm_all / $n1;

  my $sqsum_svm = 0;
  for (@svm_values) 
   {
    $sqsum_svm += ( $_ ** 2 );
    } 
  $sqsum_svm  /= $n1;
  $sqsum_svm -= ( $mean_svm ** 2 );
  my $stdev_svm = sqrt($sqsum_svm);
  my $mean_nb = $sum_nb_all / $n1;

  my $sqsum_nb = 0;
  for (@nb_values) 
   {
    $sqsum_nb += ( $_ ** 2 );
    } 
   $sqsum_nb  /= $n1;
   $sqsum_nb -= ( $mean_nb ** 2 );
   my $stdev_nb = sqrt($sqsum_nb);

  my %z_svm; my %z_nbc;

  foreach my $modelno_w2(sort  {$a<=>$b} keys %interface_model)
   {
     my $zsvm = ($svm{$modelno_w2}-$mean_svm)/$stdev_svm; 
     my $znbc = ($nb{$modelno_w2}-$mean_nb)/$stdev_nb;
    $z_svm{$modelno_w2} = sprintf "%0.5f", scale ($zsvm, -2.4,2.5);
    $z_nbc{$modelno_w2} = sprintf "%0.5f", scale ($znbc, -2.3,2.5);
   }

#######Normalize and scale###########################################################

#ZDOCK#
  my $sum_zdockscore_all= sum (@zdock_values); 
  my $n_zdockscore = scalar keys %zdock_score;
  my $mean_zdockscore = $sum_zdockscore_all / $n_zdockscore;
  my $sqsum_zdockscore = 0;
  for (@zdock_values)
   {
    $sqsum_zdockscore += ( $_ ** 2 );
    }
   $sqsum_zdockscore  /= $n_zdockscore;
   $sqsum_zdockscore -= ( $mean_zdockscore ** 2 );
   my $stdev_zdockscore = sqrt($sqsum_zdockscore);

  my %z_zdockscore;my $zzdockscore=0;
  foreach my $modelno_x1(sort  {$a<=>$b} keys %zdock_score)
   {
     $zzdockscore = ($zdock_score{$modelno_x1}-$mean_zdockscore)/$stdev_zdockscore;
     $z_zdockscore{$modelno_x1}= sprintf "%0.5f",scale ($zzdockscore,-1.47,5);
   }
#PDP#
  my $sum_pdp_all= sum (@pdp_values); 
  my $n_pdp = scalar keys %pdp_score;
  my $mean_pdp = $sum_pdp_all / $n_pdp;
  my $sqsum_pdp = 0;
  for (@pdp_values)
   {
    $sqsum_pdp += ( $_ ** 2 );
    }
   $sqsum_pdp  /= $n_pdp;
   $sqsum_pdp -= ( $mean_pdp ** 2 );
   my $stdev_pdp = sqrt($sqsum_pdp);

  my %z_pdp;my $zpdp;
  foreach my $modelno_x2(sort  {$a<=>$b} keys %pdp_score)
   {
     $zpdp = ($pdp_score{$modelno_x2}-$mean_pdp)/$stdev_pdp;
     $z_pdp{$modelno_x2}= sprintf "%0.5f",scale ($zpdp,-4.5,4);
   }

####### M/C learning input file ####################################################
 my @out_file;
 foreach my $key (sort { $a<=>$b } keys %z_nbc)
  {
   if (exists $z_zdockscore{$key} && exists $z_pdp{$key} && exists $jacc_score{$key} & exists $z_svm{$key})
   {
    if ($homo_hetero == 0)
     {
      push (@out_file,"0 1:$z_svm{$key} 2:$z_nbc{$key} 3:$jacc_score{$key} 4:$z_pdp{$key} 5:$z_zdockscore{$key}\n");
      }
    elsif ($homo_hetero == 1)
     {
      push (@out_file, "0 1:$z_svm{$key} 2:$z_nbc{$key} 3:$z_pdp{$key} 4:$z_zdockscore{$key}\n");
      }
    }
   }
 write_file($outfile_name,@out_file);
#######CLOSE#########################################################################
 sub scale 
  {
  my $in = $_[0] ;
  my $lb = $_[1]  ;
  my $ub = $_[2] ;
  my $out = ($in - $lb)/($ub-$lb)* 2.0 - 1.0;
  $out = -1.0 if ( $out < -1.0 );
  $out = 1.0 if ( $out > 1.0 );
  return $out;
  }

#####################################################################################
sub pdp
 {
 my $fpdb = $_[0];
 open (PDB, "$fpdb") || die "Cannot open $fpdb for reading.\n";
  my @pdb1=<PDB>;
  chomp(@pdb1);
 close (PDB);
 my @pdb2 = grep(/ATOM/, @pdb1);
 my $prev_chain=0;
 my $chain=0;my $flag=0;my $prev_resno=0;my $resi_name;
 my $key ; my @array1; my @array2; my %h1; my %h2; my %h3; my %h4;
 my ($x1,$x2,$x3,$x4)=(0,0,0,0); my ($y1,$y2,$y3,$y4)= (0,0,0,0); my ($z1,$z2,$z3,$z4)= (0,0,0,0);
 my @backbone=qw/C O N CA/;my $count=0;
 foreach my $wpdb2 (@pdb2)
 {
  if ( length($wpdb2) > 53 )
  {
   if ( substr($wpdb2, 0, 6) eq "ATOM  " )
   {
    my $atom_name = substr($wpdb2, 12,  4);$atom_name =~ s/^\s+//;$atom_name =~ s/\s+//;
    #$resi_name = substr($wpdb2, 17,  3);
     $chain = substr($wpdb2,21,1);
    my $resi_no = substr($wpdb2,22,4); $resi_no =~ s/^\s+//;
    my $x_coor = substr($wpdb2,30,8);
    my $y_coor = substr($wpdb2,38,8);
    my $z_coor = substr($wpdb2,46,8);
    #$flag = 1 if $chain ne $prev_chain && $prev_chain ne 0;
    if ($resi_no == $prev_resno || $prev_resno == 0)
     {
      $resi_name = substr($wpdb2, 17,  3);
      $chain = substr($wpdb2,21,1);
      if ($atom_name eq "N"){$x1 = $x_coor;$y1=$y_coor;$z1=$z_coor;}
      elsif ($atom_name eq "O"){$x2 = $x_coor;$y2=$y_coor;$z2=$z_coor;}
      elsif ($atom_name eq "CA"){$x4 = $x_coor;$y4=$y_coor;$z4=$z_coor;}
      elsif ($atom_name ne "C") 
       {
        #print "$atom_name $x_coor $y_coor $z_coor\n";
        $x3=$x3+$x_coor; $y3=$y3+$y_coor; $z3=$z3+$z_coor;
        $count++;
       }
      }
    else
     {
      my ($s1,$s2,$s3)=($x4,$y4,$z4);
      if ($count != 0) { $s1=$x3/$count;  $s2=$y3/$count; $s3=$z3/$count; }
      my @array=($x1,$y1,$z1,$x2,$y2,$z2,$s1,$s2,$s3,$x4,$y4,$z4,$resi_name, $prev_chain); 
      if ($flag==0) { $h1{$prev_resno}=[@array];}
      elsif ($flag==1) { $h2{$prev_resno}=[@array];}
      $count=0; $x3=0; $y3=0;$z3=0;
      if ($atom_name eq "N"){$x1 = $x_coor;$y1=$y_coor;$z1=$z_coor;}
     }
     $flag = 1 if $chain ne $prev_chain && $prev_chain ne 0;
     $prev_resno=$resi_no; $prev_chain = $chain;
     }
    }
   }
 my ($s1,$s2,$s3)=($x4,$y4,$z4);
 $s1=$x3/$count if $count !=0 ; $s2=$y3/$count if $count !=0; $s3=$z3/$count if $count !=0;
 my @array=($x1,$y1,$z1,$x2,$y2,$z2,$s1,$s2,$s3,$x4,$y4,$z4,$resi_name,$chain); 
 $h1{$prev_resno}=[@array];



 my $pdp1=0;
 foreach my $key_chk1(sort{$a <=> $b} keys %h1) 
  {
   my ($a1,$b1,$c1,$a2,$b2,$c2,$a3,$b3,$c3,$a4,$b4,$c4,$resi_name1,$chain1)=@{$h1{$key_chk1}};
   foreach my $key_chk2(sort{$a <=> $b} keys %h2)
    { 
      my ($d1,$e1,$f1,$d2,$e2,$f2,$d3,$e3,$f3,$d4,$e4,$f4,$resi_name2,$chain2)=@{$h2{$key_chk2}};
      my $dist_chek=sqrt( ($a4-$d4)**2 + ($b4-$e4)**2 + ($c4-$f4)**2 ); 
     
      if ( $dist_chek < 18 )
       {     
        my $NN= sqrt( ($a1-$d1)**2 + ($b1-$e1)**2 + ($c1-$f1)**2 );
        my $NO= sqrt( ($a1-$d2)**2 + ($b1-$e2)**2 + ($c1-$f2)**2 );
        my $ON= sqrt( ($a2-$d1)**2 + ($b2-$e1)**2 + ($c2-$f1)**2 );
        my $OO= sqrt( ($a2-$d2)**2 + ($b2-$e2)**2 + ($c2-$f2)**2 );

        my $NS= sqrt( ($a1-$d3)**2 + ($b1-$e3)**2 + ($c1-$f3)**2 );
        my $OS= sqrt( ($a2-$d3)**2 + ($b2-$e3)**2 + ($c2-$f3)**2 );
        my $SN= sqrt( ($a3-$d1)**2 + ($b3-$e1)**2 + ($c3-$f1)**2 );
        my $SO= sqrt( ($a3-$d2)**2 + ($b3-$e2)**2 + ($c3-$f2)**2 );

        my $SS= sqrt( ($a3-$d3)**2 + ($b3-$e3)**2 + ($c3-$f3)**2 ); 
        
        if ($NN <= 4.0){$pdp1=$pdp1+6.84};
        if ($NO <= 4.0){$pdp1=$pdp1-2.14};
        if ($ON <= 4.0){$pdp1=$pdp1-2.14};        
        if ($OO <= 4.0){$pdp1=$pdp1-1.82}; 

        if ($NS <= 5.6){my $key="NH".$resi_name2;$pdp1=$pdp1+$matrix{$key}};
        if ($OS <= 5.6){my $key="CO".$resi_name2;$pdp1=$pdp1+$matrix{$key}};
        if ($SN <= 5.6){my $key=$resi_name1."NH";$pdp1=$pdp1+$matrix{$key}};
        if ($SO <= 5.6){my $key=$resi_name1."CO";$pdp1=$pdp1+$matrix{$key}};

        if ($SS <= 6.8){my $key=$resi_name1.$resi_name2;$pdp1=$pdp1+$matrix{$key}};
        #print "$key_chk1 $key_chk2 $dist_chek $SS\n";
       }    
     }
   }
 return ($pdp1);
}
###################################################################################################
 sub pdb2seq
  {
  my %aa=qw(
	ALA 	A
	CYS 	C
	ASP 	D
	GLU 	E 
	PHE 	F
	GLY	G
	HIS	H
	ILE	I
	LYS	K
	LEU	L
	MET	M
	ASN	N
	PRO	P
	GLN	Q
	ARG	R
	SER	S
	THR	T
	VAL	V
	TRP	W
	TYR	Y);

  my @pdb_in = @_;my $residue ; my $resno;
  my $oldresno=-1;my $seq=(); #print "@pdb_in";
  foreach (@pdb_in)
   {
    if (/^ATOM/)
    {
     my $type = substr($_,13,2); 
     if ($type eq "CA")	
      {	
       my $res = substr($_, 17, 3); 
       chomp($residue=$aa{$res});$residue=~ s/^\s+//;$residue=~s/\s+$//;
       $resno=substr($_, 22, 4);
       if ($resno>$oldresno)
        {
         $seq=$seq.$residue ; 
         $oldresno=$resno;
         }
       }
     }
    }
   return ($seq);
   }

 sub blosum 
 { 
  my ($anw, $bnw) = @_;  
  my $snw = 0;  
  $snw = $nwmat1{$anw.$bnw} if ( exists $nwmat1{$anw.$bnw} );  
  return ($snw);
 }
 
 sub on_align  {  
  my ($i, $j) = @_;
  $nwseq1o = $nwseq1[$i] . $nwseq1o;
  $nwseq2o = $nwseq2[$j] . $nwseq2o;
 }
 
 sub on_shift_a 
 {  
  my $i = shift;
  $nwseq1o = $nwseq1[$i].$nwseq1o;
  $nwseq2o = "-$nwseq2o";
 };

 sub on_shift_b 
 { 
  my $j = shift;
  $nwseq1o = "-$nwseq1o";
  $nwseq2o = $nwseq2[$j].$nwseq2o;
 };

 sub seq_align 
 { 
  my $a = $_[0];my $b = $_[1];
  @nwseq1 = split (//, $a);
  @nwseq2 = split (//, $b); 
  $nwseq1o = '';$nwseq2o = '';
  $score = $matcher1->align(\@nwseq1,\@nwseq2,{align => \&on_align, shift_a => \&on_shift_a, shift_b => \&on_shift_b });
  return($nwseq1o,$nwseq2o); 
 }

 sub map_alignments
   {
    my $aln_target = $_[0];
    my $aln_template = $_[1];
    my $target_length = $_[2];
    my $length =  length ($aln_target);
    my $m1 = 0 ; my $m2 =0;my %alignlig2rec;
    for (my $i=0; $i<= $length; $i++)
     {
      my $count_m1=0; my $count_m2=0;
      if ( substr($aln_target,$i,1) ne '-') { $m1++; $count_m1 = 1;}
      if ( substr($aln_template,$i,1) ne '-') { $m2++; $count_m2 = 1;}
      if ($count_m1 == 1 && $count_m2 == 1 && $m1 <= $target_length)
       {
         $alignlig2rec{$m2}=$m1;
         }
       }
    return \%alignlig2rec ; 
   }
###################################################################################################
# find_jacc (int_pdb, con_list,%alignlig2rec1,@rec_resi_numbers )

 sub find_jacc
 {
 my @interface_file = read_file($_[0]);
 my $con=0;my %map; my $prev_chain=0;
 my $chain;my $flag=0;
 my %alignlig2rec =  %{$_[2]};
 foreach (@interface_file)
  {
   if ( length($_) > 53 && substr($_, 0, 6) eq "ATOM  ")
   {
    $chain = substr($_,21,1);
    $con++; $flag = 1 if $chain ne $prev_chain && $prev_chain ne 0;
    my $resi_no = substr($_,22,4); $resi_no =~ s/^\s+//;
    if ($flag==0)
     {
      $map{$con}=$resi_no."_R";
     }
    else
     {
      $map{$con}=$alignlig2rec{$resi_no}."_L";
      }
    $prev_chain = $chain;
   } 
 }

 my @input1 = read_file($_[1]);#(contact file)
 my %contact ;
 my @input2 = grep(/^RES  /,@input1);
 foreach (@input2)
  {
   my @split1 = split(/\s+/,$_);
   my $con1 = $split1[1];
   foreach (3..$#split1)
    {
     my $key = $map{$con1}.":".$map{$split1[$_]};
     $contact{$key}=1;
     }
   }

## Calculation of MCC and JACCKARD Index
 my $tp = 0;
 my $tn = 0;
 my $fp = 0;
 my $fn = 0;

 my %check ;
 my @rec_resi_no= @{$_[3]};
 foreach my $rn1 ( @rec_resi_no)
  {
  foreach my $rn2 ( @rec_resi_no)
   {
    my $key = $rn1."_R:".$rn2."_L";
    my $reversekey = $rn1."_L:".$rn2."_R";
    my $samekey= $rn2."_L:".$rn1."_R";
    my $tt1 = 0; my $tt2=0;
    if (!exists $check{$key} && !exists $check{$samekey})
    {
    if ( exists $contact{$key}){$tt1=1;}
    if ( exists $contact{$reversekey}){$tt2=1};  
    if ( $tt1 and $tt2 )
    {
    $tp++;
    }
    elsif ( $tt1 and !$tt2 )
    {
     $fn++;
    }
    elsif ( !$tt1 and $tt2 )
    {
     $fp++;
    }
    else
    {
     $tn++;
    }
    }
   $check{$key}=1;
   $check{$reversekey}=1;
 }
}

 my $tpr = 0.0;
 my $fpr = 0.0;
 my $acc = 0.0;
 my $spc = 0.0;
 my $ppv = 0.0;
 my $npv = 0.0;
 my $fdr = 0.0;
 my $mcc = 0.0;
 my $jac = 0.0;

 $tpr = $tp / ($tp + $fn) if ( ($tp + $fn) );
 $fpr = $fp / ($fp + $tn) if ( ($fp + $tn) );
 $acc = ($tp + $tn) / ($tp + $fn + $fp + $tn) if ( ($tp + $fn + $fp + $tn) );
 $spc = $tn / ($fp + $tn) if ( ($fp + $tn) );
 $ppv = $tp / ($tp + $fp) if ( ($tp + $fp) );
 $npv = $tn / ($tn + $fn) if ( ($tn + $fn) );
 $fdr = $fp / ($fp + $tp) if ( ($fp + $tp) );

 $mcc = ($tp*$tn - $fp*$fn)/sqrt(($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn)) if ( (($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn)) != 0.0 );
 $jac = $tp/($fp+$fn+$tp);
 return ($mcc,$jac)
}

#####################################################################################
############################         Subroutine   ialign     ########################
#####################################################################################
 
 

 #ialign($ARGV[0],$ARGV[1],$ARGV[2]);
 
 
 sub ialign 
 { 
  my $pdb_file1 = $_[0]; my $pdb_file2 = $_[1]; 
  my $out_file = $_[2];my $work_path = $_[3];

  my @pdb_set1 = (); my @pdb_set2 = ();
 
  if( defined $pdb_file1 ) { $pdb_set1[0] = {'pdb_file'=>$pdb_file1,'pchains'=>$pchains1}; }
  if( defined $pdb_file2 ) { $pdb_set2[0] = {'pdb_file'=>$pdb_file2,'pchains'=>$pchains2}; }
 
  my $switch = '';
  if( $aln_flag == 0 or $aln_flag == 1 or $aln_flag == 2 ) {$switch .= " -v $aln_flag ";}
  else { die "Error: invalid option for alignment printout format\n"; }
 
  my $out;
  if( defined $out_file ) {open $out, ">$out_file" or die "Error: could not open $out_file\n";}
  else { open $out, ">&STDOUT";}
 
  my %parsed = ();   ### save information of processed PDB files
  for (@pdb_set1) {my $name = $_->{pdb_file} . '_' . $_->{pchains}; $parsed{$name} = $_;}
  for (@pdb_set2) {my $name = $_->{pdb_file} . '_' . $_->{pchains}; $parsed{$name} = $_;}
 
  for my $name (sort keys %parsed) {
   my ($parsed_pdb_file, $pchains) = parsePDBFile( $parsed{$name}, $pdb_path, $work_path, $env_bin_path, $out );
   $parsed{$name}->{'parsed_file'} = $parsed_pdb_file;
   $parsed{$name}->{'chains'}      = $pchains;
  }
 
  for my $name (sort keys %parsed) {
   my $parsed_pdb_file = $parsed{$name}->{parsed_file};
   my $pchains         = $parsed{$name}->{chains};
 
   $parsed{$name}->{ppi} = extractProtProtInt( $env_bin_path, $parsed_pdb_file, $pchains, $out );
  }
 
  my $set1 = \@pdb_set1;my $set2 = \@pdb_set2;
  if (scalar @pdb_set2 == 0 ) { $set2 = \@pdb_set1; }  ### all-against-all
 
  my $num1 = scalar @$set1;my $num2 = scalar @$set2;my %complete = ();my $counter  = 1;
 
  for(my $m=0; $m<$num1; $m++) {
   my $name1    = $$set1[$m]->{pdb_file} . '_' . $$set1[$m]->{pchains};
   my $ppi1     = $parsed{$name1}->{ppi};
   my $num_ppi1 = scalar @$ppi1;
 
   for(my $i=0; $i<$num_ppi1; $i++) {
 
     for(my $k=0; $k<$num2; $k++) {
       my $name2     = $$set2[$k]->{pdb_file} . '_' . $$set2[$k]->{pchains};
       my $ppi2      = $parsed{$name2}->{ppi};
       my $num_ppi2  = scalar @$ppi2;
 
       for(my $j=0; $j<$num_ppi2; $j++) {
 
 	my $pair     = $$ppi1[$i]->{name} . '_' . $$ppi2[$j]->{name};
 	my $pair_rev = $$ppi2[$j]->{name} . '_' . $$ppi1[$i]->{name};
 
 	next if( exists $complete{$pair_rev} || exists $complete{$pair} );  ### avoid repeat
 	next if( $$ppi1[$i]->{name} eq $$ppi2[$j]->{name} and not defined $sel_flag );  ### against itself
 
 	my ($result, $trans_vec, $rot_mat, $parsed_out) = callISalign( $$ppi1[$i], $$ppi2[$j], $switch, $out, $tra_flag );
 	#if( defined $out_file ) { print ">>>$counter $result\n"; }
 
 
 	$counter++;
 	$complete{$pair} = 1;
 	$complete{$pair_rev} = 1;
       } ### j interface 2
     } ### k protein 2
 
   } ### i interface 1
 } ### m protein 1
 close $out;
 
 }
 
#### perform interface alignment with IS-align####################################################
 sub callISalign {
   my ($ppi1, $ppi2, $switch, $out, $tra_flag) = @_;
 
   my $name1     = $ppi1->{name};
   my $int_file1 = $ppi1->{intfile};
   my $con_file1 = $ppi1->{confile};
 
   my $name2     = $ppi2->{name};
   my $int_file2 = $ppi2->{intfile};
   my $con_file2 = $ppi2->{confile};
 
   my $pair = "$name1 vs $name2";
   print $out ">>>$name1 vs $name2\n";
 
   ### call IS-align
   my $alnout = `$env_bin_path/IS-align $switch $int_file1 $int_file2 $con_file1 $con_file2`;
 
   ### parsing output of IS-align
   my ($trans_vec, $rot_mat, $score_ln, $parsed_out) = parseISalignOut( $alnout );
 
   ### outputing results
   print $out "\n$parsed_out\n\n";
  if ( defined $score_ln )
   {
   my $result = "$pair: $score_ln";
   return ($result, $trans_vec, $rot_mat, $parsed_out);
   }
 }
 
######### parse a raw PDB file  ###################################################################
 sub parsePDBFile {
   my ($parsed_rec, $pdb_path, $work_path, $env_bin_path, $outfile) = @_;
 
   my $pdb_file   = $parsed_rec->{pdb_file};
   my $pchain_lst = $parsed_rec->{pchains};
 
   my ($base, $path, $ext) = fileparse( $pdb_file, qr/\..*/ );
   my $parsed_pdb_file = "$work_path/$base.parsed";
 
   my @pchains = ();
   unless ( -s $parsed_pdb_file ) {
     if( not (-e $pdb_file) and defined $pdb_path ) {
       $pdb_file = "$pdb_path/$pdb_file";
     }
 
     unless( -s $pdb_file ) {
       if( -s "$pdb_file.pdb" )    { $pdb_file = "$pdb_file.pdb"; }
       elsif( -s "$pdb_file.ent" ) { $pdb_file = "$pdb_file.ent"; }
       else {
 	die "Error: could not find the pdb file $pdb_file\n";
       }
     }
 
     my $out = `perl $env_bin_path/process_pdb.pl -pdb $pdb_file -o $parsed_pdb_file -nowarn`;
   }
 
   @pchains = split(//,$pchain_lst);
   if( scalar @pchains == 0 ) { @pchains = getProtChains( $parsed_pdb_file ); }
 
   print $outfile "$pdb_file:  found ", scalar @pchains, " protein chains ", join(' ',@pchains), "\n";
 
   return( $parsed_pdb_file, \@pchains );
 }
 
########### get protein chains from  Calpha atoms #################################################
 sub getProtChains {
   my $pdbfile = shift @_;
 
   open PDB, "<$pdbfile" or die "Error: could open file $pdbfile\n";
   my %pchain = ();
   while (my $line = <PDB>) {
     next if ( $line =~ /^REMARK/ );
     last if ( $line =~ /^ENDMDL/ );  #### only consider the first model
 
     #### PDB format: http://www.wwpdb.org/documentation/format23/sect9.html ####
     if ( $line =~ /^(ATOM  |HETATM)/ ) {
       (my $atomname = substr($line, 12, 4))=~ s/ //g;
       (my $altloc = substr($line, 16, 1));
       (my $resname = substr($line,17,3)) =~ s/ //g;
       (my $chain = substr($line, 21, 1));
       (my $resid = substr($line, 22, 4)) =~ s/ //g;
       (my $icode = substr($line, 26, 1)) =~ s/ //g;
 
       if( $atomname eq 'CA' && convAA( $resname ) ne 'X' ) {
 
 	if ( $chain eq '' ) { $chain = '_'; }
 
 	unless( exists $pchain{$chain} ) {
 	  $pchain{$chain} = 0;
 	}
 
 	$pchain{$chain}++;
       }
     }
 
   }
 
   close PDB;
 
   my @protein_chains = ();
   foreach my $chain (sort keys %pchain) {
     if( $pchain{$chain} > $MIN_NUM_RES ) {
       #push( @protein_chains, { 'id'=>$chain, 'len'=>$pchain{$chain} } );
       push( @protein_chains, $chain );
     }
   }
 
   return @protein_chains;
 }
 
 
######## convert residues from three-letter to one-letter and vice visa ###########################
 sub convAA {
   my $res = shift @_;
 
   my %ts=(
      'GLY'=>'G',     'ALA'=>'A',     'VAL'=>'V',     'LEU'=>'L',
      'ILE'=>'I',     'SER'=>'S',     'THR'=>'T',     'CYS'=>'C',
      'MET'=>'M',     'PRO'=>'P',     'ASP'=>'D',     'ASN'=>'N',
      'GLU'=>'E',     'GLN'=>'Q',     'LYS'=>'K',     'ARG'=>'R',
      'HIS'=>'H',     'PHE'=>'F',     'TYR'=>'Y',     'TRP'=>'W',
   );
 
   if (exists $ts{$res}) {
     return $ts{$res};
   } else {
     return "X";
   }
 }
 
 
 
#################################################################################################
 
 sub extractProtProtInt {
   my ($env_bin_path, $parsed_pdb_file, $pchains, $outfile) = @_;
 
   my ($base, $path, $ext) = fileparse( $parsed_pdb_file, qr/\..*/ );
 
   my @ppi = ();
   my $num_pchains = scalar @$pchains;
   for(my $i=0; $i<$num_pchains; $i++) {
     my $pchain_a = $$pchains[$i];
     for(my $j=$i+1; $j<$num_pchains; $j++) {
       my $pchain_b = $$pchains[$j];
       my $pair     = $pchain_a . $pchain_b;
       my $name     = $base . $pair;
       my $int_file = "$path$name\_int.pdb";
       my $con_file = "$path$name\_con.lst";
 
       my $num_int_res = 0;
       if( -s $int_file and -s $con_file ) {
 	my $output = `grep "^Total interacting residues" $con_file`;
 	if( $output =~ /Total interacting residues.* (\d+)/ ) {
 	  $num_int_res = $1;
 	}
 	else { die "Error: file $con_file seems corrupted\n"; }
       }
       else {
 	my $output = `$env_bin_path/extint -s $parsed_pdb_file -r $pchain_a -l $pchain_b -d $CONT_CUTOFF -i $int_file -c $con_file -e -ca`;
 	if( $output =~ / (\d+) interface residues extracted/ ) {
 	  $num_int_res = $1;
 	}
 	else { die "Error: extraction of interfaces failed\n"; }
       }
 
       if( $num_int_res > $MIN_INT_RES ) {
         push( @ppi, { 'name'=>$name, 'size'=>$num_int_res, 'pair'=>$pair,
 		      'intfile'=>$int_file, 'confile'=>$con_file, 'pdbfile'=>$parsed_pdb_file } );
       }
     }
   }
 
   printf $outfile "$parsed_pdb_file: found %d valid PPI(s)", scalar @ppi;
   for(my $i=0; $i<scalar @ppi; $i++) {
     print $outfile " $ppi[$i]->{pair} $ppi[$i]->{size}";
   }
   print $outfile "\n";
 
   return \@ppi;
 }

##### get translation vector and rotation matrix from IS-align output###########################
 sub parseISalignOut {
   my $alnout = shift @_;
 
   my @output = split(/\n/, $alnout);
 
   ### remove useless information for end user
   my $num = scalar @output;
   for(my $i=0; $i < $num; $i++) {
     last if ( $output[0] =~ /^Structure 1/ );
     shift @output;
   }
 
   ### find measurement of the alignment
   my ($alnlen, $rmsd, $sid, $col, $pvalue, $score);
   my $score_ln;
   my @trans_vec = ();
   my @rot_mat = ();
   for (my $i=0; $i < scalar @output; $i++)  {
     my $line = $output[$i];
     if ( $line =~ /^\w{2}-score =\s*(\d+\.\d+), P-value = \s*(\S+)/ ) {
       $score    = $1;
       $pvalue   = $2;
       $score_ln = $line;
     }
     elsif ( $line =~ /^\w{2}-score =\s*(\d+\.\d+)/ ) {
       $score    = $1;
       $score_ln = $line;
     }
     elsif ( $line =~ /^Number of aligned residues  =\s*(\d+)/ ) {
       $alnlen = $1;
     }
     elsif ( $line =~ /^Number of aligned contacts  =\s*(\d+)/ ) {
       $col = $1;
     }
     elsif ( $line =~ /^RMSD =\s*(\d+\.\d+), Seq identity  =\s*(\d+\.\d+)/ ) {
       $rmsd = $1;
       $sid  = $2;
     }
     elsif( $line =~ /Transformation matrix/ ) {
       for(my $j=0; $j<3; $j++) {
         my @fields = split(' ',$output[$i+$j+2]);
         push( @trans_vec, $fields[1] );
         push( @rot_mat, [ ($fields[2], $fields[3], $fields[4]) ] );
       }
       last;
     }
   }
 
   my $out = join("\n", @output);
   return( \@trans_vec, \@rot_mat, $score_ln, $out );
 }
###################################################################################################
__DATA__
XX PDPs ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL NH CO
XX ALA -3.56 4.17 1.69 -1.45 -3.00 -2.52 3.56 0.54 -0.07 -0.65 -0.58 4.26 5.39 -0.17 4.50 1.23 -0.21 -3.36 -0.60 -2.05 0.60 -0.44
XX ARG 4.17 -1.66 -0.12 -2.00 -1.02 4.20 -1.80 -0.04 0.29 2.80 -1.80 -0.90 -1.10 -3.40 1.33 -0.02 -0.55 -4.65 -3.18 -2.69 1.60 -1.37  
XX ASN 1.69 -0.12 -0.33 -1.09 6.52 -2.25 -2.20 1.13 0.95 -0.86 2.20 1.22 -1.33 -1.83 -0.42 -0.09 -0.99 -2.81 -1.58 0.66 0.17 -0.71
XX ASP -1.45 -2.00 -1.09 0.97 2.57 0.63 1.45 0.67 -0.71 0.01 -0.62 -2.45 1.21 6.05 1.11 -2.42 -0.92 -0.90 -1.03 -1.61 0.85 0.07 
XX CYS -3.00 -1.02 6.52 2.57 10.00 -0.63 -3.34 0.94 -4.93 -1.13 -2.02 -2.95 0.00 6.01 4.87 -1.74 5.30 2.60 -0.36 1.37 0.44 -0.72
XX GLN -2.52 4.20 -2.25 0.63 -0.63 9.72 1.25 -0.76 1.65 -1.83 -2.63 1.23 1.94 8.05 0.80 -0.75 1.01 -3.60 -2.36 -0.68 1.35 -1.28 
XX GLU 3.56 -1.80 -2.20 1.45 -3.34 1.25 1.49 0.89 -0.81 3.26 0.05 -3.01 -0.91 0.36 2.76 -1.26 -0.63 -0.53 -1.74 1.04 -0.38 0.05
XX GLY 0.54 -0.04 1.13 0.67 0.94 -0.76 0.89 -0.34 0.44 -1.86 0.35 -0.26 -0.90 -0.16 2.87 1.20 0.43 0.10 -1.63 -1.06 -1.28 0.49
XX HIS -0.07 0.29 0.95 -0.71 -4.93 1.65 -0.81 0.44 1.04 -1.82 2.20 0.42 -4.03 -2.13 1.65 1.24 0.35 1.68 -4.25 -0.20 -0.65 -1.15
XX ILE -0.65 2.80 -0.86 0.01 -1.13 -1.83 3.26 -1.86 -1.82 -3.27 -4.09 1.34 0.50 -6.05 0.25 -1.24 -.66 -5.01 -2.88 -3.78 2.55 -0.09
XX LEU -0.58 -1.80 2.20 -0.62 -2.02 -2.63 0.05 0.35 2.20 -4.09 -3.43 0.80 -1.90 -2.21 -1.40 0.29 0.72 -0.39 -1.43 -1.80 -0.70 -0.09
XX LYS 4.26 0.90 1.22 -2.45 -2.95 1.23 -3.01 -0.26 0.42 -1.34 0.80 3.24 3.64 -2.39 3.85 1.71 0.65 -3.37 -2.34 -2.01 1.02 -0.57
XX MET 5.39 -1.10 -1.33 1.21 0.00 1.94 -0.91 -0.90 -4.03 -0.50 -1.90 3.64 10.00 -4.17 -4.02 -0.55 -2.81 -5.89 -2.28 -1.32 0.74 -1.22
XX PHE -0.17 -3.04 -1.83 6.05 6.01 8.05 0.36 -0.16 -2.13 -6.05 -2.21 -2.39 -4.17 -6.42 -3.96 0.17 -1.08 -2.65 -3.51 -2.17 0.49 -0.79
XX PRO 4.50 1.33 -0.42 1.11 4.87 0.80 2.76 2.87 1.65 0.25 -1.40 3.85 -4.02 -3.96 -2.50 -0.29 -1.18 -4.56 -2.23 0.22 0.01 -0.69
XX SER 1.23 -0.02 -0.09 -2.42 -1.74 -0.75 -1.26 1.20 1.24 -1.24 0.29 1.71 -0.55 0.17 -0.29 0.45 -1.01 -0.86 -1.38 1.27 0.65 -0.94
XX THR -0.21 -0.55 -0.99 -0.92 5.30 1.01 -0.63 0.43 -0.35 1.66 0.72 0.65 -2.81 -1.08 -1.18 -1.01 0.22 2.43 -1.76 0.83 0.76 -1.22
XX TRP -3.36 -4.65 -2.81 -0.90 2.60 -3.60 -0.53 0.10 1.68 -5.01 -0.39 -3.37 -5.89 -2.65 -4.56 -0.86 2.43 8.64 -2.83 -1.87 -0.28 -0.70
XX TYR -0.60 -3.18 -1.58 -1.03 -0.36 -2.36 -1.74 -1.63 -4.25 -2.88 -1.43 -2.34 -2.28 -3.51 -2.23 -1.38 -1.76 -2.83 -2.15 -1.37 1.38 -1.50
XX VAL -2.05 -2.69 0.66 -1.61 1.37 -0.68 1.04 -1.06 -0.20 -3.78 -1.80 -2.01 -1.32 -2.71 0.22 1.27 0.83 -1.87 -1.37 -0.37 -0.22 -0.27  
XX NH 0.60 1.60 0.17 0.85 0.44 1.35 -0.38 -1.28 -0.65 2.55 -0.70 1.02 0.74 0.49 0.01 0.65 0.76 -0.28 1.38 -0.22 6.84 -2.14
XX CO -0.44 -1.37 -0.71 0.07 -0.72 -1.28 0.05 0.49 -1.15 -0.09 -0.09 -0.57 -1.22 -0.79 -0.69 -0.94 -1.22 -0.70 -1.50 -0.27 -2.14 1.82
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
A  5 -2 -1 -2 -1 -1 -1  0 -2 -1 -2 -1 -1 -3 -1  1  0 -3 -2  0 -2 -1 -1
R -2  7 -1 -2 -4  1  0 -3  0 -4 -3  3 -2 -3 -3 -1 -1 -3 -1 -3 -1  0 -1
N -1 -1  7  2 -2  0  0  0  1 -3 -4  0 -2 -4 -2  1  0 -4 -2 -3  4  0 -1
D -2 -2  2  8 -4  0  2 -1 -1 -4 -4 -1 -4 -5 -1  0 -1 -5 -3 -4  5  1 -1
C -1 -4 -2 -4 13 -3 -3 -3 -3 -2 -2 -3 -2 -2 -4 -1 -1 -5 -3 -1 -3 -3 -2
Q -1  1  0  0 -3  7  2 -2  1 -3 -2  2  0 -4 -1  0 -1 -1 -1 -3  0  4 -1
E -1  0  0  2 -3  2  6 -3  0 -4 -3  1 -2 -3 -1 -1 -1 -3 -2 -3  1  5 -1
G  0 -3  0 -1 -3 -2 -3  8 -2 -4 -4 -2 -3 -4 -2  0 -2 -3 -3 -4 -1 -2 -2
H -2  0  1 -1 -3  1  0 -2 10 -4 -3  0 -1 -1 -2 -1 -2 -3  2 -4  0  0 -1
I -1 -4 -3 -4 -2 -3 -4 -4 -4  5  2 -3  2  0 -3 -3 -1 -3 -1  4 -4 -3 -1
L -2 -3 -4 -4 -2 -2 -3 -4 -3  2  5 -3  3  1 -4 -3 -1 -2 -1  1 -4 -3 -1
K -1  3  0 -1 -3  2  1 -2  0 -3 -3  6 -2 -4 -1  0 -1 -3 -2 -3  0  1 -1
M -1 -2 -2 -4 -2  0 -2 -3 -1  2  3 -2  7  0 -3 -2 -1 -1  0  1 -3 -1 -1
F -3 -3 -4 -5 -2 -4 -3 -4 -1  0  1 -4  0  8 -4 -3 -2  1  4 -1 -4 -4 -2
P -1 -3 -2 -1 -4 -1 -1 -2 -2 -3 -4 -1 -3 -4 10 -1 -1 -4 -3 -3 -2 -1 -2
S  1 -1  1  0 -1  0 -1  0 -1 -3 -3  0 -2 -3 -1  5  2 -4 -2 -2  0  0 -1
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  2  5 -3 -2  0  0 -1  0
W -3 -3 -4 -5 -5 -1 -3 -3 -3 -3 -2 -3 -1  1 -4 -4 -3 15  2 -3 -5 -2 -3
Y -2 -1 -2 -3 -3 -1 -2 -3  2 -1 -1 -2  0  4 -3 -2 -2  2  8 -1 -3 -2 -1
V  0 -3 -3 -4 -1 -3 -3 -4 -4  4  1 -3  1 -1 -3 -2  0 -3 -1  5 -4 -3 -1
B -2 -1  4  5 -3  0  1 -1  0 -4 -4  0 -3 -4 -2  0  0 -5 -3 -4  5  2 -1
Z -1  0  0  1 -3  4  5 -2  0 -3 -3  1 -1 -4 -1  0 -1 -2 -2 -3  2  5 -1
X -1 -1 -1 -1 -2 -1 -1 -2 -1 -1 -1 -1 -1 -2 -2 -1  0 -3 -1 -1 -1 -1 -1
