# erankppi
Dimer model ranking with machine learning

Installation and requirements

Install the following Perl modules, which are available from CPAN:
* File::Temp 	http://search.cpan.org/~dagolden/File-Temp-0.2301/
* File::Slurp 	http://search.cpan.org/~uri/File-Slurp-9999.19/
* File::Copy 	http://search.cpan.org/~rjbs/perl-5.18.1/
* Compress::Zlib 	http://search.cpan.org/~pmqs/IO-Compress-2.063/
* List::Util 	http://search.cpan.org/~pevans/Scalar-List-Utils-1.35/
* Algorithm::NeedlemanWunsch 	http://search.cpan.org/~vbar/Algorithm-NeedlemanWunsch-0.04/
* Benchmark 	http://search.cpan.org/~rjbs/perl-5.18.1/
* Cwd 	http://search.cpan.org/~smueller/PathTools-3.40/
* YAML 	http://search.cpan.org/~ingy/YAML-0.88/

Additional software:
LIBSVM 	Support Vector Machines for classification and regression
iAlign  a method for the structural comparison of protein-protein 

STEP 1. Create input file (erankppi_input.txt) for machine learning or Linear regression prediction. 
For this step you'll need the following :
	- zdock output file
	- efindsiteppi prediction file for the receptor
	- ialign package

Set the environment variables : PPI_IALIGN_BIN , PPI_CREATLIG

	export PPI_IALIGN_BIN="/$DIR1/ialign/bin"
	export PPI_CREATLIG="/$DIR2/zdock3.0.2_linux_x64"

Usage : ./erankppi_input.pl

        -z <zdock output file>
        -d <any dimer pdb file (you may use 9rub-dimer.pdb from the example directory)>
        -c <receptor chain id>
        -p <efindsitePPI prediction file for the receptor>

additional options:

        -x <interface distance cutoff for ialign, default 5.0 Angstrom>
        -y <min number of interface residues, default 1>
        -t < HOMODIMER - 0, HETERODIMER -1, default 0>
        -o < output file name, default erankppi_input.txt>

Example:
../erankppi_input.pl -z 11as.zdock.out -d 9rub-dimer.pdb -c A -p 11as.sites.dat

STEP 2 
For homodimers run svm-predict on the input file generate in step 1.

/$DIR/libsvm-3.14/svm-predict erankppi_input.txt /$DIR/eRankPPI/model/model_svm output.predict.txt

For heterodimers run linear regression on the input file.

/$DIR/eRankPPI/src/regression_predict.pl erankppi_input.txt  > $output_dir/$i.predict.txt
