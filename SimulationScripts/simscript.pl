#!/usr/bin/perl

print "\n\n\n Starting simulations ..... \n\n\n";


foreach my $gamma (0, 0.1, 0.2, 0.3, 0.4, 0.5) {

    $ngenes1 = $gamma*50000;
    $ngenes2 = (1.0-$gamma)*50000;

    use constant NLOOP1 => 500;  
    $count = 0;
    while ( $count < NLOOP1)
    {
	$genelength = 1;

	# remove files from last rep
	system("/bin/rm -f simtrees.dat simtrees1.dat simtrees2.dat simtrees3.dat simtrees1sub.dat simtrees2sub.dat temp*");

	# call the python program "CoalTrees_Gamma.py" to simulate trees from parental tree 1 -- set parameters there
	system("python CoalTrees_Gamma.py > simtrees1.dat");
	$nsites1 = $ngenes1*$genelength;

	# now call the python program "CoalTrees_OneMinusGamma" to simulate trees from parental tree 2 -- set parameters there
	system("python CoalTrees_OneMinusGamma.py > simtrees2.dat");
	$nsites2 = $ngenes2*$genelength;

	# put all trees in one file
	system("tail -$ngenes1 simtrees1.dat > simtrees1sub.dat");
	system("tail -$ngenes2 simtrees2.dat > simtrees2sub.dat");
	system("cat simtrees1sub.dat simtrees2sub.dat > simtrees.dat");
	$nsites = $nsites1+$nsites2;
	$ngenes = $ngenes1+$ngenes2;

	# get the file ready for seq-gen
	#system(qq( awk '{print "[$genelength]" \$1}' simtrees.dat > simtrees3.dat));
   
	# The commands below call seq-gen to simulate sites along the gene trees
	# GTR model
	system("./seq-gen < simtrees.dat > infile -q -mGTR -r 1.0 0.2 10.0 0.75 3.2 1.6 -f 0.15 0.35 0.15 0.35 -a 5.0 -g 3 -i 0.2  -l $nsites -p $ngenes");
	# add flags -a alpha -g ncat for gamma ratevar with ncat categories and an alpha of alpha
	# JC69 model -s option removed since scaling done in CoalTrees.py
	# system("./seq-gen < simtrees.dat > infile -s 0.005 -q -mHKY -l $nsites -p $ngenes");

	# remove first line
	system("sed -i '1d' infile");
	# put sequences in order 1 to ...
	system("sort -n infile > temp1");
	system("cat topline temp1 > infile_data");

	system("./ssa > tempout1");
	system("cat tempout1 results.txt > tempout2");
	system("/bin/cp -f tempout2 results.txt");
	system("/bin/rm -f tempout1 tempout2");

	$count++;
	
    }

    system(qq(/bin/cp -f results.txt results_bl1_nc2_50K_$gamma.txt));
    system("/bin/rm -f results.txt");
    
} #end foreach gamma


print "Done.\n";


