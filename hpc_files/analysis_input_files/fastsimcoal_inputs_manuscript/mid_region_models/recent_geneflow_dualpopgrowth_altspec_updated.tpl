//Parameters for the coalescence simulation program : fastsimcoal.exe
2 samples to simulate :
//Population effective sizes (haploid number of genes)
N_POP1
N_POP2
//Haploid samples sizes
16
16
//Growth rates: negative growth implies population expansion
GR_POP1
GR_POP2
//Number of migration matrices : 0 implies no migration between demes
3
//Migration matrix 0
0	MIG01
MIG10	0
//Migration matrix 1
0	0
0	0
//Migration matrix 2
0	0   
0	0
//Historical event: time, source, sink, migrants, new deme size, new growth rate, new migration matrix
2 historical event
CHANGM 0 0 0 1 GR_POP1 1
TDIV 1 0 1 RSANC 0 2
//Number of independent chromosome
1 0
//Number of contiguous linkage blocks
1
//Per Block: Data type, No. of loci, Recombination rate to the right-side locus, plus optional parameters
FREQ 1 0 3.5e-9 OUTEXP
