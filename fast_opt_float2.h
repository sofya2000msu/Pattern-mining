/*
This it the code for fast pattern mining algorithm. It works on time-series datasets and allow the uncertainty as well as temporal contraints.
It also has item-based contraints allowing to mine only for pattern which have items from different categories/groups.
The main ideas are described in
[1] Sofya Titarenko, Valeriy Titarenko, Georgios Aivaliotis, Jan Palczewski, "Fast implementation of pattern mining algorithms with time stamp uncertainties and temporal constraints", submitted in Journal of big Data, 2019
[2]  Sofya Titarenko, Valeriy Titarenko, Georgios Aivaliotis, Jan Palczewski, "A constraint-based frequent pattern mining
algorithm and its optimisation for multicore systems", EMiT2019 proceedings
*/
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <time.h>


#define n32 32

#define nmem 50000

using namespace std;

class dataSpam{
public:
	dataSpam();
	~dataSpam();

	FILE *fin, *fout, *foi;
		

	int readSettings(); // read setting from the input file, such as input/output files, uncertainty interval value, number of openMp threads, etc..
	int readInput(); // reads the input file
	int sortByEvent();//sorting the events. See [1,2] for the details.
	
	int startProcessing();//the main function which is managing pattern mining
	int count1D();// function which counts unique events and keeps Index array to help sorting. See [1,2] for the details.
	int remove1D();// function which removes not frequent unique events. See [1,2] for the details.
	int searchIndex(int ndim);// funcion which searches for the location of sub-pattern found on the previous level. See [1,2] for the details.
	int searchIndexVec(int ndim, int *ivec);// funcion which searches for the location of sub-pattern found on the previous level. See [1,2] for the details.
		
	
	int findNewBitList(int ndim);//function which finds frequent patterns of length n
	int bit1D();//function which initialises bitmap ID lists.
	
	int printImage();// function which prints found frequent patterns.
	int bitCount(unsigned int *vec);// function which counts how many non zero bits in bitmap ID list
	int bitCountList(unsigned int *vec);// function which counts how many non zero bits in bitmap ID list


	unsigned int findBitList(unsigned int *vec, int tid);// function which finds a location of the record in 

	int printInfo(int ics);// function which prints the information about calculating time, threads, etc..
	
	bool isTrue_standard(int ncomb, unsigned int irec, int ithread);//function which does pattern search, general formulation. See [1,2] for the details.
	bool isTrue_ordered(int ncomb, unsigned int irec, int ithread);//function which does pattern search, taking prior assumption into account. See [1,2] for the details.
	bool isTrue_standard_TR(int ncomb, unsigned int irec, int ithread);//function which does pattern search, general formulation + time restriction. See [1,2] for the details.
	bool isTrue_ordered_TR(int ncomb, unsigned int irec, int ithread);//function which does pattern search, taking prior assumption into account + time restriction. See [1,2] for the details.
	
	
	int nthreads;// number of openMP threads
	int nRecords, nRecords32, nChunkRecords32;// number of records, 32 sized chunks of the records
		
	char inputFile[500], outputFile[500], outputFolder[500], prefix[100], outPrefix[500];

	int **plistOfPatterns, *indForCand;// list where each thread stores found frequent patterns; array for candidate patterns
	long long **pPrevIndex;// index related with the location of subpattern. See [1,2] for the details.
	unsigned int *iList, **pListIn;//related with location of the candidate
	bool *vIsCandOk;// stores information about if candidate pattern is frequent

	int ibeta, support1000, ncores;// stores input variables
	int nTotalEvents, nTotalEvents1D, nEvents, nEventsOriginal;// total amount of events, unique events, unique frequent events, etc..
	int n1D;
	int **pEvents, **pRecordStart;// pointers to the array of events and array which containts the starting location of a record
	float **pfStartDate, **pfEndDate;// pointers to the arrays of starting time/date and ending time/date
	int *vsort;
	int **pDistEventsPerRecord;// stores how much unique events per record
	int **pStartDistEvent, **pShiftDistEvent;//pointers to start of unique event and its end

	int *c1D, *index1D, *index1D_rev, *iorder;// arrays storing indexes of events, helps sorting

	int *vcCurrent, *vcStart, *vcEnd, *vcTest;//stores the locations of "current" index, its "starting" point and "ending" point. See in searching functin and [1] for the details.
	int **pcCurrentAll, **pcStartAll, **pcEndAll, **pcTestAll;// pointers to the "current", "starting", "ending" points for every thread.
	int *vt_top, *vt_bottom, *vt_diff;


	int maxGroup;// maximum length of the pattern
	long long maxLen;//number of frequent patterns you expect
	long long maxLen2;//related with the storing the locations of subpatterns
	long long curLen;//related with the position within the array of found frequent patterns
	long long *startGroup;//is "starting" position of the array of found frequent patterns
	int *vCombination;//the array of found frequent patterns
	long long **pCombLoc;//frequent patterns within one thread
	unsigned int *vBits, *vbt1, *vbt2, *vbt3;//arrays helping in storing and processing bitmap lists
	unsigned int **pbt;//helping in storing and processing bitmap lists
	int **pindex;//pointers which help in finding and storing the location of subpatterns
		
	int iMinSupport;// minimum support
	char *vdata;//input stream from the file

	int isPrintPatterns, isTimeRestriction;// stores input information about prining patterns/time restriction
	float fTimeRestriction, fbeta;//stores input information about prining patterns/time restriction
	int isOrdered;//stores input information version of the algorithm/FarPaM1 or FaRPaM2

	float support;// minimum support
};
