#include "fast_opt_float2.h"
#include "memory.h"

dataSpam::dataSpam(){
	nRecords = 0;
	
	pEvents = nullptr;
	pRecordStart = nullptr;
	pfEndDate = nullptr;
	pfStartDate = nullptr;
	c1D = nullptr;
	index1D = nullptr;
	vsort = nullptr;

	vBits = nullptr;
	vCombination = nullptr;

	fout = nullptr;
	fin = nullptr;
	foi = nullptr;

}

dataSpam::~dataSpam(){
	if (fout != nullptr){
		fclose(fout); fout = nullptr;
	}
	if (foi != nullptr){
		fclose(foi); foi = nullptr;
	}
	if (fin != nullptr){
		fclose(fin); fin = nullptr;
	}
	if (vdata != nullptr){
		free(vdata); vdata = nullptr;
	}
	if (vBits != nullptr){
		free(vBits); vBits = nullptr;
	}
	if (pEvents != nullptr){
		for (int i = 0; i < nthreads; i++){
			free(pEvents[i]);
		}
		free(pEvents); pEvents = nullptr;
	}
	if (pRecordStart != nullptr){
		for (int i = 0; i < nthreads; i++){
			free(pRecordStart[i]);
		}
		free(pRecordStart); pRecordStart = nullptr;
	}
	if (pfStartDate != nullptr){
		for (int i = 0; i < nthreads; i++){
			free(pfStartDate[i]);
		}
		free(pfStartDate); pfStartDate = nullptr;
	}
	if (pfEndDate != nullptr){
		for (int i = 0; i < nthreads; i++){
			free(pfEndDate[i]);
		}
		free(pfEndDate); pfEndDate = nullptr;
	}
	
	if (c1D != nullptr){
		free(c1D); c1D = nullptr;
	}
	
	if (index1D != nullptr){
		free(index1D); index1D = nullptr;
	}
	if (vsort != nullptr){
		free(vsort); vsort = nullptr;
	}
		
}


int dataSpam::readInput(){
	int kr;
	fin = fopen(inputFile, "rb");
	if (fin == nullptr) return -1;
	_fseeki64(fin, 0, SEEK_END);

	long long fsize;
	fsize = _ftelli64(fin);


	vdata = (char*)malloc(sizeof(char)*fsize);
	rewind(fin);
	fread(vdata, sizeof(char), fsize, fin);
	fclose(fin); fin = nullptr;


	nTotalEvents = 0;
	nRecords = 0;
	for (int i = 0; i < fsize; i++){
		if (vdata[i] == '\t' || vdata[i] == '\n') nTotalEvents++;
		if (vdata[i] == '\n') nRecords++;
	}

	fprintf(foi,"Records\t%i\n", nRecords);
	printf("Records: %i\n", nRecords);

	nChunkRecords32 = nRecords / n32;
	if (nRecords % n32 > 0) nChunkRecords32++;
	nRecords32 = n32 * nChunkRecords32;

	printf("Chunks (records): %i\n", nChunkRecords32);
	printf("Events: %i\n", nTotalEvents);
	fprintf(foi,"Chunks\t%i\n", nChunkRecords32);
	fprintf(foi,"Events\t%i\n", nTotalEvents);
	vCombination = (int *)malloc(size_t(sizeof(int)) * size_t(maxGroup) * size_t(maxLen));
	if (vCombination == nullptr){
		printf("Error: memory combination\n");
		return -1;
	}
	vIsCandOk = (bool *)malloc(sizeof(bool)* size_t(maxLen));
	startGroup = (long long *)malloc(sizeof(long long) * (maxGroup + 1));
	
	plistOfPatterns = (int**)malloc(sizeof(int*)*nthreads);
	pPrevIndex = (long long**)malloc(sizeof(long long*)*nthreads);
	pCombLoc = (long long**)malloc(sizeof(long long*)*nthreads);
	indForCand = (int*)malloc(sizeof(int)*nthreads);
	

	vt_top = (int *)malloc(sizeof(int)* maxGroup);
	vt_bottom = (int *)malloc(sizeof(int)* maxGroup);
	vt_diff = (int *)malloc(sizeof(int)* maxGroup);

	vcCurrent = (int *)malloc(sizeof(int) * maxGroup);
	vcStart = (int *)malloc(sizeof(int) * maxGroup);
	vcEnd = (int *)malloc(sizeof(int) * maxGroup);
	vcTest = (int *)malloc(sizeof(int)* maxGroup);

	pindex = (int **)malloc(sizeof(int*)* nthreads);
	pcCurrentAll = (int **)malloc(sizeof(int*)* nthreads);
	pcStartAll = (int **)malloc(sizeof(int*)* nthreads);
	pcEndAll = (int **)malloc(sizeof(int*)* nthreads);
	pcTestAll = (int **)malloc(sizeof(int*)* nthreads);

	for (int i = 0; i < nthreads; i++){
		pPrevIndex[i] = (long long *)malloc(size_t(sizeof(long long))* size_t(maxLen2));
		if (pPrevIndex[i] == nullptr){
			printf("Error: memory pPrev (%u)\n", size_t(sizeof(long long))* size_t(maxLen2)/(1024*1024));
			return -1;
		}
		pCombLoc[i] = (long long *)malloc(size_t(sizeof(long long))* size_t(maxLen2));
		if (pCombLoc[i] == nullptr){
			printf("Error: memory pComb (%u)\n", size_t(sizeof(long long))* size_t(maxLen2)/ (1024 * 1024));
			return -1;
		}
		plistOfPatterns[i] = (int*)malloc(sizeof(int)*maxGroup*maxGroup);
		pindex[i] = (int *)malloc(sizeof(int)* nmem);
		pcCurrentAll[i] = (int *)malloc(sizeof(int)* nmem);
		pcStartAll[i] = (int *)malloc(sizeof(int)* nmem);
		pcEndAll[i] = (int *)malloc(sizeof(int)* nmem);
		pcTestAll[i] = (int *)malloc(sizeof(int)* nmem);
	}

	vBits = (unsigned int *)malloc(size_t(sizeof(unsigned int)) * size_t(maxLen)*size_t(nChunkRecords32));
	pbt = (unsigned int **)malloc(sizeof(unsigned int*)* nthreads);
	vbt1 = (unsigned int *)malloc(sizeof(unsigned int) * nChunkRecords32);
	vbt2 = (unsigned int *)malloc(sizeof(unsigned int) * nChunkRecords32);
	vbt3 = (unsigned int *)malloc(sizeof(unsigned int) * nChunkRecords32);
	curLen = 0;

	if (vBits == nullptr) {
		printf("Error: memory Bits\n");
		return -1;
	}

	iList = (unsigned int *)malloc(sizeof(unsigned int)*(nRecords32 + 1));
	pListIn = (unsigned int **)malloc(sizeof(unsigned int*)*nthreads);
		
	pRecordStart = (int **)malloc(sizeof(int*)*nthreads);
	pEvents = (int **)malloc(sizeof(int*)*nthreads);
	pfStartDate = (float **)malloc(sizeof(float*)*nthreads);
	pfEndDate = (float **)malloc(sizeof(float*)*nthreads);

	for (int i = 0; i < nthreads; i++){
		pbt[i] = (unsigned int *)malloc(sizeof(unsigned int)* nChunkRecords32);
		pRecordStart[i] = (int *)malloc(sizeof(int)*(nRecords32 + 1));
		pEvents[i] = (int *)malloc(sizeof(int)*nTotalEvents);
		pfStartDate[i] = (float *)malloc(sizeof(float)*nTotalEvents);
		pfEndDate[i] = (float *)malloc(sizeof(float)*nTotalEvents);
		pListIn[i] = (unsigned int *)malloc(sizeof(unsigned int)*nRecords32);
	}

	int *vRecordStart, *vEvents;
	float *vStartDate, *vEndDate;

	vRecordStart = pRecordStart[0];
	vEvents = pEvents[0];
	vStartDate = pfStartDate[0];
	vEndDate = pfEndDate[0];

	vRecordStart[0] = 0;
	nRecords = 0;
	nTotalEvents = 0;
		
	kr = 0;

	nEventsOriginal = 0;

	bool nw;
	
	for (kr = 0; kr < fsize; kr++){
		vEvents[nTotalEvents] = atoi(vdata + kr);
		while (vdata[kr] != ','){
			kr++;
		}
		kr++;

		vStartDate[nTotalEvents] = (float)atof(vdata + kr) - fbeta;
		while (vdata[kr] != ','){
			kr++;
		}
		kr++;

		vEndDate[nTotalEvents] = vStartDate[nTotalEvents] + fbeta;

/*		vEndDate[nTotalEvents] = (float)atof(vdata + kr);
		while (vdata[kr] != '\t'){
			kr++;
		}
		kr++;
		
		nEventsOriginal = __max(nEventsOriginal, vEvents[nTotalEvents]);*/
		if (nEventsOriginal < vEvents[nTotalEvents]) nEventsOriginal = vEvents[nTotalEvents];

		nTotalEvents++;
		nw = false;
		
		while (!(vdata[kr] >= '0' && vdata[kr] <= '9')){
			kr++;
			if (kr >= fsize) break;
			nw = true;
		}
		if (nw){
			nRecords++;
			vRecordStart[nRecords] = nTotalEvents;
		}
		kr--;
	}
	nEventsOriginal++;

	for (int i = nRecords + 1; i <= nRecords32; i++){
		vRecordStart[i] = vRecordStart[nRecords];
	}
		

	printf("Events (original): %i\n", nEventsOriginal);

	free(vdata); vdata = nullptr;

	iMinSupport = ceil(support * nRecords);
	printf("RecordsNew: %i\n", nRecords);
	printf("MinSupport: %i\n", iMinSupport);

	fprintf(foi, "RecordsNew\t%i\n", nRecords);
	fprintf(foi, "MinSupport\t%i\n", iMinSupport);
		
	return 0;
}





int dataSpam::sortByEvent(){

	bool t;
	
	float fval1, fval2;
	int val3;

	int *vEvents, *vRecordStart;
	float *vStartDate, *vEndDate;

	vEvents = pEvents[0];
	vStartDate = pfStartDate[0];
	vEndDate = pfEndDate[0];
	vRecordStart = pRecordStart[0];

	bool iy;

	for (int i = 0; i < nRecords; i++){
		do{
			t = false;
			for (int j = vRecordStart[i]; j < vRecordStart[i + 1] - 1; j++){
				if (vEvents[j] < vEvents[j + 1])continue;
				iy = false;
				if (vEvents[j] > vEvents[j + 1]){
					iy = true;
				}else{
					if (vStartDate[j] < vStartDate[j + 1])continue;
					if (vStartDate[j] > vStartDate[j + 1]){
						iy = true;
					}
					else{
						if (vEndDate[j] < vEndDate[j + 1]) continue;
					}
				}
				//if (vEvents[j] == vEvents[j + 1] && vStartDate[j] < vStartDate[j + 1])continue;
				//|| (vEvents[j] == vEvents[j + 1] && vStartDate[j] == vStartDate[j + 1] && vEndDate[j] > vEndDate[j + 1])){
				//if (vEvents[j] > vEvents[j + 1] || (vEvents[j] == vEvents[j + 1] && vStartDate[j] > vStartDate[j + 1]) || (vEvents[j] == vEvents[j + 1] && vStartDate[j] == vStartDate[j + 1] && vEndDate[j] > vEndDate[j + 1])){
					t = true;
					fval1 = vStartDate[j];
					fval2 = vEndDate[j];
					val3 = vEvents[j];
					vStartDate[j] = vStartDate[j + 1];
					vEndDate[j] = vEndDate[j + 1];
					vEvents[j] = vEvents[j + 1];
					vStartDate[j + 1] = fval1;
					vEndDate[j + 1] = fval2;
					vEvents[j + 1] = val3;
				//}
			}
		} while (t);
		
		fval1 = vStartDate[vRecordStart[i]];

		for (int j = vRecordStart[i] + 1; j < vRecordStart[i + 1]; j++){
			fval1 = __min(vStartDate[j], fval1);
		}

		for (int j = vRecordStart[i]; j < vRecordStart[i + 1]; j++){
			vStartDate[j] -= fval1;
			vEndDate[j] -= fval1;
		}


	}
	

	
	return 0;
}


int dataSpam::count1D(){
	int ival;

	int *vEvents, *vRecordStart;
	float *vStartDate, *vEndDate;

	vEvents = pEvents[0];
	vStartDate = pfStartDate[0];
	vEndDate = pfEndDate[0];
	vRecordStart = pRecordStart[0];

	c1D = (int *)malloc(sizeof(int)* nEventsOriginal);
	index1D = (int *)malloc(sizeof(int)* nEventsOriginal);
	index1D_rev = (int *)malloc(sizeof(int)* nEventsOriginal);
	for (int i = 0; i < nEventsOriginal; i++){
		c1D[i] = 0;
	}

	for (int i = 0; i < nRecords; i++){
		for (int j = vRecordStart[i]; j < vRecordStart[i + 1]; j++){
			if (j == vRecordStart[i]){
				c1D[vEvents[j]]++;
			}else if(vEvents[j] != vEvents[j-1]){
				c1D[vEvents[j]]++;
			}
		}
	}


	nEvents = 0;

	for (int i = 0; i < nEventsOriginal; i++){
		if (c1D[i] >= iMinSupport) nEvents++;
		index1D[i] = i;
	}

	printf("Events with min support: %i\n", nEvents);
	

	bool t;
	do{
		t = false;
		for (int j = 0; j < nEventsOriginal - 1; j++){
			if (c1D[index1D[j]] < c1D[index1D[j + 1]]){
				t = true;
				ival = index1D[j];
				index1D[j] = index1D[j + 1];
				index1D[j + 1] = ival;
			}
		}
	} while (t);

	for (int i = 0; i < nEventsOriginal; i++){
		index1D_rev[index1D[i]] = i;
	}


	return 0;
}

int dataSpam::remove1D(){
	int ic, vst, minDate, u, nev;
	int mtot;
	float fval1, fval2;
	int val3;
	int *nDistEventsPerRecord, *vStartDistEvent, *vShiftDistEvent;
	int *nDistEventsPerRecord_loc, *vStartDistEvent_loc, *vShiftDistEvent_loc;
	float *vStartDate, *vEndDate;
	int *vEvents, *vRecordStart;
	float *vStartDate_loc, *vEndDate_loc;
	int *vEvents_loc, *vRecordStart_loc;

	bool t;

	mtot = 0;

	n1D = 0;
	for (int i = 0; i < nEventsOriginal; i++){
		index1D_rev[index1D[i]] = i;
		if (c1D[i] >= iMinSupport){
			n1D++;
			mtot += c1D[i];
		}
	}

	mtot += 2 * nRecords32;

	pDistEventsPerRecord = (int **)malloc(sizeof(int*)* nthreads);
	pStartDistEvent = (int **)malloc(sizeof(int*)* nthreads);
	pShiftDistEvent = (int **)malloc(sizeof(int*)* nthreads);
	
	for (int i = 0; i < nthreads; i++){
		pDistEventsPerRecord[i] = (int *)malloc(sizeof(int)* nRecords32);
		pStartDistEvent[i] = (int *)malloc(sizeof(int)* nRecords32);
	}
	for (int i = 0; i < nthreads; i++){
		pShiftDistEvent[i] = (int *)malloc(sizeof(int)* mtot);
	}

	nDistEventsPerRecord = pDistEventsPerRecord[0];
	vStartDistEvent = pStartDistEvent[0];
		
	vStartDate = pfStartDate[0];
	vEndDate = pfEndDate[0];
	vEvents = pEvents[0];
	vRecordStart = pRecordStart[0];
	vShiftDistEvent = pShiftDistEvent[0];

	u = 0;
	
	for (int i = 0; i < nRecords32; i++){
		ic = 0;
		vst = vRecordStart[i];
		t = true;
		for (int j = vRecordStart[i]; j < vRecordStart[i + 1]; j++){
			if (c1D[vEvents[j]] < iMinSupport) continue;
			
			vEvents[vst + ic] = index1D_rev[vEvents[j]];
	
			vStartDate[vst + ic] = vStartDate[j];
			vEndDate[vst + ic] = vEndDate[j];
			ic++;
		}



		do{
			t = false;
			for (int j = vRecordStart[i]; j < vRecordStart[i] + ic - 1; j++){
				if (vEvents[j] > vEvents[j + 1] || (vEvents[j] == vEvents[j + 1] && vStartDate[j] > vStartDate[j + 1]) || (vEvents[j] == vEvents[j + 1] && vStartDate[j] == vStartDate[j + 1] && vEndDate[j] > vEndDate[j + 1])){
					t = true;
					fval1 = vStartDate[j];
					fval2 = vEndDate[j];
					val3 = vEvents[j];
					vStartDate[j] = vStartDate[j + 1];
					vEndDate[j] = vEndDate[j + 1];
					vEvents[j] = vEvents[j + 1];
					vStartDate[j + 1] = fval1;
					vEndDate[j + 1] = fval2;
					vEvents[j + 1] = val3;
				}
			}
		} while (t);
		
		
		
		vStartDistEvent[i] = u;
		if (ic == 0){
			nDistEventsPerRecord[i] = 0;
			continue;
		}
		minDate = vStartDate[vst];
		for (int j = 1; j < ic; j++){
			minDate = __min(minDate, vStartDate[vst + j]);
		}
		for (int j = 0; j < ic; j++){
			vStartDate[vst + j] -= minDate;
			vEndDate[vst + j] -= minDate;
		}
		
		nev = 1;

		vShiftDistEvent[u] = 0;
		for (int j = 1; j < ic; j++){
			if (vEvents[vst + j] == vEvents[vst + j - 1]) continue;
			
			
			u++;
			vShiftDistEvent[u] = j;
			vEvents[vst + nev] = vEvents[vst + j];
			nev++;
		}
		
		u++;
		vShiftDistEvent[u] = ic;
		u++;
			

		nDistEventsPerRecord[i] = nev;
				
		
	}

	for (int k = 1; k < nthreads; k++){
		vStartDate_loc = pfStartDate[k];
		vEndDate_loc = pfEndDate[k];
		vEvents_loc = pEvents[k];
		vRecordStart_loc = pRecordStart[k];
		nDistEventsPerRecord_loc = pDistEventsPerRecord[k];
		vShiftDistEvent_loc = pShiftDistEvent[k];
		vStartDistEvent_loc = pStartDistEvent[k];
		for (int j = 0; j < nTotalEvents; j++){
			vEvents_loc[j] = vEvents[j];
			vStartDate_loc[j] = vStartDate[j];
			vEndDate_loc[j] = vEndDate[j];
		}
		for (int j = 0; j <= nRecords32; j++){
			vRecordStart_loc[j] = vRecordStart[j];
		}
		for (int j = 0; j < nRecords32; j++){
			nDistEventsPerRecord_loc[j] = nDistEventsPerRecord[j];
			vStartDistEvent_loc[j] = vStartDistEvent[j];
		}
		for (int j = 0; j < mtot; j++){
			vShiftDistEvent_loc[j] = vShiftDistEvent[j];
		}
	}
		
		
	return 0;
}

int dataSpam::bit1D(){
	unsigned int bit;
	int irec, ievent;
	if (nEvents > maxLen) return -1;
	int *vEvents, *vRecordStart;
	int *nDistEventsPerRecord;

	nDistEventsPerRecord = pDistEventsPerRecord[0];
	vEvents = pEvents[0];
	vRecordStart = pRecordStart[0];

	startGroup[0] = 0;
	startGroup[1] = nEvents;
	for (int i = 0; i < nEvents * maxGroup; i++){
		vCombination[i] = -1;
	}
	for (int i = 0; i < nEvents; i++){
		vCombination[i*maxGroup] = i;
	}
	for (int i = 0; i < nChunkRecords32*nEvents; i++){
		vBits[i] = 0;
	}
	
	for (int ic = 0; ic < nChunkRecords32; ic++){
		for (int j = 0; j < n32; j++){
			bit = 1 << j;
			irec = ic*n32 + j;
			for (int k = 0; k < nDistEventsPerRecord[irec]; k++){
				ievent = vEvents[vRecordStart[irec] + k];
				vBits[ievent * nChunkRecords32 + ic] |= bit;
			}
		}
	}

	
	curLen = nEvents;
//	printInfo(1);
	return 0;
}




bool dataSpam::isTrue_ordered(int ncomb, unsigned int irec, int ithread){
	bool t, tf;
	int n, p, nd, istartp, indvc, SD, pS, pD, oot, i1, i2;
	int *vEvents_loc;
	float *vSD, *vED, timeStart, timeEnd;
	int *vcCurrent_loc, *vcStart_loc, *vcEnd_loc, *vcTest_loc, *vindex_loc;

	int *vEvents, *vRecordStart;
	float *vStartDate, *vEndDate;
	int *nDistEventsPerRecord;
	int *vStartDistEvent, *vShiftDistEvent;
	
	n = 0;
	
	nDistEventsPerRecord = pDistEventsPerRecord[ithread];
	vEvents = pEvents[ithread];
	vRecordStart = pRecordStart[ithread];
	vStartDate = pfStartDate[ithread];
	vEndDate = pfEndDate[ithread];
	vStartDistEvent = pStartDistEvent[ithread];
	vShiftDistEvent = pShiftDistEvent[ithread];


	nd = nDistEventsPerRecord[irec];

	vEvents_loc = vEvents + vRecordStart[irec];
	vSD = vStartDate + vRecordStart[irec];
	vED = vEndDate + vRecordStart[irec];
	
	istartp = vStartDistEvent[irec];
	
	vcCurrent_loc = pcCurrentAll[ithread];
	vcStart_loc = pcStartAll[ithread];
	vcEnd_loc = pcEndAll[ithread];
	vindex_loc = pindex[ithread];

	int indPrev[100];

	for (int i = 0; i < ncomb; i++){
		indPrev[i] = -1;
	}
	
	for (int i = 0; i < ncomb; i++){
		tf = false;
		for (int j = 0; j < nd; j++){
			if (vindex_loc[i] != vEvents_loc[j])continue;
			tf = true;
			vcEnd_loc[i] = vShiftDistEvent[istartp + j + 1];
			vcStart_loc[i] = vShiftDistEvent[istartp + j];
			break;
		}
		if (!tf) return false;
	}

	

	for (int i = 0; i < ncomb; i++){
		vcCurrent_loc[i] = vcStart_loc[i];
	}

	for (int i = ncomb - 1; i > 0; i--){
		for (int j = i - 1; j >= 0; j--){
			if (vcCurrent_loc[i] == vcCurrent_loc[j]){
				indPrev[i] = j;
				break;
			}
		}
	}

	for (int i = 0; i < ncomb; i++){
		for (int j = i + 1; j < ncomb; j++){
			if (vcCurrent_loc[i] == vcCurrent_loc[j]){
				vcCurrent_loc[j]++;
				if (vcCurrent_loc[j] == vcEnd_loc[j]) return false;
			}
		}
	}

	timeStart = vSD[vcCurrent_loc[0]];
	

	for (int k = 1; k < ncomb; k++){
		t = false;
		if (indPrev[k] >= 0){
			vcCurrent_loc[k] = vcCurrent_loc[indPrev[k]] + 1;
		}
		for (; vcCurrent_loc[k] < vcEnd_loc[k]; vcCurrent_loc[k]++){
			indvc = vcCurrent_loc[k];
	
			if (vED[indvc] < timeStart) continue;
			timeStart = __max(timeStart, vSD[indvc]);
			t = true;
			break;
		}
		if (!t) return false;

	}

	return true;

	
}



bool dataSpam::isTrue_ordered_TR(int ncomb, unsigned int irec, int ithread){
	bool t, tf;
	int n, p, nd, istartp, indvc, SD, pS, pD, oot, i1, i2;
	int *vEvents_loc;
	float *vSD, *vED, timeStart, timeEnd;
	int *vcCurrent_loc, *vcStart_loc, *vcEnd_loc, *vcTest_loc, *vindex_loc;

	int *vEvents, *vRecordStart;
	float *vStartDate, *vEndDate;
	int *nDistEventsPerRecord;
	int *vStartDistEvent, *vShiftDistEvent;

	n = 0;

	nDistEventsPerRecord = pDistEventsPerRecord[ithread];
	vEvents = pEvents[ithread];
	vRecordStart = pRecordStart[ithread];
	vStartDate = pfStartDate[ithread];
	vEndDate = pfEndDate[ithread];
	vStartDistEvent = pStartDistEvent[ithread];
	vShiftDistEvent = pShiftDistEvent[ithread];


	nd = nDistEventsPerRecord[irec];

	vEvents_loc = vEvents + vRecordStart[irec];
	vSD = vStartDate + vRecordStart[irec];
	vED = vEndDate + vRecordStart[irec];

	istartp = vStartDistEvent[irec];

	vcCurrent_loc = pcCurrentAll[ithread];
	vcStart_loc = pcStartAll[ithread];
	vcEnd_loc = pcEndAll[ithread];
	vindex_loc = pindex[ithread];

	int indPrev[100];

	for (int i = 0; i < ncomb; i++){
		indPrev[i] = -1;
	}

	for (int i = 0; i < ncomb; i++){
		tf = false;
		for (int j = 0; j < nd; j++){
			if (vindex_loc[i] != vEvents_loc[j])continue;
			tf = true;
			vcEnd_loc[i] = vShiftDistEvent[istartp + j + 1];
			vcStart_loc[i] = vShiftDistEvent[istartp + j];
			break;
		}
		if (!tf) return false;
	}



	for (int i = 0; i < ncomb; i++){
		vcCurrent_loc[i] = vcStart_loc[i];
	}

	for (int i = ncomb - 1; i > 0; i--){
		for (int j = i - 1; j >= 0; j--){
			if (vcCurrent_loc[i] == vcCurrent_loc[j]){
				indPrev[i] = j;
				break;
			}
		}
	}

	for (int i = 0; i < ncomb; i++){
		for (int j = i + 1; j < ncomb; j++){
			if (vcCurrent_loc[i] == vcCurrent_loc[j]){
				vcCurrent_loc[j]++;
				if (vcCurrent_loc[j] == vcEnd_loc[j]) return false;
			}
		}
	}

	float time1, time2;

	for (; vcCurrent_loc[0] < vcEnd_loc[0]; vcCurrent_loc[0]++){
		time1 = vSD[vcCurrent_loc[0]];
		time2 = vED[vcCurrent_loc[0]];

		for (int i = 1; i < ncomb; i++){
			vcCurrent_loc[i] = vcStart_loc[i];
		}

		for (int k = 1; k < ncomb; k++){
			t = false;
			if (indPrev[k] >= 0){
				vcCurrent_loc[k] = vcCurrent_loc[indPrev[k]] + 1;
			}
			for (; vcCurrent_loc[k] < vcEnd_loc[k]; vcCurrent_loc[k]++){
				indvc = vcCurrent_loc[k];
				if (vED[indvc] < time1) continue;
				time1 = __max(time1, vSD[indvc]);
				time2 = __min(time2, vED[indvc]);
				t = true;
				break;
			}
			if (time1 > time2 + fTimeRestriction) break;
			if (!t) return false;
			if (k == ncomb - 1) return true;

		}
	}

	return false;


}


bool dataSpam::isTrue_standard(int ncomb, unsigned int irec, int ithread){
	bool t, tf;
	int n, p, nd, istartp, indvc, SD, pS, pD, oot, i1, i2;
	int *vEvents_loc;
	float *vSD, *vED, timeStart;
	int *vcCurrent_loc, *vcStart_loc, *vcEnd_loc, *vcTest_loc, *vindex_loc;
	int indPrev[50], indNext[50];
	int iu, iw, indRes[100], ic, nn;
	int valStart, valEnd;

	int *vEvents, *vRecordStart;
	float *vStartDate, *vEndDate;
	int *nDistEventsPerRecord;
	int *vStartDistEvent, *vShiftDistEvent;

	n = 0;

	nDistEventsPerRecord = pDistEventsPerRecord[ithread];
	vEvents = pEvents[ithread];
	vRecordStart = pRecordStart[ithread];
	vStartDate = pfStartDate[ithread];
	vEndDate = pfEndDate[ithread];
	vStartDistEvent = pStartDistEvent[ithread];
	vShiftDistEvent = pShiftDistEvent[ithread];


	nd = nDistEventsPerRecord[irec];

	vEvents_loc = vEvents + vRecordStart[irec];
	vSD = vStartDate + vRecordStart[irec];
	vED = vEndDate + vRecordStart[irec];

	istartp = vStartDistEvent[irec];

	vcCurrent_loc = pcCurrentAll[ithread];
	vcStart_loc = pcStartAll[ithread];
	vcEnd_loc = pcEndAll[ithread];
	vindex_loc = pindex[ithread];

	for (int i = 0; i < ncomb; i++){
		tf = false;
		for (int j = 0; j < nd; j++){
			if (vindex_loc[i] != vEvents_loc[j])continue;
			tf = true;
			vcEnd_loc[i] = vShiftDistEvent[istartp + j + 1];
			vcStart_loc[i] = vShiftDistEvent[istartp + j];
			break;
		}
		if (!tf) return false;
	}



	for (int i = 0; i < ncomb; i++){
		vcCurrent_loc[i] = vcStart_loc[i];
	}

		
	for (int i = 0; i < ncomb; i++){
		indPrev[i] = -1;
		indNext[i] = -1;
	}

	for (int i = 1; i < ncomb; i++){
		for (int j = i-1; j >= 0; j--){
			if (vcStart_loc[i] == vcStart_loc[j]){
				indPrev[i] = j;
				indNext[j] = i;
				break;
	
			}
		}
	}

	for (int i = 1; i < ncomb; i++){
		if (indPrev[i] == -1) continue;
		vcCurrent_loc[i] = vcCurrent_loc[indPrev[i]] + 1;
		if (vcCurrent_loc[i] >= vcEnd_loc[i]) return false;
	}
		
	
	do{

		timeStart = vSD[vcCurrent_loc[0]];
		oot = -1;
		for (int i = 1; i < ncomb; i++){
			indvc = vcCurrent_loc[i];
			if (timeStart > vED[indvc]){
				oot = i;
				break;
			}
			timeStart = __max(timeStart, vSD[indvc]);
		}
		if (oot == -1) return true;
				
		p = oot;

		valStart = vcStart_loc[p];
		valEnd = vcEnd_loc[p];

		iu = indPrev[p];
		ic = 0;
		nn = 0;
		while (iu >= 0){
			indRes[ic] = vcCurrent_loc[iu];
			ic++;
			iu = indPrev[iu];
		}

		do{					
			t = false; 
			vcCurrent_loc[p]++;
			if (vcCurrent_loc[p] == valEnd){
				iw = indPrev[p];
				if (iw == -1) return false;
				p = iw;
				nn++;
				t = true;
			}else{
				for (int i = nn; i < ic; i++){
					if (vcCurrent_loc[p] == indRes[i]){
						t = true;
						break;
					}
				}
			}
		} while (t);

		iu = indPrev[p];
		ic = 1;
		indRes[0] = vcCurrent_loc[p];
		while(iu >= 0){
			indRes[ic] = vcCurrent_loc[iu];
			ic++;
			iu = indPrev[iu];
		}

		iu = indNext[p];
		while (iu > 0){
			vcCurrent_loc[iu] = valStart;
			do{
				t = false;
				for (int i = 0; i < ic; i++){
					if (vcCurrent_loc[iu] == indRes[i]){
						vcCurrent_loc[iu]++;
						t = true;
						break;
					}
				}
			} while (t);
			indRes[ic] = vcCurrent_loc[iu];
			ic++;
			iu = indNext[iu];
		}
			


	} while (true);
	

	return false;
}




bool dataSpam::isTrue_standard_TR(int ncomb, unsigned int irec, int ithread){
	bool t, tf;
	int n, p, nd, istartp, indvc, SD, pS, pD, oot, i1, i2;
	int *vEvents_loc;
	float *vSD, *vED, time1, time2;
	int *vcCurrent_loc, *vcStart_loc, *vcEnd_loc, *vcTest_loc, *vindex_loc;
	int iq;
	int indPrev[100], indNext[100], indRes[100];

	int *vEvents, *vRecordStart;
	float *vStartDate, *vEndDate;
	int *nDistEventsPerRecord;
	int *vStartDistEvent, *vShiftDistEvent;

	n = 0;

	nDistEventsPerRecord = pDistEventsPerRecord[ithread];
	vEvents = pEvents[ithread];
	vRecordStart = pRecordStart[ithread];
	vStartDate = pfStartDate[ithread];
	vEndDate = pfEndDate[ithread];
	vStartDistEvent = pStartDistEvent[ithread];
	vShiftDistEvent = pShiftDistEvent[ithread];


	nd = nDistEventsPerRecord[irec];

	vEvents_loc = vEvents + vRecordStart[irec];
	vSD = vStartDate + vRecordStart[irec];
	vED = vEndDate + vRecordStart[irec];

	istartp = vStartDistEvent[irec];

	vcCurrent_loc = pcCurrentAll[ithread];
	vcStart_loc = pcStartAll[ithread];
	vcEnd_loc = pcEndAll[ithread];
	vindex_loc = pindex[ithread];

	for (int i = 0; i < ncomb; i++){
		tf = false;
		for (int j = 0; j < nd; j++){
			if (vindex_loc[i] != vEvents_loc[j])continue;
			tf = true;
			vcEnd_loc[i] = vShiftDistEvent[istartp + j + 1];
			vcStart_loc[i] = vShiftDistEvent[istartp + j];
			break;
		}
		if (!tf) return false;
	}

	

	for (int i = 0; i < ncomb; i++){
		vcCurrent_loc[i] = vcStart_loc[i];
	}
	
	for (int i = 0; i < ncomb; i++){
		indPrev[i] = -1;
		indNext[i] = -1;
	}

	for (int i = 1; i < ncomb; i++){
		for (int j = i - 1; j >= 0; j--){
			if (vcStart_loc[i] == vcStart_loc[j]){
				indPrev[i] = j;
				indNext[j] = i;
				break;

			}
		}
	}

	for (int i = 1; i < ncomb; i++){
		if (indPrev[i] == -1) continue;
		vcCurrent_loc[i] = vcCurrent_loc[indPrev[i]] + 1;
		if (vcCurrent_loc[i] >= vcEnd_loc[i]) return false;
	}


	int iu, ic, nn, valStart, valEnd, iw;
	

	do{

		time1 = vSD[vcCurrent_loc[0]];
		time2 = vED[vcCurrent_loc[0]];
		oot = -1;
		iq = 1;
		for (int i = 1; i < ncomb; i++){
			indvc = vcCurrent_loc[i];
			if (time1 > vED[indvc]){
				oot = i;
				break;
			}
			time1 = __max(time1, vSD[indvc]);
			time2 = __min(time2, vED[indvc]);
			if (time1 > time2 + fTimeRestriction){
				iq = 0;
				break;
			}
		}
		if (oot == -1 && iq == 1){
			return true;
		}

		
		
	
		
		if (iq == 0){
			vcCurrent_loc[0]++;
			if (vcCurrent_loc[0] == vcEnd_loc[0]) return false;
			for (int i = 0; i < ncomb; i++){
				for (int j = i + 1; j < ncomb; j++){
					if (vcCurrent_loc[i] == vcCurrent_loc[j]){
						vcCurrent_loc[j]++;
						if (vcCurrent_loc[j] == vcEnd_loc[j]) return false;
					}
				}
			}
			continue;
		}

		
		p = oot;

		valStart = vcStart_loc[p];
		valEnd = vcEnd_loc[p];

		iu = indPrev[p];
		ic = 0;
		nn = 0;
		while (iu >= 0){
			indRes[ic] = vcCurrent_loc[iu];
			ic++;
			iu = indPrev[iu];
		}

		do{
			t = false;
			vcCurrent_loc[p]++;
			if (vcCurrent_loc[p] == valEnd){
				iw = indPrev[p];
				if (iw == -1) return false;
				p = iw;
				nn++;
				t = true;
			}else{
				for (int i = nn; i < ic; i++){
					if (vcCurrent_loc[p] == indRes[i]){
						t = true;
						break;
					}
				}
			}
		} while (t);

		iu = indPrev[p];
		ic = 1;
		indRes[0] = vcCurrent_loc[p];
		while(iu >= 0){
			indRes[ic] = vcCurrent_loc[iu];
			ic++;
			iu = indPrev[iu];
		}

		iu = indNext[p];
		while (iu > 0){
			vcCurrent_loc[iu] = valStart;
			do{
				t = false;
				for (int i = 0; i < ic; i++){
					if (vcCurrent_loc[iu] == indRes[i]){
						vcCurrent_loc[iu]++;
						t = true;
						break;
					}
				}
			} while (t);
			indRes[ic] = vcCurrent_loc[iu];
			ic++;
			iu = indNext[iu];
		}
		

/*		do{
			p = oot;
			

			do{
				vcCurrent_loc[p]++;

				if (vcCurrent_loc[p] < vcEnd_loc[p]) break;
				vcCurrent_loc[p] = vcStart_loc[p];
				p--;
			} while (p >= 0);

			if (p < 0) return false;

			t = false;

			for (int i = 0; i < oot; i++){
				int op;
				op = vcCurrent_loc[i];
				for (int j = i + 1; j <= oot; j++){
					if (op == vcCurrent_loc[j]){
						t = true;
						oot = j;
						break;
					}
				}
				if (t)break;
			}
			if (!t && oot < ncomb - 1){
				oot = ncomb - 1;

				for (int i = 0; i < oot; i++){
					int op;
					op = vcCurrent_loc[i];
					for (int j = i + 1; j <= oot; j++){
						if (op == vcCurrent_loc[j]){
							t = true;
							oot = j;
							break;
						}
					}
					if (t)break;
				}
			}
		} while (t);
*/

	} while (true);


	return false;
}




int dataSpam::searchIndex(int ndim){
	int index_1, index_2;
	int *vindex;

	vindex = pindex[0];
	index_1 = startGroup[ndim - 1];
	index_2 = startGroup[ndim];

	/*bool u;
	for (int i = index_1; i < index_2; i++){
	u = true;
	for (int j = 0; j < ndim; j++){
	u = u && (vCombination[i*maxGroup + j] == vindex[j]);
	}
	if (!u) continue;
	return i;
	}*/

	int ind1, ind2, ind3, ind3_old;
	int ig;

	ind1 = index_1;
	ind2 = index_2;
	ind3 = -1;

	do{
		ind3_old = ind3;
		ind3 = (ind1 + ind2) / 2;
		ig = 0;
		for (int j = 0; j < ndim; j++){
			ig = vCombination[ind3*maxGroup + j] - vindex[j];
			if (ig != 0)break;
		}
		if (ig == 0) return ind3;
		if (ig > 0){
			ind2 = ind3;
		}
		else{
			ind1 = ind3;
		}
	} while (ind3_old != ind3);

	/*bool u;
	for (int i = index_1; i < index_2; i++){
		u = true;
		for (int j = 0; j < ndim; j++){
			u = u && (vCombination[i*maxGroup + j] == vindex[j]);
		}
		if (!u) continue;
		return i;
	}*/

	return -1;
}



int dataSpam::searchIndexVec(int ndim, int *ivec){
	int index_1, index_2;
	int ind1, ind2, ind3, ind3_old;
	int ig;
	
	index_1 = startGroup[ndim - 1];
	index_2 = startGroup[ndim];
	
	ind1 = index_1;
	ind2 = index_2;
	ind3 = -1;

	do{
		ind3_old = ind3;
		ind3 = (ind1 + ind2) / 2;
		ig = 0;
		for (int j = 0; j < ndim; j++){
			ig = vCombination[ind3*maxGroup + j] - ivec[j];
			if (ig != 0)break;
		}
		if (ig == 0) return ind3;
		if (ig > 0){
			ind2 = ind3;
		}
		else{
			ind1 = ind3;
		}
	} while (ind3_old != ind3);
	
	return -1;
}



int dataSpam::findNewBitList(int ndim){
	clock_t begin_time, end_time;
	
	long long nCandidates, nFound;
	long long index_1, index_2;
	long long curLen_cand;
	long long indCP[100];
	long long uBest[100], *icl;
	int ibest;
	int isOut, isq;
	bool isBetter;

	int *vCombStart, *vPrevStart;
		
	unsigned int *vu1, *vu2;
	double gt;
	
	index_1 = startGroup[ndim - 2];
	index_2 = startGroup[ndim - 1];

	
	begin_time = clock();
		

#pragma omp parallel
	{
		int tid;
		int *listOfPatterns;
		long long *iCandPrevIndex, *iPrevIndex, *vCombLoc, *iCombLoc;
		int nrec;
		long long nc;
		unsigned int *pbt_loc;
		unsigned int* vBit_loc, *vBit_loc2;
		bool t;
		tid = omp_get_thread_num();
		nthreads = omp_get_num_threads();
		listOfPatterns = plistOfPatterns[tid];
		
		vCombLoc = pCombLoc[tid];
		iPrevIndex = pPrevIndex[tid];
		pbt_loc = pbt[tid];
		nc = 0;
#pragma omp for schedule(auto)
		for (long long i = index_1; i < index_2; i++){
			
			
			
			for (int j = 0; j < ndim - 1; j++){
				listOfPatterns[j] = vCombination[i*maxGroup + j];
			}
			for (int k = 0; k < ndim - 1; k++){
				for (int j = 0; j < k; j++){
					listOfPatterns[(k + 1)*maxGroup + j] = vCombination[i*maxGroup + j];
				}
				for (int j = k; j < ndim - 2; j++){
					listOfPatterns[(k + 1)*maxGroup + j] = vCombination[i*maxGroup + j + 1];
				}
			}

			for (int j = 0; j < nEvents; j++){
				iCombLoc = vCombLoc + nc*ndim;
				iCandPrevIndex = iPrevIndex + nc*ndim;
				iCandPrevIndex[0] = i;
				t = true;

				for (int k = 0; k < ndim - 1; k++){
					listOfPatterns[(k + 1)*maxGroup + ndim - 2] = j;
				}

				for (int k = 0; k < ndim - 1; k++){
					iCandPrevIndex[k + 1] = searchIndexVec(ndim - 1, listOfPatterns + (k + 1)*maxGroup);
					if (iCandPrevIndex[k + 1] < 0){
						t = false;
						break;
					}

				}
				if (!t) continue;

				

				vBit_loc = vBits + i*nChunkRecords32;
				vBit_loc2 = vBits + iCandPrevIndex[1] * nChunkRecords32;

				for (int k = 0; k < nChunkRecords32; k++){
					pbt_loc[k] = vBit_loc[k] & vBit_loc2[k];
				}

				nrec = bitCountList(pbt_loc);
				if (nrec < iMinSupport){
					t = false;
					continue;
				}

				for (int p = 2; p < ndim; p++){
					vBit_loc = vBits + iCandPrevIndex[p] * nChunkRecords32;
					for (int k = 0; k < nChunkRecords32; k++){
						pbt_loc[k] &= vBit_loc[k];
					}
					nrec = bitCountList(pbt_loc);
					if (nrec < iMinSupport){
						t = false;
						break;
					}
				}

				if (!t) continue;
							
				for (int k = 0; k < ndim - 1; k++){
					iCombLoc[k] = listOfPatterns[k];
				}

				iCombLoc[ndim - 1] = j;
								
				nc++;

			}
		}
		indForCand[tid] = nc;
	}

	nCandidates = 0;

	for (int i = 0; i < nthreads; i++){
		nCandidates += (long long)(indForCand[i]);
	}
	

	curLen_cand = curLen + nCandidates;

	for (long long i = curLen*maxGroup; i < curLen_cand*maxGroup; i++){
		vCombination[i] = -1;
	}

	for (int i = 0; i < nthreads; i++){
		indCP[i] = 0;
	}

	
	
	vCombStart = vCombination + curLen*maxGroup;
	vPrevStart = vCombination + (curLen + nCandidates)*maxGroup;
		

	for (long long k = 0; k < nCandidates; k++){
		ibest = -1;
		for (int i = 0; i < ndim; i++){
			uBest[i] = 99999999;
		}
		for (int j = 0; j < nthreads; j++){
			if (indCP[j] >= indForCand[j]) continue;
			icl = pCombLoc[j] + indCP[j] * ndim;
			isBetter = true;
			for (int ii = 0; ii < ndim; ii++){
				if (icl[ii] > uBest[ii]){
					isBetter = false;
					break;
				}
				if (icl[ii] < uBest[ii]){
					break;
				}
			}
			if (!isBetter)continue;
			ibest = j;
			for (int ii = 0; ii < ndim; ii++){
				uBest[ii] = icl[ii];
			}
		}
		for (int ii = 0; ii < ndim; ii++){
			vCombStart[k*maxGroup + ii] = uBest[ii];
			vPrevStart[k*ndim + ii] = *(pPrevIndex[ibest] + indCP[ibest] * ndim + ii);
		}
		for (int ii = ndim; ii < maxGroup; ii++){
			vCombStart[k*maxGroup + ii] = -1;
		}
		indCP[ibest]++;

	}

	if (isTimeRestriction){

		if (isOrdered == 1){

#pragma omp parallel
			{
				int tid;
				unsigned int iCand, aa, bb;
				unsigned int bit, uu, nrec_loc, *vbt_loc, *iList_loc, *iListOut_loc;
				unsigned int *pbt_loc;
				unsigned int* vBit_loc, *vBit_loc2;
				int *vindex;
				tid = omp_get_thread_num();
				vindex = pindex[tid];
				iList_loc = pListIn[tid];
				pbt_loc = pbt[tid];

				for (int i = 0; i < maxGroup; i++){
					vindex[i] = -1;
				}

#pragma omp for schedule(auto)	
				for (long long i = 0; i < nCandidates; i++){
					vIsCandOk[i] = true;
					vbt_loc = vBits + (i + curLen)*nChunkRecords32;

					vBit_loc = vBits + vPrevStart[i*ndim + 0] * nChunkRecords32;
					vBit_loc2 = vBits + vPrevStart[i*ndim + 1] * nChunkRecords32;

					for (int k = 0; k < nChunkRecords32; k++){
						vbt_loc[k] = vBit_loc[k] & vBit_loc2[k];
					}

					for (int p = 2; p < ndim; p++){
						vBit_loc = vBits + vPrevStart[i*ndim + p] * nChunkRecords32;
						for (int k = 0; k < nChunkRecords32; k++){
							vbt_loc[k] &= vBit_loc[k];
						}
					}

					nrec_loc = findBitList(vbt_loc, tid);

					for (int p = 0; p < maxGroup; p++){
						vindex[p] = vCombStart[i*maxGroup + p];
					}
					uu = 0;
					for (int k = 0; k < nrec_loc; k++){
						iCand = iList_loc[k];
						if (isTrue_ordered_TR(ndim, iCand, tid)) continue;
						bb = iCand%n32;
						iCand = iCand / n32;
						vbt_loc[iCand] ^= (1 << bb);
						uu++;
						if (nrec_loc - uu < iMinSupport){
							vIsCandOk[i] = false;
							break;
						}
					}

				}
			}
		}
		else{

#pragma omp parallel
			{
				int tid;
				unsigned int iCand, aa, bb;
				unsigned int bit, uu, nrec_loc, *vbt_loc, *iList_loc, *iListOut_loc;
				unsigned int *pbt_loc;
				unsigned int* vBit_loc, *vBit_loc2;
				int *vindex;
				tid = omp_get_thread_num();
				vindex = pindex[tid];
				iList_loc = pListIn[tid];
				pbt_loc = pbt[tid];

				for (int i = 0; i < maxGroup; i++){
					vindex[i] = -1;
				}

#pragma omp for schedule(auto)	
				for (long long i = 0; i < nCandidates; i++){
					vIsCandOk[i] = true;
					vbt_loc = vBits + (i + curLen)*nChunkRecords32;


					vBit_loc = vBits + vPrevStart[i*ndim + 0] * nChunkRecords32;
					vBit_loc2 = vBits + vPrevStart[i*ndim + 1] * nChunkRecords32;

					for (int k = 0; k < nChunkRecords32; k++){
						vbt_loc[k] = vBit_loc[k] & vBit_loc2[k];
					}


					for (int p = 2; p < ndim; p++){
						vBit_loc = vBits + vPrevStart[i*ndim + p] * nChunkRecords32;
						for (int k = 0; k < nChunkRecords32; k++){
							vbt_loc[k] &= vBit_loc[k];
						}
					}

					nrec_loc = findBitList(vbt_loc, tid);

					for (int p = 0; p < maxGroup; p++){
						vindex[p] = vCombStart[i*maxGroup + p];
					}
					uu = 0;
					for (int k = 0; k < nrec_loc; k++){
						iCand = iList_loc[k];
						if (isTrue_standard_TR(ndim, iCand, tid)) continue;
						bb = iCand%n32;
						iCand = iCand / n32;
						vbt_loc[iCand] ^= (1 << bb);
						uu++;
						if (nrec_loc - uu < iMinSupport){
							vIsCandOk[i] = false;
							break;
						}
					}

				}
			}
		}
	}
	else{
		if (isOrdered == 1){

#pragma omp parallel
			{
				int tid;
				unsigned int iCand, aa, bb;
				unsigned int bit, uu, nrec_loc, *vbt_loc, *iList_loc, *iListOut_loc;
				unsigned int *pbt_loc;
				unsigned int* vBit_loc, *vBit_loc2;
				int *vindex;
				tid = omp_get_thread_num();
				vindex = pindex[tid];
				iList_loc = pListIn[tid];
				pbt_loc = pbt[tid];

				for (int i = 0; i < maxGroup; i++){
					vindex[i] = -1;
				}

#pragma omp for schedule(auto)	
				for (long long i = 0; i < nCandidates; i++){
					vIsCandOk[i] = true;
					vbt_loc = vBits + (i + curLen)*nChunkRecords32;

					vBit_loc = vBits + vPrevStart[i*ndim + 0] * nChunkRecords32;
					vBit_loc2 = vBits + vPrevStart[i*ndim + 1] * nChunkRecords32;

					for (int k = 0; k < nChunkRecords32; k++){
						vbt_loc[k] = vBit_loc[k] & vBit_loc2[k];
					}

					for (int p = 2; p < ndim; p++){
						vBit_loc = vBits + vPrevStart[i*ndim + p] * nChunkRecords32;
						for (int k = 0; k < nChunkRecords32; k++){
							vbt_loc[k] &= vBit_loc[k];
						}
					}

					nrec_loc = findBitList(vbt_loc, tid);

					for (int p = 0; p < maxGroup; p++){
						vindex[p] = vCombStart[i*maxGroup + p];
					}
					uu = 0;
					for (int k = 0; k < nrec_loc; k++){
						iCand = iList_loc[k];
						if (isTrue_ordered(ndim, iCand, tid)) continue;
						bb = iCand%n32;
						iCand = iCand / n32;
						vbt_loc[iCand] ^= (1 << bb);
						uu++;
						if (nrec_loc - uu < iMinSupport){
							vIsCandOk[i] = false;
							break;
						}
					}

				}
			}
		}
		else{

#pragma omp parallel
			{
				int tid;
				unsigned int iCand, aa, bb;
				unsigned int bit, uu, nrec_loc, *vbt_loc, *iList_loc, *iListOut_loc;
				unsigned int *pbt_loc;
				unsigned int* vBit_loc, *vBit_loc2;
				int *vindex;
				tid = omp_get_thread_num();
				vindex = pindex[tid];
				iList_loc = pListIn[tid];
				pbt_loc = pbt[tid];

				for (int i = 0; i < maxGroup; i++){
					vindex[i] = -1;
				}

#pragma omp for schedule(auto)	
				for (long long i = 0; i < nCandidates; i++){
					vIsCandOk[i] = true;
					vbt_loc = vBits + (i + curLen)*nChunkRecords32;


					vBit_loc = vBits + vPrevStart[i*ndim + 0] * nChunkRecords32;
					vBit_loc2 = vBits + vPrevStart[i*ndim + 1] * nChunkRecords32;

					for (int k = 0; k < nChunkRecords32; k++){
						vbt_loc[k] = vBit_loc[k] & vBit_loc2[k];
					}


					for (int p = 2; p < ndim; p++){
						vBit_loc = vBits + vPrevStart[i*ndim + p] * nChunkRecords32;
						for (int k = 0; k < nChunkRecords32; k++){
							vbt_loc[k] &= vBit_loc[k];
						}
					}

					nrec_loc = findBitList(vbt_loc, tid);

					for (int p = 0; p < maxGroup; p++){
						vindex[p] = vCombStart[i*maxGroup + p];
					}
					uu = 0;
					for (int k = 0; k < nrec_loc; k++){
						iCand = iList_loc[k];
						if (isTrue_standard(ndim, iCand, tid)) continue;
						bb = iCand%n32;
						iCand = iCand / n32;
						vbt_loc[iCand] ^= (1 << bb);
						uu++;
						if (nrec_loc - uu < iMinSupport){
							vIsCandOk[i] = false;
							break;
						}
					}

				}
			}
		}
	
	}
		

	
	isOut = 0;
	isq = curLen;
		
	for (long long i = curLen; i < curLen_cand; i++){
		if (!vIsCandOk[i-isq]){
			isOut++;
			continue;
		}
		
		vu1 = vBits + i*nChunkRecords32;
		vu2 = vBits + curLen*nChunkRecords32;

		for (int k = 0; k < nChunkRecords32; k++){
			vu2[k] = vu1[k];
		}

		for (long long k = 0; k < maxGroup; k++){
			vCombination[curLen*maxGroup + k] = vCombination[i*maxGroup + k];
		}

		curLen++;
	}
	

	startGroup[ndim] = curLen;
	nFound = startGroup[ndim] - startGroup[ndim - 1];
		
				
	end_time = clock();
		

	gt = double(end_time - begin_time) / CLOCKS_PER_SEC;
	printf("Group (%i), Candidates: %i, Size: %i, time: %f\n", ndim, nCandidates, nFound, gt);
	
	if (nCandidates == 0){
		fprintf(foi, "%i\t%i\t%6.4f\t%i\t%13.4e\t%i\t%i\t%i\t%i\t%13.4e\t%13.4e\t%13.4e\t%13.4e\n", nthreads, support1000, fbeta, isTimeRestriction, fTimeRestriction, isOrdered, ndim, nCandidates, nFound, gt, 0.0, 0.0, 0.0);
	}
	else if (nFound == 0){
		fprintf(foi, "%i\t%i\t%6.4f\t%i\t%13.4e\t%i\t%i\t%i\t%i\t%13.4e\t%13.4e\t%13.4e\t%13.4e\n", nthreads, support1000, fbeta, isTimeRestriction, fTimeRestriction, isOrdered, ndim, nCandidates, nFound, gt, float(nFound) / float(nCandidates), gt / float(nCandidates), 0.0);
	}
	else{
		fprintf(foi, "%i\t%i\t%6.4f\t%i\t%13.4e\t%i\t%i\t%i\t%i\t%13.4e\t%13.4e\t%13.4e\t%13.4e\n", nthreads, support1000, fbeta, isTimeRestriction, fTimeRestriction, isOrdered, ndim, nCandidates, nFound, gt, float(nFound) / float(nCandidates), gt / float(nCandidates), gt / float(nFound));
	}

	return 0;
}


int dataSpam::printInfo(int ics){
	char fname[500];

	//return 0;
	sprintf(fname, "%s\\%i.txt", outputFolder, ics);
	
	FILE *fp;
	fp = fopen(fname, "w");
	if (fp == nullptr) return -1;

	unsigned int a;
	int iSize, it, nrec, iu, irec;
	iSize = startGroup[ics] - startGroup[ics - 1];
	fprintf(fp, "Order: %i\n", ics);
	fprintf(fp, "Minimal support: %i\n", iMinSupport);
	fprintf(fp, "Number of combinations with minimal support: %i\n", iSize);
	it = 0;
	for (int i = startGroup[ics - 1]; i < startGroup[ics]; i++){
		fprintf(fp, "Combination #%i (global #%i), event (new [", it, i);
		for (int j = 0; j < ics; j++){
			fprintf(fp, "%4i", vCombination[i*maxGroup + j]);
			if (j < ics - 1) fprintf(fp, ", ");
		}
		fprintf(fp, "], original [");
		for (int j = 0; j < ics; j++){
			fprintf(fp, "%4i", index1D[vCombination[i*maxGroup + j]]);
			if (j < ics - 1) fprintf(fp, ", ");
		}
		nrec = bitCount(vBits + i*nChunkRecords32);
		fprintf(fp, "], total records %i): ", nrec);
		iu = 0;
		for (int ii = 0; ii < nChunkRecords32; ii++){
			a = vBits[i*nChunkRecords32 + ii];
			for (int jj = 0; jj < n32; jj++){
				irec = ii*n32 + jj;
				if ((a & 1) == 1){
					fprintf(fp, "%i", irec);
					iu++;
					if (iu < nrec) fprintf(fp, ", ");
				}
				a = a >> 1;
			}
		}
		fprintf(fp, "\n");
		it++;
	}


	fclose(fp);
	return 0;
}

int dataSpam::printImage(){
	FILE *ft;
	char fName[500];
	sprintf(fName, "%s\\%s_pattern.txt", outputFolder, outPrefix);

	ft = fopen(fName, "w");
	
	iorder = (int *)malloc(sizeof(int)*curLen);

	for (long long i = 0; i < curLen*maxGroup; i++){
		if (vCombination[i] < 0) continue;
		vCombination[i] = index1D[vCombination[i]];
	}

	for (long long i = 0; i < curLen; i++){
		iorder[i] = i;
	}


	bool t;
	int n1, n2, ss, m;
	do{
		t = false;
		for (long long i = 0; i < curLen - 1; i++){
			for (long long j = 0; j < maxGroup; j++){
				n1 = vCombination[iorder[i] * maxGroup + j];
				n2 = vCombination[iorder[i + 1] * maxGroup + j];
				ss = n1 - n2;
				if (ss != 0)break;
			}
			if (ss < 0) continue;
			m = iorder[i];
			iorder[i] = iorder[i + 1];
			iorder[i + 1] = m;
			t = true;
		}
	} while (t);


	unsigned int *bloc, bit;
	int irec;
	for (long long i = 0; i < curLen; i++){
		if(vCombination[iorder[i] * maxGroup + 2] == -1)continue;
//		if (vCombination[iorder[i] * maxGroup + 1] == -1)continue;
		for (long long j = 0; j < maxGroup; j++){
			m = vCombination[iorder[i] * maxGroup + j];
			if (m < 0) {
				break;
			}else{
				fprintf(ft, "%i -1 ", m); 
			}
			
		}
		bloc = vBits + iorder[i] * nChunkRecords32;
		fprintf(ft, "#SUP: %i #SID:", bitCount(bloc));
		for (int k = 0; k < nChunkRecords32; k++){
			bit = bloc[k];
			for (int m = 0; m < n32; m++){
				irec = k*n32 + m;

				if (((bit >> m) & 1) != 1) continue;
				fprintf(ft, " %i", irec);
			}
		}
		fprintf(ft, "\n");
	}

	
	fclose(ft);
	return 0;
}

int dataSpam::bitCount(unsigned int *vec){
	int u;
	unsigned int a;
	u = 0;
	for (int ic = 0; ic < nChunkRecords32; ic++){
		a = vec[ic];
		for (int j = 0; j < n32; j++){
			u += (1 & a);
			a = a >> 1;
		}
	}
	return u;
}



int dataSpam::bitCountList(unsigned int *vec){
	unsigned int it, v, c;
	it = 0;
	for (int ic = 0; ic < nChunkRecords32; ic++){
		v = vec[ic];
	
		v = v - ((v >> 1) & 0x55555555);                    // reuse input as temporary
		v = (v & 0x33333333) + ((v >> 2) & 0x33333333);     // temp
		c = ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24; // count

		it += c;

	}
	return it;
}

unsigned int dataSpam::findBitList(unsigned int *vec, int tid){
	unsigned int it, n, u, a;
	unsigned int *iList_loc;

	iList_loc = pListIn[tid];
	it = 0;
	n = 0;
	for (int ic = 0; ic < nChunkRecords32; ic++){
		a = vec[ic];
		for (int j = 0; j < n32; j++){
			u = a & 1;
			iList_loc[n] = it * u;
			a = a >> 1;
			n += u;
			it++;
		}
	}
	return n;
}





int dataSpam::startProcessing(){

	sprintf(outputFile, "%s\\%s_info.txt", outputFolder, outPrefix);
	foi = fopen(outputFile, "w");
	if (foi == nullptr){
		printf("Error: cannot open file \"%s\".\n", outputFile);
		return -1;
	}

	support = 0.001*float(support1000);

	
	#pragma omp parallel
	{
		if (omp_get_thread_num() == 0){
			nthreads = __min(ncores, omp_get_num_threads());
		}
	}

	omp_set_dynamic(0);     
	omp_set_num_threads(nthreads);
	

	double gt;

	clock_t global_start, global_end, begin_time, end_time;
	if (readInput() != 0) return -1;
	global_start = clock();
	if (sortByEvent() != 0) return -2;
	if (count1D() != 0) return -4;
	if (remove1D() != 0) return -5;
	if (bit1D() != 0) return -6;

	fprintf(foi, "nThr\tSupp\tBeta\tisTR\tTimeRestriction\tisOrd\tDim\tnCand\tnFound\tRunTimeDimen\tFound/Candid\tTime/Candid\tTime/Found\n");

	for (int i = 2; i <= maxGroup; i++){
		if (findNewBitList(i) != 0) return -7;
	}
	
	if (isPrintPatterns == 1){
		if (printImage() != 0) return -11;
	}
	
	global_end = clock();
	gt = double(global_end - global_start) / CLOCKS_PER_SEC;
	printf("Size\t%i\n", curLen);
	printf("MemorySize\t%llu\n", curLen*nChunkRecords32);
	printf("TimeGlobal\t%f\n", gt);
		
	fprintf(foi, "Size\t%i\n", curLen);
	fprintf(foi, "MemorySize\t%llu\n", curLen*nChunkRecords32);
	fprintf(foi, "TimeGlobal\t%f\n", gt);
	fclose(foi); foi = nullptr;
	return 0;
}


int dataSpam::readSettings(){
	char setFile[500], buff[500];
	FILE *fs;

	sprintf(setFile, "fast_set_float2.txt");
	fs = fopen(setFile, "r");
	if (fs == nullptr){
		printf("Error: settings file \"%s\"\n", setFile);
		return -1;
	}

	fscanf(fs, "%*[^\n]\n"); fgets(buff, 500, fs); sscanf(buff, "%[^[\n#]]", inputFile);
	fscanf(fs, "%*[^\n]\n"); fgets(buff, 500, fs); sscanf(buff, "%[^[\n#]]", outputFolder);
	fscanf(fs, "%*[^\n]\n"); fgets(buff, 500, fs); sscanf(buff, "%[^[\n#]]", prefix);
//	fscanf(fs, "%*[^\n]\n%i\n", &ibeta);
	fscanf(fs, "%*[^\n]\n%f\n", &fbeta);
	fscanf(fs, "%*[^\n]\n%i\n", &ncores);
	fscanf(fs, "%*[^\n]\n%i\n", &support1000);
	fscanf(fs, "%*[^\n]\n%i\n", &maxGroup);
	fscanf(fs, "%*[^\n]\n%lld\n", &maxLen);
	fscanf(fs, "%*[^\n]\n%lld\n", &maxLen2);
	fscanf(fs, "%*[^\n]\n%i\n", &isPrintPatterns);
	fscanf(fs, "%*[^\n]\n%i\n", &isTimeRestriction);
	fscanf(fs, "%*[^\n]\n%f\n", &fTimeRestriction);
	fscanf(fs, "%*[^\n]\n%i\n", &isOrdered);

	fclose(fs); fs = nullptr;

	sprintf(outPrefix, "%s_beta_%f_nthreads_%i_sup_%i_iTR_%i_TR_%f_iO_%i", prefix, fbeta, ncores, support1000, isTimeRestriction, fTimeRestriction, isOrdered);
	
	return 0;
}

int main(){
	dataSpam *ds;

	ds = new dataSpam();

	if (ds->readSettings() == 0){
		ds->startProcessing();
	}

	delete ds; ds = nullptr;

		

	return 0;
}
