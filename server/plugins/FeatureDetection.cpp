/*
	SuperCollider real time audio synthesis system
    Copyright (c) 2002 James McCartney. All rights reserved.
	http://www.audiosynth.com

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
*/

//Feature (Onset) Detection implemented by sick lincoln for sc3
//Jensen,K. & Andersen, T. H. (2003). Real-time Beat Estimation Using Feature Extraction.
//In Proceedings of the Computer Music Modeling and RetrievalSymposium, Lecture Notes in Computer Science. Springer Verlag.
//Hainsworth, S. (2003) Techniques for the Automated Analysis of Musical Audio. PhD, university of cambridge engineering dept.

//possible to make a Goto style Detector for a given band and with history of two samples-
//should do separately as PV_GotoBandTrack

//next perhaps Duxbury et al/ Mauri et al different conception of high frequency content with ratio of changes

#include "FFT_UGens.h"

struct PV_OnsetDetectionBase : public Unit
{
	float * m_prevframe;
	int m_numbins;
	int m_waiting, m_waitSamp, m_waitLen;
};

//FFT onset detector combining 4 advised features from Jensen/Andersen
struct PV_JensenAndersen : public PV_OnsetDetectionBase
{
	float m_hfc,m_hfe,m_sc,m_sf;
	int m_fourkindex;
};


//FFT onset detector combining 2 advised features from Hainsworth PhD
struct PV_HainsworthFoote : public PV_OnsetDetectionBase
{
	float m_prevNorm;
	int m_5kindex,m_30Hzindex;

};

//for time domain onset detection/RMS
struct RunningSum : public Unit {
	int msamp, mcount;
	float msum,msum2;
	//float mmeanmult;
	float* msquares;
};

// like RunningSum, but with variable size summing window. - mtmccrea
struct RunningSum2 : public Unit {
    int nsamps, maxSamps, head, tail, resetCounter, resetSamps;
    float msum, msum2;
    bool reset;
    float* msquares;
};

extern "C"
{
	void PV_OnsetDetectionBase_Ctor(PV_OnsetDetectionBase *unit);
	void PV_OnsetDetectionBase_Dtor(PV_OnsetDetectionBase *unit);

	void PV_JensenAndersen_Ctor(PV_JensenAndersen *unit);
	void PV_JensenAndersen_Dtor(PV_JensenAndersen *unit);
	void PV_JensenAndersen_next(PV_JensenAndersen *unit, int inNumSamples);

	void PV_HainsworthFoote_Ctor(PV_HainsworthFoote *unit);
	void PV_HainsworthFoote_Dtor(PV_HainsworthFoote *unit);
	void PV_HainsworthFoote_next(PV_HainsworthFoote *unit, int inNumSamples);

	void RunningSum_next_k(RunningSum *unit, int inNumSamples);
	void RunningSum_Ctor(RunningSum* unit);
	void RunningSum_Dtor(RunningSum* unit);

    void RunningSum2_Ctor(RunningSum2 *unit);
    void RunningSum2_Dtor(RunningSum2 *unit);
    void RunningSum2_next(RunningSum2 *unit, int inNumSamples);
}

#define PV_FEAT_GET_BUF_UNLOCKED \
	uint32 ibufnum = (uint32)fbufnum; \
	int bufOK = 1; \
	World *world = unit->mWorld; \
	SndBuf *buf; \
	if (ibufnum >= world->mNumSndBufs) { \
		int localBufNum = ibufnum - world->mNumSndBufs; \
		Graph *parent = unit->mParent; \
		if(localBufNum <= parent->localBufNum) { \
			buf = parent->mLocalSndBufs + localBufNum; \
		} else { \
			bufOK = 0; \
			buf = world->mSndBufs; \
			if(unit->mWorld->mVerbosity > -1){ Print("FFT Ctor error: Buffer number overrun: %i\n", ibufnum); } \
		} \
	} else { \
		buf = world->mSndBufs + ibufnum; \
	} \
	int numbins = (buf->samples - 2) >> 1; \
	if (!buf->data) { \
		if(unit->mWorld->mVerbosity > -1){ Print("FFT Ctor error: Buffer %i not initialised.\n", ibufnum); } \
		bufOK = 0; \
	}

#define PV_FEAT_GET_BUF 		\
	PV_FEAT_GET_BUF_UNLOCKED	\
	LOCK_SNDBUF(buf);


void PV_OnsetDetectionBase_Ctor(PV_OnsetDetectionBase *unit)
{
	float fbufnum = ZIN0(0);

	PV_FEAT_GET_BUF_UNLOCKED

	unit->m_numbins = numbins;
	int insize = unit->m_numbins * sizeof(float);

	if(bufOK) {
		unit->m_prevframe = (float*)RTAlloc(unit->mWorld, insize);
		memset(unit->m_prevframe, 0, insize);
	}

	unit->m_waiting=0;
	unit->m_waitSamp=0;
	unit->m_waitLen=0;
	ClearUnitOutputs(unit, 1);
}

void PV_OnsetDetectionBase_Dtor(PV_OnsetDetectionBase *unit)
{
	if(unit->m_prevframe)
		RTFree(unit->mWorld, unit->m_prevframe);
}



void PV_JensenAndersen_Ctor(PV_JensenAndersen *unit)
{
	PV_OnsetDetectionBase_Ctor(unit);

	unit->m_hfc= 0.0;
	unit->m_hfe= 0.0;
	unit->m_sf= 0.0;
	unit->m_sc= 0.0;

	unit->m_fourkindex= (int)(4000.0/(unit->mWorld->mSampleRate))*(unit->m_numbins);

	SETCALC(PV_JensenAndersen_next);
}

void PV_JensenAndersen_Dtor(PV_JensenAndersen *unit)
{
	PV_OnsetDetectionBase_Dtor(unit);
}


void PV_JensenAndersen_next(PV_JensenAndersen *unit, int inNumSamples)
{
	float outval=0.0;
	float fbufnum = ZIN0(0);

	if(unit->m_waiting==1) {
		unit->m_waitSamp+=inNumSamples;
		if(unit->m_waitSamp>=unit->m_waitLen)
			unit->m_waiting=0;
	}

	if (!(fbufnum < 0.f))
	//if buffer ready to process
	{
		PV_FEAT_GET_BUF

		SCPolarBuf *p = ToPolarApx(buf);

		//four spectral features useful for onset detection according to Jensen/Andersen

		float magsum=0.0, magsumk=0.0, magsumkk=0.0, sfsum=0.0, hfesum=0.0;

		float * q= unit->m_prevframe;

		int k4= unit->m_fourkindex;

		//ignores dc, nyquist
		for (int i=0; i<numbins; ++i) {
			float mag= ((p->bin[i]).mag);
			int k= i+1;
			float qmag= q[i];
			magsum+=mag;
			magsumk+=k*mag;
			magsumkk+=k*k*mag;
			sfsum+= fabs(mag- (qmag));
			if(i>k4) hfesum+=mag;
		}

		float binmult= 1.f/numbins;
		//normalise
		float sc= (magsumk/magsum)*binmult;
		float hfe= hfesum*binmult;
		float hfc= magsumkk*binmult*binmult*binmult;
		float sf= sfsum*binmult;

		//printf("sc %f hfe %f hfc %f sf %f \n",sc, hfe, hfc, sf);

		//if(sc<0.0) sc=0.0;
		//if(hfe<0.0) hfe=0.0;
		//if(hfc<0.0) hfc=0.0;
		//if(sf<0.0) sf=0.0;

		//ratio of current to previous frame perhaps better indicator than first derivative difference
		float scdiff= sc-(unit->m_sc);
		float hfediff= hfe-(unit->m_hfe);
		float hfcdiff= hfc-(unit->m_hfc);
		float sfdiff= sf-(unit->m_sf);

		//store as old frame values for taking difference
		unit->m_sc=sc;
		unit->m_hfe=hfe;
		unit->m_hfc=hfc;
		unit->m_sf=sf;

		//printf("sc %f hfe %f hfc %f sf %f \n",sc, hfe, hfc, sf);
		//printf("sc %f hfe %f hfc %f sf %f \n",scdiff, hfediff, hfcdiff, sfdiff);

		//does this trigger?
		//may need to take derivatives across previous frames by storing old values

		float sum = (ZIN0(1)*scdiff)+(ZIN0(2)*hfediff)+(ZIN0(3)*hfcdiff)+(ZIN0(4)*sfdiff);

		//printf("sum %f thresh %f \n",sum, ZIN0(7));

		//if over threshold, may also impose a wait here
		if(sum>ZIN0(5) && (unit->m_waiting==0)) {//printf("bang! \n");
			outval=1.0;
			unit->m_waiting=1;
			unit->m_waitSamp=inNumSamples;
			unit->m_waitLen=(int)(ZIN0(6)*(world->mSampleRate));
		}

		//take copy of this frame's magnitudes as prevframe

		for (int i=0; i<numbins; ++i)
			q[i]= p->bin[i].mag;
	}

	Fill(inNumSamples, &ZOUT0(0), outval);
}


void PV_HainsworthFoote_Ctor(PV_HainsworthFoote *unit)
{
	PV_OnsetDetectionBase_Ctor(unit);

	World *world = unit->mWorld;

	unit->m_5kindex= (int)((5000.0/(world->mSampleRate))*(unit->m_numbins));
	unit->m_30Hzindex= (int)((30.0/(world->mSampleRate))*(unit->m_numbins));

	unit->m_prevNorm= 1.0;

	//unit->m_5kindex,  unit->m_30Hzindex,
	//printf("numbins %d  sr %d \n",  unit->m_numbins, world->mSampleRate);
	//printf("test %d sr %f 5k %d 30Hz %d\n", unit->m_numbins, world->mSampleRate, unit->m_5kindex,  unit->m_30Hzindex);

	SETCALC(PV_HainsworthFoote_next);
}

void PV_HainsworthFoote_Dtor(PV_HainsworthFoote *unit)
{
	PV_OnsetDetectionBase_Dtor(unit);
}

static const float lmult= 1.442695040889; //loge(2) reciprocal

void PV_HainsworthFoote_next(PV_HainsworthFoote *unit, int inNumSamples)
{
	float outval=0.0;
	float fbufnum = ZIN0(0);

	if(unit->m_waiting==1)
	{
		unit->m_waitSamp+=inNumSamples;
		if(unit->m_waitSamp>=unit->m_waitLen) {unit->m_waiting=0;}
	}

	if (!(fbufnum < 0.f))
	//if buffer ready to process
	{
		PV_FEAT_GET_BUF

		SCPolarBuf *p = ToPolarApx(buf);

		float dnk, prevmag, mkl=0.0, footesum=0.0, norm=0.0;

		float * q= unit->m_prevframe;

		int k5= unit->m_5kindex;
		int h30= unit->m_30Hzindex;

		for (int i=0; i<numbins; ++i) {
			float mag= ((p->bin[i]).mag);
			float qmag= q[i];

			if(i>=h30 && i<k5) {
				prevmag= qmag;
				//avoid divide by zero
				if(prevmag<0.0001) prevmag=0.0001;

				//no log2 in maths library, so use log2(x)= log(x)/log(2) where log is to base e
				//could just use log and ignore scale factor but hey let's stay accurate to the source for now
				dnk= log(mag/prevmag)*lmult;

				if(dnk>0.0) mkl+=dnk;
			}

			norm+=mag*mag;
			footesum+=mag*qmag;
		}


		mkl= mkl/(k5-h30);
		//Foote measure- footediv will be zero initially
		float footediv= ((sqrt(norm))*(sqrt(unit->m_prevNorm)));
		if(footediv<0.0001f)
			footediv=0.0001f;
		float foote= 1.0- (footesum/footediv); //1.0 - similarity
		//printf("mkl %f foote %f \n",mkl, foote);

		unit->m_prevNorm= norm;
		float sum = (ZIN0(1)*mkl)+(ZIN0(2)*foote);

		//printf("sum %f thresh %f \n",sum, ZIN0(7));

		//if over threshold, may also impose a 50mS wait here
		if(sum>ZIN0(3) && (unit->m_waiting==0)) {
			outval=1.0;
			unit->m_waiting=1;
			unit->m_waitSamp=inNumSamples;
			unit->m_waitLen=(int)(ZIN0(4)*(unit->mWorld->mSampleRate));
		}

		//take copy of this frame's magnitudes as prevframe

		for (int i=0; i<numbins; ++i)
			q[i]= p->bin[i].mag;
	}

	Fill(inNumSamples, &ZOUT0(0), outval);
}


void RunningSum_Ctor( RunningSum* unit )
{
	SETCALC(RunningSum_next_k);

	unit->msamp= (int) ZIN0(1);

	//unit->mmeanmult= 1.0f/(unit->msamp);
	unit->msum=0.0f;
	unit->msum2=0.0f;
	unit->mcount=0; //unit->msamp-1;

	unit->msquares= (float*)RTAlloc(unit->mWorld, unit->msamp * sizeof(float));
	//initialise to zeroes
	for(int i=0; i<unit->msamp; ++i)
	unit->msquares[i]=0.f;
	OUT0(0) = 0.f;

}

void RunningSum_Dtor(RunningSum *unit)
{
	RTFree(unit->mWorld, unit->msquares);
}

//RMS is easy because convolution kernel can be updated just by deleting oldest sample and adding newest-
//half hanning window convolution etc requires updating values for all samples in memory on each iteration
void RunningSum_next_k( RunningSum *unit, int inNumSamples )
{
	float *in = ZIN(0);
	float *out = ZOUT(0);

	int count= unit->mcount;
	int samp= unit->msamp;

	float * data= unit->msquares;
	float sum= unit->msum;
	//avoids floating point error accumulation over time- thanks to Ross Bencina
	float sum2= unit->msum2;

	int todo=0;
	int done=0;
	while(done<inNumSamples) {
		todo= sc_min(inNumSamples-done,samp-count);

		for(int j=0;j<todo;++j) {
			sum -=data[count];
			float next= ZXP(in);
			data[count]= next;
			sum += next;
			sum2 +=next;
			ZXP(out) = sum;
			++count;
		}

		done+=todo;

		if( count == samp ) {
			count = 0;
			sum = sum2;
			sum2 = 0;
		}

	}

	unit->mcount =count;
	unit->msum =  sum;
	unit->msum2 =  sum2;
}


// RunningSum2 adds a variable window size to RunningSum - mtmccrea
void RunningSum2_Ctor( RunningSum2* unit )
{
    if ((int) ZIN0(2) == 0) {
        printf("RunningSum2 Error: maxSamps initialized to 0.\n");
        SETCALC(*ClearUnitOutputs);
        unit->mDone = true;
        return;
    }

    SETCALC(RunningSum2_next);

    unit->maxSamps  = (int) ZIN0(2);
    unit->nsamps    = sc_max(1, sc_min((int) ZIN0(1), unit->maxSamps)); // clip(1, maxSamps)
    unit->msum      = 0.0f;
    unit->msum2     = 0.0f;
    unit->resetCounter = 0;
    unit->head      = 0; // first write position
    unit->tail      = unit->maxSamps - unit->nsamps;
    unit->reset     = false;
    unit->msquares  = (float*)RTAlloc(unit->mWorld, unit->maxSamps * sizeof(float));
	unit->resetSamps = unit->nsamps;

    if (unit->msquares == nullptr) {
        SETCALC(*ClearUnitOutputs);
        ClearUnitOutputs(unit, 1);
        if (unit->mWorld->mVerbosity > -2) {
            printf("Failed to allocate memory for RunningSum2\n");
        }
        return;
    }

    //initialise to zeroes
    for (int i=0; i < unit->maxSamps; ++i)
        unit->msquares[i] = 0.f;

    ZOUT0(0) = 0.f;
}

void RunningSum2_Dtor(RunningSum2 *unit)
{
    RTFree(unit->mWorld, unit->msquares);
}

void RunningSum2_next( RunningSum2 *unit, int inNumSamples )
{
    float *in   = ZIN(0);
    float *out  = ZOUT(0);
    float *data = unit->msquares;
	
	int maxWinSize  = unit->maxSamps;
	int newWinSize  = (int) ZIN0(1);  // number of samples to average
    int prevWinSize = unit->nsamps;   // keep track of previous block's window size
	int curWinSize  = unit->nsamps;   // keep track of window size as it ramps to new size
	
	int head = unit->head; // current write index in the rolling buffer
    int tail = unit->tail; // current tail  index in the rolling buffer
	
    float sum = unit->msum;
    float sum2 = unit->msum2; // modeled after RunningSum - thanks to Ross Bencina
	
	int resetCounter = unit->resetCounter; // trigger sum<>sum2 swap
	bool reset = unit->reset;

	bool error = false; // TEMP
	
    newWinSize = sc_max(1, sc_min(newWinSize, maxWinSize));   // clamp [1, maxWinSize]
	
	if (newWinSize > prevWinSize)
		// window grows, grow reset counter to match
		unit->resetSamps = newWinSize;
	if (resetCounter < newWinSize)
		// reset count is less than newWinSize can safely set reset counter to newWinSize
		unit->resetSamps = newWinSize;
	
    float sampSlope = CALCSLOPE( (float)newWinSize, prevWinSize );
	
    for (int i = 0; i < inNumSamples; ++i) {
        if (sampSlope != 0.0f) {					// handle change in summing window size
			// step this many samples for this output samples
            int steps = prevWinSize + (int)(sampSlope * (i+1)) - curWinSize;
			
            if (steps != 0) {
                float sumChange = 0.0;

                if (steps > 0) {                    // window grows
					// printf("growing window by %d\n", step);
                    for (int j = 0; j < steps; ++j) {
						// grow window
                        tail--;
						
                        if (tail < 0)
							// wrap
							tail += maxWinSize;
						
						// add tail values back into sum
                        sumChange += data[tail];
                    }
					
                    // add all accumulated values from growing tail
                    // should be outside loop to avoid accumulating error
                    sum += sumChange;
					// sum2 += sumChange; // these will be added in the course of accumulating sum2
					curWinSize += steps;
                } else {                            // window shrinks: step < 0
                    steps = abs(steps);
					float sumPreCopy = 0.0;
					
//					printf("shrinking window by %d\n", steps);
                    for (int j=0; j < steps; ++j) {
					
						// accumulate the amount to remove from the tail
						sumChange += data[tail];
						
						// shrink the size of the sum window
						curWinSize--;
						
						if (curWinSize == resetCounter) { // window shrinks to size equal to resetCounter, trigger sum2 copy
//							/* ERROR CHECKING */
//							/* END ERROR CHECKING */
							
							sum = sum2;
							sum2 = 0.0f;
							sumChange = 0;
							resetCounter = 0;
							reset = true;
							// only update resetSamps if window size crossed resetCounter
							unit->resetSamps = newWinSize;
						}
						
						tail++;                     // move tail
						if (tail == maxWinSize)
                            tail = 0;               // wrap
                    }

					// remove sum of values accumulated by shrinking window
                    // should be outside loop to avoid accumulating error
					sum -= sumChange;
                    // sum2 -= sumChange;
					
					if (sumPreCopy != 0.0f) {
						sumPreCopy -= sumChange;
						// printf("resolves to [sum from sum2],[manual sum] %f  %f \n", sum, sumPreCopy);
					}

				}
//              curWinSize += steps;
			}
        }

        // remove the tail value from the sum
        float rmv = data[tail];
        sum -= rmv;
        // don't remove last val from sum2 if sum2 was just reset to 0
        if (!reset) {
            sum2 -= rmv;
            reset = false;
        }

        // write the new sample to the data buffer and add it to the sums
        float next= ZXP(in);
        data[head]= next;
        sum += next;
        sum2 += next;

        ZXP(out) = sum;

        // increment and wrap the head and tail indices
        head++;
        if (head == maxWinSize)
            head = 0;
        tail++;
        if (tail == maxWinSize)
            tail = 0;
		
        resetCounter++;
		
		// swap the sums once window is full to avoid
        // floating point error (tip from RunningSum)
		if (resetCounter == unit->resetSamps) {
			//		if (resetCounter == newWinSize) {
			float maxsum = sc_max(sum, sum2);
			float minsum = sc_min(sum, sum2);
			if ((minsum/maxsum) < 0.99) { // within 1% ?
				//			if (sum != sum2) {
				printf("~~~~ MISMATCH IN SWAPPING SUMS : ");
				printf("swapping sum %f and sum2 %f resetCounter %d\n", sum, sum2, resetCounter);
				error = true;
			} else {
				error = false;
			}
			
            sum = sum2;
            sum2 = 0.0f;
            resetCounter = 0;
            reset = true;
			unit->resetSamps = newWinSize; // only update resetSamps if resetCounter has been reached
        }
    }

	if (error) {
		if (fabs(sum - newWinSize) > 1) {
			printf("	MISMATCH IN SWAPPING SUMS ~~~~~~~~~~~~~~~~~~~~~~~~~");
			printf("out %f should be %d \n", sum, newWinSize);
		}
	}
	
    unit->nsamps = newWinSize;
    unit->resetCounter = resetCounter;
    unit->head  = head;
    unit->msum  = sum;
    unit->msum2 = sum2;
    unit->tail  = tail;
    unit->reset = reset;
}


//
//							/* ERROR CHECKING */
//							float maxsum = sc_max(sum, sum2);
//							float minsum = sc_min(sum, sum2);
//
//							if ((minsum/maxsum) < 0.99) { // within 1% ?
//							// if (sum != sum2) {
//
//								// printf("~~~~ MISMATCH IN SWAPPING SUMS : ");
//								// printf("swapping while shrinking sum %f and sum2 %f resetCounter %d\n", sum, sum2, resetCounter);
//								error = true;
//								sumPreCopy = sum;
//							} else {
//								error = false;
//							}
//							/* END ERROR CHECKING */

void initFeatureDetectors(InterfaceTable *it)
{
    DefineDtorUnit(PV_JensenAndersen);
    DefineDtorUnit(PV_HainsworthFoote);
    DefineDtorUnit(RunningSum);
    DefineDtorUnit(RunningSum2);
}
