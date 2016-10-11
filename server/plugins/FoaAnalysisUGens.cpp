//
//  FoaAnalysisUGens.cpp
//  SuperCollider
//
//  Created by Michael McCrea on 8/3/15.
//
//

#include "SC_Plugin.h"

const double rHalfPi = 1. / M_PI_2;

static InterfaceTable *ft;


struct Aeda : public Unit { };

// for time domain onset detection/RMS
struct VRunningSum : public Unit {
    int nsamps, maxsamps, head, tail, resetcount;
    float msum,msum2;
    bool reset;
    float* msquares;
};

extern "C" {
    
    void load(InterfaceTable *inTable);
    
    void Aeda_Ctor(Aeda *unit);
    void Aeda_next(Aeda *unit, int inNumSamples);

    void VRunningSum_Ctor(VRunningSum *unit);
    void VRunningSum_Dtor(VRunningSum *unit);
    void VRunningSum_next_a(VRunningSum *unit, int inNumSamples);
}

void Aeda_Ctor( Aeda *unit) {

    SETCALC(Aeda_next);
    Aeda_next(unit,1);
}

void Aeda_next( Aeda *unit, int inNumSamples) {

    float *Win = IN(0);
    float *Xin = IN(1);
    float *Yin = IN(2);
    float *Zin = IN(3);
    
    float *out0 = OUT(0);   // azimuth
    float *out1 = OUT(1);   // elevation
    float *out2 = OUT(2);   // directivity
    float *out3 = OUT(3);   // amp

    float w, x, y, z, w_sqrd;
    float azim, elev, dir, amp;
    float p_sqrd, v_sqrd;
    
    
    /* aed-a calculation */
    
    for (int i=0; i<inNumSamples; ++i) {
        
        w = Win[i] * sqrt(2);       // scale W
        x = Xin[i];
        y = Yin[i];
        z = Zin[i];
        w_sqrd = pow(w,2);
        
        double b_sqrd[4] = { w_sqrd, pow(x,2), pow(y,2), pow(z,2) };
        
        p_sqrd = b_sqrd[0];                             // W component
        v_sqrd = b_sqrd[1] + b_sqrd[2] + b_sqrd[3];     // summed XYZ
        
        // directivity measure
        dir = M_PI_2 - (2. * atan2(sqrt(v_sqrd), sqrt(p_sqrd)));
//        dir = 1.0 - ( fabs(dir) / M_PI_2 );             // normalized dir
        dir = 1. - ( fabs(dir) * rHalfPi );
        
        //  pv_mean = b[0] * b;
        float pv_mean[4] = { w_sqrd, x * w, y * w, z * w };
        // atan2(y,x) - see pv_mean above
        azim = atan2(pv_mean[2], pv_mean[1]);
        // atan2(z,sqrt(x^2 + y^2))
        elev = atan2( pv_mean[3], sqrt( pow( pv_mean[1],2) + pow(pv_mean[2], 2 ) ) );
        // sqrt((b**2).sum / 2);
        amp = sqrt( ( b_sqrd[0] + b_sqrd[1] + b_sqrd[2] + b_sqrd[3] ) * 0.5 );
        
        out0[i] = azim;
        out1[i] = elev;
        out2[i] = dir;
        out3[i] = amp;
    }
    
}


// ------------------------------------------------------
// VRunningSum

void VRunningSum_Ctor( VRunningSum* unit )
{
    if ((int) ZIN0(2) == 0) {
        printf("VRunningSum Error: maxsamps initialized to 0.\n");
        SETCALC(*ClearUnitOutputs);
        unit->mDone = true;
        return;
    }
    
    SETCALC(VRunningSum_next_a);
    
    unit->maxsamps  = (int) ZIN0(2);
    unit->nsamps    = sc_max(1, sc_min( (int) ZIN0(1), unit->maxsamps )); // clamp 1>maxsamp
    unit->msum      = 0.0f;
    unit->msum2     = 0.0f;
    unit->resetcount= 0; //unit->msamp-1;
    unit->head      = 0; // first write position
    unit->tail      = unit->maxsamps - unit->nsamps;
    unit->reset     = false;
    
//    printf("initnsamp: %i, maxsamps: %i\n", unit->nsamps, unit->maxsamps);
    
    unit->msquares  = (float*)RTAlloc(unit->mWorld, unit->maxsamps * sizeof(float));
    //initialise to zeroes
    for(int i=0; i<unit->maxsamps; ++i)
        unit->msquares[i] = 0.f;
    
}

void VRunningSum_Dtor(VRunningSum *unit)
{
    RTFree(unit->mWorld, unit->msquares);
}


void VRunningSum_next_a( VRunningSum *unit, int inNumSamples )
{
    float *in   = ZIN(0);
    float *out  = ZOUT(0);
    
    int prevnsamps  = unit->nsamps;
    int maxsamps    = unit->maxsamps;
    int resetcount  = unit->resetcount;
    bool reset      = unit->reset;
    
    float * data= unit->msquares;
    float sum   = unit->msum;
    float sum2  = unit->msum2;          // modeled after RunningSum - thanks to Ross Bencina
    int head    = unit->head;           // current write index in the rolling buffer
    int tail    = unit->tail;           // current tail  index in the rolling buffer
    
    int startnsamps = unit->nsamps;   // keep track of previous window size
    
    int nsamps = (int) ZIN0(1);                     // number of samples to average
    nsamps = sc_max(1, sc_min(nsamps, maxsamps));   // clamp 1>maxsamp

    
//    int startnsamps = unit->nsamps;   // keep track of previous window size
    float dsamp_slope = CALCSLOPE( (float)nsamps, startnsamps );
//    float dsamp_slope = CALCSLOPE( (float)nsamps, prevnsamps );
    
//        printf("prevnsamps %i, startnsamps %i, nsamps %i, dsamp_slope: %f\n", prevnsamps, startnsamps, nsamps, dsamp_slope);
    
//    float preSampDiff = sum2;
    
    for (int i=0; i<inNumSamples; ++i) {
        
        // handle change in summing window size
        if (dsamp_slope != 0.f) {
            int nextnsamps = startnsamps + (int)(dsamp_slope * (i+1)); // CHECK +1, remove?
            int dsamp = nextnsamps - prevnsamps;    // window size delta
            
            if (dsamp != 0) {
                float dsampDiff = 0.;
                
                if (dsamp > 0) {
                    // window grows
                    for (int j=0; j<dsamp; ++j) {
                        
                        // grow the sum range by 1, wrap in range if needed
                        tail--;
                        if (tail < 0) tail += maxsamps; // tail is last index in buffer
                        
//                        // add value overtaken by growing window
//                        sum +=  data[tail];
//                        sum2 += data[tail];
                        
                        dsampDiff += data[tail];

//                        resetcount++;
                    }
                    // add accumulated value overtaken by growing window
                    // should be outside loop to avoid accumulating error
                    sum  += dsampDiff;
                    sum2 += dsampDiff;
                    
                } else {
                    // window shrinks: dsamp < 0
                    int test = abs(dsamp);
                    for (int j=0; j < test; ++j) {
//                        // subtract the sample values at the tail
//                        sum -=  data[tail];
//                        sum2 -= data[tail];
                        
                        dsampDiff += data[tail];
                        
                        // shrink the sum range by 1, wrap if needed
                        tail++;
                        if (tail == maxsamps) tail = 0;
                        
//                        resetcount--;
//                        resetcount++; // try always incrementing operation reset count
                    }
                    // remove sum of values accumulated by shrinking window
                    // should be outside loop to avoid accumulating error
                    sum  -= dsampDiff;
                    sum2 -= dsampDiff;
                }
                
                prevnsamps += dsamp;
            }
        }
        
        // remove the tail
        float rmv = data[tail];
        sum -=  rmv;

        // only remove last val from sum2 if sum2 wasn't just reset
        if  ( !reset ) {
            sum2 -= rmv;
            reset = false;
        }

        // add and store new inval
        float next= ZXP(in);
        data[head]= next;     // write the new sample in and add it
        sum  += next;
        sum2 += next;

        ZXP(out) = sum;
        
//        // DEBUG: the true sum of the window:
//        float trueSum = 0;
//        for (int j=0; j<nsamps; ++j) {
//            int dex = tail + 1 + j;
//            dex = dex % maxsamps;
//            trueSum += data[dex];
//        }
//        ZXP(out) = trueSum;
        
        
        // increment and wrap the head and tail
        head++;
        if (head == maxsamps) head = 0;
        tail++;
        if (tail == maxsamps) tail = 0;
        
        // swap the sums once window is full to avoid floating point error (tip from RunningSum)
        resetcount++;
        
        // orig
        if (resetcount == nsamps) {
            sum  = sum2;
            sum2 = 0.;
            resetcount = 0;
            reset = true;
        }
    }
    
    unit->nsamps = nsamps;
    unit->resetcount= resetcount;
    unit->head  = head;
    unit->msum  = sum;
    unit->msum2 = sum2;
    unit->tail  = tail;
    unit->reset = reset;
}


PluginLoad(FoaAnalysisUGens) {
    
    ft = inTable; // store pointer to InterfaceTable
    
    DefineSimpleUnit(Aeda);
    
    DefineDtorUnit(VRunningSum);
}

