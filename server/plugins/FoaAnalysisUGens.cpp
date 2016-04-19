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
    unit->head      = 0;
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
    
    int nsamps = (int) ZIN0(1);                     // number of samples to average
    nsamps = sc_max(1, sc_min(nsamps, maxsamps));   // clamp 1>maxsamp

    
//    printf("prevnsamps: %i, maxsamp: %i, resetcount: %i, samp: %i, head: %i, tail: %i\n", prevnsamps, maxsamps, resetcount, nsamps, head, tail);
    //printf("sum: f%, resetcount: %i, samp: %i, head: i%, tail: i%\n", sum, resetcount, nsamps, head, tail);
    
//    int dsamp   = nsamps - prevnsamps;  // detect change of samp average size
    
    int startnsamps = unit->nsamps;   // keep track of previous average span size
    float dsamp_slope = CALCSLOPE( (float)nsamps, startnsamps );
//    printf("nsamps %i, startnsamps %i, dsamp_slope: %f\n", nsamps, startnsamps, dsamp_slope);
    
    for (int i=0; i<inNumSamples; ++i) {
        
        // handle change in summing window size
        if (dsamp_slope != 0.f) {
            // add or remove tail values if the tail moved
            int nextnsamps = startnsamps + (int)(dsamp_slope * (i+1));
            int dsamp = nextnsamps - prevnsamps;    // delta of the sample average window size
            
            if (dsamp != 0) {
                
                if (dsamp > 0) {
                    for (int j=0; j<dsamp; ++j) {
                        // grow the sum range by 1, wrap in wrange if needed
                        tail--;
                        if (tail < 0) tail += maxsamps;
                        sum +=  data[tail];
                        sum2 += data[tail];
                        
                        resetcount++;
                    }
                } else { // dsamp < 0
                    int test = abs(dsamp);
                    for (int j=0; j < test; ++j) {
                        // subtract the sample values at the tail
                        sum -=  data[tail];
                        sum2 -= data[tail];
                        // shrink the sum range by 1, wrap if needed
                        tail++;
                        if (tail == maxsamps) tail = 0;
                        
                        resetcount--;
                    }
                }
                
                prevnsamps += dsamp;
            }
        }
//        printf("nsamps: %i, prevnsamps: %i\n", nsamps, prevnsamps);
        
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
        
//        printf("sum: %f, sum2: %f\n", sum, sum2);
        
        ZXP(out) = sum;
        
        // increment and wrap the head and tail
        head++;
        if (head == maxsamps) head = 0;
        tail++;
        if (tail == maxsamps) tail = 0;
        
        // swap the sums to avoid floating point error (tip from RunningSum)
        resetcount++;
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

//{
//    float *in = ZIN(0);
//    float *out = ZOUT(0);
//
//    int prevmsamp = unit->msamp;
//    int maxsamp = unit->maxsamp;
//    int resetcount = unit->resetcount;
//
//    int samp    = (int) ZIN0(1);                // number of samples to average
//    samp = sc_max(1, sc_min(samp, maxsamp));    // clamp 1>maxsamp
//    
//    int wrpntr   = unit->wrpntr;     // the current write position in the rolling buffer
//    
////        printf("prevmsamp: %i, maxsamp: %i, resetcount: %i, samp: %i, wrpntr: %i\n", prevmsamp, maxsamp, resetcount, samp, wrpntr);
//    
//    float * data= unit->msquares;
//    float sum   = unit->msum;
//    //avoids floating point error accumulation over time- thanks to Ross Bencina
//    float sum2  = unit->msum2;
//    
//    int dsamp   = samp - prevmsamp; // detect change of samp average size
//    int head    = unit->head;       // previous head index
//    
////    printf("prevmsamp: %i, maxsamp: %i, resetcount: %i, samp: %i, wrpntr: %i, head: %i\n", prevmsamp, maxsamp, resetcount, samp, wrpntr, head);
//    
//    // add or remove values if the head moved
//    
//    if (dsamp > 0) {
//        // adjust head position and wrap in maxsamp range
//        head -= dsamp;
//        if (head < 0)
//            head += maxsamp;
//
////         printf("growing size, head: %i\n", head);
//        
//        // average span grows, add old values back
//        for (int i=0; i<dsamp; ++i) {
//            int adddex = head + i;
//            if (adddex >maxsamp) adddex -= maxsamp; // wrap negative dex around maxsamp
//            float newoldval = data[adddex];
//            sum+= newoldval;
//            sum2+= newoldval;
//        }
//    }
//    if (dsamp < 0) {
//        // adjust head position and wrap in maxsamp range
//        head += dsamp;
//        if (head >= maxsamp)
//            head -= maxsamp;
//
////        printf("shrinking size, head: %i\n", head);
//        
//        // average span shrinks
//        int iter = abs(dsamp);
//        for (int i=0; i<iter; ++i) {
//            int rmvdex = head - i;
//            if (rmvdex <0) rmvdex += maxsamp;
//            float remoldval = data[rmvdex];
//            sum -= remoldval;
//            sum2-= remoldval;
//        }
//    }
//    
//    
//    int todo    = 0;
//    int done    = 0;
//    
//    // iterate over inNumSamples
//    while(done<inNumSamples) {
//        
//        // todo is the size of a block of samples remaining to add to the sum
//        // before swapping sum with sum2, or inNumSamples if it's less than remaining buffer samples until reset
//        todo= sc_min( inNumSamples-done, samp-resetcount );
////        printf("todo: %i\n", todo);
//        
//        for(int j=0;j<todo;++j) {
//            sum -=data[head];       // subtract the sample from the head
//            float next= ZXP(in);
//            data[wrpntr]= next;     // write the new sample in and add it
//            sum += next;
//            sum2 +=next;
//            
//            ZXP(out) = sum;
//            
//            ++resetcount;
//            
//            ++wrpntr;
//            // wrap count within maxbufsize range
//            if (wrpntr == maxsamp)
//                wrpntr = 0;
//
//            ++head;
//            // wrap count within maxbufsize range
//            if (head == maxsamp)
//                head = 0;
////            printf("resetcount: %i, wrpntr: %i, head: %i, sum: %f, sum2: %f, next: %f\n",
////                   resetcount, wrpntr, head, sum, sum2, next);
//        }
////        printf("tick\n");
//       
//        done += todo;
////        printf("done: %i\n", done);
//        
//        // if count reaches the end of specified sample cycle, swap sum for sum2
//        if( resetcount >= samp ) {
//            // count = 0;
//            sum = sum2;
//            sum2 = 0;
//            resetcount = 0;
//        }
////        printf("resetcount: %i, wrpntr: %i, head: %i, sum: %f, sum2: %f\n",
////               resetcount, wrpntr, head, sum, sum2);
//
//    }
//
////    printf("out\n");
//    
//    unit->msamp = samp;
//    unit->resetcount= resetcount;
//    unit->wrpntr  = wrpntr;
//    unit->msum  = sum;
//    unit->msum2 = sum2;
//    unit->head  = head;
//}


PluginLoad(FoaAnalysisUGens) {
    
    ft = inTable; // store pointer to InterfaceTable
    
    DefineSimpleUnit(Aeda);
    
    DefineDtorUnit(VRunningSum);
}

