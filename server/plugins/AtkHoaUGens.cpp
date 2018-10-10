//
//  AtkHoaUGens.cpp
//  SuperCollider
//
//  Created by Michael McCrea on 10/3/18.
//

#include "AtkHoaUGens.hpp"
#include "SC_PlugIn.hpp"

#include "boost/math/special_functions.hpp"
#include "boost/math/special_functions/spherical_harmonic.hpp"
#include "boost/math/special_functions/legendre.hpp"


// InterfaceTable contains pointers to functions in the host (server).
static InterfaceTable *ft;

// declare struct to hold unit generator state
struct HoaSphCoeff : public SCUnit{
    
    // Constructor usually does 3 things.
    // 1. set the calculation function.
    // 2. initialize the unit generator state variables.
    // 3. calculate one sample of output.
public:
    HoaSphCoeff() {
       
        // 2. initialize the unit generator state variables.
        mDeg = in0(0); // scalar, non-modulatable, could be ints
        mIdx = in0(1);
        
        mPrevTheta = in0(2);
        mPrevPhi = in0(3);
		
		// calc mPrevOut now in case calc function is kk or ss (could move to those cases)
        mPrevOut = boost::math::spherical_harmonic_r<float>(mDeg, mIdx, mPrevTheta, mPrevPhi);

        // 1. set the calculation function.
        switch (inRate(2)) {
            case calc_ScalarRate:
                switch (inRate(3)) {
                    case calc_BufRate:
                        Print("sk calc\n");
                        set_calc_function<HoaSphCoeff,&HoaSphCoeff::next_sk>();
                        next_sk(1);
						break;
                    case calc_FullRate:
                        Print("sa calc\n");
                        set_calc_function<HoaSphCoeff,&HoaSphCoeff::next_sa>();
                        next_sa(1);
						break;
					case calc_ScalarRate:
						Print("ss calc\n");
						set_calc_function<HoaSphCoeff,&HoaSphCoeff::next_ss>();
						next_ss(1);
						break;
                }
				break;
            case calc_BufRate:
                switch (inRate(3)) {
                    case calc_BufRate:
                        Print("kk calc\n");
                        set_calc_function<HoaSphCoeff,&HoaSphCoeff::next_kk>();
                        next_kk(1);
						break;
                    case calc_FullRate:
                        Print("ka calc\n");
                        set_calc_function<HoaSphCoeff,&HoaSphCoeff::next_ka>();
                        next_ka(1);
						break;
					case calc_ScalarRate:
						Print("ks calc\n");
						set_calc_function<HoaSphCoeff,&HoaSphCoeff::next_ks>();
						next_ks(1);
						break;
                }
				break;
            case calc_FullRate:
                switch (inRate(3)) {
                    case calc_BufRate:
                        Print("ak calc\n");
                        set_calc_function<HoaSphCoeff,&HoaSphCoeff::next_ak>();
                        next_ak(1);
						break;
                    case calc_FullRate:
                        Print("aa calc\n");
                        set_calc_function<HoaSphCoeff,&HoaSphCoeff::next_aa>();
                        next_aa(1);
						break;
					case calc_ScalarRate:
						Print("as calc\n");
						set_calc_function<HoaSphCoeff,&HoaSphCoeff::next_as>();
						next_as(1);
						break;
                }
				break;
        }
		Print("\n");
    }
    
private:
    int mDeg;
    int mIdx;
    float mPrevTheta;
    float mPrevPhi;
    float mPrevOut;
    
    //////////////////////////////////////////////////////////////////
    
    bool inputsChanged(float theta, float phi)
    {
        bool changed = (theta != mPrevTheta || phi != mPrevPhi);
//        if (changed) {
            mPrevTheta = theta;
            mPrevPhi = phi;
//        }
        return changed;
    }
    
    template <typename SignalT, typename SignalP>
    void generateOutput (int inNumSamples, SignalT theta, SignalP phi) {
        float newTheta;
        float newPhi;
        
        float *outBuf = out(0);     // pointer to the output buffer
        
        // inNumSamples will be blocksize (e.g. 64) for ar, 1 for kr.
        for (int i=0; i < inNumSamples; ++i)
        {
            newTheta = theta.consume();
            newPhi = phi.consume();
            
            // check if inputs have changed, if not, output previous calc
            if (!inputsChanged(newTheta, newPhi)) {
				printf("no input change %d\n", i);
                outBuf[i] = mPrevOut;
            } else {
                // unsigned n, int m, T theta, T phi
                float sphHarmR = boost::math::spherical_harmonic_r<float>(mDeg, mIdx, newTheta, newPhi);
                
                outBuf[i] = sphHarmR;
                mPrevOut = sphHarmR;
            }
        }
        
//        mPrevTheta = newTheta;
//        mPrevPhi = newPhi;
    }
    
//    auto getSignal(int index, float prev) {
//        if (isScalarRateIn(index)) {
//            return makeScalar(in0(index))
//        }
//        if (isControlRateIn(index)) {
//            return makeSlope<float>(in0(index), prev)
//        }
//        if (isAudioRateIn(index)) {
//            return makeSignal(index)
//        }
//    }
    
//    void next_test(int inNumSamples)
//    {
//        auto theta {makeScalar(in0(2))};
//        auto phi = {};
//
//        // determine type arg 1 signal
//        if (isScalarRateIn(2)) {
//            theta = makeScalar(in0(2));
//        }
//        if (isControlRateIn(2)) {
//            theta = makeSlope<float>(in0(2), mPrevTheta)
//        }
//        if (isAudioRateIn(2)) {
//            theta = makeSignal(2);
//        }
//
//        generateOutput<decltype(theta), decltype(phi)>(inNumSamples, theta, phi);
//    }
    
//    // calculation function: ar theta, ar phi
//    void build_calc_func(int inNumSamples)
//    {
//        AudioSignal<float> theta = makeSignal(2);
//        AudioSignal<float> phi = makeSignal(3);
//
//#define DefineOutputGenFunc(type1, type2, inNumSamples, theta, phi) \
//    HoaSphCoeff::generateOutput<type1, type2>(inNumSamples, theta, phi);
//
//
//        DefineOutputGenFunc(AudioSignal<float>, AudioSignal<float>, inNumSamples, theta, phi);
//    }
    
//    // slimmer version?
//    void next_aa(int inNumSamples, thetaSig, phiSig)
//    {
//        auto theta = makeSignal(2);
//        auto phi = makeSignal(3);
//
//        generateOutput<decltype(theta), decltype(theta)>(inNumSamples, theta, phi);
//    }
    
    template <typename FloatType>
    void updateAudioSignal(AudioSignal<FloatType> audioSignal, FloatType * pointer)
    {
        audioSignal.pointer = pointer;
    }
    
    template <typename FloatType>
    void updateSlopeSignal(SlopeSignal<FloatType> slopeSignal, FloatType next, FloatType prev)
    {
        slopeSignal.slope = calcSlope<FloatType>(next, prev);
        slopeSignal.value = prev;
    }
    
    
    
    /* calculation functions: next_thetaRatephiRate */
    
    void next_aa(int inNumSamples)
    {
        AudioSignal<float> theta = makeSignal(2);
        AudioSignal<float> phi = makeSignal(3);
        
        generateOutput<AudioSignal<float>, AudioSignal<float>>(inNumSamples, theta, phi);
    }
    
    void next_kk(int inNumSamples)
    {
        SlopeSignal<float> theta = makeSlope<float>(in0(2), mPrevTheta);
        SlopeSignal<float> phi = makeSlope<float>(in0(3), mPrevPhi);
        
        generateOutput<SlopeSignal<float>, SlopeSignal<float>>(inNumSamples, theta, phi);
    }

    void next_ka(int inNumSamples)
    {
        SlopeSignal<float> theta = makeSlope<float>(in0(2), mPrevTheta);
        AudioSignal<float> phi = makeSignal(3);
        
        generateOutput<SlopeSignal<float>, AudioSignal<float>>(inNumSamples, theta, phi);
    }

    void next_ak(int inNumSamples)
    {
        AudioSignal<float> theta = makeSignal(2);
        SlopeSignal<float> phi = makeSlope<float>(in0(3), mPrevPhi);

        generateOutput<AudioSignal<float>, SlopeSignal<float>>(inNumSamples, theta, phi);
    }
    
    void next_as(int inNumSamples)
    {
        AudioSignal<float> theta = makeSignal(2);
        ScalarSignal<float> phi = makeScalar<float>(in0(3));
        
        generateOutput<AudioSignal<float>, ScalarSignal<float>>(inNumSamples, theta, phi);
    }
    
    void next_sa(int inNumSamples)
    {
        ScalarSignal<float> theta = makeScalar<float>(in0(2));
        AudioSignal<float> phi = makeSignal(3);
        
        generateOutput<ScalarSignal<float>, AudioSignal<float>>(inNumSamples, theta, phi);
    }
    
    void next_ks(int inNumSamples)
    {
        SlopeSignal<float> theta = makeSlope<float>(in0(2), mPrevTheta);
        ScalarSignal<float> phi = makeScalar<float>(in0(3));
        
        generateOutput<SlopeSignal<float>, ScalarSignal<float>>(inNumSamples, theta, phi);
    }
    
    void next_sk(int inNumSamples)
    {
        ScalarSignal<float> theta = makeScalar<float>(in0(2));
        SlopeSignal<float> phi = makeSlope<float>(in0(3), mPrevPhi);
        
        generateOutput<ScalarSignal<float>, SlopeSignal<float>>(inNumSamples, theta, phi);
    }
    
    // special case: output value is constant
    void next_ss (int inNumSamples) {
        float *outBuf = out(0);     // pointer to the output buffer
        
        // inNumSamples will be blocksize (e.g. 64) for ar, 1 for kr.
        for (int i=0; i < inNumSamples; ++i)
        {
            outBuf[i] = mPrevOut;
        }
    }

    
//    // calculation function for an audio rate frequency argument
//    void next_a(int inNumSamples)
//    {
//
//        float *outBuf = out(0);     // pointer to the output buffer
//        const float *freq = in(0);  // pointer to the input buffer
//
//        // get phase and freqmul constant from struct and store it in a local variable.
//        // The optimizer will cause them to be loaded it into a register.
////        float freqmul = mFreqMul;
////        double phase = mPhase;
//        int aDeg = mDeg.consume(); // ambisonic degree
//        int aIdx = 0; // ambisonic index
//        float theta = in0(2); //mTheta;
//        float phi =  in0(3); //mPhi;
//
//        // perform a loop for the number of samples in the control period.
//        // If this unit is audio rate then inNumSamples will be 64 or whatever
//        // the block size is. If this unit is control rate then inNumSamples will
//        // be 1.
//        for (int i=0; i < inNumSamples; ++i)
//        {
//            // unsigned n, int m, T theta, T phi
//            float sphHarmR = boost::math::spherical_harmonic_r(aDeg, aIdx, theta, phi);
//
//            // out must be written last for in place operation
////            float z = phase;
////            phase += freq[i] * freqmul;
//
//            // these if statements wrap the phase a +1 or -1.
////            if (phase >= 1.f) phase -= 2.f;
////            else if (phase <= -1.f) phase += 2.f;
//
//            // write the output
//            //            outBuf[i] = z;
//            outBuf[i] = sphHarmR;
//        }
//
//        // store the phase back to the struct
////        mPhase = phase;
//    }
//
//    //////////////////////////////////////////////////////////////////
//
//    // calculation function for a control rate frequency argument
//    void next_k(int inNumSamples)
//    {
////        // get the pointer to the output buffer
////        float *outBuf = out(0);
////
////        // freq is control rate, so calculate it once.
////        float freq = in0(0) * mFreqMul;
////
////        // get phase from struct and store it in a local variable.
////        // The optimizer will cause it to be loaded it into a register.
////        double phase = mPhase;
////
////        // since the frequency is not changing then we can simplify the loops
////        // by separating the cases of positive or negative frequencies.
////        // This will make them run faster because there is less code inside the loop.
////        if (freq >= 0.f) {
////            // positive frequencies
////            for (int i=0; i < inNumSamples; ++i)
////            {
////                outBuf[i] = phase;
////                phase += freq;
////                if (phase >= 1.f) phase -= 2.f;
////            }
////        } else {
////            // negative frequencies
////            for (int i=0; i < inNumSamples; ++i)
////            {
////                outBuf[i] = phase;
////                phase += freq;
////                if (phase <= -1.f) phase += 2.f;
////            }
////        }
////
////        // store the phase back to the struct
////        mPhase = phase;
//    }
};

// the entry point is called by the host when the plug-in is loaded
PluginLoad(AtkHoaUGens)
{
    // InterfaceTable *inTable implicitly given as argument to the load function
    ft = inTable; // store pointer to InterfaceTable
    
    // registerUnit takes the place of the Define*Unit functions.
    // It automatically checks for the presence of a destructor function.
    // However, it does not seem to be possible to disable buffer aliasing with the C++ header.
    registerUnit<HoaSphCoeff>(ft, "HoaSphCoeff");
}
