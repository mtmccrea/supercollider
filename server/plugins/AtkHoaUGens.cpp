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
struct HoaSphericalCoeff : public SCUnit{
	
	
	// Constructor usually does 3 things.
	// 1. set the calculation function.
	// 2. initialize the unit generator state variables.
	// 3. calculate one sample of output.
public:
	HoaSphericalCoeff() {
		// 1. set the calculation function.
		if (isAudioRateIn(2) || isAudioRateIn(3)) {
			// if the theta or phi argument are audio rate
			set_calc_function<HoaSphericalCoeff,&HoaSphericalCoeff::next_a>();
		} else {
			// if thene frequency argument is control rate (or a scalar).
			set_calc_function<HoaSphericalCoeff,&HoaSphericalCoeff::next_k>();
		}
		
		// 2. initialize the unit generator state variables.
		// initialize a constant for multiplying the frequency
		mFreqMul = 2.0 * sampleDur();
		// get initial phase of oscillator
		mPhase = in0(1);
		
		mDeg = in0(0);
		mIdx = in0(1);
		mTheta = in0(2);
		mPhi = in0(3);
		
		// 3. calculate one sample of output.
		if (isAudioRateIn(0)) {
			next_a(1);
		} else {
			next_k(1);
		}
		
	}
	
private:
	double mPhase;  // phase of the oscillator, from -1 to 1.
	float mFreqMul; // a constant for multiplying frequency
	int mDeg;       // ambisonic degree (l)
	int mIdx;       // ambisonic index  (m)
	double mTheta;
	double mPhi;
	
	//////////////////////////////////////////////////////////////////
	
	// The calculation function executes once per control period
	// which is typically 64 samples.
	
	// calculation function for an audio rate frequency argument
	void next_a(int inNumSamples)
	{
		
		float *outBuf = out(0);     // pointer to the output buffer
		const float *freq = in(0);  // pointer to the input buffer
		
		// get phase and freqmul constant from struct and store it in a local variable.
		// The optimizer will cause them to be loaded it into a register.
		float freqmul = mFreqMul;
		double phase = mPhase;
		int aDeg = mDeg; // ambisonic degree
		int aIdx = mIdx; // ambisonic index
		double theta = in0(2); //mTheta;
		double phi =  in0(3); //mPhi;
		
		// perform a loop for the number of samples in the control period.
		// If this unit is audio rate then inNumSamples will be 64 or whatever
		// the block size is. If this unit is control rate then inNumSamples will
		// be 1.
		for (int i=0; i < inNumSamples; ++i)
		{
			// unsigned n, int m, T theta, T phi
			float shr = boost::math::spherical_harmonic_r(aDeg, aIdx, theta, phi);
			
			// out must be written last for in place operation
			float z = phase;
			phase += freq[i] * freqmul;
			
			// these if statements wrap the phase a +1 or -1.
			if (phase >= 1.f) phase -= 2.f;
			else if (phase <= -1.f) phase += 2.f;
			
			// write the output
//			outBuf[i] = z;
			outBuf[i] = shr;
		}
		
		// store the phase back to the struct
		mPhase = phase;
	}
	
	//////////////////////////////////////////////////////////////////
	
	// calculation function for a control rate frequency argument
	void next_k(int inNumSamples)
	{
		// get the pointer to the output buffer
		float *outBuf = out(0);
		
		// freq is control rate, so calculate it once.
		float freq = in0(0) * mFreqMul;
		
		// get phase from struct and store it in a local variable.
		// The optimizer will cause it to be loaded it into a register.
		double phase = mPhase;
		
		// since the frequency is not changing then we can simplify the loops
		// by separating the cases of positive or negative frequencies.
		// This will make them run faster because there is less code inside the loop.
		if (freq >= 0.f) {
			// positive frequencies
			for (int i=0; i < inNumSamples; ++i)
			{
				outBuf[i] = phase;
				phase += freq;
				if (phase >= 1.f) phase -= 2.f;
			}
		} else {
			// negative frequencies
			for (int i=0; i < inNumSamples; ++i)
			{
				outBuf[i] = phase;
				phase += freq;
				if (phase <= -1.f) phase += 2.f;
			}
		}
		
		// store the phase back to the struct
		mPhase = phase;
	}
};

// the entry point is called by the host when the plug-in is loaded
PluginLoad(AtkHoaUGens)
{
	// InterfaceTable *inTable implicitly given as argument to the load function
	ft = inTable; // store pointer to InterfaceTable
	
	// registerUnit takes the place of the Define*Unit functions. It automatically checks for the presence of a
	// destructor function.
	// However, it does not seem to be possible to disable buffer aliasing with the C++ header.
	registerUnit<HoaSphericalCoeff>(ft, "HoaSphericalCoeff");
}
