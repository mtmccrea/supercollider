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
//		mPrevTheta = 0.0;
//		mPrevPhi = in0(3);
//		mCnt = 0;
//		mPrevBlockOut = 0.0;
//		set_calc_function<HoaSphCoeff,&HoaSphCoeff::next_ks>();
		
		// 2. initialize the unit generator state variables.
		mDeg = in0(0); // scalar, non-modulatable, could be ints
		mIdx = in0(1);

		mPrevTheta = in0(2);
		mPrevPhi = in0(3);

		// calc mPrevBlockOut now in case calc function is kk or ss (could move to those cases)
		mPrevBlockOut = boost::math::spherical_harmonic_r<float>(mDeg, mIdx, mPrevTheta, mPrevPhi);

		// 1. set the calculation function.
		switch (inRate(2)) {
			case calc_ScalarRate:
				switch (inRate(3)) {
					case calc_BufRate:
						Print("sk calc\n");
						set_calc_function<HoaSphCoeff,&HoaSphCoeff::next_sk>();
						// next_sk(1);
						break;
					case calc_FullRate:
						Print("sa calc\n");
						set_calc_function<HoaSphCoeff,&HoaSphCoeff::next_sa>();
						// next_sa(1);
						break;
					case calc_ScalarRate:
						Print("ss calc\n");
						set_calc_function<HoaSphCoeff,&HoaSphCoeff::next_ss>();
						// next_ss(1);
						break;
				}
				break;
			case calc_BufRate:
				switch (inRate(3)) {
					case calc_BufRate:
						Print("kk calc\n");
						set_calc_function<HoaSphCoeff,&HoaSphCoeff::next_kk>();
						// next_kk(1);
						break;
					case calc_FullRate:
						Print("ka calc\n");
						set_calc_function<HoaSphCoeff,&HoaSphCoeff::next_ka>();
						// next_ka(1);
						break;
					case calc_ScalarRate:
						Print("ks calc\n");
						set_calc_function<HoaSphCoeff,&HoaSphCoeff::next_ks>();
						// next_ks(1);
						break;
				}
				break;
			case calc_FullRate:
				switch (inRate(3)) {
					case calc_BufRate:
						Print("ak calc\n");
						set_calc_function<HoaSphCoeff,&HoaSphCoeff::next_ak>();
						// next_ak(1);
						break;
					case calc_FullRate:
						Print("aa calc\n");
						set_calc_function<HoaSphCoeff,&HoaSphCoeff::next_aa>();
						// next_aa(1);
						break;
					case calc_ScalarRate:
						Print("as calc\n");
						set_calc_function<HoaSphCoeff,&HoaSphCoeff::next_as>();
						// next_as(1);
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
	float mPrevBlockOut;
	// int mCnt;
	
	//////////////////////////////////////////////////////////////////
	
	bool inputsChanged(float theta, float phi)
	{
		bool tChanged = (theta != mPrevTheta);
		bool pChanged = (phi != mPrevPhi);
		
		if (tChanged) mPrevTheta = theta;
		if (pChanged) mPrevPhi = phi;
		
		return tChanged || pChanged;
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
				outBuf[i] = mPrevBlockOut;
			} else {

				// unsigned n, int m, T theta, T phi
				float sphHarmR = boost::math::spherical_harmonic_r<float>(mDeg, mIdx, newTheta, newPhi);

				outBuf[i] = sphHarmR;
				mPrevBlockOut = sphHarmR;
			}
		}
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
			outBuf[i] = mPrevBlockOut;
		}
	}
	
	// currently not used:
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

/////// SCRATCH /////////

//		template <typename SignalT, typename SignalP>
//		void generateOutput (int inNumSamples, SignalT theta, SignalP phi) {
//			float newTheta;
//			float newPhi;
//
//			float *outBuf = out(0);     // pointer to the output buffer
//
////			printf("mPrevBlockOut: %f\n", mPrevBlockOut);
//			printf("SlopeSignal: value(%f), slope(%f)\n", theta.value, theta.slope);
//			printf("samples:");
//
//			// inNumSamples will be blocksize (e.g. 64) for ar, 1 for kr.
//			for (int i=0; i < inNumSamples; ++i)
//			{
//				printf("%d", i);
//				newTheta = theta.consume();
//				newPhi = phi.consume();
//
//				// check if inputs have changed, if not, output previous calc
//				if (!inputsChanged(newTheta, newPhi)) {
//					// printf("no input change %d\n", i);
//					outBuf[i] = mPrevBlockOut;
//				} else {
//					outBuf[i] = newTheta;
//				}
//			}
//			// set prevOut as last output sample from this control period
//			mPrevBlockOut = newTheta;
////			printf("\n");
//		}


//    void build_calc_func(int inNumSamples)
//    {
//        AudioSignal<float> theta = makeSignal(2);
//        AudioSignal<float> phi = makeSignal(3);
//
//#define DefineOutputGenFunc(type1, type2, inNumSamples, theta, phi) \
//    HoaSphCoeff::generateOutput<type1, type2>(inNumSamples, theta, phi);
//
//        DefineOutputGenFunc(AudioSignal<float>, AudioSignal<float>, inNumSamples, theta, phi);
//    }


//	void next_ks(int inNumSamples)
//	{
//
//		printf("%d mPrevTheta: %f, mPrevOut: %f, in0(0): %f\n", mCnt, mPrevTheta, mPrevBlockOut, in0(0));
//		SlopeSignal<float> theta = makeSlope<float>(in0(0), mPrevTheta);
//		ScalarSignal<float> phi = makeScalar<float>(in0(1));
//
//		generateOutput<SlopeSignal<float>, ScalarSignal<float>>(inNumSamples, theta, phi);
//		mCnt++;
//	}

//    // slimmer version?
//    void next_aa(int inNumSamples, thetaSig, phiSig)
//    {
//        auto theta = makeSignal(2);
//        auto phi   = makeSignal(3);
//
//        generateOutput<decltype(theta), decltype(theta)>(inNumSamples, theta, phi);
//    }
