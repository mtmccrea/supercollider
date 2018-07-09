#include "SC_PlugIn.h"
#include "SIMD_Unit.hpp"
#include <limits.h>
#include <cstdio>
#include "function_attributes.h"

#include <boost/align/is_aligned.hpp>

#define BOOST_MATH_DISABLE_FLOAT128 1
#include "boost/math/special_functions.hpp"
#include "boost/math/special_functions/spherical_harmonic.hpp"
#include "boost/math/special_functions/legendre.hpp"

#include <vector>


static InterfaceTable *ft;

//////////////////////////////////////////////////////////////////////////////////////////////////
// This version has a dynamic number of arguments

// Inherits from SCUnit, not Unit.
// Also note that the constructor, destructor, and calc functions are methods of the SCUnit.
struct NAryOpUGen : public SCUnit {

public:
	
	// Constructor function
	NAryOpUGen (*primFunc) {
		Print("NAryOpUGen_Ctor A\n");
		m_numInputs = numInputs();
		m_primFuncToCall = primFunc;
		
		// store indices of ar inputs and kr inputs
		// and store input rates
		for (int i =0; i < m_numInputs; i++) {
			int rate;
			
			if (isAudioRateIn(i)) {
				m_arIdxList.push_back(i);
				m_arCount++;
				rate = 0;
			} else if (isControlRateIn(i)) {
				m_krIdxList.push_back(i);
				m_interpSlopes.push_back(0.0); // initialize kr interp slopes to 0
				m_krCount++;
				rate = 1;
			} else {
				rate = 2; // TODO: scalar and demand rates
			}
			
			m_inRates.push_back(rate); // add this rate to m_inRates
		}
		
		// debug: getting proper num inputs and rates?
		printf("NAryOpUGen numInputs %d\n", m_numInputs);
		printf("NAryOpUGen rates [");
		for (int i =0; i < m_numInputs; i++) {
			printf(" %d ", m_inRates[i]);
		}
		printf("]\n");
		
		// set calc function.
		if(bufferSize() == 1) {
			set_calc_function<NAryOpUGen, &NAryOpUGen::next_k>();
		} else {
			if(m_arCount > 0) {
				set_calc_function<NAryOpUGen, &NAryOpUGen::next_a>();
			} else {
				set_calc_function<NAryOpUGen, &NAryOpUGen::next_k>();
			}
		}
		
		next_k(1);    // TODO: apppropriate to call next_k here to calc first sample?
	}
	
	// If you want a destructor, you would declare "~NAryOpUGen() { ... }" here.
	
private:
	// declare state variables here.
	int m_numInputs, m_arCount=0, m_krCount=0;
	std::vector<int> m_arIdxList;
	std::vector<int> m_krIdxList;
	std::vector<int> m_inRates;
	std::vector<double> m_prevInputs; // TODO: will inputs always be doubles? floats? ints?
	std::vector<double> m_interpSlopes; // interpSlopes will be of size m_krCount, storing slopes for input(m_krIdxList[i])
	double m_prevOutput = 0.0;
	NAryOpFunc * m_primFuncToCall;
	
	
	inline int checkInputsChanged (int frame)
	{
		// iterate over all the inputs, check if their values changed
		for (int i =0; i < m_numInputs; i++) {
			double newVal;
			// get rate for this input
			if (m_inRates[i] == 0) {         // kr
				newVal = in0(i);             // frame == 0 for kr
			} else if (m_inRates[i] == 1) {  // ar
				newVal = in(i)[frame];       // TODO: need proper index in the loop
			} else {
				// shouldn't get here
				// TODO: scalar, demand
				newVal = 0.0;
			}
			
			if (m_prevInputs[i] != newVal)
				return i;  // changedAt index i, return now, signaling new calc
		}
		
		return -1;	// inputs didn't change
	}
	
	inline int krInputsChanged ()
	{
		// iterate over all the inputs, check if their values changed
		for (int i=0; i<m_krCount; i++) {
			
			int inputIdx = m_krIdxList[i];
			double newVal = in0(inputIdx);
			
			if (m_prevInputs[inputIdx] != newVal)
				return i;  // changedAt index i of krIndxList, return now, signaling new calc
		}
		
		return -1;	// inputs didn't change
	}
	
	inline int arInputsChanged (int frame)
	{
		// iterate over all the inputs, check if their values changed
		for (int i=0; i<m_arCount; i++) {
			
			int inputIdx = m_arIdxList[i];
			double newVal = in(inputIdx)[frame];
			
			if (m_prevInputs[inputIdx] != newVal)
				return i;  // changedAt index i of arIndxList, return now, signaling new calc
		}
		
		return -1;	// inputs didn't change
	}
		
	/* Calc functions */
	
	// calc func for audio-rate output, all audio-rate inputs
	void next_a(int inNumSamples) {
		float* outbuf = out(0);
		
		// TODO: use consume()
		
		for (int i = 0; i < inNumSamples; i++) {
			int changedAt = arInputsChanged(i);       // input changed at this index of the m_arIdxList (not input index)
			
			if (changedAt > -1) {                     // at least one input has changed, need to call the op with new inputs
				for (int j=changedAt; j < m_arCount; j++) {
					int inputIdx = m_arIdxList[j];
					m_prevInputs[inputIdx] = in(inputIdx)[i];       // new inputs become previous
				}
				
//				sc_fold(inbuf[i], curlo, curhi, range, range2);
				auto z = primToCall(&m_prevInputs...);  // TODO: how to unpack?, need template for each number of args?
				outbuf[i] = (double)z;                // TODO: double?
				m_prevOutput = (double)z;
			} else {                                  /// inputs haven't changed, no need to call the op
				outbuf[i] = m_prevOutput;             // TODO: unnecessary to re-write the same value here?
			}
		}
	}
	
	// calc func for control-rate output, all control-rate inputs
	void next_k(int inNumSamples) {
		float outbuf = out0(0);  // TODO: output pointer for kr? "reference to first sample of output signal" (?)
		
		for (int i = 0; i < inNumSamples; i++) {     // inNumSamples = 1 for kr
			int changedAt = krInputsChanged();       // input changed at this index of the m_arIdxList (not input index)
			
			// TODO: use consume()
			
			changedAt = checkInputsChanged(0);
			
			if (changedAt > -1) {   /// at least one input has changed, need to call the op with new inputs
				for (int j=changedAt; j < m_krCount; j++) {
					int inputIdx = m_krIdxList[j];
					m_prevInputs[inputIdx] = in0(inputIdx);       // new inputs become previous
				}
				
				auto z = primToCall(&m_prevInputs...);  // TODO: how to unroll?, need template for each number of args?
				outbuf[i] = (double)z;                // TODO: is this properly writing to the "reference" to the out signal? TODO: double?
				m_prevOutput = (double)z;
			} else {
				// inputs haven't changed, no need to call the op
				outbuf[i] = m_prevOutput;             // TODO: unnecessary to re-write the same value here?
			}
		}
	}
	
	// calc func for audio-rate output, mixed control and audio-rate inputs
	void next_ka(int inNumSamples) {
		float* outbuf = out(0);
		int krChangedAt; // index of first changed kr arg
		int arChangedAt; // index of first changed kr arg
		bool krChanged;
		bool recalc = false;
		
		// TODO: use consume()
		
		// zero out interp slopes
		// NOTE: interp slope for a-rate is always 0
		for (int i=0; i<m_numInputs; i++) {
			m_interpSlopes[i] = 0.0;
		}
		
		// check if kr inputs have changed, if so, calc slope
		krChangedAt = krInputsChanged();
		krChanged = krChangedAt > -1;
		if (krChanged) {
			for (int i=krChangedAt; i<m_krCount; i++) {
				int inputIdx = m_krIdxList[i];
				double next = in0(inputIdx);
				double prev = m_prevInputs[inputIdx];
				// TODO: we already know that krChangedAt is a changed value,
				// the following check is really only needed to check for subsequent input values
				if (next != prev) {
					// TODO: worth keeping track of only those that have changed?
					m_interpSlopes[i] = calcSlope(next, prev);
				}
			}
		}
		
		for (int i = 0; i < inNumSamples; i++) {
			if (krChanged) {
				// inputs are changing, need to call op
				//				increment kr by slope
				for (int j=0; j<m_krCount; j++) {
					int inputIdx = m_krIdxList[j];
					m_prevInputs[inputIdx] += m_interpSlopes[j];
				}
				recalc = true;
			}
			
			if (arChangedAt > -1) {
				for (int j=0; j<m_arCount; j++) {
					int inputIdx = m_arIdxList[j];
					m_prevInputs[inputIdx] += in(inputIdx)[i];
				}
				recalc = true;
			}
			
			if (recalc) {
				auto z = primToCall(&m_prevInputs...);  // TODO: how to unroll?, need template for each number of args?
				outbuf[i] = (double)z;                // TODO: double?
				m_prevOutput = (double)z;
			} else {
				outbuf[i] = m_prevOutput;  			  // TODO: unnecessary?
			}
			
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////

PluginLoad(Math)
{
	ft = inTable;
	
	registerUnit<NAryOpUGen>(ft, "NAryOpUGen");
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//
//// Inherits from SCUnit, not Unit.
//// Also note that the constructor, destructor, and calc functions are methods of the SCUnit.
//struct FoldTest2 : public SCUnit {
//public:
//    // Constructor function
//    FoldTest2() {
//        // set calc function.
//        // int numInputs = numInputs();
//		//Print("FoldTest2_Ctor A\n");
//
//        if(bufferSize() == 1) {
//    			// _aa? Well, yes - that calc func doesn't interpolate
//    			// and interpolation is not needed for kr (1 sample/block)
//    		set_calc_function<FoldTest2, &FoldTest2::next_aa>();
//    	} else {
//    		if(isAudioRateIn(1)) {
//    			if(isAudioRateIn(2))
//    				set_calc_function<FoldTest2, &FoldTest2::next_aa>();
//    			else
//    				set_calc_function<FoldTest2, &FoldTest2::next_ak>();
//    		} else {
//    			if(isAudioRateIn(2))
//    				set_calc_function<FoldTest2, &FoldTest2::next_ka>();
//    			else
//    				set_calc_function<FoldTest2, &FoldTest2::next_kk>();
//    		}
//    	}
//
//        m_lo = in0(1);
//        m_hi = in0(2);
//
//        // FoldTest_next_kk(unit, 1);
//		// TODO: apppropriate to call next_kk here to calc first sample?
//        next_kk(1);
//    }
//
//    // If you want a destructor, you would declare "~FoldTest2() { ... }" here.
//
//private:
//    // You can declare state variables here.
//    float m_lo, m_hi;
//
//    /* Calc functions */
//
//    // void FoldTest_next_aa(FoldTest* unit, int inNumSamples)
//    // {
//    //     float *out = ZOUT(0);
//    //     float *in  = ZIN(0);
//    //     float *lo = ZIN(1);
//    //     float *hi = ZIN(2);
//    //
//    //     LOOP1(inNumSamples,
//    //         float curhi = ZXP(hi);
//    //         float curlo = ZXP(lo);
//    //         float range = curhi - curlo;
//    //         float range2 = range * 2.0;
//    //         ZXP(out) = sc_fold(ZXP(in), curlo, curhi, range, range2);
//    //     );
//    // }
//
//    void next_aa(int inNumSamples) {
//        // ins are const float*, not float*.
//        const float* inbuf = in(0); /// get input signal pointer at index
//    	const float* lo = in(1);
//    	const float* hi = in(2);
//
//        float* outbuf = out(0);
//
//        // TODO: use consume()
//
//        for (int i = 0; i < inNumSamples; i++) {
//            float curhi = hi[i];
//    		float curlo = lo[i];
//    		float range = curhi - curlo;
//    		float range2 = range * 2.0;
//            // ZXP(out) = sc_fold(ZXP(in), curlo, curhi, range, range2);
//            outbuf[i] = sc_fold(inbuf[i], curlo, curhi, range, range2);
//        }
//    }
//
//    // void FoldTest_next_kk(FoldTest* unit, int inNumSamples)
//    // {
//    // 	float *out = ZOUT(0);
//    // 	float *in  = ZIN(0);
//    // 	float next_lo = ZIN0(1);
//    // 	float next_hi = ZIN0(2);
//    // 	float lo = unit->m_lo;
//    // 	float lo_slope = CALCSLOPE(next_lo, lo);
//    // 	float hi = unit->m_hi;
//    // 	float hi_slope = CALCSLOPE(next_hi, hi);
//    //
//    // 	LOOP1(inNumSamples,
//    // 		float range = hi - lo;
//    // 		float range2 = range * 2.f;
//    // 		ZXP(out) = sc_fold(ZXP(in), lo, hi, range, range2);
//    //
//    // 		lo += lo_slope;
//    // 		hi += hi_slope;
//    // 	);
//    // 	unit->m_lo = lo;
//    // 	unit->m_hi = hi;
//    // }
//
//    void next_kk(int inNumSamples)
//    {
//    	const float* inbuf  = in(0);
//    	const float next_lo = in0(1); /// get first sample of input signal at index
//    	const float next_hi = in0(2);
//
//        float* outbuf = out(0);
//
//    	float lo_slope = calcSlope(next_lo, m_lo);
//    	float hi_slope = calcSlope(next_hi, m_hi);
//
//    	for (int i = 0; i < inNumSamples; i++) {
//    		float range = m_hi - m_lo;
//    		float range2 = range * 2.f;
//
//    		outbuf[i] = sc_fold(inbuf[i], m_lo, m_hi, range, range2);
//
//            // TODO: becase this isn't using LOOP, should this
//            // happen before setting out[i] or/and before calculating range?
//    		m_lo += lo_slope;
//    		m_hi += hi_slope;
//		};
//    }
//
//    // void FoldTest_next_ka(FoldTest* unit, int inNumSamples)
//    // {
//    // 	float *out = ZOUT(0);
//    // 	float *in  = ZIN(0);
//    // 	float next_lo = ZIN0(1);
//    // 	float *hi = ZIN(2);
//    // 	float lo = unit->m_lo;
//    // 	float lo_slope = CALCSLOPE(next_lo, lo);
//    //
//    // 	LOOP1(inNumSamples,
//    // 		float curhi = ZXP(hi);
//    // 		float range = curhi - lo;
//    // 		float range2 = range * 2.f;
//    // 		ZXP(out) = sc_fold(ZXP(in), lo, curhi, range, range2);
//    // 		lo += lo_slope;
//    // 	);
//    // 	unit->m_lo = lo;
//    // }
//
//    void next_ka(int inNumSamples)
//    {
//        const float* inbuf = in(0);
//        const float next_lo = in0(1); /// get first sample of input signal at index
//        const float* hi = in(2); /// get input signal pointer at index
//
//        float* outbuf = out(0);
//
//        float lo_slope = calcSlope(next_lo, m_lo);
//
//        for (int i = 0; i < inNumSamples; i++) {
//            float curhi = hi[i];
//            float range = curhi - m_lo;
//            float range2 = range * 2.f;
//
//            outbuf[i] = sc_fold(inbuf[i], m_lo, curhi, range, range2);
//
//            m_lo += lo_slope; // only need to set previous lo, hi is AR so isn't interpolated
//		};
//    }
//
//    // void FoldTest_next_ak(FoldTest* unit, int inNumSamples)
//    // {
//    // 	float *out = ZOUT(0);
//    // 	float *in  = ZIN(0);
//    // 	float *lo = ZIN(1);
//    // 	float next_hi = ZIN0(2);
//    // 	float hi = unit->m_hi;
//    // 	float hi_slope = CALCSLOPE(next_hi, hi);
//    //
//    // 	LOOP1(inNumSamples,
//    // 		float curlo = ZXP(lo);
//    // 		float range = hi - curlo;
//    // 		float range2 = range * 2.f;
//    // 		ZXP(out) = sc_fold(ZXP(in), curlo, hi, range, range2);
//    // 		hi += hi_slope;
//    // 	);
//    // 	unit->m_hi = hi;
//    // }
//
//    void next_ak(int inNumSamples)
//    {
//        const float* inbuf = in(0);
//        const float* lo = in(1);
//        float next_hi = in0(2);
//
//        float* outbuf = out(0);
//
//        float hi_slope = calcSlope(next_hi, m_hi);
//
//        for (int i = 0; i < inNumSamples; i++) {
//            float curlo = lo[i];
//            float range = m_hi - curlo;
//            float range2 = range * 2.f;
//
//            outbuf[i] = sc_fold(inbuf[i], curlo, m_hi, range, range2);
//
//            m_hi += hi_slope;
//		};
//    }
//
//};
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//PluginLoad(Math)
//{
//	ft = inTable;
//
//	registerUnit<FoldTest2>(ft, "FoldTest2");
//}


// #include "SC_PlugIn.h"
// #include "SIMD_Unit.hpp"
// #include <limits.h>
// #include <cstdio>
// #include "function_attributes.h"
//
// #include <boost/align/is_aligned.hpp>
//
// #define BOOST_MATH_DISABLE_FLOAT128 1
// #include "boost/math/special_functions.hpp"
// #include "boost/math/special_functions/spherical_harmonic.hpp"
// #include "boost/math/special_functions/legendre.hpp"
//
//
// static InterfaceTable *ft;
//
// //////////////////////////////////////////////////////////////////////////////////////////////////
//
// struct FoldTest : public Unit
// {
// 	float m_lo, m_hi, m_range;
// };
//
// // SphericalHarmonicReal
// struct SHR : public Unit
// {
// 	float m_m, m_theta, m_phi; // TODO: m_m needs to be an INT
// };
//
// //////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// void FoldTest_next_kk(FoldTest *unit, int inNumSamples);
// void FoldTest_next_ak(FoldTest *unit, int inNumSamples);
// void FoldTest_next_ka(FoldTest *unit, int inNumSamples);
// void FoldTest_next_aa(FoldTest *unit, int inNumSamples);
// void FoldTest_Ctor(FoldTest* unit);
//
// void SHR_next_kk(SHR *unit, int inNumSamples);
// void SHR_Ctor(SHR *unit);
//
// void SHR_next_kk(SHR *unit, int inNumSamples)
// {
// 	float *out = ZOUT(0);
// 	float *in  = ZIN(0);
//
// 	float next_m = ZIN0(1);
// 	float next_theta = ZIN0(2);
// 	float next_phi = ZIN0(3);
//
// 	// don't interpolate m
// 	float theta = unit->m_theta;
// 	float theta_slope = CALCSLOPE(next_theta, theta);
// 	float phi = unit->m_phi;
// 	float phi_slope = CALCSLOPE(next_phi, phi);
//
// 	LOOP1(inNumSamples,
// 		  //		ZXP(out) = sc_fold(ZXP(in), lo, hi, range, range2);
// 		  ZXP(out) = boost::math::spherical_harmonic_r(ZXP(in), next_m, theta, phi);
//
// 		  theta += theta_slope;
// 		  phi += phi_slope;
// 		  );
// 	unit->m_theta = theta;
// 	unit->m_phi = phi;
// }
//
// void SHR_Ctor(SHR* unit)
// {
// 	SETCALC(FoldTest_next_kk);
//
// 	unit->m_m = ZIN0(1);
// 	unit->m_theta = ZIN0(2);
// 	unit->m_phi = ZIN0(3);
//
// 	SHR_next_kk(unit, 1);
// }
//
//
// //void LFSaw_next_a(LFSaw *unit, int inNumSamples)
// //{
// //	float *out = ZOUT(0);
// //	float *freq = ZIN(0);
// //
// //	float freqmul = unit->mFreqMul;
// //	double phase = unit->mPhase;
// //	LOOP1(inNumSamples,
// //		  float z = phase; // out must be written last for in place operation
// //		  phase += ZXP(freq) * freqmul;
// //		  if (phase >= 1.f) phase -= 2.f;
// //		  else if (phase <= -1.f) phase += 2.f;
// //		  ZXP(out) = z;
// //		  );
// //
// //	unit->mPhase = phase;
// //}
// //void LFSaw_next_k(LFSaw *unit, int inNumSamples)
// //{
// //	float *out = ZOUT(0);
// //	float freq = ZIN0(0) * unit->mFreqMul;
// //
// //	double phase = unit->mPhase;
// //	if (freq >= 0.f) {
// //		LOOP1(inNumSamples,
// //			  ZXP(out) = phase;
// //			  phase += freq;
// //			  if (phase >= 1.f) phase -= 2.f;
// //			  );
// //	} else {
// //		LOOP1(inNumSamples,
// //			  ZXP(out) = phase;
// //			  phase += freq;
// //			  if (phase <= -1.f) phase += 2.f;
// //			  );
// //	}
// //
// //	unit->mPhase = phase;
// //}
// //
// //void LFSaw_Ctor(LFSaw* unit)
// //{
// //	if (INRATE(0) == calc_FullRate)
// //		SETCALC(LFSaw_next_a);
// //	else
// //		SETCALC(LFSaw_next_k);
// //
// //	unit->mFreqMul = 2.0 * unit->mRate->mSampleDur;
// //	unit->mPhase = ZIN0(1);
// //
// //	LFSaw_next_k(unit, 1);
// //}
//
// //////////////////////////////////////////////////////////////////////////////////////////////////
//
// void FoldTest_next_kk(FoldTest* unit, int inNumSamples)
// {
// 	float *out = ZOUT(0);
// 	float *in  = ZIN(0);
// 	float next_lo = ZIN0(1);
// 	float next_hi = ZIN0(2);
// 	float lo = unit->m_lo;
// 	float lo_slope = CALCSLOPE(next_lo, lo);
// 	float hi = unit->m_hi;
// 	float hi_slope = CALCSLOPE(next_hi, hi);
//
// 	LOOP1(inNumSamples,
// 		float range = hi - lo;
// 		float range2 = range * 2.f;
// 		ZXP(out) = sc_fold(ZXP(in), lo, hi, range, range2);
//
// 		lo += lo_slope;
// 		hi += hi_slope;
// 	);
// 	unit->m_lo = lo;
// 	unit->m_hi = hi;
// }
//
// void FoldTest_next_ka(FoldTest* unit, int inNumSamples)
// {
// 	float *out = ZOUT(0);
// 	float *in  = ZIN(0);
// 	float next_lo = ZIN0(1);
// 	float *hi = ZIN(2);
// 	float lo = unit->m_lo;
// 	float lo_slope = CALCSLOPE(next_lo, lo);
//
// 	LOOP1(inNumSamples,
// 		float curhi = ZXP(hi);
// 		float range = curhi - lo;
// 		float range2 = range * 2.f;
// 		ZXP(out) = sc_fold(ZXP(in), lo, curhi, range, range2);
// 		lo += lo_slope;
// 	);
// 	unit->m_lo = lo;
// }
//
// void FoldTest_next_ak(FoldTest* unit, int inNumSamples)
// {
// 	float *out = ZOUT(0);
// 	float *in  = ZIN(0);
// 	float *lo = ZIN(1);
// 	float next_hi = ZIN0(2);
// 	float hi = unit->m_hi;
// 	float hi_slope = CALCSLOPE(next_hi, hi);
//
// 	LOOP1(inNumSamples,
// 		float curlo = ZXP(lo);
// 		float range = hi - curlo;
// 		float range2 = range * 2.f;
// 		ZXP(out) = sc_fold(ZXP(in), curlo, hi, range, range2);
// 		hi += hi_slope;
// 	);
// 	unit->m_hi = hi;
// }
//
// void FoldTest_next_aa(FoldTest* unit, int inNumSamples)
// {
// 	float *out = ZOUT(0);
// 	float *in  = ZIN(0);
// 	float *lo = ZIN(1);
// 	float *hi = ZIN(2);
//
// 	LOOP1(inNumSamples,
// 		float curhi = ZXP(hi);
// 		float curlo = ZXP(lo);
// 		float range = curhi - curlo;
// 		float range2 = range * 2.0;
// 		ZXP(out) = sc_fold(ZXP(in), curlo, curhi, range, range2);
// 	);
// }
//
// void FoldTest_Ctor(FoldTest* unit)
// {
// 	if(BUFLENGTH == 1) {
// 			// _aa? Well, yes - that calc func doesn't interpolate
// 			// and interpolation is not needed for kr (1 sample/block)
// 		SETCALC(FoldTest_next_aa);
// 	} else {
// 		if(INRATE(1) == calc_FullRate) {
// 			if(INRATE(2) == calc_FullRate)
// 				SETCALC(FoldTest_next_aa);
// 			else
// 				SETCALC(FoldTest_next_ak);
// 		} else {
// 			if(INRATE(2) == calc_FullRate)
// 				SETCALC(FoldTest_next_ka);
// 			else
// 				SETCALC(FoldTest_next_kk);
// 		}
// 	}
//
// 	unit->m_lo = ZIN0(1);
// 	unit->m_hi = ZIN0(2);
//
// 	FoldTest_next_kk(unit, 1);
// }
//
// //////////////////////////////////////////////////////////////////////////////////////////////////
//
// PluginLoad(Math)
// {
// 	ft = inTable;
//
// 	DefineSimpleUnit(SHR);
// 	DefineSimpleUnit(FoldTest);
//
// }
