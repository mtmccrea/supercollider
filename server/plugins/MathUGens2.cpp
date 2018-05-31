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


static InterfaceTable *ft;

//////////////////////////////////////////////////////////////////////////////////////////////////

// Inherits from SCUnit, not Unit.
// Also note that the constructor, destructor, and calc functions are methods of the SCUnit.
struct FoldTest2 : public SCUnit {
public:
    // Constructor function
    FoldTest2() {
        // set calc function.
        // int numInputs = numInputs();

        if(bufferSize() == 1) {
    			// _aa? Well, yes - that calc func doesn't interpolate
    			// and interpolation is not needed for kr (1 sample/block)
    		set_calc_function<FoldTest2, &FoldTest2::next_aa>();
    	} else {
    		if(isAudioRateIn(1)) {
    			if(isAudioRateIn(2))
    				set_calc_function<FoldTest2, &FoldTest2::next_aa>();
    			else
    				set_calc_function<FoldTest2, &FoldTest2::next_ak>();
    		} else {
    			if(isAudioRateIn(2))
    				set_calc_function<FoldTest2, &FoldTest2::next_ka>();
    			else
    				set_calc_function<FoldTest2, &FoldTest2::next_kk>();
    		}
    	}

        m_lo = in0(1);
        m_hi = in0(2);

        // FoldTest_next_kk(unit, 1);
        next(1);
    }

    // If you want a destructor, you would declare "~FoldTest2() { ... }" here.

private:
    // You can declare state variables here.
    float m_lo, m_hi;

    /* Calc functions */

    // void FoldTest_next_aa(FoldTest* unit, int inNumSamples)
    // {
    //     float *out = ZOUT(0);
    //     float *in  = ZIN(0);
    //     float *lo = ZIN(1);
    //     float *hi = ZIN(2);
    //
    //     LOOP1(inNumSamples,
    //         float curhi = ZXP(hi);
    //         float curlo = ZXP(lo);
    //         float range = curhi - curlo;
    //         float range2 = range * 2.0;
    //         ZXP(out) = sc_fold(ZXP(in), curlo, curhi, range, range2);
    //     );
    // }

    void next_aa(int inNumSamples) {
        // ins are const float*, not float*.
        const float *in  = in(0); /// get input signal pointer at index
    	const float *lo = in(1);
    	const float *hi = in(2);

        float* outbuf = out(0);

        // TODO: use consume()

        for (int i = 0; i < inNumSamples; i++) {
            float curhi = hi[i];
    		float curlo = lo[i];
    		float range = curhi - curlo;
    		float range2 = range * 2.0;
            // ZXP(out) = sc_fold(ZXP(in), curlo, curhi, range, range2);
            outbuf[i] = sc_fold(in[i], curlo, curhi, range, range2);
        }
    }

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

    void next_kk(int inNumSamples)
    {
    	const float *in  = in(0);
    	const float next_lo = in0(1); /// get first sample of input signal at index
    	const float next_hi = in0(2);

        float *out = out(0);

    	float lo_slope = calcSlope(next_lo, m_lo);
    	float hi_slope = calcSlope(next_hi, m_hi);

    	for (int i = 0; i < inNumSamples; i++) {
    		float range = hi - lo;
    		float range2 = range * 2.f;

    		out[i] = sc_fold(in[i], lo, hi, range, range2);

            // TODO: becase this isn't using LOOP, should this
            // happen before setting out[i] or/and before calculating range?
    		lo += lo_slope;
    		hi += hi_slope;
    	);

        m_lo = lo;
    	m_hi = hi;
    }

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

    void next_ka(int inNumSamples)
    {
        const float *in  = in(0);
        const float next_lo = in0(1); /// get first sample of input signal at index
        const float *hi = in(2); /// get input signal pointer at index

        float *out = out(0);

        float lo_slope = calcSlope(next_lo, m_lo);

        for (int i = 0; i < inNumSamples; i++) {
            float curhi = hi[i];
            float range = curhi - lo;
            float range2 = range * 2.f;

            out[i] = sc_fold(in[i], lo, hi, range, range2);

            lo += lo_slope;
        );

        m_lo = lo; // only need to set previous lo, hi is AR so isn't interpolated
    }

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

    void next_ak(int inNumSamples)
    {
        const float *in = in(0);
        const float *lo = in(1);
        float next_hi = in0(2);

        float *out = out(0);

        float hi_slope = calcSlope(next_hi, m_hi);

        for (int i = 0; i < inNumSamples; i++) {
            float curlo = lo[i];
            float range = hi - curlo;
            float range2 = range * 2.f;

            out[i] = sc_fold(in[i], curlo, hi, range, range2);

            hi += hi_slope;
        );

        m_hi = hi;
    }

};


//////////////////////////////////////////////////////////////////////////////////////////////////

PluginLoad(Math)
{
	ft = inTable;

	registerUnit<BoringMixer2>(ft, "FoldTest2");
}
