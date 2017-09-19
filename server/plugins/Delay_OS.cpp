// #include <SC_PlugIn.h>
// #include <SC_InterfaceTable.h>
// #include <cstdio>
//
// #include "Oversamp.h"

#include "SC_PlugIn.h"
// #include "SC_InterfaceTable.h"
#include "Oversamp.h"

#include <cstdio>



static InterfaceTable *ft;

struct DelayUnit_OS : public Unit
{
	// from Delay Unit struct
	float *m_dlybuf;

	float m_dsamp, m_fdelaylen;
	float m_delaytime, m_maxdelaytime;
	long m_iwrphase, m_idelaylen, m_mask;
	long m_numoutput;
};

struct DelayC_OS8 : public DelayUnit_OS
{
	// from Oversamp8 struct
	float m_upmem[8];
	float m_domem[584];
};

struct DelayC_OS4 : public DelayUnit_OS
{
	float m_upmem[8];
	float m_domem[302];
};

struct DelayC_OS3 : public DelayUnit_OS
{
	float m_upmem[8];
	float m_domem[238];
};

struct DelayL_OS8 : public DelayUnit_OS
{
	float m_upmem[8];
	float m_domem[584];
};

struct DelayL_OS4 : public DelayUnit_OS
{
	float m_upmem[8];
	float m_domem[302];
};

struct DelayL_OS3 : public DelayUnit_OS
{
	float m_upmem[8];
	float m_domem[238];
};

extern "C"  {

	void load(InterfaceTable *inTable);

	void DelayUnit_OS_AllocDelayLine(DelayUnit_OS *unit, int osSize);
	void DelayUnit_OS_Reset(DelayUnit_OS *unit, int osSize);

//	void DelayUnit_OS8_Dtor(DelayUnit_OS *unit);

	void DelayC_OS8_Ctor(DelayC_OS8 *unit);
	void DelayC_OS8_Dtor(DelayC_OS8 *unit);
	void DelayC_OS8_next(DelayC_OS8 *unit, int inNumSamples);
	void DelayC_OS8_next_z(DelayC_OS8 *unit, int inNumSamples);


	void DelayC_OS4_Ctor(DelayC_OS4 *unit);
	void DelayC_OS4_Dtor(DelayC_OS4 *unit);
	void DelayC_OS4_next(DelayC_OS4 *unit, int inNumSamples);
	void DelayC_OS4_next_z(DelayC_OS4 *unit, int inNumSamples);

	void DelayC_OS3_Ctor(DelayC_OS3 *unit);
	void DelayC_OS3_Dtor(DelayC_OS3 *unit);
	void DelayC_OS3_next(DelayC_OS3 *unit, int inNumSamples);
	void DelayC_OS3_next_z(DelayC_OS3 *unit, int inNumSamples);


	void DelayL_OS8_Ctor(DelayL_OS8 *unit);
	void DelayL_OS8_next(DelayL_OS8 *unit, int inNumSamples);
	void DelayL_OS8_next_z(DelayL_OS8 *unit, int inNumSamples);

	void DelayL_OS4_Ctor(DelayL_OS4 *unit);
	void DelayL_OS4_next(DelayL_OS4 *unit, int inNumSamples);
	void DelayL_OS4_next_z(DelayL_OS4 *unit, int inNumSamples);

	void DelayL_OS3_Ctor(DelayL_OS3 *unit);
	void DelayL_OS3_next(DelayL_OS3 *unit, int inNumSamples);
	void DelayL_OS3_next_z(DelayL_OS3 *unit, int inNumSamples);
}

/* Functions used by all Delay_OS UGens here */
void DelayUnit_OS_AllocDelayLine(DelayUnit_OS *unit, int osSize)
{
	long delaybufsize = (long)ceil(unit->m_maxdelaytime * SAMPLERATE * osSize + 1.f);
	delaybufsize = delaybufsize + BUFLENGTH;
	delaybufsize = NEXTPOWEROFTWO(delaybufsize);  // round up to next power of two
	unit->m_fdelaylen = unit->m_idelaylen = delaybufsize;

	RTFree(unit->mWorld, unit->m_dlybuf);
	unit->m_dlybuf = (float*)RTAlloc(unit->mWorld, delaybufsize * sizeof(float));
	//bzero(unit->m_dlybuf, delaybufsize * sizeof(float)); // the less-efficient approach
	unit->m_mask = delaybufsize - 1;
}

float CalcDelay_OS(DelayUnit_OS *unit, float delaytime, int osSize);
float CalcDelay_OS(DelayUnit_OS *unit, float delaytime, int osSize)
{
	float next_dsamp = delaytime * SAMPLERATE * osSize;
	return sc_clip(next_dsamp, 1.f, unit->m_fdelaylen);
}

void DelayUnit_OS_Reset(DelayUnit_OS *unit, int osSize)
{
	unit->m_maxdelaytime = ZIN0(1);
	unit->m_delaytime = ZIN0(2);
	unit->m_dlybuf = 0;
	DelayUnit_OS_AllocDelayLine(unit, osSize);
	unit->m_dsamp = CalcDelay_OS(unit, unit->m_delaytime, osSize);
	unit->m_numoutput = 0;
	unit->m_iwrphase = 0;
}



//////////////////////////
//////////////////////////
/* Ctors & Dtors */
//////////////////////////
//////////////////////////

void DelayC_OS8_Ctor(DelayC_OS8 *unit)
{
	OVERSAMPLE8_INIT;
	DelayUnit_OS_Reset(unit, 8);
	SETCALC(DelayC_OS8_next_z);
	ZOUT0(0) = 0.f;
}

void DelayC_OS4_Ctor(DelayC_OS4 *unit)
{
	OVERSAMPLE4_INIT;
	DelayUnit_OS_Reset(unit, 4);
	SETCALC(DelayC_OS4_next_z);
	ZOUT0(0) = 0.f;
}

void DelayC_OS3_Ctor(DelayC_OS3 *unit)
{
	OVERSAMPLE3_INIT;
	DelayUnit_OS_Reset(unit, 3);
	SETCALC(DelayC_OS3_next_z);
	ZOUT0(0) = 0.f;
}

void DelayL_OS8_Ctor(DelayL_OS8 *unit)
{
	OVERSAMPLE8_INIT;
	DelayUnit_OS_Reset(unit, 8);
	SETCALC(DelayL_OS8_next_z);
	ZOUT0(0) = 0.f;
}

void DelayL_OS4_Ctor(DelayL_OS4 *unit)
{
	OVERSAMPLE4_INIT;
	DelayUnit_OS_Reset(unit, 4);
	SETCALC(DelayL_OS4_next_z);
	ZOUT0(0) = 0.f;
}

void DelayL_OS3_Ctor(DelayL_OS3 *unit)
{
	OVERSAMPLE3_INIT;
	DelayUnit_OS_Reset(unit, 3);
	SETCALC(DelayL_OS3_next_z);
	ZOUT0(0) = 0.f;
}

void DelayC_OS8_Dtor(DelayC_OS8 *unit)
{
	RTFree(unit->mWorld, unit->m_dlybuf);
}
void DelayC_OS4_Dtor(DelayC_OS4 *unit)
{
	RTFree(unit->mWorld, unit->m_dlybuf);
}
void DelayC_OS3_Dtor(DelayC_OS3 *unit)
{
	RTFree(unit->mWorld, unit->m_dlybuf);
}
void DelayL_OS8_Dtor(DelayL_OS8 *unit)
{
	RTFree(unit->mWorld, unit->m_dlybuf);
}
void DelayL_OS4_Dtor(DelayL_OS4 *unit)
{
	RTFree(unit->mWorld, unit->m_dlybuf);
}
void DelayL_OS3_Dtor(DelayL_OS3 *unit)
{
	RTFree(unit->mWorld, unit->m_dlybuf);
}


////////////////////////
////////////////////////
/*   NEXT Functions   */
////////////////////////
////////////////////////

void DelayC_OS8_next(DelayC_OS8 *unit, int inNumSamples)
{

	float *out = ZOUT(0);
	float *in = ZIN(0);
	float delaytime = ZIN0(2);

	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	long mask = unit->m_mask;
	int upbufdx = 0;

	UPSAMPLE8;

	if (delaytime == unit->m_delaytime) {	// if the delay time hasn't changed
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		LOOP1(inNumSamples * 8,
			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  long irdphase1 = iwrphase - idsamp;
			  long irdphase2 = irdphase1 - 1;
			  long irdphase3 = irdphase1 - 2;
			  long irdphase0 = irdphase1 + 1;
			  float d0 = dlybuf[irdphase0 & mask];
			  float d1 = dlybuf[irdphase1 & mask];
			  float d2 = dlybuf[irdphase2 & mask];
			  float d3 = dlybuf[irdphase3 & mask];
			  domemoff[upbufdx] = cubicinterp(frac, d0, d1, d2, d3);
			  iwrphase++;
			  upbufdx++;
			  );
	} else {

		float next_dsamp = CalcDelay_OS(unit, delaytime, 8);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp) * 0.125;
		LOOP1(inNumSamples * 8,
			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  dsamp += dsamp_slope;
			  long idsamp = (long)dsamp;
			  float frac = dsamp - idsamp;
			  long irdphase1 = iwrphase - idsamp;
			  long irdphase2 = irdphase1 - 1;
			  long irdphase3 = irdphase1 - 2;
			  long irdphase0 = irdphase1 + 1;
			  float d0 = dlybuf[irdphase0 & mask];
			  float d1 = dlybuf[irdphase1 & mask];
			  float d2 = dlybuf[irdphase2 & mask];
			  float d3 = dlybuf[irdphase3 & mask];
			  domemoff[upbufdx] = cubicinterp(frac, d0, d1, d2, d3);
			  iwrphase++;
			  upbufdx++;
			  );
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
	}

	DOWNSAMPLE8;

	unit->m_iwrphase = iwrphase;

}


void DelayC_OS8_next_z(DelayC_OS8 *unit, int inNumSamples)
{

	float *out = ZOUT(0);
	float *in = ZIN(0);
	float delaytime = ZIN0(2);

	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	long mask = unit->m_mask;
	float d0, d1, d2, d3;
	int upbufdx = 0;

	UPSAMPLE8;

	if (delaytime == unit->m_delaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		LOOP1(inNumSamples * 8,
			  long irdphase1 = iwrphase - idsamp;
			  long irdphase2 = irdphase1 - 1;
			  long irdphase3 = irdphase1 - 2;
			  long irdphase0 = irdphase1 + 1;

			  dlybuf[iwrphase & mask] = domemoff[upbufdx];

			  if (irdphase0 < 0) {
				  domemoff[upbufdx] = 0.f;
			  } else {
				  if (irdphase1 < 0) {
					  d1 = d2 = d3 = 0.f;
					  d0 = dlybuf[irdphase0 & mask];
				  } else if (irdphase2 < 0) {
					  d1 = d2 = d3 = 0.f;
					  d0 = dlybuf[irdphase0 & mask];
					  d1 = dlybuf[irdphase1 & mask];
				  } else if (irdphase3 < 0) {
					  d3 = 0.f;
					  d0 = dlybuf[irdphase0 & mask];
					  d1 = dlybuf[irdphase1 & mask];
					  d2 = dlybuf[irdphase2 & mask];
				  } else {
					  d0 = dlybuf[irdphase0 & mask];
					  d1 = dlybuf[irdphase1 & mask];
					  d2 = dlybuf[irdphase2 & mask];
					  d3 = dlybuf[irdphase3 & mask];
				  }
				  domemoff[upbufdx] = cubicinterp(frac, d0, d1, d2, d3);
			  }
			  iwrphase++;
			  upbufdx++;
			  );
	} else {

		float next_dsamp = CalcDelay_OS(unit, delaytime, 8);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp) * 0.125;

		LOOP1(inNumSamples * 8,
			  dsamp += dsamp_slope;
			  long idsamp = (long)dsamp;
			  float frac = dsamp - idsamp;
			  long irdphase1 = iwrphase - idsamp;
			  long irdphase2 = irdphase1 - 1;
			  long irdphase3 = irdphase1 - 2;
			  long irdphase0 = irdphase1 + 1;

			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  if (irdphase0 < 0) {
				  domemoff[upbufdx] = 0.f;
			  } else {
				  if (irdphase1 < 0) {
					  d1 = d2 = d3 = 0.f;
					  d0 = dlybuf[irdphase0 & mask];
				  } else if (irdphase2 < 0) {
					  d1 = d2 = d3 = 0.f;
					  d0 = dlybuf[irdphase0 & mask];
					  d1 = dlybuf[irdphase1 & mask];
				  } else if (irdphase3 < 0) {
					  d3 = 0.f;
					  d0 = dlybuf[irdphase0 & mask];
					  d1 = dlybuf[irdphase1 & mask];
					  d2 = dlybuf[irdphase2 & mask];
				  } else {
					  d0 = dlybuf[irdphase0 & mask];
					  d1 = dlybuf[irdphase1 & mask];
					  d2 = dlybuf[irdphase2 & mask];
					  d3 = dlybuf[irdphase3 & mask];
				  }
				  domemoff[upbufdx] = cubicinterp(frac, d0, d1, d2, d3);
			  }
			  iwrphase++;
			  upbufdx++;
			  );
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
	}

	unit->m_iwrphase = iwrphase;
	unit->m_numoutput += inNumSamples * 8;							/* x8 added */

	DOWNSAMPLE8;

	if (unit->m_numoutput >= unit->m_idelaylen) {
		SETCALC(DelayC_OS8_next);
	}
	//	else {
	//		Print("Keep going wiht the z\n");
	//	}
}


void DelayC_OS4_next(DelayC_OS4 *unit, int inNumSamples)
{

	float *out = ZOUT(0);
	float *in = ZIN(0);
	float delaytime = ZIN0(2);

	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	long mask = unit->m_mask;
	int upbufdx = 0;

	UPSAMPLE4;

	if (delaytime == unit->m_delaytime) {	// if the delay time hasn't changed
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		LOOP1(inNumSamples * 4,
			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  long irdphase1 = iwrphase - idsamp;
			  long irdphase2 = irdphase1 - 1;
			  long irdphase3 = irdphase1 - 2;
			  long irdphase0 = irdphase1 + 1;
			  float d0 = dlybuf[irdphase0 & mask];
			  float d1 = dlybuf[irdphase1 & mask];
			  float d2 = dlybuf[irdphase2 & mask];
			  float d3 = dlybuf[irdphase3 & mask];
			  domemoff[upbufdx] = cubicinterp(frac, d0, d1, d2, d3);
			  iwrphase++;
			  upbufdx++;
			  );
	} else {

		float next_dsamp = CalcDelay_OS(unit, delaytime, 4);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp) * 0.25;
		LOOP1(inNumSamples * 4,
			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  dsamp += dsamp_slope;
			  long idsamp = (long)dsamp;
			  float frac = dsamp - idsamp;
			  long irdphase1 = iwrphase - idsamp;
			  long irdphase2 = irdphase1 - 1;
			  long irdphase3 = irdphase1 - 2;
			  long irdphase0 = irdphase1 + 1;
			  float d0 = dlybuf[irdphase0 & mask];
			  float d1 = dlybuf[irdphase1 & mask];
			  float d2 = dlybuf[irdphase2 & mask];
			  float d3 = dlybuf[irdphase3 & mask];
			  domemoff[upbufdx] = cubicinterp(frac, d0, d1, d2, d3);
			  iwrphase++;
			  upbufdx++;
			  );
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
	}

	DOWNSAMPLE4;

	unit->m_iwrphase = iwrphase;

}


void DelayC_OS4_next_z(DelayC_OS4 *unit, int inNumSamples)
{

	float *out = ZOUT(0);
	float *in = ZIN(0);
	float delaytime = ZIN0(2);

	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	long mask = unit->m_mask;
	float d0, d1, d2, d3;
	int upbufdx = 0;

	UPSAMPLE4;

	if (delaytime == unit->m_delaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		LOOP1(inNumSamples * 4,
			  long irdphase1 = iwrphase - idsamp;
			  long irdphase2 = irdphase1 - 1;
			  long irdphase3 = irdphase1 - 2;
			  long irdphase0 = irdphase1 + 1;

			  dlybuf[iwrphase & mask] = domemoff[upbufdx];

			  if (irdphase0 < 0) {
				  domemoff[upbufdx] = 0.f;
			  } else {
				  if (irdphase1 < 0) {
					  d1 = d2 = d3 = 0.f;
					  d0 = dlybuf[irdphase0 & mask];
				  } else if (irdphase2 < 0) {
					  d1 = d2 = d3 = 0.f;
					  d0 = dlybuf[irdphase0 & mask];
					  d1 = dlybuf[irdphase1 & mask];
				  } else if (irdphase3 < 0) {
					  d3 = 0.f;
					  d0 = dlybuf[irdphase0 & mask];
					  d1 = dlybuf[irdphase1 & mask];
					  d2 = dlybuf[irdphase2 & mask];
				  } else {
					  d0 = dlybuf[irdphase0 & mask];
					  d1 = dlybuf[irdphase1 & mask];
					  d2 = dlybuf[irdphase2 & mask];
					  d3 = dlybuf[irdphase3 & mask];
				  }
				  domemoff[upbufdx] = cubicinterp(frac, d0, d1, d2, d3);
			  }
			  iwrphase++;
			  upbufdx++;
			  );
	} else {

		float next_dsamp = CalcDelay_OS(unit, delaytime, 4);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp) * 0.25;

		LOOP1(inNumSamples * 4,
			  dsamp += dsamp_slope;
			  long idsamp = (long)dsamp;
			  float frac = dsamp - idsamp;
			  long irdphase1 = iwrphase - idsamp;
			  long irdphase2 = irdphase1 - 1;
			  long irdphase3 = irdphase1 - 2;
			  long irdphase0 = irdphase1 + 1;

			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  if (irdphase0 < 0) {
				  domemoff[upbufdx] = 0.f;
			  } else {
				  if (irdphase1 < 0) {
					  d1 = d2 = d3 = 0.f;
					  d0 = dlybuf[irdphase0 & mask];
				  } else if (irdphase2 < 0) {
					  d1 = d2 = d3 = 0.f;
					  d0 = dlybuf[irdphase0 & mask];
					  d1 = dlybuf[irdphase1 & mask];
				  } else if (irdphase3 < 0) {
					  d3 = 0.f;
					  d0 = dlybuf[irdphase0 & mask];
					  d1 = dlybuf[irdphase1 & mask];
					  d2 = dlybuf[irdphase2 & mask];
				  } else {
					  d0 = dlybuf[irdphase0 & mask];
					  d1 = dlybuf[irdphase1 & mask];
					  d2 = dlybuf[irdphase2 & mask];
					  d3 = dlybuf[irdphase3 & mask];
				  }
				  domemoff[upbufdx] = cubicinterp(frac, d0, d1, d2, d3);
			  }
			  upbufdx++;
			  iwrphase++;
			  );
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
	}

	unit->m_iwrphase = iwrphase;
	unit->m_numoutput += inNumSamples * 4;

	DOWNSAMPLE4;

	if (unit->m_numoutput >= unit->m_idelaylen) {
		SETCALC(DelayC_OS4_next);
	}
}


void DelayC_OS3_next(DelayC_OS3 *unit, int inNumSamples)
{

	float *out = ZOUT(0);
	float *in = ZIN(0);
	float delaytime = ZIN0(2);

	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	long mask = unit->m_mask;
	int upbufdx = 0;

	UPSAMPLE3;

	if (delaytime == unit->m_delaytime) {	// if the delay time hasn't changed
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		LOOP1(inNumSamples * 3,
			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  long irdphase1 = iwrphase - idsamp;
			  long irdphase2 = irdphase1 - 1;
			  long irdphase3 = irdphase1 - 2;
			  long irdphase0 = irdphase1 + 1;
			  float d0 = dlybuf[irdphase0 & mask];
			  float d1 = dlybuf[irdphase1 & mask];
			  float d2 = dlybuf[irdphase2 & mask];
			  float d3 = dlybuf[irdphase3 & mask];
			  domemoff[upbufdx] = cubicinterp(frac, d0, d1, d2, d3);
			  iwrphase++;
			  upbufdx++;
			  );
	} else {

		float next_dsamp = CalcDelay_OS(unit, delaytime, 3);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp) * 0.3333333333;
		LOOP1(inNumSamples * 3,
			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  dsamp += dsamp_slope;
			  long idsamp = (long)dsamp;
			  float frac = dsamp - idsamp;
			  long irdphase1 = iwrphase - idsamp;
			  long irdphase2 = irdphase1 - 1;
			  long irdphase3 = irdphase1 - 2;
			  long irdphase0 = irdphase1 + 1;
			  float d0 = dlybuf[irdphase0 & mask];
			  float d1 = dlybuf[irdphase1 & mask];
			  float d2 = dlybuf[irdphase2 & mask];
			  float d3 = dlybuf[irdphase3 & mask];
			  domemoff[upbufdx] = cubicinterp(frac, d0, d1, d2, d3);
			  iwrphase++;
			  upbufdx++;
			  );
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
	}

	DOWNSAMPLE3;

	unit->m_iwrphase = iwrphase;

}


void DelayC_OS3_next_z(DelayC_OS3 *unit, int inNumSamples)
{

	float *out = ZOUT(0);
	float *in = ZIN(0);
	float delaytime = ZIN0(2);

	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	long mask = unit->m_mask;
	float d0, d1, d2, d3;
	int upbufdx = 0;

	UPSAMPLE3;

	if (delaytime == unit->m_delaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		LOOP1(inNumSamples * 3,
			  long irdphase1 = iwrphase - idsamp;
			  long irdphase2 = irdphase1 - 1;
			  long irdphase3 = irdphase1 - 2;
			  long irdphase0 = irdphase1 + 1;

			  dlybuf[iwrphase & mask] = domemoff[upbufdx];

			  if (irdphase0 < 0) {
				  domemoff[upbufdx] = 0.f;
			  } else {
				  if (irdphase1 < 0) {
					  d1 = d2 = d3 = 0.f;
					  d0 = dlybuf[irdphase0 & mask];
				  } else if (irdphase2 < 0) {
					  d1 = d2 = d3 = 0.f;
					  d0 = dlybuf[irdphase0 & mask];
					  d1 = dlybuf[irdphase1 & mask];
				  } else if (irdphase3 < 0) {
					  d3 = 0.f;
					  d0 = dlybuf[irdphase0 & mask];
					  d1 = dlybuf[irdphase1 & mask];
					  d2 = dlybuf[irdphase2 & mask];
				  } else {
					  d0 = dlybuf[irdphase0 & mask];
					  d1 = dlybuf[irdphase1 & mask];
					  d2 = dlybuf[irdphase2 & mask];
					  d3 = dlybuf[irdphase3 & mask];
				  }
				  domemoff[upbufdx] = cubicinterp(frac, d0, d1, d2, d3);
			  }
			  iwrphase++;
			  upbufdx++;
			  );
	} else {

		float next_dsamp = CalcDelay_OS(unit, delaytime, 3);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp) * 0.33333333;

		LOOP1(inNumSamples * 3,
			  dsamp += dsamp_slope;
			  long idsamp = (long)dsamp;
			  float frac = dsamp - idsamp;
			  long irdphase1 = iwrphase - idsamp;
			  long irdphase2 = irdphase1 - 1;
			  long irdphase3 = irdphase1 - 2;
			  long irdphase0 = irdphase1 + 1;

			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  if (irdphase0 < 0) {
				  domemoff[upbufdx] = 0.f;
			  } else {
				  if (irdphase1 < 0) {
					  d1 = d2 = d3 = 0.f;
					  d0 = dlybuf[irdphase0 & mask];
				  } else if (irdphase2 < 0) {
					  d1 = d2 = d3 = 0.f;
					  d0 = dlybuf[irdphase0 & mask];
					  d1 = dlybuf[irdphase1 & mask];
				  } else if (irdphase3 < 0) {
					  d3 = 0.f;
					  d0 = dlybuf[irdphase0 & mask];
					  d1 = dlybuf[irdphase1 & mask];
					  d2 = dlybuf[irdphase2 & mask];
				  } else {
					  d0 = dlybuf[irdphase0 & mask];
					  d1 = dlybuf[irdphase1 & mask];
					  d2 = dlybuf[irdphase2 & mask];
					  d3 = dlybuf[irdphase3 & mask];
				  }
				  domemoff[upbufdx] = cubicinterp(frac, d0, d1, d2, d3);
			  }
			  upbufdx++;
			  iwrphase++;
			  );
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
	}

	unit->m_iwrphase = iwrphase;
	unit->m_numoutput += inNumSamples * 3;

	DOWNSAMPLE3;

	if (unit->m_numoutput >= unit->m_idelaylen) {
		SETCALC(DelayC_OS3_next);
	}
}



/* DelayL_OS *//////////////////////////////////////////

void DelayL_OS8_next(DelayL_OS8 *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float *in = ZIN(0);
	float delaytime = ZIN0(2);

	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	long mask = unit->m_mask;
	int upbufdx = 0;

	UPSAMPLE8;

	if (delaytime == unit->m_delaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		LOOP1(inNumSamples * 8,
			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  //dlybuf[iwrphase & mask] = ZXP(in);
			  long irdphase = iwrphase - idsamp;
			  long irdphaseb = irdphase - 1;
			  float d1 = dlybuf[irdphase & mask];
			  float d2 = dlybuf[irdphaseb & mask];
			  domemoff[upbufdx] = lininterp(frac, d1, d2);
			  //ZXP(out) = lininterp(frac, d1, d2);
			  iwrphase++;
			  upbufdx++;
			  );
	} else {

		float next_dsamp = CalcDelay_OS(unit, delaytime, 8);
		//float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp) * 0.125;

		LOOP1(inNumSamples * 8,
			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  //dlybuf[iwrphase & mask] = ZXP(in);
			  dsamp += dsamp_slope;
			  long idsamp = (long)dsamp;
			  float frac = dsamp - idsamp;
			  long irdphase = iwrphase - idsamp;
			  long irdphaseb = irdphase - 1;
			  float d1 = dlybuf[irdphase & mask];
			  float d2 = dlybuf[irdphaseb & mask];
			  domemoff[upbufdx] = lininterp(frac, d1, d2);
			  //ZXP(out) = lininterp(frac, d1, d2);
			  iwrphase++;
			  upbufdx++;
			  );
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
	}

	DOWNSAMPLE8;

	unit->m_iwrphase = iwrphase;

}


void DelayL_OS8_next_z(DelayL_OS8 *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float *in = ZIN(0);
	float delaytime = ZIN0(2);

	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	long mask = unit->m_mask;
	int upbufdx = 0;

	UPSAMPLE8;

	if (delaytime == unit->m_delaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		LOOP1(inNumSamples * 8,
			  long irdphase = iwrphase - idsamp;
			  long irdphaseb = irdphase - 1;

			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  //dlybuf[iwrphase & mask] = ZXP(in);
			  if (irdphase < 0) {
				  domemoff[upbufdx] = 0.f;
				  //ZXP(out) = 0.f;
			  } else if (irdphaseb < 0) {
				  float d1 = dlybuf[irdphase & mask];
				  //postbuf("A %d d1 %g fr %g v %g in %g fb %g\n", irdphase, d1, frac, value, zin, feedbk);
				  domemoff[upbufdx] = d1 - frac * d1;
				  //ZXP(out) = d1 - frac * d1;
			  } else {
				  float d1 = dlybuf[irdphase & mask];
				  float d2 = dlybuf[irdphaseb & mask];
				  //postbuf("B %d d1 %g d2 %g fr %g v %g in %g fb %g\n", irdphase, d1, d2, frac, value, zin, feedbk);
				  domemoff[upbufdx] = lininterp(frac, d1, d2);
				  //ZXP(out) = lininterp(frac, d1, d2);
			  }
			  iwrphase++;
			  upbufdx++;
			  );
	} else {

		float next_dsamp = CalcDelay_OS(unit, delaytime, 8);
		//float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp) * 0.125;

		LOOP1(inNumSamples * 8,
			  dsamp += dsamp_slope;
			  long idsamp = (long)dsamp;
			  float frac = dsamp - idsamp;
			  long irdphase = iwrphase - idsamp;
			  long irdphaseb = irdphase - 1;

			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  //dlybuf[iwrphase & mask] = ZXP(in);
			  if (irdphase < 0) {
				  domemoff[upbufdx] = 0.f;
				  //ZXP(out) = 0.f;
			  } else if (irdphaseb < 0) {
				  float d1 = dlybuf[irdphase & mask];
				  domemoff[upbufdx] = d1 - frac * d1;
				  //ZXP(out) = d1 - frac * d1;
			  } else {
				  float d1 = dlybuf[irdphase & mask];
				  float d2 = dlybuf[irdphaseb & mask];
				  domemoff[upbufdx] = lininterp(frac, d1, d2);
				  //ZXP(out) = lininterp(frac, d1, d2);
			  }
			  iwrphase++;
			  upbufdx++;
			  );
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
	}

	unit->m_iwrphase = iwrphase;
	unit->m_numoutput += inNumSamples * 8;

	DOWNSAMPLE8;

	/* DOES THIS CHECK WORK?
		unit->m_numoutput is incremented one line above with inNumSamples
		while unit->m_idelaylen is set to inNumSamples * 8 in DelayUnit_OS_AllocDelayLine
	*/
	if (unit->m_numoutput >= unit->m_idelaylen) {
		SETCALC(DelayL_OS8_next);
	}
}


void DelayL_OS4_next(DelayL_OS4 *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float *in = ZIN(0);
	float delaytime = ZIN0(2);

	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	long mask = unit->m_mask;
	int upbufdx = 0;

	UPSAMPLE4;

	if (delaytime == unit->m_delaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		LOOP1(inNumSamples * 4,
			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  //dlybuf[iwrphase & mask] = ZXP(in);
			  long irdphase = iwrphase - idsamp;
			  long irdphaseb = irdphase - 1;
			  float d1 = dlybuf[irdphase & mask];
			  float d2 = dlybuf[irdphaseb & mask];
			  domemoff[upbufdx] = lininterp(frac, d1, d2);
			  //ZXP(out) = lininterp(frac, d1, d2);
			  iwrphase++;
			  upbufdx++;
			  );
	} else {

		float next_dsamp = CalcDelay_OS(unit, delaytime, 4);
		//float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp) * 0.25;

		LOOP1(inNumSamples * 4,
			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  //dlybuf[iwrphase & mask] = ZXP(in);
			  dsamp += dsamp_slope;
			  long idsamp = (long)dsamp;
			  float frac = dsamp - idsamp;
			  long irdphase = iwrphase - idsamp;
			  long irdphaseb = irdphase - 1;
			  float d1 = dlybuf[irdphase & mask];
			  float d2 = dlybuf[irdphaseb & mask];
			  domemoff[upbufdx] = lininterp(frac, d1, d2);
			  //ZXP(out) = lininterp(frac, d1, d2);
			  iwrphase++;
			  upbufdx++;
			  );
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
	}

	DOWNSAMPLE4;

	unit->m_iwrphase = iwrphase;

}


void DelayL_OS4_next_z(DelayL_OS4 *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float *in = ZIN(0);
	float delaytime = ZIN0(2);

	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	long mask = unit->m_mask;
	int upbufdx = 0;

	UPSAMPLE4;

	if (delaytime == unit->m_delaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		LOOP1(inNumSamples * 4,
			  long irdphase = iwrphase - idsamp;
			  long irdphaseb = irdphase - 1;

			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  //dlybuf[iwrphase & mask] = ZXP(in);
			  if (irdphase < 0) {
				  domemoff[upbufdx] = 0.f;
				  //ZXP(out) = 0.f;
			  } else if (irdphaseb < 0) {
				  float d1 = dlybuf[irdphase & mask];
				  //postbuf("A %d d1 %g fr %g v %g in %g fb %g\n", irdphase, d1, frac, value, zin, feedbk);
				  domemoff[upbufdx] = d1 - frac * d1;
				  //ZXP(out) = d1 - frac * d1;
			  } else {
				  float d1 = dlybuf[irdphase & mask];
				  float d2 = dlybuf[irdphaseb & mask];
				  //postbuf("B %d d1 %g d2 %g fr %g v %g in %g fb %g\n", irdphase, d1, d2, frac, value, zin, feedbk);
				  domemoff[upbufdx] = lininterp(frac, d1, d2);
				  //ZXP(out) = lininterp(frac, d1, d2);
			  }
			  iwrphase++;
			  upbufdx++;
			  );
	} else {

		float next_dsamp = CalcDelay_OS(unit, delaytime, 4);
		//float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp) * 0.25;

		LOOP1(inNumSamples * 4,
			  dsamp += dsamp_slope;
			  long idsamp = (long)dsamp;
			  float frac = dsamp - idsamp;
			  long irdphase = iwrphase - idsamp;
			  long irdphaseb = irdphase - 1;

			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  //dlybuf[iwrphase & mask] = ZXP(in);
			  if (irdphase < 0) {
				  domemoff[upbufdx] = 0.f;
				  //ZXP(out) = 0.f;
			  } else if (irdphaseb < 0) {
				  float d1 = dlybuf[irdphase & mask];
				  domemoff[upbufdx] = d1 - frac * d1;
				  //ZXP(out) = d1 - frac * d1;
			  } else {
				  float d1 = dlybuf[irdphase & mask];
				  float d2 = dlybuf[irdphaseb & mask];
				  domemoff[upbufdx] = lininterp(frac, d1, d2);
				  //ZXP(out) = lininterp(frac, d1, d2);
			  }
			  iwrphase++;
			  upbufdx++;
			  );
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
	}

	unit->m_iwrphase = iwrphase;
	unit->m_numoutput += inNumSamples * 4;

	DOWNSAMPLE4;

	/* DOES THIS CHECK WORK?
	 unit->m_numoutput is incremented one line above with inNumSamples
	 while unit->m_idelaylen is set to inNumSamples * 8 in DelayUnit_OS_AllocDelayLine
	 */
	if (unit->m_numoutput >= unit->m_idelaylen) {
		SETCALC(DelayL_OS4_next);
	}
}

void DelayL_OS3_next(DelayL_OS3 *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float *in = ZIN(0);
	float delaytime = ZIN0(2);

	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	long mask = unit->m_mask;
	int upbufdx = 0;

	UPSAMPLE3;

	if (delaytime == unit->m_delaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		LOOP1(inNumSamples * 3,
			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  //dlybuf[iwrphase & mask] = ZXP(in);
			  long irdphase = iwrphase - idsamp;
			  long irdphaseb = irdphase - 1;
			  float d1 = dlybuf[irdphase & mask];
			  float d2 = dlybuf[irdphaseb & mask];
			  domemoff[upbufdx] = lininterp(frac, d1, d2);
			  //ZXP(out) = lininterp(frac, d1, d2);
			  iwrphase++;
			  upbufdx++;
			  );
	} else {

		float next_dsamp = CalcDelay_OS(unit, delaytime, 3);
		//float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp) * 0.3333333333333;

		LOOP1(inNumSamples * 3,
			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  //dlybuf[iwrphase & mask] = ZXP(in);
			  dsamp += dsamp_slope;
			  long idsamp = (long)dsamp;
			  float frac = dsamp - idsamp;
			  long irdphase = iwrphase - idsamp;
			  long irdphaseb = irdphase - 1;
			  float d1 = dlybuf[irdphase & mask];
			  float d2 = dlybuf[irdphaseb & mask];
			  domemoff[upbufdx] = lininterp(frac, d1, d2);
			  //ZXP(out) = lininterp(frac, d1, d2);
			  iwrphase++;
			  upbufdx++;
			  );
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
	}

	DOWNSAMPLE3;

	unit->m_iwrphase = iwrphase;

}


void DelayL_OS3_next_z(DelayL_OS3 *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float *in = ZIN(0);
	float delaytime = ZIN0(2);

	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	long mask = unit->m_mask;
	int upbufdx = 0;

	UPSAMPLE3;

	if (delaytime == unit->m_delaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		LOOP1(inNumSamples * 3,
			  long irdphase = iwrphase - idsamp;
			  long irdphaseb = irdphase - 1;

			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  //dlybuf[iwrphase & mask] = ZXP(in);
			  if (irdphase < 0) {
				  domemoff[upbufdx] = 0.f;
				  //ZXP(out) = 0.f;
			  } else if (irdphaseb < 0) {
				  float d1 = dlybuf[irdphase & mask];
				  //postbuf("A %d d1 %g fr %g v %g in %g fb %g\n", irdphase, d1, frac, value, zin, feedbk);
				  domemoff[upbufdx] = d1 - frac * d1;
				  //ZXP(out) = d1 - frac * d1;
			  } else {
				  float d1 = dlybuf[irdphase & mask];
				  float d2 = dlybuf[irdphaseb & mask];
				  //postbuf("B %d d1 %g d2 %g fr %g v %g in %g fb %g\n", irdphase, d1, d2, frac, value, zin, feedbk);
				  domemoff[upbufdx] = lininterp(frac, d1, d2);
				  //ZXP(out) = lininterp(frac, d1, d2);
			  }
			  iwrphase++;
			  upbufdx++;
			  );
	} else {

		float next_dsamp = CalcDelay_OS(unit, delaytime, 3);
		//float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp) * 0.3333333333333;

		LOOP1(inNumSamples * 3,
			  dsamp += dsamp_slope;
			  long idsamp = (long)dsamp;
			  float frac = dsamp - idsamp;
			  long irdphase = iwrphase - idsamp;
			  long irdphaseb = irdphase - 1;

			  dlybuf[iwrphase & mask] = domemoff[upbufdx];
			  //dlybuf[iwrphase & mask] = ZXP(in);
			  if (irdphase < 0) {
				  domemoff[upbufdx] = 0.f;
				  //ZXP(out) = 0.f;
			  } else if (irdphaseb < 0) {
				  float d1 = dlybuf[irdphase & mask];
				  domemoff[upbufdx] = d1 - frac * d1;
				  //ZXP(out) = d1 - frac * d1;
			  } else {
				  float d1 = dlybuf[irdphase & mask];
				  float d2 = dlybuf[irdphaseb & mask];
				  domemoff[upbufdx] = lininterp(frac, d1, d2);
				  //ZXP(out) = lininterp(frac, d1, d2);
			  }
			  iwrphase++;
			  upbufdx++;
			  );
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
	}

	unit->m_iwrphase = iwrphase;
	unit->m_numoutput += inNumSamples * 3;

	DOWNSAMPLE3;

	/* DOES THIS CHECK WORK?
	 unit->m_numoutput is incremented one line above with inNumSamples
	 while unit->m_idelaylen is set to inNumSamples * 8 in DelayUnit_OS_AllocDelayLine
	 */
	if (unit->m_numoutput >= unit->m_idelaylen) {
		SETCALC(DelayL_OS3_next);
	}
}




////////////////////////
////////////////////////
/*  LOAD  */
////////////////////////
////////////////////////
// >= sc 3.9
PluginLoad(Delay_OS) {

    ft = inTable; // store pointer to InterfaceTable

		DefineDtorCantAliasUnit(DelayC_OS8);
		DefineDtorCantAliasUnit(DelayC_OS4);
		DefineDtorCantAliasUnit(DelayC_OS3);
		DefineDtorCantAliasUnit(DelayL_OS8);
		DefineDtorCantAliasUnit(DelayL_OS4);
		DefineDtorCantAliasUnit(DelayL_OS3);
}

// // <= sc 3.8
// void load(InterfaceTable *inTable)
// {
// 	ft = inTable;
//
// //#define DefineMyDelayUnit(name) \
// //	(*ft->fDefineUnit)(#name, sizeof(name), (UnitCtorFunc)&name##_Ctor, \
// //		(UnitDtorFunc)&DelayUnit_OS8_Dtor, 0);
//
// 	DefineDtorCantAliasUnit(DelayC_OS8);
// 	DefineDtorCantAliasUnit(DelayC_OS4);
// 	DefineDtorCantAliasUnit(DelayC_OS3);
// 	DefineDtorCantAliasUnit(DelayL_OS8);
// 	DefineDtorCantAliasUnit(DelayL_OS4);
// 	DefineDtorCantAliasUnit(DelayL_OS3);
// }
