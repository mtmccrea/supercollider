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

#include "SC_PlugIn.h"
#include "function_attributes.h"
#include <limits>
#include <string.h>

static InterfaceTable *ft;

struct BufUnit : public Unit
{
	SndBuf *m_buf;
	float m_fbufnum;
};

struct TableLookup : public BufUnit
{
	double m_cpstoinc, m_radtoinc;
	int32 mTableSize;
	int32 m_lomask;
};

struct DegreeToKey : public BufUnit
{
	int32 mPrevIndex;
	float mPrevKey;
	int32 mOctave;
};

struct Select : public Unit
{
};

struct TWindex : public Unit
{
	int32 m_prevIndex;
	float m_trig;
};

struct Index : public BufUnit
{
};

struct IndexL : public BufUnit
{
};

struct WrapIndex : public BufUnit
{
};

struct FoldIndex : public BufUnit
{
};

struct IndexInBetween : public BufUnit
{
};

struct DetectIndex : public BufUnit
{
	float mPrev;
	float mPrevIn;
};

struct Shaper : public BufUnit
{
	float mOffset;
	float mPrevIn;
};

struct FSinOsc : public Unit
{
	double m_b1, m_y1, m_y2, m_freq;
};

struct PSinGrain : public Unit
{
	double m_b1, m_y1, m_y2;
	double m_level, m_slope, m_curve;
	int32 mCounter;
};

struct Osc : public TableLookup
{
	int32 m_phase;
	float m_phasein;
};

struct SinOsc : public TableLookup
{
	int32 m_phase;
	float m_phasein;
};

struct SinOscFB : public TableLookup
{
	int32 m_phase;
	float m_prevout, m_feedback;
};

struct OscN : public TableLookup
{
	int32 m_phase;
	float m_phasein;
};

struct COsc : public TableLookup
{
	int32 m_phase1, m_phase2;
};

struct VOsc : public Unit
{
	double m_cpstoinc, m_radtoinc;
	int32 mTableSize;
	int32 m_lomask;
	int32 m_phase, m_phaseoffset;
	float m_phasein, m_bufpos;
};

struct VOsc3 : public Unit
{
	double m_cpstoinc;
	int32 mTableSize;
	int32 m_lomask;
	int32 m_phase1, m_phase2, m_phase3;
	float m_bufpos;
};

struct Formant : public Unit
{
	int32 m_phase1, m_phase2, m_phase3;
	double m_cpstoinc;
};



struct Blip : public Unit
{
	int32 m_phase, m_numharm, m_N;
	float m_freqin, m_scale;
	double m_cpstoinc;
};

struct Saw : public Unit
{
	int32 m_phase, m_N;
	float m_freqin, m_scale, m_y1;
	double m_cpstoinc;
};

struct Pulse : public Unit
{
	int32 m_phase, m_phaseoff, m_N;
	float m_freqin, m_scale, m_y1;
	double m_cpstoinc;
};

struct Klang : public Unit
{
	float *m_coefs;
	int32 m_numpartials;
};

struct Klank : public Unit
{
	float *m_coefs;
	float *m_buf;
	float m_x1, m_x2;
	int32 m_numpartials;
};



#define xlomask8  0x000003FC
#define xlomask9  0x000007FC
#define xlomask10 0x00000FFC
#define xlomask11 0x00001FFC
#define xlomask12 0x00003FFC
#define xlomask13 0x00007FFC

#define xlomask8i  0x000007F8
#define xlomask9i  0x00000FF8
#define xlomask10i 0x00001FF8
#define xlomask11i 0x00003FF8
#define xlomask12i 0x00007FF8
#define xlomask13i 0x0000FFF8

#define onecyc13 0x20000000



//////////////////////////////////////////////////////////////////////////////////////////////////

extern "C"
{

void DegreeToKey_Ctor(DegreeToKey *unit);
void DegreeToKey_next_1(DegreeToKey *unit, int inNumSamples);
void DegreeToKey_next_k(DegreeToKey *unit, int inNumSamples);
void DegreeToKey_next_a(DegreeToKey *unit, int inNumSamples);

void Select_Ctor(Select *unit);
void Select_next_1(Select *unit, int inNumSamples);
void Select_next_k(Select *unit, int inNumSamples);
void Select_next_a(Select *unit, int inNumSamples);

void TWindex_Ctor(TWindex *unit);
void TWindex_next_k(TWindex *unit, int inNumSamples);
void TWindex_next_ak(TWindex *unit, int inNumSamples);

void Index_Ctor(Index *unit);
void Index_next_1(Index *unit, int inNumSamples);
void Index_next_k(Index *unit, int inNumSamples);
void Index_next_a(Index *unit, int inNumSamples);

void IndexL_Ctor(IndexL *unit);
void IndexL_next_1(IndexL *unit, int inNumSamples);
void IndexL_next_k(IndexL *unit, int inNumSamples);
void IndexL_next_a(IndexL *unit, int inNumSamples);

void FoldIndex_Ctor(FoldIndex *unit);
void FoldIndex_next_1(FoldIndex *unit, int inNumSamples);
void FoldIndex_next_k(FoldIndex *unit, int inNumSamples);
void FoldIndex_next_a(FoldIndex *unit, int inNumSamples);

void WrapIndex_Ctor(WrapIndex *unit);
void WrapIndex_next_1(WrapIndex *unit, int inNumSamples);
void WrapIndex_next_k(WrapIndex *unit, int inNumSamples);
void WrapIndex_next_a(WrapIndex *unit, int inNumSamples);

void Shaper_Ctor(Shaper *unit);
void Shaper_next_1(Shaper *unit, int inNumSamples);
void Shaper_next_k(Shaper *unit, int inNumSamples);
void Shaper_next_a(Shaper *unit, int inNumSamples);


void DetectIndex_Ctor(DetectIndex *unit);
void DetectIndex_next_1(DetectIndex *unit, int inNumSamples);
void DetectIndex_next_k(DetectIndex *unit, int inNumSamples);
void DetectIndex_next_a(DetectIndex *unit, int inNumSamples);

void IndexInBetween_Ctor(IndexInBetween *unit);
void IndexInBetween_next_1(IndexInBetween *unit, int inNumSamples);
void IndexInBetween_next_k(IndexInBetween *unit, int inNumSamples);
void IndexInBetween_next_a(IndexInBetween *unit, int inNumSamples);


void FSinOsc_Ctor(FSinOsc *unit);
void FSinOsc_next(FSinOsc *unit, int inNumSamples);
void FSinOsc_next_i(FSinOsc *unit, int inNumSamples);

void PSinGrain_Ctor(PSinGrain *unit);
void PSinGrain_next(PSinGrain *unit, int inNumSamples);

void SinOsc_Ctor(SinOsc *unit);
void SinOsc_next_ikk(SinOsc *unit, int inNumSamples);
void SinOsc_next_ika(SinOsc *unit, int inNumSamples);
void SinOsc_next_iak(SinOsc *unit, int inNumSamples);
void SinOsc_next_iaa(SinOsc *unit, int inNumSamples);

void Osc_Ctor(Osc *unit);
void Osc_next_ikk(Osc *unit, int inNumSamples);
void Osc_next_ika(Osc *unit, int inNumSamples);
void Osc_next_iak(Osc *unit, int inNumSamples);
void Osc_next_iaa(Osc *unit, int inNumSamples);

void OscN_Ctor(OscN *unit);
void OscN_next_nkk(OscN *unit, int inNumSamples);
void OscN_next_nka(OscN *unit, int inNumSamples);
void OscN_next_nak(OscN *unit, int inNumSamples);
void OscN_next_naa(OscN *unit, int inNumSamples);

void COsc_Ctor(COsc *unit);
void COsc_next(COsc *unit, int inNumSamples);

void VOsc_Ctor(VOsc *unit);
void VOsc_next_ikk(VOsc *unit, int inNumSamples);
void VOsc_next_ika(VOsc *unit, int inNumSamples);

void VOsc3_Ctor(VOsc3 *unit);
void VOsc3_next_ik(VOsc3 *unit, int inNumSamples);

void Formant_Ctor(Formant *unit);
void Formant_next(Formant *unit, int inNumSamples);

void Blip_Ctor(Blip *unit);
void Blip_next(Blip *unit, int inNumSamples);

void Saw_Ctor(Saw *unit);
void Saw_next(Saw *unit, int inNumSamples);

void Pulse_Ctor(Pulse *unit);
void Pulse_next(Pulse *unit, int inNumSamples);

void Klang_Dtor(Klang *unit);
void Klang_Ctor(Klang *unit);
void Klang_next(Klang *unit, int inNumSamples);

void Klank_Dtor(Klank *unit);
void Klank_Ctor(Klank *unit);
void Klank_next(Klank *unit, int inNumSamples);

}

//////////////////////////////////////////////////////////////////////////////////////////////////

force_inline bool UnitGetTable(BufUnit * unit, int inNumSamples, const SndBuf * & buf,
							   const float * & bufData, int & tableSize)
{
	float fbufnum = ZIN0(0);
	if (fbufnum != unit->m_fbufnum) {
		uint32 bufnum = (uint32)fbufnum;
		World *world = unit->mWorld;
		if (bufnum >= world->mNumSndBufs) {
			uint32 localBufNum = bufnum - world->mNumSndBufs;
			Graph *parent = unit->mParent;
			if(localBufNum <= parent->localBufNum)
				unit->m_buf = parent->mLocalSndBufs + localBufNum;
			else {
				bufnum = 0;
				unit->m_buf = world->mSndBufs + bufnum;
			}
		} else
			unit->m_buf = world->mSndBufs + bufnum;

		unit->m_fbufnum = fbufnum;
	}
	buf = unit->m_buf;
	if(!buf) {
		ClearUnitOutputs(unit, inNumSamples);
		return false;
	}

	bufData = buf->data;
	if (!bufData) {
		ClearUnitOutputs(unit, inNumSamples);
		return false;
	}
	tableSize = buf->samples;
	return true;
}

#define GET_TABLE \
	const SndBuf * buf; \
	const float * bufData; \
	int tableSize; \
	do { \
		bool tableValid = UnitGetTable(unit, inNumSamples, buf, bufData, tableSize); \
		if (!tableValid) return; \
	} while (0);

static inline bool verify_wavetable(Unit * unit, const char * name, int tableSize, int inNumSamples)
{
	// phase computation is not precise for large wavetables.
	if (tableSize > 131072) {
		if (unit->mWorld->mVerbosity >= -1)
			Print("Warning: wave table too big (%s)\n", name);
		ClearUnitOutputs(unit, inNumSamples);
		return false;
	}

	if (!ISPOWEROFTWO(tableSize)) {
		if (unit->mWorld->mVerbosity >= -1)
			Print("Warning: size of wavetable not a power of two (%s)\n", name);
		ClearUnitOutputs(unit, inNumSamples);
		return false;
	}

	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
void TableLookup_SetTable(TableLookup* unit, int32 inSize, float* inTable)
{
	unit->mTable0 = inTable;
	unit->mTable1 = inTable + 1;
	unit->mTableSize = inSize;
	unit->mMaxIndex = unit->mTableSize - 1;
	unit->mFMaxIndex = unit->mMaxIndex;
	unit->m_radtoinc = unit->mTableSize * (rtwopi * 65536.);
	unit->m_cpstoinc = unit->mTableSize * SAMPLEDUR * 65536.;
	//Print("TableLookup_SetTable unit->m_radtoinc %g %d %g\n", m_radtoinc, unit->mTableSize, rtwopi);
}
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////

void DegreeToKey_Ctor(DegreeToKey *unit)
{
	unit->m_fbufnum = std::numeric_limits<float>::quiet_NaN();
	if (BUFLENGTH == 1) {
		SETCALC(DegreeToKey_next_1);
//	} else if (INRATE(0) == calc_FullRate) {
	} else if (INRATE(1) == calc_FullRate) { // mtm fix: check index input, not bufnum (bufnum wouldn't be fullRate)
		SETCALC(DegreeToKey_next_a);
	} else {
		SETCALC(DegreeToKey_next_k);
	}
	unit->mOctave = (int32)ZIN0(2);
	unit->mPrevIndex = std::numeric_limits<int32>::max();
	unit->mPrevKey = 0.;
	printf("[DegreeToKey] init sample:\n\t");
	DegreeToKey_next_1(unit, 1);
	printf("[DegreeToKey] first sample:\n\t");
}

void DegreeToKey_next_1(DegreeToKey *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		int32 maxindex = tableSize - 1;

	int32 key, oct;

	int32 octave = unit->mOctave;
	float val;

	int32 index = (int32)floor(ZIN0(1));
	if (index == unit->mPrevIndex) {
		val = unit->mPrevKey;
	} else if (index < 0) {
		unit->mPrevIndex = index;
		key = tableSize + index % tableSize;
		oct = (index + 1) / tableSize - 1;
		val = unit->mPrevKey = table[key] + octave * oct;
	} else if (index > maxindex) {
		unit->mPrevIndex = index;
		key = index % tableSize;
		oct = index / tableSize;
		val = unit->mPrevKey = table[key] + octave * oct;
	} else {
		unit->mPrevIndex = index;
		val = unit->mPrevKey = table[index];
	}
	printf("\t[DegreeToKey_next_1] %f\n", val);
	ZOUT0(0) = val;
}

void DegreeToKey_next_k(DegreeToKey *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		int32 maxindex = tableSize - 1;

	float *out = ZOUT(0);

	int32 key, oct;
	float octave = unit->mOctave;
	float val;

	int32 index = (int32)floor(ZIN0(1));
	if (index == unit->mPrevIndex) {
		val = unit->mPrevKey;
	} else if (index < 0) {
		unit->mPrevIndex = index;
		key = tableSize + index % tableSize;
		oct = (index + 1) / tableSize - 1;
		val = unit->mPrevKey = table[key] + octave * oct;
	} else if (index > maxindex) {
		unit->mPrevIndex = index;
		key = index % tableSize;
		oct = index / tableSize;
		val = unit->mPrevKey = table[key] + octave * oct;
	} else {
		unit->mPrevIndex = index;
		val = unit->mPrevKey = table[index];
	}
	printf("\t[DegreeToKey_next_k] %f\n", val);
	LOOP1(inNumSamples,
		ZXP(out) = val;
	);
}


void DegreeToKey_next_a(DegreeToKey *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		int32 maxindex = tableSize - 1;

	float *out = ZOUT(0);
	float *in = ZIN(1);
	int32 previndex = unit->mPrevIndex;
	float prevkey = unit->mPrevKey;
	int32 key, oct;

	float octave = unit->mOctave;

	LOOP1(inNumSamples,
		int32 index = (int32)floor(ZXP(in));
		if (index == previndex) {
			printf("\t[DegreeToKey_next_a (index == previndex)] %f\n", prevkey);
			ZXP(out) = prevkey;
		} else if (index < 0) {
			previndex = index;
			key = tableSize + index % tableSize;
			oct = (index + 1) / tableSize - 1;
			printf("\t[DegreeToKey_next_a (index < 0)] %f\n", table[key] + octave * oct);
			ZXP(out) = prevkey = table[key] + octave * oct;
		} else if (index > maxindex) {
			previndex = index;
			key = index % tableSize;
			oct = index / tableSize;
			printf("\t[DegreeToKey_next_a (index > maxindex)] %f\n", table[key] + octave * oct);
			ZXP(out) = prevkey = table[key] + octave * oct;
		} else {
			previndex = index;
			printf("\t[DegreeToKey_next_a else] %f\n", table[index]);
			ZXP(out) = prevkey = table[index];
		}
	);
	unit->mPrevIndex = previndex;
	unit->mPrevKey = prevkey;
}

////////////////////////////////////////////////////////////////////////////////////

void Select_Ctor(Select *unit)
{
	if (BUFLENGTH == 1) {
		SETCALC(Select_next_1);
	} else if (INRATE(0) == calc_FullRate) {
		SETCALC(Select_next_a);
	} else {
		SETCALC(Select_next_k);
	}
	printf("[Select] init sample:\n\t");
	Select_next_1(unit, 1);
	printf("[Select] first sample:\n\t");
}

void Select_next_1(Select *unit, int inNumSamples)
{
	int32 maxindex = unit->mNumInputs - 1;
	int32 index = (int32)ZIN0(0) + 1;
	index = sc_clip(index, 1, maxindex);
float val = ZIN0(index);
printf("[Select_next_1] \t%f\n", val);
ZOUT0(0) = val;
//	ZOUT0(0) = ZIN0(index);
}

void Select_next_k(Select *unit, int inNumSamples)
{
	int32 maxindex = unit->mNumInputs - 1;
	int32 index = (int32)ZIN0(0) + 1;
	index = sc_clip(index, 1, maxindex);

	float *out = OUT(0);
	float *in = IN(index);

	Copy(inNumSamples, out, in);
}

void Select_next_a(Select *unit, int inNumSamples)
{
	int32 maxindex = unit->mNumInputs - 1;

	float *out = ZOUT(0);
	float *in0 = ZIN(0);
	float **in = unit->mInBuf;

	for (int i=0; i<inNumSamples; ++i) {
		int32 index = (int32)ZXP(in0) + 1;
		index = sc_clip(index, 1, maxindex);
float val = in[index][i];
printf("[Select_next_a] \t%f\n", val);
ZXP(out) = val;
//		ZXP(out) = in[index][i];
	}
}
////////////////////////////////////////////////////////////////////////////////////

void TWindex_Ctor(TWindex *unit)
{
	if (INRATE(0) == calc_FullRate) {
		SETCALC(TWindex_next_ak); // todo : ar
	} else {
		SETCALC(TWindex_next_k);
	}
	
	// Generate a random value here for the initialization sample
	// and to manage the trigger state of the first output frame
	int maxindex = unit->mNumInputs;
	int32 index = maxindex;
	float sum = 0.f;
	float maxSum = 0.f;
	float normalize = ZIN0(1);

	if(normalize == 1) {
		for (int32 k=2; k<maxindex; ++k) { maxSum += ZIN0(k); }
	} else {
		maxSum = 1.f;
	}
	RGen& rgen = *unit->mParent->mRGen;
	float max = maxSum * rgen.frand();

	for (int32 k=2; k<maxindex; ++k) {
		sum += ZIN0(k);
		if(sum >= max) {
			index = k - 2;
			break;
		}
	}

	printf("[TWindex] init sample %d:\n", index);
	OUT0(0) = index;

	// ensure a first-sample trigger doesn't cause another rand val,
	// so initialization sample will equal first output sample
	// (from user's perspective, the first sample is random either way)
	unit->m_trig = 1.f;
	unit->m_prevIndex = index;

//	unit->m_prevIndex = 0; // mtm
//	unit->m_trig = -1.f;   // mtm: make it trigger the first time << mtm: this won't ensure it triggers the first time, if input <= 0
	TWindex_next_k(unit, 1); // don't call this
	printf("[TWindex] first sample:\n\t");
}


void TWindex_next_k(TWindex *unit, int inNumSamples)
{
	float trig = ZIN0(0);
	float *out = ZOUT(0);
	int maxindex = unit->mNumInputs;
	int32 index = maxindex;
	
	if (trig > 0.f && unit->m_trig <= 0.f) {
		float maxSum = 0.f;
		float normalize = ZIN0(1); // switch normalisation on or off
		
		if(normalize == 1) {
			for (int32 k=2; k<maxindex; ++k) {
				maxSum += ZIN0(k);
			}
		} else {
			maxSum = 1.f;
		}
		RGen& rgen = *unit->mParent->mRGen;
		float max = maxSum * rgen.frand();
		float sum = 0.f;
		for (int32 k=2; k<maxindex; ++k) {
			sum += ZIN0(k);
			if(sum >= max) {
				index = k - 2;
				break;
			}
		}

		unit->m_prevIndex = index;
	} else {
		index = unit->m_prevIndex;
	}

	LOOP1(inNumSamples,
printf("[TWindex_next_k] \t%d\n", index); //mtm
			ZXP(out) = index;
	)
	unit->m_trig = trig;
}

void TWindex_next_ak(TWindex *unit, int inNumSamples)
{
	int maxindex = unit->mNumInputs;
	int32 index = maxindex;

//	float sum = 0.f; // mtm
	float maxSum = -1.f;
	float normalize = ZIN0(1);//switch normalisation on or off
	float *trig = ZIN(0);
	float *out = ZOUT(0);
	float curtrig;
	RGen& rgen = *unit->mParent->mRGen;

	LOOP1(inNumSamples,
		curtrig = ZXP(trig);
		if (curtrig > 0.f && unit->m_trig <= 0.f) {
			if (maxSum < 0.f) { // maxSum not yet set
				if (normalize == 1) {
					for (int32 k=2; k<maxindex; ++k) {
						maxSum += ZIN0(k);
					}
				} else {
					maxSum = 1.f;
				}
			}
			float sum = 0.f;
			float max = maxSum * rgen.frand();
			for (int32 k=2; k<maxindex; ++k) {
				sum += ZIN0(k);
				if(sum >= max) {
					index = k - 2;
					break;
				}
			}

			unit->m_prevIndex = index;
		} else {
			index = unit->m_prevIndex;
		}

printf("[TWindex_next_ak] \t%d\n", index); //mtm
		ZXP(out) = index;
		unit->m_trig = curtrig;
	)
}

////////////////////////////////////////////////////////////////////////////////////

void Index_Ctor(Index *unit)
{
	unit->m_fbufnum = std::numeric_limits<float>::quiet_NaN();
	if (BUFLENGTH == 1) {
		SETCALC(Index_next_1);
	} else if (INRATE(1) == calc_FullRate) {
		SETCALC(Index_next_a);
	} else {
		SETCALC(Index_next_k);
	}
	printf("[Index] init sample:\n\t");
	Index_next_1(unit, 1);
	printf("[Index] first sample:\n\t");
}

void Index_next_1(Index *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		int32 maxindex = tableSize - 1;

	int32 index = (int32)ZIN0(1);
	index = sc_clip(index, 0, maxindex);
	printf("[Index_next_1] \t%f\n", table[index]); //mtm
	ZOUT0(0) = table[index];
}

void Index_next_k(Index *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		int32 maxindex = tableSize - 1;

	float *out = ZOUT(0);

	int32 index = (int32)ZIN0(1);

	index = sc_clip(index, 0, maxindex);
	float val = table[index];
	printf("[Index_next_k] \t%f\n", val); //mtm
	LOOP1(inNumSamples,
		ZXP(out) = val;
	);
}


void Index_next_a(Index *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		int32 maxindex = tableSize - 1;

	float *out = ZOUT(0);
	float *in = ZIN(1);

	LOOP1(inNumSamples,
		int32 index = (int32)ZXP(in);
		index = sc_clip(index, 0, maxindex);
		printf("[Index_next_a] \t%f\n", table[index]); //mtm
		ZXP(out) = table[index];
	);
}


////////////////////////////////////////////////////////////////////////////////////

void IndexL_Ctor(IndexL *unit)
{
	unit->m_fbufnum = std::numeric_limits<float>::quiet_NaN();
	if (BUFLENGTH == 1) {
		SETCALC(IndexL_next_1);
	} else if (INRATE(1) == calc_FullRate) {
		SETCALC(IndexL_next_a);
	} else {
		SETCALC(IndexL_next_k);
	}
	printf("[IndexL] init sample:\n\t");
	IndexL_next_1(unit, 1);
	printf("[IndexL] first sample:\n\t");
}

void IndexL_next_1(IndexL *unit, int inNumSamples)
{
	// get table
	GET_TABLE
	const float *table = bufData;
	int32 maxindex = tableSize - 1;

	float findex = ZIN0(1);
	float frac = sc_frac(findex);

	int32 index = (int32)findex;
	index = sc_clip(index, 0, maxindex);

	float a = table[index];
	float b = table[sc_clip(index + 1, 0, maxindex)];
	printf("[IndexL_next_1] \t%f\n", lininterp(frac, a, b)); //mtm
	ZOUT0(0) = lininterp(frac, a, b);
}

void IndexL_next_k(IndexL *unit, int inNumSamples)
{
	// get table
	GET_TABLE
	const float *table = bufData;
	int32 maxindex = tableSize - 1;

	float *out = ZOUT(0);

	float findex = ZIN0(1);
	float frac = sc_frac(findex);

	int32 index = (int32)findex;
	index = sc_clip(index, 0, maxindex);


	float a = table[index];
	float b = table[sc_clip(index + 1, 0, maxindex)];
	float val = lininterp(frac, a, b);

	LOOP1(inNumSamples,
		  printf("[IndexL_next_k] \t%f\n", val); //mtm
		ZXP(out) = val;
	);
}


void IndexL_next_a(IndexL *unit, int inNumSamples)
{
	// get table
	GET_TABLE
	const float *table = bufData;
	int32 maxindex = tableSize - 1;

	float *out = ZOUT(0);
	float *in = ZIN(1);

	LOOP1(inNumSamples,
		float findex = ZXP(in);
		float frac = sc_frac(findex);
		int32 i1 = sc_clip((int32)findex, 0, maxindex);
		int32 i2 = sc_clip(i1 + 1, 0, maxindex);
		float a = table[i1];
		float b = table[i2];
		  printf("[IndexL_next_a] \t%f\n", lininterp(frac, a, b)); //mtm
		ZXP(out) =  lininterp(frac, a, b);
	);

}


////////////////////////////////////////////////////////////////////////////////////


void FoldIndex_Ctor(FoldIndex *unit)
{
	unit->m_fbufnum = std::numeric_limits<float>::quiet_NaN();
	if (BUFLENGTH == 1) {
		SETCALC(FoldIndex_next_1);
	} else if (INRATE(1) == calc_FullRate) {
		SETCALC(FoldIndex_next_a);
	} else {
		SETCALC(FoldIndex_next_k);
	}
	printf("[FoldIndex] init sample:\n\t");
	FoldIndex_next_1(unit, 1);
	printf("[FoldIndex] first sample:\n\t");
}

void FoldIndex_next_1(FoldIndex *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		int32 maxindex = tableSize - 1;

	int32 index = (int32)ZIN0(1);
	index = sc_fold(index, 0, maxindex);
	printf("[FoldIndex_next_1] \t%f\n", table[index]); //mtm
	ZOUT0(0) = table[index];
}

void FoldIndex_next_k(FoldIndex *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		int32 maxindex = tableSize - 1;

	int32 index = (int32)ZIN0(1);
	float *out = ZOUT(0);

	index = sc_fold(index, 0, maxindex);
	float val = table[index];
	LOOP1(inNumSamples,
		  printf("[FoldIndex_next_k] \t%f\n", val); //mtm
		ZXP(out) = val;
	);
}


void FoldIndex_next_a(FoldIndex *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		int32 maxindex = tableSize - 1;

	float *out = ZOUT(0);
	float *in = ZIN(1);

	LOOP1(inNumSamples,
		int32 index = (int32)ZXP(in);
		index = sc_fold(index, 0, maxindex);
		  printf("[FoldIndex_next_a] \t%f\n", table[index]); //mtm
		ZXP(out) = table[index];
	);
}


////////////////////////////////////////////////////////////////////////////////////

void WrapIndex_Ctor(WrapIndex *unit)
{
	unit->m_fbufnum = std::numeric_limits<float>::quiet_NaN();
	if (BUFLENGTH == 1) {
		SETCALC(WrapIndex_next_1);
	} else if (INRATE(1) == calc_FullRate) {
		SETCALC(WrapIndex_next_a);
	} else {
		SETCALC(WrapIndex_next_k);
	}
	printf("[WrapIndex] init sample:\n\t");
	WrapIndex_next_1(unit, 1);
	printf("[WrapIndex] first sample:\n\t");
}

void WrapIndex_next_1(WrapIndex *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		int32 maxindex = tableSize - 1;

	int32 index = (int32)floor(ZIN0(1));
	index = sc_wrap(index, 0, maxindex);
	printf("[WrapIndex_next_1] \t%f\n", table[index]); //mtm
	ZOUT0(0) = table[index];
}

void WrapIndex_next_k(WrapIndex *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		int32 maxindex = tableSize - 1;

	float *out = ZOUT(0);

	int32 index = (int32)ZIN0(1);
	index = sc_wrap(index, 0, maxindex);
	float val = table[index];
	LOOP1(inNumSamples,
		  printf("[WrapIndex_next_k] \t%f\n", val); //mtm
		ZXP(out) = val;
	);
}


void WrapIndex_next_a(WrapIndex *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		int32 maxindex = tableSize - 1;

	float *out = ZOUT(0);
	float *in = ZIN(1);

	LOOP1(inNumSamples,
		int32 index = (int32)ZXP(in);
		index = sc_wrap(index, 0, maxindex);
		  printf("[WrapIndex_next_a] \t%f\n", table[index]); //mtm
		ZXP(out) = table[index];
	);
}

////////////////////////////////////////////////////////////////////////////////////

static float IndexInBetween_FindIndex(const float* table, float in, int32 maxindex)
{
	for(int32 i = 0; i <= maxindex; i++) {
		if(table[i] > in) {
			if(i == 0) {
				return 0.f;
			} else {
				return ((in - table[i - 1]) / (table[i] - table[i - 1]) + i - 1);
			}
		}
	}
	return (float)maxindex;
}

void IndexInBetween_Ctor(IndexInBetween *unit)
{
	unit->m_fbufnum = std::numeric_limits<float>::quiet_NaN();
	if (BUFLENGTH == 1) {
		SETCALC(IndexInBetween_next_1);
	} else if (INRATE(1) == calc_FullRate) {
		SETCALC(IndexInBetween_next_a);
	} else {
		SETCALC(IndexInBetween_next_k);
	}
	printf("[IndexInBetween] init sample:\n\t");
	IndexInBetween_next_1(unit, 1);
	printf("[IndexInBetween] first sample:\n\t");
}

void IndexInBetween_next_1(IndexInBetween *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		int32 maxindex = tableSize - 1;

	float in = ZIN0(1);
	  float val = IndexInBetween_FindIndex(table, in, maxindex); // mtm
	  printf("[IndexInBetween_next_1] \t%f\n", val); //mtm
	ZOUT0(0) = val; // mtm
//	ZOUT0(0) = IndexInBetween_FindIndex(table, in, maxindex);
}

void IndexInBetween_next_k(IndexInBetween *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		int32 maxindex = tableSize - 1;

	float *out = ZOUT(0);
	float in = ZIN0(1);

	float val = IndexInBetween_FindIndex(table, in, maxindex);
	LOOP1(inNumSamples,
		    printf("[IndexInBetween_next_k] \t%f\n", val); //mtm
		ZXP(out) = val;
	);
}


void IndexInBetween_next_a(IndexInBetween *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		int32 maxindex = tableSize - 1;

	float *out = ZOUT(0);
	float *in = ZIN(1);

	LOOP1(inNumSamples,
		  float val = IndexInBetween_FindIndex(table, ZXP(in), maxindex); // mtm
		  printf("[IndexInBetween_next_a] \t%f\n", val); //mtm
		  ZXP(out) = val; // mtm
//		ZXP(out) = IndexInBetween_FindIndex(table, ZXP(in), maxindex);
	);
}


////////////////////////////////////////////////////////////////////////////////////

static int32 DetectIndex_FindIndex(const float* table, float in, int32 maxindex)
{
	int32 index;
	for(index = 0; index <= maxindex; index+=1) {
		if(table[index] == in) {
			return index;
		}
	}
	return -1;
}

void DetectIndex_Ctor(DetectIndex *unit)
{
	unit->m_fbufnum = std::numeric_limits<float>::quiet_NaN();
	if (BUFLENGTH == 1) {
		SETCALC(DetectIndex_next_1);
	} else if (INRATE(1) == calc_FullRate) {
		SETCALC(DetectIndex_next_a);
	} else {
		SETCALC(DetectIndex_next_k);
	}
	unit->mPrev = -1.f;
	// ensure in != unit->mPrevIn on first frame
	unit->mPrevIn = std::numeric_limits<float>::quiet_NaN(); //mtm
	
	printf("[DetectIndex] init sample:\n\t");
	DetectIndex_next_1(unit, 1);
	printf("[DetectIndex] first sample:\n\t");
}

void DetectIndex_next_1(DetectIndex *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		int32 maxindex = tableSize - 1;

	float in = ZIN0(1);
	int32 index;
	if(in == unit->mPrevIn) {
		index = (int32)unit->mPrev;
	} else {
		index = DetectIndex_FindIndex(table, in, maxindex);
		unit->mPrev = index;
		unit->mPrevIn = in;
	}
	    printf("[DetectIndex_next_1] \t%f\n", (float)index); //mtm
	ZOUT0(0) = (float)index;
}

void DetectIndex_next_k(DetectIndex *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		int32 maxindex = tableSize - 1;

	float *out = ZOUT(0);
	float in = ZIN0(1);
	int32 index;
	float val;
	if(in == unit->mPrevIn) {
		index = (int32)unit->mPrev;
	} else {
		index = DetectIndex_FindIndex(table, in, maxindex);
		unit->mPrev = index;
		unit->mPrevIn = in;
	};
	val = (float)index;
	LOOP1(inNumSamples,
		  printf("[DetectIndex_next_k] \t%f\n", val); //mtm
		ZXP(out) = val;
	);
}


void DetectIndex_next_a(DetectIndex *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		int32 maxindex = tableSize - 1;

	float *out = ZOUT(0);
	float *in = ZIN(1);
	float prev = unit->mPrevIn;
	int32 prevIndex = (int32)unit->mPrev;
	float inval;

	LOOP1(inNumSamples,
		inval = ZXP(in);
		if(inval != prev) {
			prevIndex = DetectIndex_FindIndex(table, inval, maxindex);
		}
		prev = inval;
		    printf("[DetectIndex_next_a] \t%f\n", (float)prevIndex); //mtm
		ZXP(out) = (float)prevIndex;
	);

	unit->mPrev = prevIndex;
	unit->mPrevIn = inval;
}


////////////////////////////////////////////////////////////////////////////////////

void Shaper_Ctor(Shaper *unit)
{
	unit->m_fbufnum = std::numeric_limits<float>::quiet_NaN();
	if (BUFLENGTH == 1) {
		SETCALC(Shaper_next_1);
	} else if (INRATE(1) == calc_FullRate) {
		SETCALC(Shaper_next_a);
	} else {
		SETCALC(Shaper_next_k);
	}
//	unit->mPrevIn = ZIN0(0); // mtm
	unit->mPrevIn = ZIN0(1); // mtm -- FIXed
	printf("[Shaper] init sample:\n\t"); //mtm
	Shaper_next_1(unit, 1);
	printf("[Shaper] first sample:\n\t");  //mtm
}

float force_inline ShaperPerform(const float * table0, const float * table1, float in, float offset, float fmaxindex)
{
	float findex = offset + in * offset;
	findex = sc_clip(findex, 0.f, fmaxindex);
	int32 index = (int32)findex;
	float pfrac = findex - (index - 1);
	index <<= 3;
	float val1 = *(const float*)((const char*)table0 + index);
	float val2 = *(const float*)((const char*)table1 + index);
	float val = val1 + val2 * pfrac;
	return val;
}

void Shaper_next_1(Shaper *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table0 = bufData;
		const float *table1 = table0 + 1;
		float fmaxindex = (float)(tableSize>>1) - 0.001;
		float offset = tableSize * 0.25;
	float val = ShaperPerform(table0, table1, ZIN0(1), offset, fmaxindex); // mtm
	printf("\t[Shaper_next_1] %f\n", val); //mtm
	ZOUT0(0) = val; //mtm
//	ZOUT0(0) = ShaperPerform(table0, table1, ZIN0(1), offset, fmaxindex); //mtm
}

void Shaper_next_k(Shaper *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table0 = bufData;
		const float *table1 = table0 + 1;
		float fmaxindex = (float)(tableSize>>1) - 0.001;
		float offset = tableSize * 0.25;

	float *out = ZOUT(0);
	float fin = ZIN0(1);

	if (fin == unit->mPrevIn) {
		LOOP1(inNumSamples,
			ZXP(out) = ShaperPerform(table0, table1, fin, offset, fmaxindex);
		);
	} else {
//		float phaseinc = (fin - unit->mPrevIn) * offset; //mtm
//		unit->mPrevIn = fin; //mtm
		
		float prevIn = unit->mPrevIn; //mtm
		float phaseinc = CALCSLOPE(fin, prevIn); //mtm
		fin = prevIn; //mtm

		LOOP1(inNumSamples,
			  float val = ShaperPerform(table0, table1, fin, offset, fmaxindex); // mtm
			  printf("\t[Shaper_next_k] %f\n", val); //mtm
			  ZXP(out) = val; //mtm
//			ZXP(out) = ShaperPerform(table0, table1, fin, offset, fmaxindex); // mtm
			fin += phaseinc;
		);
		unit->mPrevIn = fin; //mtm
	}
}

void Shaper_next_a(Shaper *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table0 = bufData;
		const float *table1 = table0 + 1;
		float fmaxindex = (float)(tableSize>>1) - 0.001;
		float offset = tableSize * 0.25;

	float *out = ZOUT(0);
	const float *in = ZIN(1);

	LOOP1(inNumSamples,
		float fin = ZXP(in);
		  float val = ShaperPerform(table0, table1, fin, offset, fmaxindex); // mtm
		  printf("\t[Shaper_next_a] %f\n", val); //mtm
		  ZXP(out) = val; //mtm
//		ZXP(out) = ShaperPerform(table0, table1, fin, offset, fmaxindex);
	);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

void FSinOsc_Ctor(FSinOsc *unit)
{
	double b1, y1, y2;
	if (INRATE(0) == calc_ScalarRate)
		SETCALC(FSinOsc_next_i);
	else
		SETCALC(FSinOsc_next);
	unit->m_freq = ZIN0(0);
	float iphase = ZIN0(1);
	float w = unit->m_freq * unit->mRate->mRadiansPerSample;
//	unit->m_b1 = 2. * cos(w); //mtm
	unit->m_b1 = b1 = 2. * cos(w);

//	unit->m_y1 = sin(iphase);
//	unit->m_y2 = sin(iphase - w);
//	printf("[FSinOsc] init sample:\n\t%f\n", unit->m_y1);//mtm
//	ZOUT0(0) = unit->m_y1; // this is y(-1) not (y0) mtm
//	printf("[FSinOsc] first sample:\n\t"); //mtm

	unit->m_y1 = y1 = sin(iphase - w);
	unit->m_y2 = y2 = sin(iphase - 2 * w);

	float outn = b1 * y1 - y2; //mtm
		printf("[FSinOsc] init sample:\n\t%f\n", outn);//mtm
	ZOUT0(0) = outn;//mtm
		printf("[FSinOsc] first sample:\n\t"); //mtm
}

void FSinOsc_next(FSinOsc *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	double freq = ZIN0(0);
	double b1;
	if (freq != unit->m_freq) {
		unit->m_freq = freq;
		double w = freq * unit->mRate->mRadiansPerSample;
		unit->m_b1 = b1 = 2.f * cos(w);
	} else {
		b1 = unit->m_b1;
	}
	double y0;
	double y1 = unit->m_y1;
	double y2 = unit->m_y2;
	//Print("y %g %g  b1 %g\n", y1, y2, b1);
	//Print("%d %d\n", unit->mRate->mFilterLoops, unit->mRate->mFilterRemain);
	LOOP(unit->mRate->mFilterLoops,
		ZXP(out) = y0 = b1 * y1 - y2;
		 printf("[FSinOsc_next_mFilterLoops] %f\n", y0);//mtm
		ZXP(out) = y2 = b1 * y0 - y1;
		  printf("[FSinOsc_next_mFilterLoops] %f\n", y2);//mtm
		ZXP(out) = y1 = b1 * y2 - y0;
		  printf("[FSinOsc_next_mFilterLoops] %f\n", y1);//mtm
	);
	LOOP(unit->mRate->mFilterRemain,
		ZXP(out) = y0 = b1 * y1 - y2;
		 printf("[FSinOsc_next_mFilterRemain] %f\n", y0);//mtm
		y2 = y1;
		y1 = y0;
	);
	//Print("y %g %g  b1 %g\n", y1, y2, b1);
	unit->m_y1 = y1;
	unit->m_y2 = y2;
}

void FSinOsc_next_i(FSinOsc *unit, int inNumSamples)
{
	float * __restrict__ out = ZOUT(0);
	double b1 = unit->m_b1;

	double y0;
	double y1 = unit->m_y1;
	double y2 = unit->m_y2;
	//Print("y %g %g  b1 %g\n", y1, y2, b1);
	//Print("%d %d\n", unit->mRate->mFilterLoops, unit->mRate->mFilterRemain);
	LOOP(unit->mRate->mFilterLoops,
		y0 = b1 * y1 - y2;
		y2 = b1 * y0 - y1;
		y1 = b1 * y2 - y0;
		 printf("[FSinOsc_next_i_mFilterLoops] %f\n", y0);//mtm
		 printf("[FSinOsc_next_i_mFilterLoops] %f\n", y2);//mtm
		 printf("[FSinOsc_next_i_mFilterLoops] %f\n", y1);//mtm
		ZXP(out) = y0;
		ZXP(out) = y2;
		ZXP(out) = y1;
	);
	LOOP(unit->mRate->mFilterRemain,
		ZXP(out) = y0 = b1 * y1 - y2;
		 printf("[FSinOsc_next_i_mFilterRemain] %f\n", y0);//mtm
		y2 = y1;
		y1 = y0;
	);
	//Print("y %g %g  b1 %g\n", y1, y2, b1);
	unit->m_y1 = y1;
	unit->m_y2 = y2;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

void PSinGrain_Ctor(PSinGrain *unit)
{
	SETCALC(PSinGrain_next);
	float freq = ZIN0(0);
	float dur  = ZIN0(1);
	float amp  = ZIN0(2);

	float w = freq * unit->mRate->mRadiansPerSample;
	float sdur = SAMPLERATE * dur;

	float rdur = 1.f / sdur;
	float rdur2 = rdur * rdur;

	unit->m_level = 0.f;
	unit->m_slope = 4.0 * (rdur - rdur2);	// ampslope
	unit->m_curve = -8.0 * rdur2;			// ampcurve
	unit->mCounter = (int32)(sdur + .5);

	/* calc feedback param and initial conditions */
//	unit->m_b1 = 2. * cos(w);
//	unit->m_y1 = 0.f;
//	unit->m_y2 = -sin(w) * amp;
//	ZOUT0(0) = 0.f;
	
	double b1, y1, y2;
	unit->m_b1 = b1 = 2. * cos(w);
	unit->m_y1 = y1 = -sin(w) * amp;;
	unit->m_y2 = y2 = -sin(w + w) * amp;
	
	float outn = b1 * y1 - y2; //mtm
	printf("[PSinGrain] init sample:\n\t%f\n", outn);//mtm
	ZOUT0(0) = outn;//mtm
	printf("[PSinGrain] first sample:\n\t"); //mtm
}


void PSinGrain_next(PSinGrain *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float y0;
	float y1 = unit->m_y1;
	float y2 = unit->m_y2;
	float b1 = unit->m_b1;
	float level = unit->m_level;
	float slope = unit->m_slope;
	float curve = unit->m_curve;
	int32 counter = unit->mCounter;
	int32 remain = inNumSamples;
	int32 nsmps;
	do {
		if (counter<=0) {
			nsmps = remain;
			remain = 0;
				printf("[PSinGrain_next_clear] %f\n", 0.f);//mtm
			LOOP(nsmps, ZXP(out) = 0.f;); // can't use Clear bcs might not be aligned
		} else {
			nsmps = sc_min(remain, counter);
			remain -= nsmps;
			counter -= nsmps;
			if (nsmps == inNumSamples) {
				nsmps = unit->mRate->mFilterLoops;
				LOOP(nsmps,
					y0 = b1 * y1 - y2;
					 printf("[PSinGrain_next_mFilterLoops] %f\n", y0 * level);//mtm
					ZXP(out) = y0 * level;
					level += slope;
					slope += curve;
					y2 = b1 * y0 - y1;
					 printf("[PSinGrain_next_mFilterLoops] %f\n", y2 * level);//mtm
					ZXP(out) = y2 * level;
					level += slope;
					slope += curve;
					y1 = b1 * y2 - y0;
					 printf("[PSinGrain_next_mFilterLoops] %f\n", y1 * level);//mtm
					ZXP(out) = y1 * level;
					level += slope;
					slope += curve;
				);
				nsmps = unit->mRate->mFilterRemain;
				LOOP(nsmps,
					y0 = b1 * y1 - y2;
					y2 = y1;
					y1 = y0;
					 printf("[PSinGrain_next_mFilterRemain] %f\n", y0 * level);//mtm
					ZXP(out) = y0 * level;
					level += slope;
					slope += curve;
				);
			} else {
				LOOP(nsmps,
					y0 = b1 * y1 - y2;
					y2 = y1;
					y1 = y0;
					 printf("[PSinGrain_next_else2] %f\n", y0 * level);//mtm
					ZXP(out) = y0 * level;
					level += slope;
					slope += curve;
				);
			}
			if (counter == 0) {
				NodeEnd(&unit->mParent->mNode);
			}
		}
	} while (remain>0);
	unit->mCounter = counter;
	unit->m_level = level;
	unit->m_slope = slope;
	unit->m_y1 = y1;
	unit->m_y2 = y2;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename OscType, int FreqInputIndex>
force_inline void Osc_ikk_perform(OscType *unit, const float * table0, const float * table1, int inNumSamples)
{
	float *out = ZOUT(0);
	float freqin = ZIN0(FreqInputIndex);
	float phasein = ZIN0(FreqInputIndex + 1);

	int32 phase = unit->m_phase;
	int32 lomask = unit->m_lomask;

	int32 freq = (int32)(unit->m_cpstoinc * freqin);
	int32 phaseinc = freq + (int32)(CALCSLOPE(phasein, unit->m_phasein) * unit->m_radtoinc);
	unit->m_phasein = phasein;

	LOOP1(inNumSamples,
		  float val = lookupi1(table0, table1, phase, lomask);//mtm
		  printf("[Osc_ikk_perform] %f\n", val);//mtm
		  ZXP(out) = val;//mtm
//		ZXP(out) = lookupi1(table0, table1, phase, lomask);
		phase += phaseinc;
	);
	unit->m_phase = phase;
}

void SinOsc_next_ikk(SinOsc *unit, int inNumSamples)
{
	float *table0 = ft->mSineWavetable;
	float *table1 = table0 + 1;

	Osc_ikk_perform<SinOsc, 0>(unit, table0, table1, inNumSamples);
}


template <typename OscType, int FreqInputIndex>
force_inline void Osc_ika_perform(OscType *unit, const float * table0, const float * table1, int inNumSamples)
{
	float *out = ZOUT(0);
	float freqin = ZIN0(FreqInputIndex);
	float *phasein = ZIN(FreqInputIndex + 1);

	int32 phase = unit->m_phase;
	int32 lomask = unit->m_lomask;

	int32 freq = (int32)(unit->m_cpstoinc * freqin);
	float radtoinc = unit->m_radtoinc;
	LOOP1(inNumSamples,
		int32 phaseoffset = phase + (int32)(radtoinc * ZXP(phasein));
		  float val = lookupi1(table0, table1, phaseoffset, lomask);//mtm
		  printf("[Osc_ika_perform] %f\n", val);//mtm
		  ZXP(out) = val;//mtm
//		ZXP(out) = lookupi1(table0, table1, phaseoffset, lomask);//mtm
		phase += freq;
	);
	unit->m_phase = phase;
}

void SinOsc_next_ika(SinOsc *unit, int inNumSamples)
{
	const float *table0 = ft->mSineWavetable;
	const float *table1 = table0 + 1;
	Osc_ika_perform<SinOsc, 0>(unit, table0, table1, inNumSamples);
}


template <typename OscType, int FreqInputIndex>
force_inline void Osc_iaa_perform(OscType * unit, const float * table0, const float * table1, int inNumSamples)
{
	float *out = ZOUT(0);
	float *freqin = ZIN(FreqInputIndex);
	float *phasein = ZIN(FreqInputIndex + 1);

	int32 phase = unit->m_phase;
	int32 lomask = unit->m_lomask;

	float cpstoinc = unit->m_cpstoinc;
	float radtoinc = unit->m_radtoinc;
	LOOP1(inNumSamples,
		float phaseIn = ZXP(phasein);
		float freqIn  = ZXP(freqin);
		int32 phaseoffset = phase + (int32)(radtoinc * phaseIn);
		float z = lookupi1(table0, table1, phaseoffset, lomask);
		phase += (int32)(cpstoinc * freqIn);
		  printf("[Osc_iaa_perform] %f\n", z);//mtm
		ZXP(out) = z;
	);
	unit->m_phase = phase;
}

void SinOsc_next_iaa(SinOsc *unit, int inNumSamples)
{
	const float *table0 = ft->mSineWavetable;
	const float *table1 = table0 + 1;

	Osc_iaa_perform<SinOsc, 0>(unit, table0, table1, inNumSamples);
}


template <typename OscType, int FreqInputIndex>
force_inline void Osc_iak_perform(OscType *unit, const float * table0, const float * table1, int inNumSamples)
{
	float *out = ZOUT(0);
	float *freqin = ZIN(FreqInputIndex);
	float phasein = ZIN0(FreqInputIndex + 1);

	int32 phase = unit->m_phase;
	int32 lomask = unit->m_lomask;

	float cpstoinc = unit->m_cpstoinc;
	float radtoinc = unit->m_radtoinc;
	float phasemod = unit->m_phasein;

	if (phasein != phasemod) {
		float phaseslope = CALCSLOPE(phasein, phasemod);

		LOOP1(inNumSamples,
			int32 pphase = phase + (int32)(radtoinc * phasemod);
			phasemod += phaseslope;
			float z = lookupi1(table0, table1, pphase, lomask);
			phase += (int32)(cpstoinc * ZXP(freqin));
			  printf("[Osc_iak_perform-phasemod] %f\n", z);//mtm
			ZXP(out) = z;
		);
	} else {
		LOOP1(inNumSamples,
			int32 pphase = phase + (int32)(radtoinc * phasemod);
			float z = lookupi1(table0, table1, pphase, lomask);
			  printf("[Osc_iak_perform-staticphase] %f\n", z);//mtm
			phase += (int32)(cpstoinc * ZXP(freqin));
			ZXP(out) = z;
		);
	}
	unit->m_phase = phase;
	unit->m_phasein = phasein;
}

void SinOsc_next_iak(SinOsc *unit, int inNumSamples)
{
	float *table0 = ft->mSineWavetable;
	float *table1 = table0 + 1;

	Osc_iak_perform<SinOsc, 0>(unit, table0, table1, inNumSamples);
}


template <typename OscType, int FreqInputIndex>
force_inline void Osc_iai_perform(OscType *unit, const float * table0, const float * table1, int inNumSamples)
{
	float *out = ZOUT(0);
	float *freqin = ZIN(FreqInputIndex);

	int32 phase = unit->m_phase;
	int32 lomask = unit->m_lomask;

	float cpstoinc = unit->m_cpstoinc;
	float radtoinc = unit->m_radtoinc;
	float phasemod = unit->m_phasein;

	LOOP1(inNumSamples,
		  int32 pphase = phase + (int32)(radtoinc * phasemod);
		  float z = lookupi1(table0, table1, pphase, lomask);
		  phase += (int32)(cpstoinc * ZXP(freqin));
			  printf("[Osc_iai_perform] %f\n", z);//mtm
		  ZXP(out) = z;
	);
	unit->m_phase = phase;
}

void SinOsc_next_iai(SinOsc *unit, int inNumSamples)
{
	float *table0 = ft->mSineWavetable;
	float *table1 = table0 + 1;
	Osc_iai_perform<SinOsc, 0>(unit, table0, table1, inNumSamples);
}


void SinOsc_Ctor(SinOsc *unit)
{
	int tableSize2 = ft->mSineSize;
	unit->m_phasein = ZIN0(1);
	unit->m_radtoinc = tableSize2 * (rtwopi * 65536.);
	unit->m_cpstoinc = tableSize2 * SAMPLEDUR * 65536.;
	unit->m_lomask = (tableSize2 - 1) << 3;
	
	int32 initPhase;//mtm
	printf("[SinOsc] init sample:\n\t");//mtm
	
	if (INRATE(0) == calc_FullRate) {
		if (INRATE(1) == calc_FullRate)
			SETCALC(SinOsc_next_iaa);
		else if (INRATE(1) == calc_BufRate)
			SETCALC(SinOsc_next_iak);
		else
			SETCALC(SinOsc_next_iai);

//		unit->m_phase = 0;//mtm
		unit->m_phase = initPhase = 0;//mtm
		SinOsc_next_iaa(unit, 1);//mtm
	} else {
		if (INRATE(1) == calc_FullRate) {
			//Print("next_ika\n");
			SETCALC(SinOsc_next_ika);
//			unit->m_phase = 0;//mtm
			unit->m_phase = initPhase = 0;//mtm
			SinOsc_next_iaa(unit, 1);//mtm
		} else {
			SETCALC(SinOsc_next_ikk);
//			unit->m_phase = (int32)(unit->m_phasein * unit->m_radtoinc);//mtm
			unit->m_phase = initPhase = (int32)(unit->m_phasein * unit->m_radtoinc);//mtm
			SinOsc_next_ikk(unit, 1);//mtm
		}
	}
	
//	int32 initPhase;
//	unit->m_phase = initPhase = (int32)(unit->m_phasein * unit->m_radtoinc);//mtm
//	int32 initPhase = unit->m_phase;
//	SinOsc_next_ikk(unit, 1);
	unit->m_phase = initPhase; // restore initial state mtm
	printf("[SinOsc] first sample:\n\t");
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////!!!

void SinOscFB_next_kk(SinOscFB *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freqin = ZIN0(0);

	float feedback = unit->m_feedback;
	float nextFeedback = ZIN0(1) * unit->m_radtoinc;

	float *table0 = ft->mSineWavetable;
	float *table1 = table0 + 1;

	int32 phase = unit->m_phase;
	int32 lomask = unit->m_lomask;
	float prevout = unit->m_prevout;
	float feedback_slope = CALCSLOPE(nextFeedback, feedback);
	int32 freq = (int32)(unit->m_cpstoinc * freqin);

	LooP(inNumSamples) {
		prevout = lookupi1(table0, table1, phase + (int32)(feedback * prevout), lomask);
		  printf("[SinOscFB_next_kk] %f\n", prevout);//mtm
		ZXP(out) = prevout;
		phase += freq;
		feedback += feedback_slope;
	}
	unit->m_phase = phase;
	unit->m_prevout = prevout;
	unit->m_feedback = feedback;
}

void SinOscFB_Ctor(SinOscFB *unit)
{
	//Print("next_ik\n");
	SETCALC(SinOscFB_next_kk);

	int tableSize2 = ft->mSineSize;
	unit->m_lomask = (tableSize2 - 1) << 3;
	unit->m_radtoinc = tableSize2 * (rtwopi * 65536.);
	unit->m_cpstoinc = tableSize2 * SAMPLEDUR * 65536.;
	unit->m_prevout = 0.;
	unit->m_feedback = ZIN0(1) * unit->m_radtoinc;

//	float *table0 = ft->mSineWavetable;//mtm
//	float *table1 = table0 + 1;//mtm
//	int32 lomask = unit->m_lomask;//mtm
//	int32 prevphase = -(int32)(unit->m_cpstoinc * ZIN0(0));//mtm
//	unit->m_prevout = lookupi1(table0, table1, prevphase, lomask);//mtm

	unit->m_phase = 0;
	printf("[SinOscFB] init sample:\n\t");
	SinOscFB_next_kk(unit, 1);
	printf("[SinOscFB] first sample:\n\t");
	unit->m_phase = 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

void Osc_Ctor(Osc *unit)
{
	unit->mTableSize = -1;

	float fbufnum = ZIN0(0);
	uint32 bufnum = (uint32)fbufnum;
	World *world = unit->mWorld;

	SndBuf *buf;
	if (bufnum >= world->mNumSndBufs) {
		int localBufNum = bufnum - world->mNumSndBufs;
		Graph *parent = unit->mParent;
		if(localBufNum <= parent->localBufNum) {
			buf = unit->m_buf = parent->mLocalSndBufs + localBufNum;
		} else {
			buf = unit->m_buf = world->mSndBufs;
		}
	} else {
		buf = unit->m_buf = world->mSndBufs + bufnum;
	}

	int tableSize = buf->samples;
	int tableSize2 = tableSize >> 1;
	unit->m_radtoinc = tableSize2 * (rtwopi * 65536.); // Osc, OscN, PMOsc

	unit->m_phasein = ZIN0(2);
	int32 initphase;//mtm
	printf("[Osc] init sample:\n\t");//mtm
	
	if (INRATE(1) == calc_FullRate) {
		if (INRATE(2) == calc_FullRate) {
			//Print("next_iaa\n");
			SETCALC(Osc_next_iaa);
//			unit->m_phase = 0;//mtm
		} else {
			//Print("next_iak\n");
			SETCALC(Osc_next_iak);
//			unit->m_phase = 0; //mtm
		}
		unit->m_phase = initphase = 0;
		Osc_next_iaa(unit, 1);
	} else {
		if (INRATE(2) == calc_FullRate) {
			//Print("next_ika\n");
			SETCALC(Osc_next_ika);
			unit->m_phase = initphase = 0;
			Osc_next_iaa(unit, 1);
		} else {
			//Print("next_ikk\n");
			SETCALC(Osc_next_ikk);
			unit->m_phase = initphase = (int32)(unit->m_phasein * unit->m_radtoinc);
			Osc_next_ikk(unit, 1);
		}
	}

//	Osc_next_ikk(unit, 1);//mtm
	printf("[Osc] first sample:\n\t");//mtm
	unit->m_phase = initphase;
}

force_inline bool Osc_get_table(Osc *unit, const float *& table0, const float *& table1, int inNumSamples)
{
	const SndBuf * buf;
	const float * bufData;
	int tableSize;
	bool tableValid = UnitGetTable(unit, inNumSamples, buf, bufData, tableSize);
	if (!tableValid)
		return false;

	table0 = bufData;
	table1 = table0 + 1;
	if (tableSize != unit->mTableSize) {
		unit->mTableSize = tableSize;
		int tableSize2 = tableSize >> 1;
		unit->m_lomask = (tableSize2 - 1) << 3; // Osc, OscN, COsc, COsc, COsc2, OscX4, OscX2
		unit->m_radtoinc = tableSize2 * (rtwopi * 65536.); // Osc, OscN, PMOsc
		// Osc, OscN, PMOsc, COsc, COsc2, OscX4, OscX2
		unit->m_cpstoinc = tableSize2 * SAMPLEDUR * 65536.;
	}
	if (!verify_wavetable(unit, "Osc", tableSize, inNumSamples)) return false;

	return true;
}

void Osc_next_ikk(Osc *unit, int inNumSamples)
{
	const float * table0;
	const float * table1;
	bool tableValid = Osc_get_table(unit, table0, table1, inNumSamples);
	if (!tableValid) return;

	Osc_ikk_perform<Osc, 1>(unit, table0, table1, inNumSamples);
}

void Osc_next_ika(Osc *unit, int inNumSamples)
{
	const float * table0;
	const float * table1;
	bool tableValid = Osc_get_table(unit, table0, table1, inNumSamples);
	if (!tableValid) return;

	Osc_ika_perform<Osc, 1>(unit, table0, table1, inNumSamples);
}

void Osc_next_iaa(Osc *unit, int inNumSamples)
{
	const float * table0;
	const float * table1;
	bool tableValid = Osc_get_table(unit, table0, table1, inNumSamples);
	if (!tableValid) return;

	Osc_iaa_perform<Osc, 1>(unit, table0, table1, inNumSamples);
}


void Osc_next_iak(Osc *unit, int inNumSamples)
{
	const float * table0;
	const float * table1;
	bool tableValid = Osc_get_table(unit, table0, table1, inNumSamples);
	if (!tableValid) return;

	Osc_iak_perform<Osc, 1>(unit, table0, table1, inNumSamples);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void OscN_Ctor(OscN *unit)
{
	unit->mTableSize = -1;

	float fbufnum = ZIN0(0);
	uint32 bufnum = (uint32)fbufnum;
	World *world = unit->mWorld;
	SndBuf *buf;
	if (bufnum >= world->mNumSndBufs) {
		int localBufNum = bufnum - world->mNumSndBufs;
		Graph *parent = unit->mParent;
		if(localBufNum <= parent->localBufNum) {
			buf = unit->m_buf = parent->mLocalSndBufs + localBufNum;
		} else {
			buf = unit->m_buf = world->mSndBufs;
		}
	} else {
		buf = unit->m_buf = world->mSndBufs + bufnum;
	}

	int tableSize = buf->samples;
	unit->m_radtoinc = tableSize * (rtwopi * 65536.);
	
	unit->m_phasein = ZIN0(2);
	int32 initphase;
	if (INRATE(1) == calc_FullRate) {
		if (INRATE(2) == calc_FullRate) {
			SETCALC(OscN_next_naa);
		} else {
			SETCALC(OscN_next_nak);
		}
		unit->m_phase = initphase = 0;//mtm
		OscN_next_naa(unit, 1);//mtm
	} else {
		if (INRATE(2) == calc_FullRate) {
			SETCALC(OscN_next_nka);
			unit->m_phase = initphase = 0;//mtm
			OscN_next_naa(unit, 1);//mtm
		} else {
			SETCALC(OscN_next_nkk);
			unit->m_phase = initphase = (int32)(unit->m_phasein * unit->m_radtoinc);//mtm
			OscN_next_nkk(unit, 1);//mtm
		}
	}

	printf("[OscN] init sample:\n\t");
//	OscN_next_nkk(unit, 1);
	printf("[OscN] first sample:\n\t");
	unit->m_phase = initphase;//mtm
}


void OscN_next_nkk(OscN *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		if (tableSize != unit->mTableSize) {
			unit->mTableSize = tableSize;
			unit->m_lomask = (tableSize - 1) << 2;
			unit->m_radtoinc = tableSize * (rtwopi * 65536.);
			unit->m_cpstoinc = tableSize * SAMPLEDUR * 65536.;
		}

	if (!verify_wavetable(unit, "OscN", tableSize, inNumSamples)) return;

	float *out = ZOUT(0);
	float freqin = ZIN0(1);
	float phasein = ZIN0(2);

	int32 phase = unit->m_phase;
	int32 lomask = unit->m_lomask;

	int32 freq = (int32)(unit->m_cpstoinc * freqin);
	int32 phaseinc = freq + (int32)(CALCSLOPE(phasein, unit->m_phasein) * unit->m_radtoinc);
	unit->m_phasein = phasein;

	LOOP1(inNumSamples,
		  float val = *(float*)((char*)table + ((phase >> xlobits) & lomask));
		  printf("[OscN_next_nkk] %f\n", val);//mtm
		  ZXP(out) = val;
//		ZXP(out) = *(float*)((char*)table + ((phase >> xlobits) & lomask));//mtm
		phase += phaseinc;
	);
	unit->m_phase = phase;
}



void OscN_next_nka(OscN *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		if (tableSize != unit->mTableSize) {
			unit->mTableSize = tableSize;
			unit->m_lomask = (tableSize - 1) << 2;
			unit->m_radtoinc = tableSize * (rtwopi * 65536.);
			unit->m_cpstoinc = tableSize * SAMPLEDUR * 65536.;
		}

	if (!verify_wavetable(unit, "OscN", tableSize, inNumSamples)) return;

	float *out = ZOUT(0);
	float freqin = ZIN0(1);
	float *phasein = ZIN(2);

	int32 phase = unit->m_phase;
	int32 lomask = unit->m_lomask;

	int32 freq = (int32)(unit->m_cpstoinc * freqin);
	float radtoinc = unit->m_radtoinc;
	LOOP1(inNumSamples,
		int32 pphase = phase + (int32)(radtoinc * ZXP(phasein));
		  printf("[OscN_next_nka] %f\n", *(float*)((char*)table + ((pphase >> xlobits) & lomask)));//mtm
		ZXP(out) = *(float*)((char*)table + ((pphase >> xlobits) & lomask));
		phase += freq;
	);
	unit->m_phase = phase;
}

void OscN_next_naa(OscN *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		if (tableSize != unit->mTableSize) {
			unit->mTableSize = tableSize;
			unit->m_lomask = (tableSize - 1) << 2;
			unit->m_radtoinc = tableSize * (rtwopi * 65536.);
			unit->m_cpstoinc = tableSize * SAMPLEDUR * 65536.;
		}

	if (!verify_wavetable(unit, "OscN", tableSize, inNumSamples)) return;

	float *out = ZOUT(0);
	float *freqin = ZIN(1);
	float *phasein = ZIN(2);

	int32 phase = unit->m_phase;
	int32 lomask = unit->m_lomask;

	float cpstoinc = unit->m_cpstoinc;
	float radtoinc = unit->m_radtoinc;
	LOOP1(inNumSamples,
		int32 pphase = phase + (int32)(radtoinc * ZXP(phasein));
		float z = *(float*)((char*)table + ((pphase >> xlobits) & lomask));
		phase += (int32)(cpstoinc * ZXP(freqin));
		  printf("[OscN_next_naa] %f\n", z);//mtm
		ZXP(out) = z;
	);
	unit->m_phase = phase;
}


void OscN_next_nak(OscN *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table = bufData;
		if (tableSize != unit->mTableSize) {
			unit->mTableSize = tableSize;
			unit->m_lomask = (tableSize - 1) << 2;
			unit->m_radtoinc = tableSize * (rtwopi * 65536.);
			unit->m_cpstoinc = tableSize * SAMPLEDUR * 65536.;
		}

	if (!verify_wavetable(unit, "OscN", tableSize, inNumSamples)) return;

	float *out = ZOUT(0);
	float *freqin = ZIN(1);
	float phasein = ZIN0(2);

	int32 phase = unit->m_phase;
	int32 lomask = unit->m_lomask;

	float cpstoinc = unit->m_cpstoinc;
	float radtoinc = unit->m_radtoinc;
	float phasemod = unit->m_phasein;
	float phaseslope = CALCSLOPE(phasein, phasemod);

	LOOP1(inNumSamples,
		int32 pphase = phase + (int32)(radtoinc * phasemod);
		phasemod += phaseslope;
		float z = *(float*)((char*)table + ((pphase >> xlobits) & lomask));
		phase += (int32)(cpstoinc * ZXP(freqin));
		    printf("[OscN_next_nak] %f\n", z);//mtm
		ZXP(out) = z;
	);
	unit->m_phase = phase;
	unit->m_phasein = phasein;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

void COsc_Ctor(COsc *unit)
{
	unit->m_fbufnum = std::numeric_limits<float>::quiet_NaN();
	SETCALC(COsc_next);
	unit->m_phase1 = 0;
	unit->m_phase2 = 0;
	unit->mTableSize = -1;
	printf("[COsc] init sample:\n\t");
	COsc_next(unit, 1);
	printf("[COsc] first sample:\n\t");
	unit->m_phase1 = 0;//mtm
	unit->m_phase2 = 0;//mtm
}


void COsc_next(COsc *unit, int inNumSamples)
{
	// get table
	GET_TABLE
		const float *table0 = bufData;
		const float *table1 = table0 + 1;
		if (tableSize != unit->mTableSize) {
			unit->mTableSize = tableSize;
			int tableSize2 = tableSize >> 1;
			unit->m_lomask = (tableSize2 - 1) << 3; // Osc, OscN, COsc, COsc, COsc2, OscX4, OscX2
			// Osc, OscN, PMOsc, COsc, COsc2, OscX4, OscX2
			unit->m_cpstoinc = tableSize2 * SAMPLEDUR * 65536.;
		}
	if (!verify_wavetable(unit, "COsc", tableSize, inNumSamples)) return;

	float *out = ZOUT(0);
	float freqin = ZIN0(1);
	float beats = ZIN0(2) * 0.5f;

	int32 phase1 = unit->m_phase1;
	int32 phase2 = unit->m_phase2;
	int32 lomask = unit->m_lomask;

	int32 cfreq = (int32)(unit->m_cpstoinc * freqin);
	int32 beatf = (int32)(unit->m_cpstoinc * beats);
	int32 freq1 = cfreq + beatf;
	int32 freq2 = cfreq - beatf;
	LOOP1(inNumSamples,
		float a = lookupi1(table0, table1, phase1, lomask);
		float b = lookupi1(table0, table1, phase2, lomask);
		   printf("[COsc_next] %f\n", a + b);//mtm
		ZXP(out) = a + b;
		phase1 += freq1;
		phase2 += freq2;
	);
	unit->m_phase1 = phase1;
	unit->m_phase2 = phase2;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

static inline const SndBuf * VOscGetBuf(int & bufnum, World * world, Unit * unit)
{
	if (bufnum < 0)
		bufnum = 0;

	const SndBuf * bufs;
//	 printf("[VOsc] bufnum %d\n", bufnum);//mtm
	if (bufnum+1 >= world->mNumSndBufs) {
		int localBufNum = bufnum - world->mNumSndBufs;
		Graph *parent = unit->mParent;
		if(localBufNum <= parent->localBufNum) {
			bufs = parent->mLocalSndBufs + localBufNum;
		} else {
			bufnum = 0;
			bufs = world->mSndBufs + bufnum;
		}
	} else {
		if (bufnum >= world->mNumSndBufs)
			bufnum = 0;
		bufs = world->mSndBufs + sc_max(0, bufnum);
	}
	return bufs;
}

void VOsc_Ctor(VOsc *unit)
{
	float nextbufpos = ZIN0(0);
	unit->m_bufpos = nextbufpos;
	int bufnum = sc_floor(nextbufpos);
	World *world = unit->mWorld;

	const SndBuf *bufs = VOscGetBuf(bufnum, world, unit);

	int tableSize = bufs[0].samples;

	unit->mTableSize = tableSize;
	int tableSize2 = tableSize >> 1;
	unit->m_lomask = (tableSize2 - 1) << 3;
	unit->m_radtoinc = tableSize2 * (rtwopi * 65536.);
	unit->m_cpstoinc = tableSize2 * SAMPLEDUR * 65536.;

	unit->m_phasein = ZIN0(2);
	unit->m_phaseoffset = (int32)(unit->m_phasein * unit->m_radtoinc);
	
	double initphase; //mtm
	printf("[VOsc] init sample:\n\t");//mtm
	if (INRATE(2) == calc_FullRate) {
		SETCALC(VOsc_next_ika);
		unit->m_phase = initphase = 0;//mtm
		VOsc_next_ika(unit, 1);//mtm
	} else {
		SETCALC(VOsc_next_ikk);
		unit->m_phase = initphase =unit->m_phaseoffset;//mtm
		VOsc_next_ikk(unit, 1);//mtm
	}
	
	unit->m_phase = initphase;
//	VOsc_next_ikk(unit, 1);//mtm
	printf("[VOsc] first sample:\n\t");//mtm
}

void VOsc_next_ikk(VOsc *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float nextbufpos = ZIN0(0);
	float freqin = ZIN0(1);
	float phasein = ZIN0(2);

	float prevbufpos = unit->m_bufpos;
	float bufdiff = nextbufpos - prevbufpos;

	int32 phase = unit->m_phase;
	int32 lomask = unit->m_lomask;

	int32 freq = (int32)(unit->m_cpstoinc * freqin);
	int32 phaseinc = freq + (int32)(CALCSLOPE(phasein, unit->m_phasein) * unit->m_radtoinc);
	unit->m_phasein = phasein;
	int tableSize = unit->mTableSize;
	float cur = prevbufpos;
	World *world = unit->mWorld;

	if (bufdiff == 0.f) {
		float level = cur - sc_floor(cur);
		int32 bufnum = (int)sc_floor(cur);
		
//		printf("[VOsc] level %f, bufnum %d\n", level, bufnum);//mtm
		const SndBuf *bufs = VOscGetBuf(bufnum, world, unit);
		if (!verify_wavetable(unit, "VOsc", tableSize, inNumSamples)) return;//mtm
		
		const float *table0  = bufs[0].data;
		const float *table2  = bufs[1].data;
		if (!table0 || !table2 || tableSize != bufs[0].samples|| tableSize != bufs[1].samples) {
			ClearUnitOutputs(unit, inNumSamples);
			return;
		}

		const float *table1 = table0 + 1;
		const float *table3 = table2 + 1;

		LOOP1(inNumSamples,
			float pfrac = PhaseFrac1(phase);
			uint32 index = ((phase >> xlobits1) & lomask);
			float val0 = *(float*)((char*)table0 + index);
			float val1 = *(float*)((char*)table1 + index);
			float val2 = *(float*)((char*)table2 + index);
			float val3 = *(float*)((char*)table3 + index);
			float a = val0 + val1 * pfrac;
			float b = val2 + val3 * pfrac;
			  printf("[VOsc_next_ikk] %f\n", a + level * (b - a));//mtm
			ZXP(out) = a + level * (b - a);
			phase += phaseinc;
		);
	} else {
		int nsmps;
		int donesmps = 0;
		int remain = inNumSamples;
		while (remain) {
			float level = cur - sc_floor(cur);

			float cut;
			if (bufdiff > 0.) {
				cut = sc_min(nextbufpos, sc_floor(cur+1.f));
			} else {
				cut = sc_max(nextbufpos, sc_ceil(cur-1.f));
			}

			float sweepdiff = cut - cur;
			if (cut == nextbufpos) nsmps = remain;
			else {
				float sweep = (float)inNumSamples / bufdiff;
				nsmps = (int)sc_floor(sweep * sweepdiff + 0.5f) - donesmps;
				nsmps = sc_clip(nsmps, 1, remain);
			}

			float slope = sweepdiff / (float)nsmps;

			int32 bufnum = (int32)sc_floor(cur);
			const SndBuf *bufs = VOscGetBuf(bufnum, world, unit);
			if (!verify_wavetable(unit, "VOsc", tableSize, inNumSamples)) return;

			const float *table0  = bufs[0].data;
			const float *table2  = bufs[1].data;
			if (!table0 || !table2 || tableSize != bufs[0].samples|| tableSize != bufs[1].samples) {
				ClearUnitOutputs(unit, inNumSamples);
				return;
			}

			const float *table1 = table0 + 1;
			const float *table3 = table2 + 1;

			LOOP(nsmps,
				float pfrac = PhaseFrac1(phase);
				uint32 index = ((phase >> xlobits1) & lomask);
				float val0 = *(float*)((char*)table0 + index);
				float val1 = *(float*)((char*)table1 + index);
				float val2 = *(float*)((char*)table2 + index);
				float val3 = *(float*)((char*)table3 + index);
				float a = val0 + val1 * pfrac;
				float b = val2 + val3 * pfrac;
				 printf("[VOsc_next_ikk2] %f\n", a + level * (b - a));//mtm
				ZXP(out) = a + level * (b - a);
				phase += phaseinc;
				level += slope;
			);
			donesmps += nsmps;
			remain -= nsmps;
			cur = cut;
		}
	}
	unit->m_bufpos = nextbufpos;
	unit->m_phase = phase;
}

void VOsc_next_ika(VOsc *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float nextbufpos = ZIN0(0);
	float freqin = ZIN0(1);
	float *phasein = ZIN(2);

	float prevbufpos = unit->m_bufpos;
	float bufdiff = nextbufpos - prevbufpos;

	int32 phase = unit->m_phase;
	int32 lomask = unit->m_lomask;

	int32 freq = (int32)(unit->m_cpstoinc * freqin);
	int32 phaseinc = freq;
	int tableSize = unit->mTableSize;
	float cur = prevbufpos;
	World *world = unit->mWorld;

	if (bufdiff == 0.f) {
		float level = cur - sc_floor(cur);
		int32 bufnum = (int)sc_floor(cur);

		const SndBuf *bufs = VOscGetBuf(bufnum, world, unit);
		if (!verify_wavetable(unit, "VOsc", tableSize, inNumSamples)) return;

		const float *table0  = bufs[0].data;
		const float *table2  = bufs[1].data;
		if (!table0 || !table2 || tableSize != bufs[0].samples|| tableSize != bufs[1].samples) {
			ClearUnitOutputs(unit, inNumSamples);
			return;
		}

		const float *table1 = table0 + 1;
		const float *table3 = table2 + 1;

		LOOP1(inNumSamples,
			int32 pphase = phase + (int32)(ZXP(phasein) * unit->m_radtoinc);
			float pfrac = PhaseFrac1(pphase);
			uint32 index = ((pphase >> xlobits1) & lomask);
			float val0 = *(float*)((char*)table0 + index);
			float val1 = *(float*)((char*)table1 + index);
			float val2 = *(float*)((char*)table2 + index);
			float val3 = *(float*)((char*)table3 + index);
			float a = val0 + val1 * pfrac;
			float b = val2 + val3 * pfrac;
			  printf("[VOsc_next_ika] %f\n", a + level * (b - a));//mtm
			ZXP(out) = a + level * (b - a);
			phase += phaseinc;
		);
	} else {
		int nsmps;
		int donesmps = 0;
		int remain = inNumSamples;
		while (remain) {
			float level = cur - sc_floor(cur);

			float cut;
			if (bufdiff > 0.) {
				cut = sc_min(nextbufpos, sc_floor(cur+1.f));
			} else {
				cut = sc_max(nextbufpos, sc_ceil(cur-1.f));
			}

			float sweepdiff = cut - cur;
			if (cut == nextbufpos) nsmps = remain;
			else {
				float sweep = (float)inNumSamples / bufdiff;
				nsmps = (int)sc_floor(sweep * sweepdiff + 0.5f) - donesmps;
				nsmps = sc_clip(nsmps, 1, remain);
			}

			float slope = sweepdiff / (float)nsmps;

			int32 bufnum = (int32)sc_floor(cur);
			const SndBuf *bufs = VOscGetBuf(bufnum, world, unit);
			if (!verify_wavetable(unit, "VOsc", tableSize, inNumSamples)) return;

			const float *table0  = bufs[0].data;
			const float *table2  = bufs[1].data;
			if (!table0 || !table2 || tableSize != bufs[0].samples|| tableSize != bufs[1].samples) {
				ClearUnitOutputs(unit, inNumSamples);
				return;
			}

			const float *table1 = table0 + 1;
			const float *table3 = table2 + 1;

			LOOP(nsmps,
				int32 pphase = phase + (int32)(ZXP(phasein) * unit->m_radtoinc);
				float pfrac = PhaseFrac1(pphase);
				uint32 index = ((pphase >> xlobits1) & lomask);
				float val0 = *(float*)((char*)table0 + index);
				float val1 = *(float*)((char*)table1 + index);
				float val2 = *(float*)((char*)table2 + index);
				float val3 = *(float*)((char*)table3 + index);
				float a = val0 + val1 * pfrac;
				float b = val2 + val3 * pfrac;
				 printf("[VOsc_next_ika2] %f\n", a + level * (b - a));//mtm
				ZXP(out) = a + level * (b - a);
				phase += phaseinc;
				level += slope;
			);
			donesmps += nsmps;
			remain -= nsmps;
			cur = cut;
		}
	}
	unit->m_bufpos = nextbufpos;
	unit->m_phase = phase;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

void VOsc3_Ctor(VOsc3 *unit)
{
	SETCALC(VOsc3_next_ik);

	float nextbufpos = ZIN0(0);
	unit->m_bufpos = nextbufpos;
	int32 bufnum = (int32)sc_floor(nextbufpos);
	World *world = unit->mWorld;

	const SndBuf *bufs = VOscGetBuf(bufnum, world, unit);
	int tableSize = bufs[0].samples;

	unit->mTableSize = tableSize;
	int tableSize2 = tableSize >> 1;
	unit->m_lomask = (tableSize2 - 1) << 3;
	unit->m_cpstoinc = tableSize2 * SAMPLEDUR * 65536.;

	unit->m_phase1 = 0;
	unit->m_phase2 = 0;
	unit->m_phase3 = 0;

	printf("[VOsc3] init sample:\n\t");//mtm
	VOsc3_next_ik(unit, 1);
	printf("[VOsc3] first sample:\n\t");//mtm
	
	unit->m_phase1 = 0;//mtm
	unit->m_phase2 = 0;//mtm
	unit->m_phase3 = 0;//mtm
}

void VOsc3_next_ik(VOsc3 *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float nextbufpos = ZIN0(0);
	float freq1in = ZIN0(1);
	float freq2in = ZIN0(2);
	float freq3in = ZIN0(3);

	float prevbufpos = unit->m_bufpos;
	float bufdiff = nextbufpos - prevbufpos;

	int32 phase1 = unit->m_phase1;
	int32 phase2 = unit->m_phase2;
	int32 phase3 = unit->m_phase3;

	int32 freq1 = (int32)(unit->m_cpstoinc * freq1in);
	int32 freq2 = (int32)(unit->m_cpstoinc * freq2in);
	int32 freq3 = (int32)(unit->m_cpstoinc * freq3in);

	int32 lomask = unit->m_lomask;
	int tableSize = unit->mTableSize;
	float cur = prevbufpos;
	World *world = unit->mWorld;

	if (bufdiff == 0.f) {
		float level = cur - (int)cur;

		int bufnum = (int)cur;

		const SndBuf *bufs = VOscGetBuf(bufnum, world, unit);
		if (!verify_wavetable(unit, "VOsc3", tableSize, inNumSamples)) return;

		const float *table0  = bufs[0].data;
		const float *table2  = bufs[1].data;
		if (!table0 || !table2 || tableSize != bufs[0].samples|| tableSize != bufs[1].samples) {
			ClearUnitOutputs(unit, inNumSamples);
			return;
		}

		const float *table1 = table0 + 1;
		const float *table3 = table2 + 1;

		LOOP1(inNumSamples,

			float pfrac1 = PhaseFrac1(phase1);
			float pfrac2 = PhaseFrac1(phase2);
			float pfrac3 = PhaseFrac1(phase3);

			int index1 = ((phase1 >> xlobits1) & lomask);
			int index2 = ((phase2 >> xlobits1) & lomask);
			int index3 = ((phase3 >> xlobits1) & lomask);

			phase1 += freq1;
			phase2 += freq2;
			phase3 += freq3;

			float val10 = *(float*)((char*)table0 + index1);
			float val11 = *(float*)((char*)table1 + index1);
			float val12 = *(float*)((char*)table2 + index1);
			float val13 = *(float*)((char*)table3 + index1);
			float a = val10 + val11 * pfrac1;
			float b = val12 + val13 * pfrac1;

			float val20 = *(float*)((char*)table0 + index2);
			float val21 = *(float*)((char*)table1 + index2);
			float val22 = *(float*)((char*)table2 + index2);
			float val23 = *(float*)((char*)table3 + index2);
			a += val20 + val21 * pfrac2;
			b += val22 + val23 * pfrac2;

			float val30 = *(float*)((char*)table0 + index3);
			float val31 = *(float*)((char*)table1 + index3);
			float val32 = *(float*)((char*)table2 + index3);
			float val33 = *(float*)((char*)table3 + index3);
			a += val30 + val31 * pfrac3;
			b += val32 + val33 * pfrac3;
			  printf("[VOsc3_next_ik1] %f\n", a + level * (b - a));//mtm
			ZXP(out) = a + level * (b - a);
		);
	} else {
		int nsmps;
		int donesmps = 0;
		int remain = inNumSamples;
		do {
			float level = cur - sc_trunc(cur);

			float cut;
			if (bufdiff >= 0.)
				cut = sc_min(nextbufpos, sc_trunc(cur+1.f));
			else
				cut = sc_max(nextbufpos, sc_ceil(cur-1.f));

			float sweepdiff = cut - cur;
			if (cut == nextbufpos) nsmps = remain;
			else {
				float sweep = (float)inNumSamples / bufdiff;
				nsmps = sc_floor(sweep * sweepdiff + 0.5f) - donesmps;
				nsmps = sc_clip(nsmps, 1, remain);
			}

			float slope = sweepdiff / (float)nsmps;

			int bufnum = (int)cur;

			const SndBuf *bufs = VOscGetBuf(bufnum, world, unit);
			if (!verify_wavetable(unit, "VOsc3", tableSize, inNumSamples)) return;

			const float *table0  = bufs[0].data;
			const float *table2  = bufs[1].data;
			if (!table0 || !table2 || tableSize != bufs[0].samples|| tableSize != bufs[1].samples) {
				ClearUnitOutputs(unit, inNumSamples);
				return;
			}

			const float *table1 = table0 + 1;
			const float *table3 = table2 + 1;

			LOOP(nsmps,

				float pfrac1 = PhaseFrac1(phase1);
				float pfrac2 = PhaseFrac1(phase2);
				float pfrac3 = PhaseFrac1(phase3);

				int index1 = ((phase1 >> xlobits1) & lomask);
				int index2 = ((phase2 >> xlobits1) & lomask);
				int index3 = ((phase3 >> xlobits1) & lomask);

				phase1 += freq1;
				phase2 += freq2;
				phase3 += freq3;

				float val10 = *(float*)((char*)table0 + index1);
				float val11 = *(float*)((char*)table1 + index1);
				float val12 = *(float*)((char*)table2 + index1);
				float val13 = *(float*)((char*)table3 + index1);
				float a = val10 + val11 * pfrac1;
				float b = val12 + val13 * pfrac1;

				float val20 = *(float*)((char*)table0 + index2);
				float val21 = *(float*)((char*)table1 + index2);
				float val22 = *(float*)((char*)table2 + index2);
				float val23 = *(float*)((char*)table3 + index2);
				a += val20 + val21 * pfrac2;
				b += val22 + val23 * pfrac2;

				float val30 = *(float*)((char*)table0 + index3);
				float val31 = *(float*)((char*)table1 + index3);
				float val32 = *(float*)((char*)table2 + index3);
				float val33 = *(float*)((char*)table3 + index3);
				a += val30 + val31 * pfrac3;
				b += val32 + val33 * pfrac3;
				   printf("[VOsc3_next_ik2] %f\n", a + level * (b - a));//mtm
				ZXP(out) = a + level * (b - a);
				level += slope;
			);
			donesmps += nsmps;
			remain -= nsmps;
			cur = cut;
		} while (remain);
	}
	unit->m_bufpos = nextbufpos;
	unit->m_phase1 = phase1;
	unit->m_phase2 = phase2;
	unit->m_phase3 = phase3;
}

//////////////////////////////////////////////////////////////////////////////////////////

void Formant_Ctor(Formant *unit)
{
	SETCALC(Formant_next);
	unit->m_cpstoinc = ft->mSineSize * SAMPLEDUR * 65536.;
	unit->m_phase1 = 0;
	unit->m_phase2 = 0;
	unit->m_phase3 = 0;
	printf("[Formant] init sample:\n\t");//mtm
	Formant_next(unit, 1);
	printf("[Formant] first sample:\n\t");//mtm
	unit->m_phase1 = 0;//mtm
	unit->m_phase2 = 0;//mtm
	unit->m_phase3 = 0;//mtm
}

#define tqcyc13 0x18000000

void Formant_next(Formant* unit, int inNumSamples) {
    float* out = ZOUT(0);
    float freq1in = ZIN0(0);
    float freq2in = ZIN0(1);
    float freq3in = ZIN0(2);

    int32 phase1 = unit->m_phase1;
    int32 phase2 = unit->m_phase2;
    int32 phase3 = unit->m_phase3;
    float cpstoinc = unit->m_cpstoinc;
    int32 freq1 = (int32)(cpstoinc * freq1in);
    int32 freq2 = (int32)(cpstoinc * freq2in);
    int32 freq3 = (int32)(cpstoinc * freq3in);
    float* sine = ft->mSine;
    int32 formfreq = sc_max(freq1, freq3);
    LOOP1(
        inNumSamples,
        if (phase3 < onecyc13) {
            ZXP(out) = (*(float*)((char*)sine + (((phase3 + tqcyc13) >> xlobits) & xlomask13)) + 1.f)
                * *(float*)((char*)sine + ((phase2 >> xlobits) & xlomask13));
            phase3 += formfreq;
        } else { ZXP(out) = 0.f; } phase1 += freq1;
        phase2 += freq2; if (phase1 > onecyc13) {
            phase1 -= onecyc13;
            phase2 = phase1 * freq2 / freq1;
            phase3 = phase1 * freq3 / freq1;
        });

    unit->m_phase1 = phase1;
    unit->m_phase2 = phase2;
    unit->m_phase3 = phase3;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

inline float lookup13(float* table, int32 pphase) {
    float pfrac = PhaseFrac(pphase);
    float* tbl = (float*)((char*)table + ((pphase >> xlobits) & xlomask13));
    return lininterp(pfrac, tbl[0], tbl[1]);
}

void Blip_Ctor(Blip* unit) {
    SETCALC(Blip_next);
    unit->m_freqin = ZIN0(0);
    unit->m_numharm = (int32)ZIN0(1);

    unit->m_cpstoinc = ft->mSineSize * SAMPLEDUR * 65536. * 0.5;
    int32 N = unit->m_numharm;
    int32 maxN = (int32)((SAMPLERATE * 0.5) / unit->m_freqin);
    if (N > maxN)
        N = maxN;
    if (N < 1)
        N = 1;
    unit->m_N = N;
    unit->m_scale = 0.5 / N;
    unit->m_phase = 0;

    Blip_next(unit, 1);
}

void Blip_next(Blip* unit, int inNumSamples) {
    float* out = ZOUT(0);
    float freqin = ZIN0(0);
    int numharm = (int32)ZIN0(1);

    int32 phase = unit->m_phase;

    float* numtbl = ft->mSine;
    float* dentbl = ft->mCosecant;

    int32 freq, N, prevN;
    float scale, prevscale;
    bool crossfade;
    if (numharm != unit->m_numharm || freqin != unit->m_freqin) {
        N = numharm;
        int32 maxN = (int32)((SAMPLERATE * 0.5) / freqin);
        if (N > maxN) {
            float maxfreqin;
            N = maxN;
            maxfreqin = sc_max(unit->m_freqin, freqin);
            freq = (int32)(unit->m_cpstoinc * maxfreqin);
        } else {
            if (N < 1) {
                N = 1;
            }
            freq = (int32)(unit->m_cpstoinc * freqin);
        }
        crossfade = N != unit->m_N;
        prevN = unit->m_N;
        prevscale = unit->m_scale;
        unit->m_N = N;
        unit->m_scale = scale = 0.5 / N;
    } else {
        N = unit->m_N;
        freq = (int32)(unit->m_cpstoinc * freqin);
        scale = unit->m_scale;
        crossfade = false;
    }
    int32 N2 = 2 * N + 1;

    if (crossfade) {
        int32 prevN2 = 2 * prevN + 1;
        float xfade_slope = unit->mRate->mSlopeFactor;
        float xfade = 0.f;
        LOOP1(
            inNumSamples, float* tbl = (float*)((char*)dentbl + ((phase >> xlobits) & xlomask13)); float t0 = tbl[0];
            float t1 = tbl[1]; if (t0 == kBadValue || t1 == kBadValue) {
                tbl = (float*)((char*)numtbl + ((phase >> xlobits) & xlomask13));
                t0 = tbl[0];
                t1 = tbl[1];
                float pfrac = PhaseFrac(phase);
                float denom = t0 + (t1 - t0) * pfrac;
                if (std::abs(denom) < 0.0005f) {
                    ZXP(out) = 1.f;
                } else {
                    int32 rphase = phase * prevN2;
                    pfrac = PhaseFrac(rphase);
                    tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                    float numer = lininterp(pfrac, tbl[0], tbl[1]);
                    float n1 = (numer / denom - 1.f) * prevscale;

                    rphase = phase * N2;
                    pfrac = PhaseFrac(rphase);
                    tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                    numer = lininterp(pfrac, tbl[0], tbl[1]);
                    float n2 = (numer / denom - 1.f) * scale;

                    ZXP(out) = lininterp(xfade, n1, n2);
                }
            } else {
                float pfrac = PhaseFrac(phase);
                float denom = t0 + (t1 - t0) * pfrac;

                int32 rphase = phase * prevN2;
                pfrac = PhaseFrac(rphase);
                float* tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                float numer = lininterp(pfrac, tbl[0], tbl[1]);
                float n1 = (numer * denom - 1.f) * prevscale;

                rphase = phase * N2;
                pfrac = PhaseFrac(rphase);
                tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                numer = lininterp(pfrac, tbl[0], tbl[1]);
                float n2 = (numer * denom - 1.f) * scale;

                ZXP(out) = lininterp(xfade, n1, n2);
            } phase += freq;
            xfade += xfade_slope;);
    } else {
        // hmm, if freq is above sr/4 then revert to sine table osc w/ no interpolation ?
        // why bother, it isn't a common choice for a fundamental.
        LOOP1(
            inNumSamples, float* tbl = (float*)((char*)dentbl + ((phase >> xlobits) & xlomask13)); float t0 = tbl[0];
            float t1 = tbl[1]; if (t0 == kBadValue || t1 == kBadValue) {
                tbl = (float*)((char*)numtbl + ((phase >> xlobits) & xlomask13));
                t0 = tbl[0];
                t1 = tbl[1];
                float pfrac = PhaseFrac(phase);
                float denom = t0 + (t1 - t0) * pfrac;
                if (std::abs(denom) < 0.0005f) {
                    ZXP(out) = 1.f;
                } else {
                    int32 rphase = phase * N2;
                    pfrac = PhaseFrac(rphase);
                    tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                    float numer = lininterp(pfrac, tbl[0], tbl[1]);
                    ZXP(out) = (numer / denom - 1.f) * scale;
                }
            } else {
                float pfrac = PhaseFrac(phase);
                float denom = t0 + (t1 - t0) * pfrac;
                int32 rphase = phase * N2;
                pfrac = PhaseFrac(rphase);
                tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                float numer = lininterp(pfrac, tbl[0], tbl[1]);
                ZXP(out) = (numer * denom - 1.f) * scale;
            } phase += freq;);
    }

    unit->m_phase = phase;
    unit->m_freqin = freqin;
    unit->m_numharm = numharm;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////


void Saw_Ctor(Saw* unit) {
    SETCALC(Saw_next);
    unit->m_freqin = ZIN0(0);

    unit->m_cpstoinc = ft->mSineSize * SAMPLEDUR * 65536. * 0.5;
    unit->m_N = (int32)((SAMPLERATE * 0.5) / unit->m_freqin);
    unit->m_scale = 0.5 / unit->m_N;
    unit->m_phase = 0;
    unit->m_y1 = -0.46f;

    ZOUT0(0) = 0.f;
}

void Saw_next(Saw* unit, int inNumSamples) {
    float* out = ZOUT(0);
    float freqin = ZIN0(0);

    int32 phase = unit->m_phase;
    float y1 = unit->m_y1;

    float* numtbl = ft->mSine;
    float* dentbl = ft->mCosecant;

    int32 freq, N, prevN;
    float scale, prevscale;
    bool crossfade;
    if (freqin != unit->m_freqin) {
        N = (int32)((SAMPLERATE * 0.5) / freqin);
        if (N != unit->m_N) {
            float maxfreqin;
            maxfreqin = sc_max(unit->m_freqin, freqin);
            freq = (int32)(unit->m_cpstoinc * maxfreqin);
            crossfade = true;
        } else {
            freq = (int32)(unit->m_cpstoinc * freqin);
            crossfade = false;
        }
        prevN = unit->m_N;
        prevscale = unit->m_scale;
        unit->m_N = N;
        unit->m_scale = scale = 0.5 / N;
    } else {
        N = unit->m_N;
        freq = (int32)(unit->m_cpstoinc * freqin);
        scale = unit->m_scale;
        crossfade = false;
    }
    int32 N2 = 2 * N + 1;

    if (crossfade) {
        int32 prevN2 = 2 * prevN + 1;
        float xfade_slope = unit->mRate->mSlopeFactor;
        float xfade = 0.f;
        LOOP1(
            inNumSamples, float* tbl = (float*)((char*)dentbl + ((phase >> xlobits) & xlomask13)); float t0 = tbl[0];
            float t1 = tbl[1]; if (t0 == kBadValue || t1 == kBadValue) {
                tbl = (float*)((char*)numtbl + ((phase >> xlobits) & xlomask13));
                t0 = tbl[0];
                t1 = tbl[1];
                float pfrac = PhaseFrac(phase);
                float denom = t0 + (t1 - t0) * pfrac;
                if (std::abs(denom) < 0.0005f) {
                    ZXP(out) = y1 = 1.f + 0.999f * y1;
                } else {
                    int32 rphase = phase * prevN2;
                    pfrac = PhaseFrac(rphase);
                    tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                    float numer = lininterp(pfrac, tbl[0], tbl[1]);
                    float n1 = (numer / denom - 1.f) * prevscale;

                    rphase = phase * N2;
                    pfrac = PhaseFrac(rphase);
                    tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                    numer = lininterp(pfrac, tbl[0], tbl[1]);
                    float n2 = (numer / denom - 1.f) * scale;

                    ZXP(out) = y1 = n1 + xfade * (n2 - n1) + 0.999f * y1;
                }
            } else {
                float pfrac = PhaseFrac(phase);
                float denom = t0 + (t1 - t0) * pfrac;

                int32 rphase = phase * prevN2;
                pfrac = PhaseFrac(rphase);
                tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                float numer = lininterp(pfrac, tbl[0], tbl[1]);
                float n1 = (numer * denom - 1.f) * prevscale;

                rphase = phase * N2;
                pfrac = PhaseFrac(rphase);
                tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                numer = lininterp(pfrac, tbl[0], tbl[1]);
                float n2 = (numer * denom - 1.f) * scale;

                ZXP(out) = y1 = n1 + xfade * (n2 - n1) + 0.999f * y1;
            } phase += freq;
            xfade += xfade_slope;);
    } else {
        // hmm, if freq is above sr/4 then revert to sine table osc ?
        // why bother, it isn't a common choice for a fundamental.
        LOOP1(
            inNumSamples, float* tbl = (float*)((char*)dentbl + ((phase >> xlobits) & xlomask13)); float t0 = tbl[0];
            float t1 = tbl[1]; if (t0 == kBadValue || t1 == kBadValue) {
                tbl = (float*)((char*)numtbl + ((phase >> xlobits) & xlomask13));
                t0 = tbl[0];
                t1 = tbl[1];
                float pfrac = PhaseFrac(phase);
                float denom = t0 + (t1 - t0) * pfrac;
                if (std::abs(denom) < 0.0005f) {
                    ZXP(out) = y1 = 1.f + 0.999f * y1;
                } else {
                    int32 rphase = phase * N2;
                    pfrac = PhaseFrac(rphase);
                    tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                    float numer = lininterp(pfrac, tbl[0], tbl[1]);
                    ZXP(out) = y1 = (numer / denom - 1.f) * scale + 0.999f * y1;
                }
            } else {
                float pfrac = PhaseFrac(phase);
                float denom = t0 + (t1 - t0) * pfrac;
                int32 rphase = phase * N2;
                pfrac = PhaseFrac(rphase);
                float* tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                float numer = lininterp(pfrac, tbl[0], tbl[1]);
                ZXP(out) = y1 = (numer * denom - 1.f) * scale + 0.999f * y1;
            } phase += freq;);
    }

    unit->m_y1 = y1;
    unit->m_phase = phase;
    unit->m_freqin = freqin;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////


void Pulse_Ctor(Pulse* unit) {
    SETCALC(Pulse_next);
    unit->m_freqin = ZIN0(0);

    unit->m_cpstoinc = ft->mSineSize * SAMPLEDUR * 65536. * 0.5;
    unit->m_N = (int32)((SAMPLERATE * 0.5) / unit->m_freqin);
    unit->m_scale = 0.5 / unit->m_N;
    unit->m_phase = 0;
    unit->m_phaseoff = 0;
    unit->m_y1 = 0.f;
    ZOUT0(0) = 0.f;
}

void Pulse_next(Pulse* unit, int inNumSamples) {
    float* out = ZOUT(0);
    float freqin = ZIN0(0);
    float duty = ZIN0(1);

    int32 phase = unit->m_phase;
    float y1 = unit->m_y1;

    float* numtbl = ft->mSine;
    float* dentbl = ft->mCosecant;

    int32 freq, N, prevN;
    float scale, prevscale;
    bool crossfade;

    if (freqin != unit->m_freqin) {
        N = (int32)((SAMPLERATE * 0.5) / freqin);
        if (N != unit->m_N) {
            float maxfreqin;
            maxfreqin = sc_max(unit->m_freqin, freqin);
            freq = (int32)(unit->m_cpstoinc * maxfreqin);
            crossfade = true;
        } else {
            freq = (int32)(unit->m_cpstoinc * freqin);
            crossfade = false;
        }
        prevN = unit->m_N;
        prevscale = unit->m_scale;
        unit->m_N = N;
        unit->m_scale = scale = 0.5 / N;
    } else {
        N = unit->m_N;
        freq = (int32)(unit->m_cpstoinc * freqin);
        scale = unit->m_scale;
        crossfade = false;
    }
    int32 N2 = 2 * N + 1;

    int32 phaseoff = unit->m_phaseoff;
    int32 next_phaseoff = (int32)(duty * (1L << 28));
    int32 phaseoff_slope = (int32)((next_phaseoff - phaseoff) * unit->mRate->mSlopeFactor);
    unit->m_phaseoff = next_phaseoff;
    float rscale = 1.f / scale + 1.f;
    float pul1, pul2;

    if (crossfade) {
        int32 prevN2 = 2 * prevN + 1;
        float xfade_slope = unit->mRate->mSlopeFactor;
        float xfade = 0.f;
        LOOP1(
            inNumSamples, float* tbl = (float*)((char*)dentbl + ((phase >> xlobits) & xlomask13)); float t0 = tbl[0];
            float t1 = tbl[1]; if (t0 == kBadValue || t1 == kBadValue) {
                tbl = (float*)((char*)numtbl + ((phase >> xlobits) & xlomask13));
                t0 = tbl[0];
                t1 = tbl[1];
                float pfrac = PhaseFrac(phase);
                float denom = t0 + (t1 - t0) * pfrac;
                if (std::abs(denom) < 0.0005f) {
                    pul1 = 1.f;
                } else {
                    int32 rphase = phase * prevN2;
                    pfrac = PhaseFrac(rphase);
                    tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                    float numer = lininterp(pfrac, tbl[0], tbl[1]);
                    float n1 = (numer / denom - 1.f) * prevscale;

                    rphase = phase * N2;
                    pfrac = PhaseFrac(rphase);
                    tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                    numer = lininterp(pfrac, tbl[0], tbl[1]);
                    float n2 = (numer / denom - 1.f) * scale;

                    pul1 = lininterp(xfade, n1, n2);
                }
            } else {
                float pfrac = PhaseFrac(phase);
                float denom = lininterp(pfrac, t0, t1);

                int32 rphase = phase * prevN2;
                pfrac = PhaseFrac(rphase);
                tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                float numer = lininterp(pfrac, tbl[0], tbl[1]);
                float n1 = (numer * denom - 1.f) * prevscale;

                rphase = phase * N2;
                pfrac = PhaseFrac(rphase);
                tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                numer = lininterp(pfrac, tbl[0], tbl[1]);
                float n2 = (numer * denom - 1.f) * scale;

                pul1 = lininterp(xfade, n1, n2);
            }

            int32 phase2 = phase + phaseoff;
            tbl = (float*)((char*)dentbl + ((phase2 >> xlobits) & xlomask13)); t0 = tbl[0]; t1 = tbl[1];
            if (t0 == kBadValue || t1 == kBadValue) {
                tbl = (float*)((char*)numtbl + ((phase2 >> xlobits) & xlomask13));
                t0 = tbl[0];
                t1 = tbl[1];
                float pfrac = PhaseFrac(phase2);
                float denom = t0 + (t1 - t0) * pfrac;
                if (std::abs(denom) < 0.0005f) {
                    pul2 = 1.f;
                } else {
                    int32 rphase = phase2 * prevN2;
                    pfrac = PhaseFrac(rphase);
                    tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                    float numer = lininterp(pfrac, tbl[0], tbl[1]);
                    float n1 = (numer / denom - 1.f) * prevscale;

                    rphase = phase2 * N2;
                    pfrac = PhaseFrac(rphase);
                    tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                    numer = lininterp(pfrac, tbl[0], tbl[1]);
                    float n2 = (numer / denom - 1.f) * scale;

                    pul2 = lininterp(xfade, n1, n2);
                }
            } else {
                float pfrac = PhaseFrac(phase2);
                float denom = t0 + (t1 - t0) * pfrac;

                int32 rphase = phase2 * prevN2;
                pfrac = PhaseFrac(rphase);
                tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                float numer = lininterp(pfrac, tbl[0], tbl[1]);
                float n1 = (numer * denom - 1.f) * prevscale;

                rphase = phase2 * N2;
                pfrac = PhaseFrac(rphase);
                tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                numer = lininterp(pfrac, tbl[0], tbl[1]);
                float n2 = (numer * denom - 1.f) * scale;

                pul2 = lininterp(xfade, n1, n2);
            }

            ZXP(out) = y1 = pul1 - pul2 + 0.999f * y1;
            phase += freq; phaseoff += phaseoff_slope; xfade += xfade_slope;);
    } else {
        LOOP1(
            inNumSamples, float* tbl = (float*)((char*)dentbl + ((phase >> xlobits) & xlomask13)); float t0 = tbl[0];
            float t1 = tbl[1]; if (t0 == kBadValue || t1 == kBadValue) {
                tbl = (float*)((char*)numtbl + ((phase >> xlobits) & xlomask13));
                t0 = tbl[0];
                t1 = tbl[1];
                float pfrac = PhaseFrac(phase);
                float denom = t0 + (t1 - t0) * pfrac;
                if (std::abs(denom) < 0.0005f) {
                    pul1 = rscale;
                } else {
                    int32 rphase = phase * N2;
                    pfrac = PhaseFrac(rphase);
                    tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                    float numer = lininterp(pfrac, tbl[0], tbl[1]);
                    pul1 = numer / denom;
                }
            } else {
                float pfrac = PhaseFrac(phase);
                float denom = t0 + (t1 - t0) * pfrac;
                int32 rphase = phase * N2;
                pfrac = PhaseFrac(rphase);
                float* tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                float numer = lininterp(pfrac, tbl[0], tbl[1]);
                pul1 = (numer * denom);
            }

            int32 phase2 = phase + phaseoff;
            tbl = (float*)((char*)dentbl + ((phase2 >> xlobits) & xlomask13)); t0 = tbl[0]; t1 = tbl[1];
            if (t0 == kBadValue || t1 == kBadValue) {
                tbl = (float*)((char*)numtbl + ((phase2 >> xlobits) & xlomask13));
                t0 = tbl[0];
                t1 = tbl[1];
                float pfrac = PhaseFrac(phase2);
                float denom = t0 + (t1 - t0) * pfrac;
                if (std::abs(denom) < 0.0005f) {
                    pul2 = rscale;
                } else {
                    int32 rphase = phase2 * N2;
                    pfrac = PhaseFrac(rphase);
                    tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                    float numer = lininterp(pfrac, tbl[0], tbl[1]);
                    pul2 = numer / denom;
                }
            } else {
                float pfrac = PhaseFrac(phase2);
                float denom = t0 + (t1 - t0) * pfrac;
                int32 rphase = phase2 * N2;
                pfrac = PhaseFrac(rphase);
                float* tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
                float numer = lininterp(pfrac, tbl[0], tbl[1]);

                pul2 = (numer * denom);
            }

            ZXP(out) = y1 = (pul1 - pul2) * scale + 0.999f * y1;
            phase += freq; phaseoff += phaseoff_slope;);
    }

    unit->m_y1 = y1;
    unit->m_phase = phase;
    unit->m_freqin = freqin;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

static float Klang_SetCoefs(Klang* unit) {
    float freqscale = ZIN0(0) * unit->mRate->mRadiansPerSample;
    float freqoffset = ZIN0(1) * unit->mRate->mRadiansPerSample;

    float outf = 0.;
    float* coefs = unit->m_coefs - 1;

    for (int i = 0, j = 2; i < unit->m_numpartials; ++i, j += 3) {
        float w = ZIN0(j) * freqscale + freqoffset;
        float level = ZIN0(j + 1);
        float phase = ZIN0(j + 2);

        if (phase != 0.f) {
            outf += * ++coefs = level * sin(phase); // y1
            *++coefs = level * sin(phase - w); // y2
        } else {
            outf += * ++coefs = 0.f; // y1
            *++coefs = level * -sin(w); // y2
        }
        *++coefs = 2. * cos(w); // b1
    }
    return outf;
}

void Klang_Ctor(Klang* unit) {
    SETCALC(Klang_next);
    unit->m_numpartials = (unit->mNumInputs - 2) / 3;
    int numcoefs = unit->m_numpartials * 3;
    unit->m_coefs = (float*)RTAlloc(unit->mWorld, numcoefs * sizeof(float));
    ClearUnitIfMemFailed(unit->m_coefs);
    ZOUT0(0) = Klang_SetCoefs(unit);
}

void Klang_Dtor(Klang* unit) { RTFree(unit->mWorld, unit->m_coefs); }

void Klang_next(Klang* unit, int inNumSamples) {
    float* out0 = ZOUT(0);

    float* out;
    float y0_0, y1_0, y2_0, b1_0;
    float y0_1, y1_1, y2_1, b1_1;
    float y0_2, y1_2, y2_2, b1_2;
    float y0_3, y1_3, y2_3, b1_3;
    float outf;

    float* coefs = unit->m_coefs - 1;
    int32 numpartials = unit->m_numpartials;

    switch (numpartials & 3) {
    case 3:
        y1_0 = *++coefs;
        y2_0 = *++coefs;
        b1_0 = *++coefs;
        y1_1 = *++coefs;
        y2_1 = *++coefs;
        b1_1 = *++coefs;
        y1_2 = *++coefs;
        y2_2 = *++coefs;
        b1_2 = *++coefs;

        out = out0;
        LOOP(unit->mRate->mFilterLoops, outf = y0_0 = b1_0 * y1_0 - y2_0; outf += y0_1 = b1_1 * y1_1 - y2_1;
             outf += y0_2 = b1_2 * y1_2 - y2_2; ZXP(out) = outf;

             outf = y2_0 = b1_0 * y0_0 - y1_0; outf += y2_1 = b1_1 * y0_1 - y1_1; outf += y2_2 = b1_2 * y0_2 - y1_2;
             ZXP(out) = outf;

             outf = y1_0 = b1_0 * y2_0 - y0_0; outf += y1_1 = b1_1 * y2_1 - y0_1; outf += y1_2 = b1_2 * y2_2 - y0_2;
             ZXP(out) = outf;);
        LOOP(unit->mRate->mFilterRemain, outf = y0_0 = b1_0 * y1_0 - y2_0; outf += y0_1 = b1_1 * y1_1 - y2_1;
             outf += y0_2 = b1_2 * y1_2 - y2_2; y2_0 = y1_0; y1_0 = y0_0; y2_1 = y1_1; y1_1 = y0_1; y2_2 = y1_2;
             y1_2 = y0_2; ZXP(out) = outf;);
        coefs -= 9;
        *++coefs = y1_0;
        *++coefs = y2_0;
        ++coefs;
        *++coefs = y1_1;
        *++coefs = y2_1;
        ++coefs;
        *++coefs = y1_2;
        *++coefs = y2_2;
        ++coefs;
        break;
    case 2:
        y1_0 = *++coefs;
        y2_0 = *++coefs;
        b1_0 = *++coefs;
        y1_1 = *++coefs;
        y2_1 = *++coefs;
        b1_1 = *++coefs;

        out = out0;
        LOOP(unit->mRate->mFilterLoops, outf = y0_0 = b1_0 * y1_0 - y2_0; outf += y0_1 = b1_1 * y1_1 - y2_1;
             ZXP(out) = outf;

             outf = y2_0 = b1_0 * y0_0 - y1_0; outf += y2_1 = b1_1 * y0_1 - y1_1; ZXP(out) = outf;

             outf = y1_0 = b1_0 * y2_0 - y0_0; outf += y1_1 = b1_1 * y2_1 - y0_1; ZXP(out) = outf;);
        LOOP(unit->mRate->mFilterRemain, outf = y0_0 = b1_0 * y1_0 - y2_0; outf += y0_1 = b1_1 * y1_1 - y2_1;
             y2_0 = y1_0; y1_0 = y0_0; y2_1 = y1_1; y1_1 = y0_1; ZXP(out) = outf;);

        coefs -= 6;
        *++coefs = y1_0;
        *++coefs = y2_0;
        ++coefs;
        *++coefs = y1_1;
        *++coefs = y2_1;
        ++coefs;
        break;
    case 1:
        y1_0 = *++coefs;
        y2_0 = *++coefs;
        b1_0 = *++coefs;

        out = out0;
        LOOP(unit->mRate->mFilterLoops, ZXP(out) = y0_0 = b1_0 * y1_0 - y2_0;

             ZXP(out) = y2_0 = b1_0 * y0_0 - y1_0;

             ZXP(out) = y1_0 = b1_0 * y2_0 - y0_0;);
        LOOP(unit->mRate->mFilterRemain, ZXP(out) = y0_0 = b1_0 * y1_0 - y2_0; y2_0 = y1_0; y1_0 = y0_0;);

        coefs -= 3;
        *++coefs = y1_0;
        *++coefs = y2_0;
        ++coefs;
        break;
    case 0:
        out = out0;
        ZClear(inNumSamples, out);
        break;
    }

    int32 imax = numpartials >> 2;

    for (int i = 0; i < imax; ++i) {
        y1_0 = *++coefs;
        y2_0 = *++coefs;
        b1_0 = *++coefs;
        y1_1 = *++coefs;
        y2_1 = *++coefs;
        b1_1 = *++coefs;
        y1_2 = *++coefs;
        y2_2 = *++coefs;
        b1_2 = *++coefs;
        y1_3 = *++coefs;
        y2_3 = *++coefs;
        b1_3 = *++coefs;

        out = out0;
        LOOP(unit->mRate->mFilterLoops, outf = y0_0 = b1_0 * y1_0 - y2_0; outf += y0_1 = b1_1 * y1_1 - y2_1;
             outf += y0_2 = b1_2 * y1_2 - y2_2; outf += y0_3 = b1_3 * y1_3 - y2_3; ZXP(out) += outf;

             outf = y2_0 = b1_0 * y0_0 - y1_0; outf += y2_1 = b1_1 * y0_1 - y1_1; outf += y2_2 = b1_2 * y0_2 - y1_2;
             outf += y2_3 = b1_3 * y0_3 - y1_3; ZXP(out) += outf;

             outf = y1_0 = b1_0 * y2_0 - y0_0; outf += y1_1 = b1_1 * y2_1 - y0_1; outf += y1_2 = b1_2 * y2_2 - y0_2;
             outf += y1_3 = b1_3 * y2_3 - y0_3; ZXP(out) += outf;);
        LOOP(unit->mRate->mFilterRemain, outf = y0_0 = b1_0 * y1_0 - y2_0; outf += y0_1 = b1_1 * y1_1 - y2_1;
             outf += y0_2 = b1_2 * y1_2 - y2_2; outf += y0_3 = b1_3 * y1_3 - y2_3; y2_0 = y1_0; y1_0 = y0_0;
             y2_1 = y1_1; y1_1 = y0_1; y2_2 = y1_2; y1_2 = y0_2; y2_3 = y1_3; y1_3 = y0_3; ZXP(out) += outf;);
        coefs -= 12;
        *++coefs = y1_0;
        *++coefs = y2_0;
        ++coefs;
        *++coefs = y1_1;
        *++coefs = y2_1;
        ++coefs;
        *++coefs = y1_2;
        *++coefs = y2_2;
        ++coefs;
        *++coefs = y1_3;
        *++coefs = y2_3;
        ++coefs;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

static void Klank_SetCoefs(Klank* unit) {
    int numpartials = (unit->mNumInputs - 4) / 3;
    unit->m_numpartials = numpartials;

    int numcoefs = ((unit->m_numpartials + 3) & ~3) * 5;
    unit->m_coefs = (float*)RTAlloc(unit->mWorld, (numcoefs + unit->mWorld->mBufLength) * sizeof(float));
    ClearUnitIfMemFailed(unit->m_coefs);

    unit->m_buf = unit->m_coefs + numcoefs;

    float freqscale = ZIN0(1) * unit->mRate->mRadiansPerSample;
    float freqoffset = ZIN0(2) * unit->mRate->mRadiansPerSample;
    float decayscale = ZIN0(3);

    float* coefs = unit->m_coefs;

    float sampleRate = SAMPLERATE;

    for (int i = 0, j = 4; i < numpartials; ++i, j += 3) {
        float w = ZIN0(j) * freqscale + freqoffset;
        float level = ZIN0(j + 1);
        float time = ZIN0(j + 2) * decayscale;

        float R = time == 0.f ? 0.f : exp(log001 / (time * sampleRate));
        float twoR = 2.f * R;
        float R2 = R * R;
        float cost = (twoR * cos(w)) / (1.f + R2);

        int k = 20 * (i >> 2) + (i & 3);
        coefs[k + 0] = 0.f; // y1
        coefs[k + 4] = 0.f; // y2
        coefs[k + 8] = twoR * cost; // b1
        coefs[k + 12] = -R2; // b2
        coefs[k + 16] = level * 0.25; // a0
        // Print("coefs %d  %g %g %g\n", i, twoR * cost, -R2, ampf * 0.25);
    }
}

void Klank_Ctor(Klank* unit) {
    SETCALC(Klank_next);
    unit->m_x1 = unit->m_x2 = 0.f;
    Klank_SetCoefs(unit);
    ZOUT0(0) = 0.f;
}

void Klank_Dtor(Klank* unit) { RTFree(unit->mWorld, unit->m_coefs); }

void Klank_next(Klank* unit, int inNumSamples) {
    float* out0 = ZOUT(0);
    float* in0 = ZIN(0);

    float *in, *out;
    float inf;
    float y0_0, y1_0, y2_0, a0_0, b1_0, b2_0;
    float y0_1, y1_1, y2_1, a0_1, b1_1, b2_1;
    float y0_2, y1_2, y2_2, a0_2, b1_2, b2_2;
    float y0_3, y1_3, y2_3, a0_3, b1_3, b2_3;

    int32 numpartials = unit->m_numpartials;
    int32 imax = numpartials >> 2;

    float* coefs = unit->m_coefs + imax * 20;

    switch (numpartials & 3) {
    case 3:
        y1_0 = coefs[0];
        y2_0 = coefs[4];
        b1_0 = coefs[8];
        b2_0 = coefs[12];
        a0_0 = coefs[16];
        y1_1 = coefs[1];
        y2_1 = coefs[5];
        b1_1 = coefs[9];
        b2_1 = coefs[13];
        a0_1 = coefs[17];
        y1_2 = coefs[2];
        y2_2 = coefs[6];
        b1_2 = coefs[10];
        b2_2 = coefs[14];
        a0_2 = coefs[18];

        in = in0;
        out = unit->m_buf - 1;
        LooP(unit->mRate->mFilterLoops) {
            inf = *++in;
            y0_0 = inf + b1_0 * y1_0 + b2_0 * y2_0;
            y0_1 = inf + b1_1 * y1_1 + b2_1 * y2_1;
            y0_2 = inf + b1_2 * y1_2 + b2_2 * y2_2;
            *++out = a0_0 * y0_0 + a0_1 * y0_1 + a0_2 * y0_2;

            inf = *++in;
            y2_0 = inf + b1_0 * y0_0 + b2_0 * y1_0;
            y2_1 = inf + b1_1 * y0_1 + b2_1 * y1_1;
            y2_2 = inf + b1_2 * y0_2 + b2_2 * y1_2;
            *++out = a0_0 * y2_0 + a0_1 * y2_1 + a0_2 * y2_2;

            inf = *++in;
            y1_0 = inf + b1_0 * y2_0 + b2_0 * y0_0;
            y1_1 = inf + b1_1 * y2_1 + b2_1 * y0_1;
            y1_2 = inf + b1_2 * y2_2 + b2_2 * y0_2;
            *++out = a0_0 * y1_0 + a0_1 * y1_1 + a0_2 * y1_2;
        }
        LooP(unit->mRate->mFilterRemain) {
            inf = *++in;
            y0_0 = inf + b1_0 * y1_0 + b2_0 * y2_0;
            y0_1 = inf + b1_1 * y1_1 + b2_1 * y2_1;
            y0_2 = inf + b1_2 * y1_2 + b2_2 * y2_2;
            *++out = a0_0 * y0_0 + a0_1 * y0_1 + a0_2 * y0_2;
            y2_0 = y1_0;
            y1_0 = y0_0;
            y2_1 = y1_1;
            y1_1 = y0_1;
            y2_2 = y1_2;
            y1_2 = y0_2;
        }
        coefs[0] = zapgremlins(y1_0);
        coefs[4] = zapgremlins(y2_0);
        coefs[1] = zapgremlins(y1_1);
        coefs[5] = zapgremlins(y2_1);
        coefs[2] = zapgremlins(y1_2);
        coefs[6] = zapgremlins(y2_2);
        break;
    case 2:
        y1_0 = coefs[0];
        y2_0 = coefs[4];
        b1_0 = coefs[8];
        b2_0 = coefs[12];
        a0_0 = coefs[16];
        y1_1 = coefs[1];
        y2_1 = coefs[5];
        b1_1 = coefs[9];
        b2_1 = coefs[13];
        a0_1 = coefs[17];

        in = in0;
        out = unit->m_buf - 1;
        LooP(unit->mRate->mFilterLoops) {
            inf = *++in;
            y0_0 = inf + b1_0 * y1_0 + b2_0 * y2_0;
            y0_1 = inf + b1_1 * y1_1 + b2_1 * y2_1;
            *++out = a0_0 * y0_0 + a0_1 * y0_1;

            inf = *++in;
            y2_0 = inf + b1_0 * y0_0 + b2_0 * y1_0;
            y2_1 = inf + b1_1 * y0_1 + b2_1 * y1_1;
            *++out = a0_0 * y2_0 + a0_1 * y2_1;

            inf = *++in;
            y1_0 = inf + b1_0 * y2_0 + b2_0 * y0_0;
            y1_1 = inf + b1_1 * y2_1 + b2_1 * y0_1;
            *++out = a0_0 * y1_0 + a0_1 * y1_1;
        }
        LooP(unit->mRate->mFilterRemain) {
            inf = *++in;
            y0_0 = inf + b1_0 * y1_0 + b2_0 * y2_0;
            y0_1 = inf + b1_1 * y1_1 + b2_1 * y2_1;
            *++out = a0_0 * y0_0 + a0_1 * y0_1;
            y2_0 = y1_0;
            y1_0 = y0_0;
            y2_1 = y1_1;
            y1_1 = y0_1;
        }
        coefs[0] = zapgremlins(y1_0);
        coefs[4] = zapgremlins(y2_0);
        coefs[1] = zapgremlins(y1_1);
        coefs[5] = zapgremlins(y2_1);
        break;
    case 1:
        y1_0 = coefs[0];
        y2_0 = coefs[4];
        b1_0 = coefs[8];
        b2_0 = coefs[12];
        a0_0 = coefs[16];

        // Print("rcoefs %g %g %g %g %g\n", y1_0, y2_0, b1_0, b2_0, a0_0);
        in = in0;
        out = unit->m_buf - 1;
        LooP(unit->mRate->mFilterLoops) {
            inf = *++in;
            y0_0 = inf + b1_0 * y1_0 + b2_0 * y2_0;
            *++out = a0_0 * y0_0;

            inf = *++in;
            y2_0 = inf + b1_0 * y0_0 + b2_0 * y1_0;
            *++out = a0_0 * y2_0;

            inf = *++in;
            y1_0 = inf + b1_0 * y2_0 + b2_0 * y0_0;
            *++out = a0_0 * y1_0;
            // Print("out %g %g %g\n", y0_0, y2_0, y1_0);
        }
        LooP(unit->mRate->mFilterRemain) {
            inf = *++in;
            y0_0 = inf + b1_0 * y1_0 + b2_0 * y2_0;
            *++out = a0_0 * y0_0;
            y2_0 = y1_0;
            y1_0 = y0_0;
            // Print("out %g\n", y0_0);
        }
        /*
        coefs[0] = y1_0;	coefs[4] = y2_0;
        */
        coefs[0] = zapgremlins(y1_0);
        coefs[4] = zapgremlins(y2_0);
        break;
    case 0:
        out = unit->m_buf - 1;
        LooP(inNumSamples) { *++out = 0.f; }
        break;
    }

    coefs = unit->m_coefs;

    for (int i = 0; i < imax; ++i) {
        y1_0 = coefs[0];
        y2_0 = coefs[4];
        b1_0 = coefs[8];
        b2_0 = coefs[12];
        a0_0 = coefs[16];
        y1_1 = coefs[1];
        y2_1 = coefs[5];
        b1_1 = coefs[9];
        b2_1 = coefs[13];
        a0_1 = coefs[17];
        y1_2 = coefs[2];
        y2_2 = coefs[6];
        b1_2 = coefs[10];
        b2_2 = coefs[14];
        a0_2 = coefs[18];
        y1_3 = coefs[3];
        y2_3 = coefs[7];
        b1_3 = coefs[11];
        b2_3 = coefs[15];
        a0_3 = coefs[19];

        in = in0;
        out = unit->m_buf - 1;
        LooP(unit->mRate->mFilterLoops) {
            inf = *++in;
            y0_0 = inf + b1_0 * y1_0 + b2_0 * y2_0;
            y0_1 = inf + b1_1 * y1_1 + b2_1 * y2_1;
            y0_2 = inf + b1_2 * y1_2 + b2_2 * y2_2;
            y0_3 = inf + b1_3 * y1_3 + b2_3 * y2_3;
            *++out += a0_0 * y0_0 + a0_1 * y0_1 + a0_2 * y0_2 + a0_3 * y0_3;

            inf = *++in;
            y2_0 = inf + b1_0 * y0_0 + b2_0 * y1_0;
            y2_1 = inf + b1_1 * y0_1 + b2_1 * y1_1;
            y2_2 = inf + b1_2 * y0_2 + b2_2 * y1_2;
            y2_3 = inf + b1_3 * y0_3 + b2_3 * y1_3;
            *++out += a0_0 * y2_0 + a0_1 * y2_1 + a0_2 * y2_2 + a0_3 * y2_3;

            inf = *++in;
            y1_0 = inf + b1_0 * y2_0 + b2_0 * y0_0;
            y1_1 = inf + b1_1 * y2_1 + b2_1 * y0_1;
            y1_2 = inf + b1_2 * y2_2 + b2_2 * y0_2;
            y1_3 = inf + b1_3 * y2_3 + b2_3 * y0_3;
            *++out += a0_0 * y1_0 + a0_1 * y1_1 + a0_2 * y1_2 + a0_3 * y1_3;
        }
        LooP(unit->mRate->mFilterRemain) {
            inf = *++in;
            y0_0 = inf + b1_0 * y1_0 + b2_0 * y2_0;
            y0_1 = inf + b1_1 * y1_1 + b2_1 * y2_1;
            y0_2 = inf + b1_2 * y1_2 + b2_2 * y2_2;
            y0_3 = inf + b1_3 * y1_3 + b2_3 * y2_3;
            *++out += a0_0 * y0_0 + a0_1 * y0_1 + a0_2 * y0_2 + a0_3 * y0_3;
            y2_0 = y1_0;
            y1_0 = y0_0;
            y2_1 = y1_1;
            y1_1 = y0_1;
            y2_2 = y1_2;
            y1_2 = y0_2;
            y2_3 = y1_3;
            y1_3 = y0_3;
        }
        coefs[0] = zapgremlins(y1_0);
        coefs[4] = zapgremlins(y2_0);
        coefs[1] = zapgremlins(y1_1);
        coefs[5] = zapgremlins(y2_1);
        coefs[2] = zapgremlins(y1_2);
        coefs[6] = zapgremlins(y2_2);
        coefs[3] = zapgremlins(y1_3);
        coefs[7] = zapgremlins(y2_3);
        coefs += 20;
    }

    float x0;
    float x1 = unit->m_x1;
    float x2 = unit->m_x2;

    in = unit->m_buf - 1;
    out = out0;
    LooP(unit->mRate->mFilterLoops) {
        x0 = *++in;
        *++out = x0 - x2;
        x2 = *++in;
        *++out = x2 - x1;
        x1 = *++in;
        *++out = x1 - x0;
    }
    LooP(unit->mRate->mFilterRemain) {
        x0 = *++in;
        *++out = x0 - x2;
        x2 = x1;
        x1 = x0;
    }

    unit->m_x1 = x1;
    unit->m_x2 = x2;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

static void normalize_samples(int size, float* data, float peak) {
    float maxamp = 0.f;
    for (int i = 0; i < size; ++i) {
        float absamp = std::abs(data[i]);
        if (absamp > maxamp)
            maxamp = absamp;
    }
    if (maxamp != 0.f && maxamp != peak) {
        float ampfac = peak / maxamp;
        for (int i = 0; i < size; ++i) {
            data[i] *= ampfac;
        }
    }
}

static void normalize_wsamples(int size, float* data, float peak) {
    float maxamp = 0.f;
    for (int i = 0; i < size; i += 2) {
        float absamp = std::abs(data[i] + data[i + 1]);
        if (absamp > maxamp)
            maxamp = absamp;
    }
    if (maxamp != 0.f && maxamp != peak) {
        float ampfac = peak / maxamp;
        for (int i = 0; i < size; ++i) {
            data[i] *= ampfac;
        }
    }
}

static void add_partial(int size, float* data, double partial, double amp, double phase) {
    if (amp == 0.0)
        return;
    double w = (partial * 2.0 * 3.1415926535897932384626433832795) / (double)size;
    for (int i = 0; i < size; ++i) {
        data[i] += amp * sin(phase);
        phase += w;
    }
}

static void add_wpartial(int size, float* data, double partial, double amp, double phase) {
    if (amp == 0.0)
        return;
    int size2 = size >> 1;
    double w = (partial * 2.0 * 3.1415926535897932384626433832795) / (double)size2;
    double cur = amp * sin(phase);
    phase += w;
    for (int i = 0; i < size; i += 2) {
        double next = amp * sin(phase);
        data[i] += 2 * cur - next;
        data[i + 1] += next - cur;
        cur = next;
        phase += w;
    }
}

static void add_chebyshev(int size, float* data, double partial, double amp) {
    if (amp == 0.0)
        return;
    double w = 2.0 / (double)size;
    double phase = -1.0;
    double offset = -amp * cos(partial * pi2);
    for (int i = 0; i < size; ++i) {
        data[i] += amp * cos(partial * acos(phase)) + offset;
        phase += w;
    }
}

static void add_wchebyshev(int size, float* data, double partial, double amp) {
    if (amp == 0.0)
        return;
    int size2 = size >> 1;
    double w = 2.0 / (double)size2;
    double phase = -1.0;
    double offset = -amp * cos(partial * pi2);
    double cur = amp * cos(partial * acos(phase)) + offset;
    phase += w;
    for (int i = 0; i < size; i += 2) {
        double next = amp * cos(partial * acos(phase)) + offset;
        data[i] += 2 * cur - next;
        data[i + 1] += next - cur;
        cur = next;
        phase += w;
    }
}

static void cantorFill(int size, float* data) // long offset, double amp)
{
	float *out = ZOUT(0);
	float freq1in = ZIN0(0);
	float freq2in = ZIN0(1);
	float freq3in = ZIN0(2);

	int32 phase1 = unit->m_phase1;
	int32 phase2 = unit->m_phase2;
	int32 phase3 = unit->m_phase3;
	float cpstoinc = unit->m_cpstoinc;
	int32 freq1 = (int32)(cpstoinc * freq1in);
	int32 freq2 = (int32)(cpstoinc * freq2in);
	int32 freq3 = (int32)(cpstoinc * freq3in);
	float* sine = ft->mSine;
	int32 formfreq = sc_max(freq1, freq3);
	LOOP1(inNumSamples,
		if (phase3 < onecyc13) {
			 float val = (*(float*)((char*)sine + (((phase3 + tqcyc13) >> xlobits) & xlomask13)) + 1.f)
			*  *(float*)((char*)sine + ((phase2 >> xlobits) & xlomask13)); //mtm
			 printf("[Formant_next] %f\n", val);//mtm
			 ZXP(out) = val;//mtm
//			ZXP(out) = (*(float*)((char*)sine + (((phase3 + tqcyc13) >> xlobits) & xlomask13)) + 1.f)
//			         *  *(float*)((char*)sine + ((phase2 >> xlobits) & xlomask13));
			phase3 += formfreq;
		} else {
			printf("[Formant_next else] %f\n", 0.f);//mtm
			ZXP(out) = 0.f;
		}
		phase1 += freq1;
		phase2 += freq2;
		if (phase1 > onecyc13) {
			phase1 -= onecyc13;
			phase2 = phase1 * freq2 / freq1;
			phase3 = phase1 * freq3 / freq1;
		}
	);

	unit->m_phase1 = phase1;
	unit->m_phase2 = phase2;
	unit->m_phase3 = phase3;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

inline float lookup13(float* table, int32 pphase)
{
	float pfrac = PhaseFrac(pphase);
	float* tbl = (float*)((char*)table + ((pphase >> xlobits) & xlomask13));
	return lininterp(pfrac, tbl[0], tbl[1]);
}

void Blip_Ctor(Blip *unit)
{
	SETCALC(Blip_next);
	unit->m_freqin = ZIN0(0);
	unit->m_numharm = (int32)ZIN0(1);

	unit->m_cpstoinc = ft->mSineSize * SAMPLEDUR * 65536. * 0.5;
	int32 N = unit->m_numharm;
	int32 maxN = (int32)((SAMPLERATE * 0.5) / unit->m_freqin);
	if (N > maxN) N = maxN;
	if (N < 1) N = 1;
	unit->m_N = N;
	unit->m_scale = 0.5/N;
	unit->m_phase = 0;

	printf("[Blip] init sample:\n\t");
	Blip_next(unit, 1);
	printf("[Blip] first sample:\n\t");
	
	unit->m_N = N;//mtm
	unit->m_scale = 0.5/N;//mtm
	unit->m_phase = 0;//mtm
}

void Blip_next(Blip *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freqin = ZIN0(0);
	int numharm = (int32)ZIN0(1);

	int32 phase = unit->m_phase;

	float* numtbl = ft->mSine;
	float* dentbl = ft->mCosecant;

	int32 freq, N, prevN;
	float scale, prevscale;
	bool crossfade;
	if (numharm != unit->m_numharm || freqin != unit->m_freqin) {
		N = numharm;
		int32 maxN = (int32)((SAMPLERATE * 0.5) / freqin);
		if (N > maxN) {
			float maxfreqin;
			N = maxN;
			maxfreqin = sc_max(unit->m_freqin, freqin);
			freq = (int32)(unit->m_cpstoinc * maxfreqin);
		} else {
			if (N < 1) { N = 1; }
			freq = (int32)(unit->m_cpstoinc * freqin);
		}
		crossfade = N != unit->m_N;
		prevN = unit->m_N;
		prevscale = unit->m_scale;
		unit->m_N = N;
		unit->m_scale = scale = 0.5/N;
	} else {
		N = unit->m_N;
		freq = (int32)(unit->m_cpstoinc * freqin);
		scale = unit->m_scale;
		crossfade = false;
	}
	int32 N2 = 2*N+1;

	if (crossfade) {
		int32 prevN2 = 2 * prevN + 1;
		float xfade_slope = unit->mRate->mSlopeFactor;
		float xfade = 0.f;
		LOOP1(inNumSamples,
			float* tbl = (float*)((char*)dentbl + ((phase >> xlobits) & xlomask13));
			float t0 = tbl[0];
			float t1 = tbl[1];
			if (t0 == kBadValue || t1 == kBadValue) {
				tbl = (float*)((char*)numtbl + ((phase >> xlobits) & xlomask13));
				t0 = tbl[0];
				t1 = tbl[1];
				float pfrac = PhaseFrac(phase);
				float denom = t0 + (t1 - t0) * pfrac;
				if (std::abs(denom) < 0.0005f) {
					   printf("[Blip_next 1] %f\n", 1.f);//mtm
					ZXP(out) = 1.f;
				} else {
					int32 rphase = phase * prevN2;
					pfrac = PhaseFrac(rphase);
					tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
					float numer = lininterp(pfrac, tbl[0], tbl[1]);
					float n1 = (numer / denom - 1.f) * prevscale;

					rphase = phase * N2;
					pfrac = PhaseFrac(rphase);
					tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
					numer = lininterp(pfrac, tbl[0], tbl[1]);
					float n2 = (numer / denom - 1.f) * scale;
					   printf("[Blip_next 2] %f\n", lininterp(xfade, n1, n2));//mtm
					ZXP(out) = lininterp(xfade, n1, n2);
				}
			} else {
				float pfrac = PhaseFrac(phase);
				float denom = t0 + (t1 - t0) * pfrac;

				int32 rphase = phase * prevN2;
				pfrac = PhaseFrac(rphase);
				float* tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
				float numer = lininterp(pfrac, tbl[0], tbl[1]);
				float n1 = (numer * denom - 1.f) * prevscale;

				rphase = phase * N2;
				pfrac = PhaseFrac(rphase);
				tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
				numer = lininterp(pfrac, tbl[0], tbl[1]);
				float n2 = (numer * denom - 1.f) * scale;
				   printf("[Blip_next 3] %f\n", lininterp(xfade, n1, n2));//mtm
				ZXP(out) = lininterp(xfade, n1, n2);
			}
			phase += freq;
			xfade += xfade_slope;
		);
	} else {
		// hmm, if freq is above sr/4 then revert to sine table osc w/ no interpolation ?
		// why bother, it isn't a common choice for a fundamental.
		LOOP1(inNumSamples,
			float* tbl = (float*)((char*)dentbl + ((phase >> xlobits) & xlomask13));
			float t0 = tbl[0];
			float t1 = tbl[1];
			if (t0 == kBadValue || t1 == kBadValue) {
				tbl = (float*)((char*)numtbl + ((phase >> xlobits) & xlomask13));
				t0 = tbl[0];
				t1 = tbl[1];
				float pfrac = PhaseFrac(phase);
				float denom = t0 + (t1 - t0) * pfrac;
				if (std::abs(denom) < 0.0005f) {
					  printf("[Blip_next else 4] %f\n", 1.f);//mtm
					ZXP(out) = 1.f;
				} else {
					int32 rphase = phase * N2;
					pfrac = PhaseFrac(rphase);
					tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
					float numer = lininterp(pfrac, tbl[0], tbl[1]);
					  printf("[Blip_next else 5] %f\n", (numer / denom - 1.f) * scale);//mtm
					ZXP(out) = (numer / denom - 1.f) * scale;
				}
			} else {
				float pfrac = PhaseFrac(phase);
				float denom = t0 + (t1 - t0) * pfrac;
				int32 rphase = phase * N2;
				pfrac = PhaseFrac(rphase);
				tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
				float numer = lininterp(pfrac, tbl[0], tbl[1]);
				 printf("[Blip_next else 6] %f\n", (numer * denom - 1.f) * scale);//mtm
				ZXP(out) = (numer * denom - 1.f) * scale;
			}
			phase += freq;
		);
	}

	unit->m_phase = phase;
	unit->m_freqin = freqin;
	unit->m_numharm = numharm;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////


void Saw_Ctor(Saw *unit)
{
	SETCALC(Saw_next);
	unit->m_freqin = ZIN0(0);

	unit->m_cpstoinc = ft->mSineSize * SAMPLEDUR * 65536. * 0.5;
	unit->m_N = (int32)((SAMPLERATE * 0.5) / unit->m_freqin);
	unit->m_scale = 0.5/unit->m_N;
	unit->m_phase = 0;
	unit->m_y1 = -0.46f;

	printf("[Saw] init sample:\n\t");
//	ZOUT0(0) = 0.f; // <<<<<<<<<<<<< needs to be set properly
	Saw_next(unit, 1);//mtm
	printf("[Saw] first sample:\n\t");
	
	unit->m_scale = 0.5/unit->m_N;//mtm
	unit->m_phase = 0;//mtm
	unit->m_y1 = -0.46f;//mtm
}

void Saw_next(Saw *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freqin = ZIN0(0);

	int32 phase = unit->m_phase;
	float y1 = unit->m_y1;

	float* numtbl = ft->mSine;
	float* dentbl = ft->mCosecant;

	int32 freq, N, prevN;
	float scale, prevscale;
	bool crossfade;
	if (freqin != unit->m_freqin) {
		N = (int32)((SAMPLERATE * 0.5) / freqin);
		if (N != unit->m_N) {
			float maxfreqin;
			maxfreqin = sc_max(unit->m_freqin, freqin);
			freq = (int32)(unit->m_cpstoinc * maxfreqin);
			crossfade = true;
		} else {
			freq = (int32)(unit->m_cpstoinc * freqin);
			crossfade = false;
		}
		prevN = unit->m_N;
		prevscale = unit->m_scale;
		unit->m_N = N;
		unit->m_scale = scale = 0.5/N;
	} else {
		N = unit->m_N;
		freq = (int32)(unit->m_cpstoinc * freqin);
		scale = unit->m_scale;
		crossfade = false;
	}
	int32 N2 = 2*N+1;

	if (crossfade) {
		int32 prevN2 = 2 * prevN + 1;
		float xfade_slope = unit->mRate->mSlopeFactor;
		float xfade = 0.f;
		LOOP1(inNumSamples,
			float* tbl = (float*)((char*)dentbl + ((phase >> xlobits) & xlomask13));
			float t0 = tbl[0];
			float t1 = tbl[1];
			if (t0 == kBadValue || t1 == kBadValue) {
				tbl = (float*)((char*)numtbl + ((phase >> xlobits) & xlomask13));
				t0 = tbl[0];
				t1 = tbl[1];
				float pfrac = PhaseFrac(phase);
				float denom = t0 + (t1 - t0) * pfrac;
				if (std::abs(denom) < 0.0005f) {
					 printf("[Saw_next 1] %f\n", 1.f + 0.999f * y1);//mtm
					ZXP(out) = y1 = 1.f + 0.999f * y1;
				} else {
					int32 rphase = phase * prevN2;
					pfrac = PhaseFrac(rphase);
					tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
					float numer = lininterp(pfrac, tbl[0], tbl[1]);
					float n1 = (numer / denom - 1.f) * prevscale;

					rphase = phase * N2;
					pfrac = PhaseFrac(rphase);
					tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
					numer = lininterp(pfrac, tbl[0], tbl[1]);
					float n2 = (numer / denom - 1.f) * scale;
					 printf("[Saw_next 2] %f\n", n1 + xfade * (n2 - n1) + 0.999f * y1);//mtm
					ZXP(out) = y1 = n1 + xfade * (n2 - n1) + 0.999f * y1;
				}

			} else {
				float pfrac = PhaseFrac(phase);
				float denom = t0 + (t1 - t0) * pfrac;

				int32 rphase = phase * prevN2;
				pfrac = PhaseFrac(rphase);
				tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
				float numer = lininterp(pfrac, tbl[0], tbl[1]);
				float n1 = (numer * denom - 1.f) * prevscale;

				rphase = phase * N2;
				pfrac = PhaseFrac(rphase);
				tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
				numer = lininterp(pfrac, tbl[0], tbl[1]);
				float n2 = (numer * denom - 1.f) * scale;
 printf("[Saw_next 3] %f\n", n1 + xfade * (n2 - n1) + 0.999f * y1);//mtm
				ZXP(out) = y1 = n1 + xfade * (n2 - n1) + 0.999f * y1;
			}
			phase += freq;
			xfade += xfade_slope;
		);
	} else {
		// hmm, if freq is above sr/4 then revert to sine table osc ?
		// why bother, it isn't a common choice for a fundamental.
		LOOP1(inNumSamples,
			float* tbl = (float*)((char*)dentbl + ((phase >> xlobits) & xlomask13));
			float t0 = tbl[0];
			float t1 = tbl[1];
			if (t0 == kBadValue || t1 == kBadValue) {
				tbl = (float*)((char*)numtbl + ((phase >> xlobits) & xlomask13));
				t0 = tbl[0];
				t1 = tbl[1];
				float pfrac = PhaseFrac(phase);
				float denom = t0 + (t1 - t0) * pfrac;
				if (std::abs(denom) < 0.0005f) {
					 printf("[Saw_next 4] %f\n", 1.f + 0.999f * y1);//mtm
					ZXP(out) = y1 = 1.f + 0.999f * y1;
				} else {
					int32 rphase = phase * N2;
					pfrac = PhaseFrac(rphase);
					tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
					float numer = lininterp(pfrac, tbl[0], tbl[1]);
					 printf("[Saw_next 5] %f\n", (numer / denom - 1.f) * scale + 0.999f * y1);//mtm
					ZXP(out) = y1 = (numer / denom - 1.f) * scale + 0.999f * y1;
				}
			} else {
				float pfrac = PhaseFrac(phase);
				float denom = t0 + (t1 - t0) * pfrac;
				int32 rphase = phase * N2;
				pfrac = PhaseFrac(rphase);
				float* tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
				float numer = lininterp(pfrac, tbl[0], tbl[1]);
				 printf("[Saw_next 6] %f\n", (numer * denom - 1.f) * scale + 0.999f * y1);//mtm
				ZXP(out) = y1 = (numer * denom - 1.f) * scale + 0.999f * y1;
			}
			phase += freq;
		);
	}

	unit->m_y1 = y1;
	unit->m_phase = phase;
	unit->m_freqin = freqin;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////


void Pulse_Ctor(Pulse *unit)
{
	SETCALC(Pulse_next);
	unit->m_freqin = ZIN0(0);

	unit->m_cpstoinc = ft->mSineSize * SAMPLEDUR * 65536. * 0.5;
	unit->m_N = (int32)((SAMPLERATE * 0.5) / unit->m_freqin);
	unit->m_scale = 0.5/unit->m_N;
//	unit->m_phase = 0;
//	unit->m_phaseoff = 0;
//	unit->m_y1 = 0.f;//mtm
//	ZOUT0(0) = 0.f;//mtm
	
	printf("[Pulse] mscale - before: %f\n\t", unit->m_scale);
	int32 initphaseoff;
	unit->m_phase = 0.f;
	unit->m_phaseoff = initphaseoff = (int32)(ZIN0(1) * (1L << 28));
//	unit->m_y1 = 0.f;//mtm
	unit->m_y1 = -0.476432;//mtm

	printf("[Pulse] init sample:\n\t");
	Pulse_next(unit, 1);
	printf("[Pulse] first sample:\n\t");

	unit->m_phase = 0.f;
//	unit->m_phaseoff = 0;//mtm
	unit->m_phaseoff = initphaseoff;
//	unit->m_y1 = 0.f;//mtm
	unit->m_y1 = -0.476432;//mtm
//	unit->m_scale = 0.5/unit->m_N;//mtm
	printf("[Pulse] mscale - after: %f\n\t", unit->m_scale);
}

void Pulse_next(Pulse *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freqin = ZIN0(0);
	float duty = ZIN0(1);

	int32 phase = unit->m_phase;
	float y1 = unit->m_y1;

	float* numtbl = ft->mSine;
	float* dentbl = ft->mCosecant;

	int32 freq, N, prevN;
	float scale, prevscale;
	bool crossfade;

	if (freqin != unit->m_freqin) {
		N = (int32)((SAMPLERATE * 0.5) / freqin);
		if (N != unit->m_N) {
			float maxfreqin;
			maxfreqin = sc_max(unit->m_freqin, freqin);
			freq = (int32)(unit->m_cpstoinc * maxfreqin);
			crossfade = true;
		} else {
			freq = (int32)(unit->m_cpstoinc * freqin);
			crossfade = false;
		}
		prevN = unit->m_N;
		prevscale = unit->m_scale;
		unit->m_N = N;
		unit->m_scale = scale = 0.5/N;
	} else {
		N = unit->m_N;
		freq = (int32)(unit->m_cpstoinc * freqin);
		scale = unit->m_scale;
		crossfade = false;
	}
	int32 N2 = 2*N+1;

	int32 phaseoff = unit->m_phaseoff;
	int32 next_phaseoff = (int32)(duty * (1L << 28));
	int32 phaseoff_slope = (int32)((next_phaseoff - phaseoff) * unit->mRate->mSlopeFactor);
	unit->m_phaseoff = next_phaseoff;
	
	float rscale = 1.f / scale + 1.f;
	float pul1, pul2;

	if (crossfade) {
		int32 prevN2 = 2 * prevN + 1;
		float xfade_slope = unit->mRate->mSlopeFactor;
		float xfade = 0.f;
		LOOP1(inNumSamples,
			float* tbl = (float*)((char*)dentbl + ((phase >> xlobits) & xlomask13));
			float t0 = tbl[0];
			float t1 = tbl[1];
			if (t0 == kBadValue || t1 == kBadValue) {
				tbl = (float*)((char*)numtbl + ((phase >> xlobits) & xlomask13));
				t0 = tbl[0];
				t1 = tbl[1];
				float pfrac = PhaseFrac(phase);
				float denom = t0 + (t1 - t0) * pfrac;
				if (std::abs(denom) < 0.0005f) {
					pul1 = 1.f;
				} else {
					int32 rphase = phase * prevN2;
					pfrac = PhaseFrac(rphase);
					tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
					float numer = lininterp(pfrac, tbl[0], tbl[1]);
					float n1 = (numer / denom - 1.f) * prevscale;

					rphase = phase * N2;
					pfrac = PhaseFrac(rphase);
					tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
					numer = lininterp(pfrac, tbl[0], tbl[1]);
					float n2 = (numer / denom - 1.f) * scale;

					pul1 = lininterp(xfade, n1, n2);
				}
			} else {
				float pfrac = PhaseFrac(phase);
				float denom = lininterp(pfrac, t0, t1);

				int32 rphase = phase * prevN2;
				pfrac = PhaseFrac(rphase);
				tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
				float numer = lininterp(pfrac, tbl[0], tbl[1]);
				float n1 = (numer * denom - 1.f) * prevscale;

				rphase = phase * N2;
				pfrac = PhaseFrac(rphase);
				tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
				numer = lininterp(pfrac, tbl[0], tbl[1]);
				float n2 = (numer * denom - 1.f) * scale;

				pul1 = lininterp(xfade, n1, n2);
			}

			int32 phase2 = phase + phaseoff;
			tbl = (float*)((char*)dentbl + ((phase2 >> xlobits) & xlomask13));
			t0 = tbl[0];
			t1 = tbl[1];
			if (t0 == kBadValue || t1 == kBadValue) {
				tbl = (float*)((char*)numtbl + ((phase2 >> xlobits) & xlomask13));
				t0 = tbl[0];
				t1 = tbl[1];
				float pfrac = PhaseFrac(phase2);
				float denom = t0 + (t1 - t0) * pfrac;
				if (std::abs(denom) < 0.0005f) {
					pul2 = 1.f;
				} else {
					int32 rphase = phase2 * prevN2;
					pfrac = PhaseFrac(rphase);
					tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
					float numer = lininterp(pfrac, tbl[0], tbl[1]);
					float n1 = (numer / denom - 1.f) * prevscale;

					rphase = phase2 * N2;
					pfrac = PhaseFrac(rphase);
					tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
					numer = lininterp(pfrac, tbl[0], tbl[1]);
					float n2 = (numer / denom - 1.f) * scale;

					pul2 = lininterp(xfade, n1, n2);
				}
			} else {
				float pfrac = PhaseFrac(phase2);
				float denom = t0 + (t1 - t0) * pfrac;

				int32 rphase = phase2 * prevN2;
				pfrac = PhaseFrac(rphase);
				tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
				float numer = lininterp(pfrac, tbl[0], tbl[1]);
				float n1 = (numer * denom - 1.f) * prevscale;

				rphase = phase2 * N2;
				pfrac = PhaseFrac(rphase);
				tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
				numer = lininterp(pfrac, tbl[0], tbl[1]);
				float n2 = (numer * denom - 1.f) * scale;

				pul2 = lininterp(xfade, n1, n2);
			}
			  float val = pul1 - pul2 + 0.999f * y1;
			  printf("[Pulse_next 5] %f\n", val); //mtm
			  ZXP(out) = y1 = val;//mtm
//			ZXP(out) = y1 = pul1 - pul2 + 0.999f * y1;
			phase += freq;
			phaseoff += phaseoff_slope;
			xfade += xfade_slope;
		);
	} else {
		LOOP1(inNumSamples,
			float* tbl = (float*)((char*)dentbl + ((phase >> xlobits) & xlomask13));
			float t0 = tbl[0];
			float t1 = tbl[1];
			if (t0 == kBadValue || t1 == kBadValue) {
				tbl = (float*)((char*)numtbl + ((phase >> xlobits) & xlomask13));
				t0 = tbl[0];
				t1 = tbl[1];
				float pfrac = PhaseFrac(phase);
				float denom = t0 + (t1 - t0) * pfrac;
				if (std::abs(denom) < 0.0005f) {
					pul1 = rscale;
				} else {
					int32 rphase = phase * N2;
					pfrac = PhaseFrac(rphase);
					tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
					float numer = lininterp(pfrac, tbl[0], tbl[1]);
					pul1 = numer / denom;
				}
			} else {
				float pfrac = PhaseFrac(phase);
				float denom = t0 + (t1 - t0) * pfrac;
				int32 rphase = phase * N2;
				pfrac = PhaseFrac(rphase);
				float* tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
				float numer = lininterp(pfrac, tbl[0], tbl[1]);
				pul1 = (numer * denom);
			}

			int32 phase2 = phase + phaseoff;
			  printf("\tphase: %d  offset: %d\n\t", phase, phaseoff);//mtm
			tbl = (float*)((char*)dentbl + ((phase2 >> xlobits) & xlomask13));
			t0 = tbl[0];
			t1 = tbl[1];
			if (t0 == kBadValue || t1 == kBadValue) {
				tbl = (float*)((char*)numtbl + ((phase2 >> xlobits) & xlomask13));
				t0 = tbl[0];
				t1 = tbl[1];
				float pfrac = PhaseFrac(phase2);
				float denom = t0 + (t1 - t0) * pfrac;
				if (std::abs(denom) < 0.0005f) {
					pul2 = rscale;
				} else {
					int32 rphase = phase2 * N2;
					pfrac = PhaseFrac(rphase);
					tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
					float numer = lininterp(pfrac, tbl[0], tbl[1]);
					pul2 = numer / denom;
				}
			} else {
				float pfrac = PhaseFrac(phase2);
				float denom = t0 + (t1 - t0) * pfrac;
				int32 rphase = phase2 * N2;
				pfrac = PhaseFrac(rphase);
				float* tbl = (float*)((char*)numtbl + ((rphase >> xlobits) & xlomask13));
				float numer = lininterp(pfrac, tbl[0], tbl[1]);

				pul2 = (numer * denom);
			}
			  float val = (pul1 - pul2) * scale + 0.999f * y1;//mtm
			  printf("[Pulse_next 6] %f\n", val); //mtm
			  ZXP(out) = y1 = val;//mtm
//			ZXP(out) = y1 = (pul1 - pul2) * scale + 0.999f * y1;//mtm
			phase += freq;
			phaseoff += phaseoff_slope;
		);
	}

	unit->m_y1 = y1;
	unit->m_phase = phase;
	unit->m_freqin = freqin;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

static float Klang_SetCoefs(Klang *unit)
{
	unit->m_numpartials = (unit->mNumInputs - 2)/3;

	int numcoefs = unit->m_numpartials * 3;
	unit->m_coefs = (float*)RTAlloc(unit->mWorld, numcoefs * sizeof(float));

	if (!unit->m_coefs) {
		Print("Klang: RT memory allocation failed\n");
		SETCALC(ClearUnitOutputs);
		return 0.f;
	}

	float freqscale = ZIN0(0) * unit->mRate->mRadiansPerSample;
	float freqoffset = ZIN0(1) * unit->mRate->mRadiansPerSample;

	float outf = 0.;
	float* coefs = unit->m_coefs - 1;

	for (int i=0,j=2; i<unit->m_numpartials; ++i,j+=3) {
		float w = ZIN0(j) * freqscale + freqoffset;
		float level = ZIN0(j+1);
		float phase = ZIN0(j+2);

		if (phase != 0.f) {
			outf += *++coefs = level * sin(phase);		// y1
			        *++coefs = level * sin(phase - w);	// y2
		} else {
			outf += *++coefs = 0.f;						// y1
			        *++coefs = level * -sin(w);	// y2
		}
		*++coefs = 2. * cos(w);		// b1
	}
	return outf;
}

void Klang_Ctor(Klang *unit)
{
	SETCALC(Klang_next);
	printf("[Klang] init sample:\n\t%f\n", Klang_SetCoefs(unit));
	ZOUT0(0) = Klang_SetCoefs(unit);
	printf("[Klang] first sample:\n\t");
}

void Klang_Dtor(Klang *unit)
{
	RTFree(unit->mWorld, unit->m_coefs);
}

void Klang_next(Klang *unit, int inNumSamples)
{
	float *out0 = ZOUT(0);

	float* out;
	float y0_0, y1_0, y2_0, b1_0;
	float y0_1, y1_1, y2_1, b1_1;
	float y0_2, y1_2, y2_2, b1_2;
	float y0_3, y1_3, y2_3, b1_3;
	float outf;

	float* coefs = unit->m_coefs - 1;
	int32 numpartials = unit->m_numpartials;

	switch (numpartials & 3) {
		case 3 :
			y1_0 = *++coefs;	y2_0 = *++coefs;	b1_0 = *++coefs;
			y1_1 = *++coefs;	y2_1 = *++coefs;	b1_1 = *++coefs;
			y1_2 = *++coefs;	y2_2 = *++coefs;	b1_2 = *++coefs;

			out = out0;
			LOOP(unit->mRate->mFilterLoops,
				outf  = y0_0 = b1_0 * y1_0 - y2_0;
				outf += y0_1 = b1_1 * y1_1 - y2_1;
				outf += y0_2 = b1_2 * y1_2 - y2_2;
				ZXP(out) = outf;

				outf  = y2_0 = b1_0 * y0_0 - y1_0;
				outf += y2_1 = b1_1 * y0_1 - y1_1;
				outf += y2_2 = b1_2 * y0_2 - y1_2;
				ZXP(out) = outf;

				outf  = y1_0 = b1_0 * y2_0 - y0_0;
				outf += y1_1 = b1_1 * y2_1 - y0_1;
				outf += y1_2 = b1_2 * y2_2 - y0_2;
				ZXP(out) = outf;
			);
			LOOP(unit->mRate->mFilterRemain,
				outf  = y0_0 = b1_0 * y1_0 - y2_0;
				outf += y0_1 = b1_1 * y1_1 - y2_1;
				outf += y0_2 = b1_2 * y1_2 - y2_2;
				y2_0 = y1_0;	y1_0 = y0_0;
				y2_1 = y1_1;	y1_1 = y0_1;
				y2_2 = y1_2;	y1_2 = y0_2;
				ZXP(out) = outf;
			);
			coefs -= 9;
			*++coefs = y1_0;	*++coefs = y2_0;	++coefs;
			*++coefs = y1_1;	*++coefs = y2_1;	++coefs;
			*++coefs = y1_2;	*++coefs = y2_2;	++coefs;
			break;
		case 2 :
			y1_0 = *++coefs;	y2_0 = *++coefs;	b1_0 = *++coefs;
			y1_1 = *++coefs;	y2_1 = *++coefs;	b1_1 = *++coefs;

			out = out0;
			LOOP(unit->mRate->mFilterLoops,
				outf  = y0_0 = b1_0 * y1_0 - y2_0;
				outf += y0_1 = b1_1 * y1_1 - y2_1;
				ZXP(out) = outf;

				outf  = y2_0 = b1_0 * y0_0 - y1_0;
				outf += y2_1 = b1_1 * y0_1 - y1_1;
				ZXP(out) = outf;

				outf  = y1_0 = b1_0 * y2_0 - y0_0;
				outf += y1_1 = b1_1 * y2_1 - y0_1;
				ZXP(out) = outf;
			);
			LOOP(unit->mRate->mFilterRemain,
				outf  = y0_0 = b1_0 * y1_0 - y2_0;
				outf += y0_1 = b1_1 * y1_1 - y2_1;
				y2_0 = y1_0;	y1_0 = y0_0;
				y2_1 = y1_1;	y1_1 = y0_1;
				ZXP(out) = outf;
			);

			coefs -= 6;
			*++coefs = y1_0;	*++coefs = y2_0;	++coefs;
			*++coefs = y1_1;	*++coefs = y2_1;	++coefs;
			break;
		case 1 :
			y1_0 = *++coefs;	y2_0 = *++coefs;	b1_0 = *++coefs;

			out = out0;
			LOOP(unit->mRate->mFilterLoops,
				ZXP(out) = y0_0 = b1_0 * y1_0 - y2_0;

				ZXP(out) = y2_0 = b1_0 * y0_0 - y1_0;

				ZXP(out) = y1_0 = b1_0 * y2_0 - y0_0;
			);
			LOOP(unit->mRate->mFilterRemain,
				ZXP(out) = y0_0 = b1_0 * y1_0 - y2_0;
				y2_0 = y1_0;	y1_0 = y0_0;
			);

			coefs -= 3;
			*++coefs = y1_0;	*++coefs = y2_0;	++coefs;
			break;
		case 0 :
			out = out0;
			ZClear(inNumSamples, out);
			break;
	}

	int32 imax = numpartials >> 2;

	for (int i=0; i<imax; ++i) {
		y1_0 = *++coefs;	y2_0 = *++coefs;	b1_0 = *++coefs;
		y1_1 = *++coefs;	y2_1 = *++coefs;	b1_1 = *++coefs;
		y1_2 = *++coefs;	y2_2 = *++coefs;	b1_2 = *++coefs;
		y1_3 = *++coefs;	y2_3 = *++coefs;	b1_3 = *++coefs;

		out = out0;
		LOOP(unit->mRate->mFilterLoops,
			outf  = y0_0 = b1_0 * y1_0 - y2_0;
			outf += y0_1 = b1_1 * y1_1 - y2_1;
			outf += y0_2 = b1_2 * y1_2 - y2_2;
			outf += y0_3 = b1_3 * y1_3 - y2_3;
			ZXP(out) += outf;

			outf  = y2_0 = b1_0 * y0_0 - y1_0;
			outf += y2_1 = b1_1 * y0_1 - y1_1;
			outf += y2_2 = b1_2 * y0_2 - y1_2;
			outf += y2_3 = b1_3 * y0_3 - y1_3;
			ZXP(out) += outf;

			outf  = y1_0 = b1_0 * y2_0 - y0_0;
			outf += y1_1 = b1_1 * y2_1 - y0_1;
			outf += y1_2 = b1_2 * y2_2 - y0_2;
			outf += y1_3 = b1_3 * y2_3 - y0_3;
			ZXP(out) += outf;
		);
		LOOP(unit->mRate->mFilterRemain,
			outf  = y0_0 = b1_0 * y1_0 - y2_0;
			outf += y0_1 = b1_1 * y1_1 - y2_1;
			outf += y0_2 = b1_2 * y1_2 - y2_2;
			outf += y0_3 = b1_3 * y1_3 - y2_3;
			y2_0 = y1_0;	y1_0 = y0_0;
			y2_1 = y1_1;	y1_1 = y0_1;
			y2_2 = y1_2;	y1_2 = y0_2;
			y2_3 = y1_3;	y1_3 = y0_3;
			ZXP(out) += outf;
		);
		coefs -= 12;
		*++coefs = y1_0;	*++coefs = y2_0;	++coefs;
		*++coefs = y1_1;	*++coefs = y2_1;	++coefs;
		*++coefs = y1_2;	*++coefs = y2_2;	++coefs;
		*++coefs = y1_3;	*++coefs = y2_3;	++coefs;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

static void Klank_SetCoefs(Klank *unit)
{
	int numpartials = (unit->mNumInputs - 4) / 3;
	unit->m_numpartials = numpartials;

	int numcoefs = ((unit->m_numpartials + 3) & ~3) * 5;
	unit->m_coefs = (float*)RTAlloc(unit->mWorld, (numcoefs + unit->mWorld->mBufLength) * sizeof(float));
	if (!unit->m_coefs) {
		Print("Klang: RT memory allocation failed\n");
		SETCALC(ClearUnitOutputs);
		return;
	}

	unit->m_buf = unit->m_coefs + numcoefs;

	float freqscale = ZIN0(1) * unit->mRate->mRadiansPerSample;
	float freqoffset = ZIN0(2) * unit->mRate->mRadiansPerSample;
	float decayscale = ZIN0(3);

	float* coefs = unit->m_coefs;

	float sampleRate = SAMPLERATE;

	for (int i=0,j=4; i<numpartials; ++i,j+=3) {
		float w = ZIN0(j) * freqscale + freqoffset;
		float level = ZIN0(j+1);
		float time = ZIN0(j+2) * decayscale;

		float R = time == 0.f ? 0.f : exp(log001/(time * sampleRate));
		float twoR = 2.f * R;
		float R2 = R * R;
		float cost = (twoR * cos(w)) / (1.f + R2);

		int k = 20 * (i>>2) + (i & 3);
		coefs[k+0] = 0.f;					// y1
		coefs[k+4] = 0.f;					// y2
		coefs[k+8] = twoR * cost;			// b1
		coefs[k+12] = -R2;					// b2
		coefs[k+16] = level * 0.25;		// a0
		//Print("coefs %d  %g %g %g\n", i, twoR * cost, -R2, ampf * 0.25);
	}
}

void Klank_Ctor(Klank *unit)
{
	SETCALC(Klank_next);
	unit->m_x1 = unit->m_x2 = 0.f;
	Klank_SetCoefs(unit);
	printf("[Klank] init sample:\n\t%f\n", 0.f);
	ZOUT0(0) = 0.f;
	printf("[Klank] first sample:\n\t");
}

void Klank_Dtor(Klank *unit)
{
	RTFree(unit->mWorld, unit->m_coefs);
}

void Klank_next(Klank *unit, int inNumSamples)
{
	float *out0 = ZOUT(0);
	float *in0 = ZIN(0);

	float *in, *out;
	float inf;
	float y0_0, y1_0, y2_0, a0_0, b1_0, b2_0;
	float y0_1, y1_1, y2_1, a0_1, b1_1, b2_1;
	float y0_2, y1_2, y2_2, a0_2, b1_2, b2_2;
	float y0_3, y1_3, y2_3, a0_3, b1_3, b2_3;

	int32 numpartials = unit->m_numpartials;
	int32 imax = numpartials >> 2;

	float* coefs = unit->m_coefs + imax * 20;

	switch (numpartials & 3) {
		case 3 :
			y1_0 = coefs[0];	y2_0 = coefs[4];	b1_0 = coefs[8];	b2_0 = coefs[12];	a0_0 = coefs[16];
			y1_1 = coefs[1];	y2_1 = coefs[5];	b1_1 = coefs[9];	b2_1 = coefs[13];	a0_1 = coefs[17];
			y1_2 = coefs[2];	y2_2 = coefs[6];	b1_2 = coefs[10];	b2_2 = coefs[14];	a0_2 = coefs[18];

			in = in0;
			out = unit->m_buf - 1;
			LooP(unit->mRate->mFilterLoops) {
				inf = *++in;
				y0_0 = inf + b1_0 * y1_0 + b2_0 * y2_0;
				y0_1 = inf + b1_1 * y1_1 + b2_1 * y2_1;
				y0_2 = inf + b1_2 * y1_2 + b2_2 * y2_2;
				*++out = a0_0 * y0_0 + a0_1 * y0_1 + a0_2 * y0_2;

				inf = *++in;
				y2_0 = inf + b1_0 * y0_0 + b2_0 * y1_0;
				y2_1 = inf + b1_1 * y0_1 + b2_1 * y1_1;
				y2_2 = inf + b1_2 * y0_2 + b2_2 * y1_2;
				*++out = a0_0 * y2_0 + a0_1 * y2_1 + a0_2 * y2_2;

				inf = *++in;
				y1_0 = inf + b1_0 * y2_0 + b2_0 * y0_0;
				y1_1 = inf + b1_1 * y2_1 + b2_1 * y0_1;
				y1_2 = inf + b1_2 * y2_2 + b2_2 * y0_2;
				*++out = a0_0 * y1_0 + a0_1 * y1_1 + a0_2 * y1_2;
			}
			LooP(unit->mRate->mFilterRemain) {
				inf = *++in;
				y0_0 = inf + b1_0 * y1_0 + b2_0 * y2_0;
				y0_1 = inf + b1_1 * y1_1 + b2_1 * y2_1;
				y0_2 = inf + b1_2 * y1_2 + b2_2 * y2_2;
				*++out = a0_0 * y0_0 + a0_1 * y0_1 + a0_2 * y0_2;
				y2_0 = y1_0;	y1_0 = y0_0;
				y2_1 = y1_1;	y1_1 = y0_1;
				y2_2 = y1_2;	y1_2 = y0_2;
			}
			coefs[0] = zapgremlins(y1_0);	coefs[4] = zapgremlins(y2_0);
			coefs[1] = zapgremlins(y1_1);	coefs[5] = zapgremlins(y2_1);
			coefs[2] = zapgremlins(y1_2);	coefs[6] = zapgremlins(y2_2);
			break;
		case 2 :
			y1_0 = coefs[0];	y2_0 = coefs[4];	b1_0 = coefs[8];	b2_0 = coefs[12];	a0_0 = coefs[16];
			y1_1 = coefs[1];	y2_1 = coefs[5];	b1_1 = coefs[9];	b2_1 = coefs[13];	a0_1 = coefs[17];

			in = in0;
			out = unit->m_buf - 1;
			LooP(unit->mRate->mFilterLoops) {
				inf = *++in;
				y0_0 = inf + b1_0 * y1_0 + b2_0 * y2_0;
				y0_1 = inf + b1_1 * y1_1 + b2_1 * y2_1;
				*++out = a0_0 * y0_0 + a0_1 * y0_1;

				inf = *++in;
				y2_0 = inf + b1_0 * y0_0 + b2_0 * y1_0;
				y2_1 = inf + b1_1 * y0_1 + b2_1 * y1_1;
				*++out = a0_0 * y2_0 + a0_1 * y2_1;

				inf = *++in;
				y1_0 = inf + b1_0 * y2_0 + b2_0 * y0_0;
				y1_1 = inf + b1_1 * y2_1 + b2_1 * y0_1;
				*++out = a0_0 * y1_0 + a0_1 * y1_1;
			}
			LooP(unit->mRate->mFilterRemain) {
				inf = *++in;
				y0_0 = inf + b1_0 * y1_0 + b2_0 * y2_0;
				y0_1 = inf + b1_1 * y1_1 + b2_1 * y2_1;
				*++out = a0_0 * y0_0 + a0_1 * y0_1;
				y2_0 = y1_0;	y1_0 = y0_0;
				y2_1 = y1_1;	y1_1 = y0_1;
			}
			coefs[0] = zapgremlins(y1_0);	coefs[4] = zapgremlins(y2_0);
			coefs[1] = zapgremlins(y1_1);	coefs[5] = zapgremlins(y2_1);
			break;
		case 1 :
			y1_0 = coefs[0];	y2_0 = coefs[4];	b1_0 = coefs[8];	b2_0 = coefs[12];	a0_0 = coefs[16];

			//Print("rcoefs %g %g %g %g %g\n", y1_0, y2_0, b1_0, b2_0, a0_0);
			in = in0;
			out = unit->m_buf - 1;
			LooP(unit->mRate->mFilterLoops) {
				inf = *++in;
				y0_0 = inf + b1_0 * y1_0 + b2_0 * y2_0;
				*++out = a0_0 * y0_0;

				inf = *++in;
				y2_0 = inf + b1_0 * y0_0 + b2_0 * y1_0;
				*++out = a0_0 * y2_0;

				inf = *++in;
				y1_0 = inf + b1_0 * y2_0 + b2_0 * y0_0;
				*++out = a0_0 * y1_0;
				//Print("out %g %g %g\n", y0_0, y2_0, y1_0);
			}
			LooP(unit->mRate->mFilterRemain) {
				inf = *++in;
				y0_0 = inf + b1_0 * y1_0 + b2_0 * y2_0;
				*++out = a0_0 * y0_0;
				y2_0 = y1_0;	y1_0 = y0_0;
				//Print("out %g\n", y0_0);
			}
			/*
			coefs[0] = y1_0;	coefs[4] = y2_0;
			*/
			coefs[0] = zapgremlins(y1_0);	coefs[4] = zapgremlins(y2_0);
			break;
		case 0 :
			out = unit->m_buf - 1;
			LooP(inNumSamples) { *++out = 0.f; }
			break;
	}

	coefs = unit->m_coefs;

	for (int i=0; i<imax; ++i) {
		y1_0 = coefs[0];	y2_0 = coefs[4];	b1_0 = coefs[8];	b2_0 = coefs[12];	a0_0 = coefs[16];
		y1_1 = coefs[1];	y2_1 = coefs[5];	b1_1 = coefs[9];	b2_1 = coefs[13];	a0_1 = coefs[17];
		y1_2 = coefs[2];	y2_2 = coefs[6];	b1_2 = coefs[10];	b2_2 = coefs[14];	a0_2 = coefs[18];
		y1_3 = coefs[3];	y2_3 = coefs[7];	b1_3 = coefs[11];	b2_3 = coefs[15];	a0_3 = coefs[19];

		in = in0;
		out = unit->m_buf - 1;
		LooP(unit->mRate->mFilterLoops) {
			inf = *++in;
			y0_0 = inf + b1_0 * y1_0 + b2_0 * y2_0;
			y0_1 = inf + b1_1 * y1_1 + b2_1 * y2_1;
			y0_2 = inf + b1_2 * y1_2 + b2_2 * y2_2;
			y0_3 = inf + b1_3 * y1_3 + b2_3 * y2_3;
			*++out += a0_0 * y0_0 + a0_1 * y0_1 + a0_2 * y0_2 + a0_3 * y0_3;

			inf = *++in;
			y2_0 = inf + b1_0 * y0_0 + b2_0 * y1_0;
			y2_1 = inf + b1_1 * y0_1 + b2_1 * y1_1;
			y2_2 = inf + b1_2 * y0_2 + b2_2 * y1_2;
			y2_3 = inf + b1_3 * y0_3 + b2_3 * y1_3;
			*++out += a0_0 * y2_0 + a0_1 * y2_1 + a0_2 * y2_2 + a0_3 * y2_3;

			inf = *++in;
			y1_0 = inf + b1_0 * y2_0 + b2_0 * y0_0;
			y1_1 = inf + b1_1 * y2_1 + b2_1 * y0_1;
			y1_2 = inf + b1_2 * y2_2 + b2_2 * y0_2;
			y1_3 = inf + b1_3 * y2_3 + b2_3 * y0_3;
			*++out += a0_0 * y1_0 + a0_1 * y1_1 + a0_2 * y1_2 + a0_3 * y1_3;
		}
		LooP(unit->mRate->mFilterRemain) {
			inf = *++in;
			y0_0 = inf + b1_0 * y1_0 + b2_0 * y2_0;
			y0_1 = inf + b1_1 * y1_1 + b2_1 * y2_1;
			y0_2 = inf + b1_2 * y1_2 + b2_2 * y2_2;
			y0_3 = inf + b1_3 * y1_3 + b2_3 * y2_3;
			*++out += a0_0 * y0_0 + a0_1 * y0_1 + a0_2 * y0_2 + a0_3 * y0_3;
			y2_0 = y1_0;	y1_0 = y0_0;
			y2_1 = y1_1;	y1_1 = y0_1;
			y2_2 = y1_2;	y1_2 = y0_2;
			y2_3 = y1_3;	y1_3 = y0_3;
		}
		coefs[0] = zapgremlins(y1_0);	coefs[4] = zapgremlins(y2_0);
		coefs[1] = zapgremlins(y1_1);	coefs[5] = zapgremlins(y2_1);
		coefs[2] = zapgremlins(y1_2);	coefs[6] = zapgremlins(y2_2);
		coefs[3] = zapgremlins(y1_3);	coefs[7] = zapgremlins(y2_3);
		coefs += 20;
	}

	float x0;
	float x1 = unit->m_x1;
	float x2 = unit->m_x2;

	in = unit->m_buf - 1;
	out = out0;
	LooP(unit->mRate->mFilterLoops) {
		x0 = *++in;
		*++out = x0 - x2;
		x2 = *++in;
		*++out = x2 - x1;
		x1 = *++in;
		*++out = x1 - x0;
	}
	LooP(unit->mRate->mFilterRemain) {
		x0 = *++in;
		*++out = x0 - x2;
		x2 = x1;
		x1 = x0;
	}

	unit->m_x1 = x1;
	unit->m_x2 = x2;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

static void normalize_samples(int size, float* data, float peak)
{
	float maxamp = 0.f;
	for (int i=0; i<size; ++i) {
		float absamp = std::abs(data[i]);
		if (absamp > maxamp) maxamp = absamp;
	}
	if (maxamp != 0.f && maxamp != peak) {
		float ampfac = peak / maxamp;
		for (int i=0; i<size; ++i) {
			data[i] *= ampfac;
		}
	}
}

static void normalize_wsamples(int size, float* data, float peak)
{
	float maxamp = 0.f;
	for (int i=0; i<size; i+=2) {
		float absamp = std::abs(data[i] + data[i+1]);
		if (absamp > maxamp) maxamp = absamp;
	}
	if (maxamp != 0.f && maxamp != peak) {
		float ampfac = peak / maxamp;
		for (int i=0; i<size; ++i) {
			data[i] *= ampfac;
		}
	}
}

static void add_partial(int size, float *data, double partial, double amp, double phase)
{
	if (amp == 0.0) return;
	double w = (partial * 2.0 * 3.1415926535897932384626433832795) / (double)size;
	for (int i=0; i<size; ++i) {
		data[i] += amp * sin(phase);
		phase += w;
	}
}

static void add_wpartial(int size, float *data, double partial, double amp, double phase)
{
	if (amp == 0.0) return;
	int size2 = size >> 1;
	double w = (partial * 2.0 * 3.1415926535897932384626433832795) / (double)size2;
	double cur = amp * sin(phase);
	phase += w;
	for (int i=0; i<size; i+=2) {
		double next = amp * sin(phase);
		data[i] += 2 * cur - next;
		data[i+1] += next - cur;
		cur = next;
		phase += w;
	}
}

static void add_chebyshev(int size, float *data, double partial, double amp)
{
	if (amp == 0.0) return;
	double w = 2.0 / (double)size;
	double phase = -1.0;
	double offset = -amp * cos (partial * pi2);
	for (int i=0; i<size; ++i) {
		data[i] += amp * cos (partial * acos (phase)) + offset;
		phase += w;
	}
}

static void add_wchebyshev(int size, float *data, double partial, double amp)
{
	if (amp == 0.0) return;
	int size2 = size >> 1;
	double w = 2.0 / (double)size2;
	double phase = -1.0;
	double offset = -amp * cos (partial * pi2);
	double cur = amp * cos (partial * acos (phase)) + offset;
	phase += w;
	for (int i=0; i<size; i+=2) {
		double next = amp * cos (partial * acos (phase)) + offset;
		data[i] += 2 * cur - next;
		data[i+1] += next - cur;
		cur = next;
		phase += w;
	}
}

static void cantorFill(int size, float *data) // long offset, double amp)
{
//	if (amp == 0.0) return;
	for (int i=0; i<(size); ++i) {
		int j = i;
		float flag = 1.f;
		while ((j > 0) && (flag == 1.f) ) {
			if (j % 3 == 1) { flag = 0.f; }
			j = j / 3;
		}
		if(flag) { data[i] += 1.f; }
	}
}


enum {
	flag_Normalize = 1,
	flag_Wavetable = 2,
	flag_Clear = 4
};

void ChebyFill(World *world, struct SndBuf *buf, struct sc_msg_iter *msg)
{
	if (buf->channels != 1) return;

	int flags = msg->geti();

	int size = buf->samples;
	int byteSize = size * sizeof(float);
	float *data = (float*)malloc(byteSize);

	if (flags & flag_Clear) Fill(size, data, 0.);
	else memcpy(data, buf->data, byteSize);

	for (int partial=1; msg->remain(); partial++) {
		double amp = msg->getf();
		if (flags & flag_Wavetable) add_wchebyshev(size, data, partial, amp);
		else add_chebyshev(size, data, partial, amp);
	}

	if (flags & flag_Normalize) {
		if (flags & flag_Wavetable) normalize_wsamples(size, data, 1.);
		else normalize_samples(size, data, 1.);
	}

	memcpy(buf->data, data, byteSize);
	free(data);
}


void SineFill1(World *world, struct SndBuf *buf, struct sc_msg_iter *msg)
{
	if (buf->channels != 1) return;

	int flags = msg->geti();

	int size = buf->samples;
	int byteSize = size * sizeof(float);
	float *data = (float*)malloc(byteSize);

	if (flags & flag_Clear) Fill(size, data, 0.);
	else memcpy(data, buf->data, byteSize);

	for (int partial=1; msg->remain(); partial++) {
		double amp = msg->getf();
		if (flags & flag_Wavetable) add_wpartial(size, data, partial, amp, 0.);
		else add_partial(size, data, partial, amp, 0.);
	}

	if (flags & flag_Normalize) {
		if (flags & flag_Wavetable) normalize_wsamples(size, data, 1.);
		else normalize_samples(size, data, 1.);
	}

	memcpy(buf->data, data, byteSize);
	free(data);
}

void SineFill2(World *world, struct SndBuf *buf, struct sc_msg_iter *msg)
{
	if (buf->channels != 1) return;

	int flags = msg->geti();

	int size = buf->samples;
	int byteSize = size * sizeof(float);
	float *data = (float*)malloc(byteSize);

	if (flags & flag_Clear) Fill(size, data, 0.);
	else memcpy(data, buf->data, byteSize);

	while (msg->remain()) {
		double partial = msg->getf();
		double amp = msg->getf();
		if (flags & flag_Wavetable) add_wpartial(size, data, partial, amp, 0.);
		else add_partial(size, data, partial, amp, 0.);
	}
	if (flags & flag_Normalize) {
		if (flags & flag_Wavetable) normalize_wsamples(size, data, 1.);
		else normalize_samples(size, data, 1.);
	}

	memcpy(buf->data, data, byteSize);
	free(data);
}

void SineFill3(World *world, struct SndBuf *buf, struct sc_msg_iter *msg)
{
	if (buf->channels != 1) return;

	int flags = msg->geti();

	int size = buf->samples;
	int byteSize = size * sizeof(float);
	float *data = (float*)malloc(byteSize);

	if (flags & flag_Clear) Fill(size, data, 0.);
	else memcpy(data, buf->data, byteSize);

	while (msg->remain()) {
		double partial = msg->getf();
		double amp = msg->getf();
		double phase = msg->getf();
		if (flags & flag_Wavetable) add_wpartial(size, data, partial, amp, phase);
		else add_partial(size, data, partial, amp, phase);
	}
	if (flags & flag_Normalize) {
		if (flags & flag_Wavetable) normalize_wsamples(size, data, 1.);
		else normalize_samples(size, data, 1.);
	}

	memcpy(buf->data, data, byteSize);
	free(data);
}

void NormalizeBuf(World *world, struct SndBuf *buf, struct sc_msg_iter *msg)
{
	float newmax;
	if(msg->remain() != 0){
		newmax = msg->getf();
	}else{
		newmax = 1.f;
	}
	float *data = buf->data;
	int size = buf->samples;
	normalize_samples(size, data, newmax);
}

void NormalizeWaveBuf(World *world, struct SndBuf *buf, struct sc_msg_iter *msg)
{
	float newmax;
	if(msg->remain() != 0){
		newmax = msg->getf();
	}else{
		newmax = 1.f;
	}
	float *data = buf->data;
	int size = buf->samples;
	normalize_wsamples(size, data, newmax);
}

void CopyBuf(World *world, struct SndBuf *buf, struct sc_msg_iter *msg)
{
	int frames1 = buf->frames;
	int channels1 = buf->channels;

	int toPos = msg->geti();
	uint32 bufnum2 = msg->geti();
	int fromPos = msg->geti();
	int length = msg->geti();

	if (bufnum2 >= world->mNumSndBufs) bufnum2 = 0;
	SndBuf* buf2 = world->mSndBufs + bufnum2;
	int frames2 = buf2->frames;
	int channels2 = buf2->channels;

	if (channels1 != channels2) return;

	fromPos = sc_clip(fromPos, 0, frames2-1);
	toPos = sc_clip(toPos, 0, frames1-1);

	int maxLength = sc_min(frames2 - fromPos, frames1 - toPos);
	if (length < 0) {
		length = maxLength;
	} else {
		length = sc_min(length, maxLength);
	}

	if (length <= 0) return;

	int numbytes = length * sizeof(float) * channels1;
	float *data1 = buf->data + toPos * channels1;
	float *data2 = buf2->data + fromPos * channels2;

	if ((((char*)data1 + numbytes) > (char*)data2) || (((char*)data2 + numbytes) > (char*)data1)) {
		memmove(data1, data2, numbytes);
	} else {
		memcpy(data1, data2, numbytes);
	}
}

void CantorFill(World *world, struct SndBuf *buf, struct sc_msg_iter *msg)
{
	float *data = buf->data;
	int size = buf->samples;
//	double offset = msg->getf();
//	double amp = msg->getf();
//	long offs = (long) offset;
	cantorFill(size, data);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////

PluginLoad(Osc)
{
	ft = inTable;

	DefineSimpleUnit(DegreeToKey);
	DefineSimpleUnit(Select);
	DefineSimpleUnit(TWindex);
	DefineSimpleUnit(Index);
	DefineSimpleUnit(IndexL);
	DefineSimpleUnit(FoldIndex);
	DefineSimpleUnit(WrapIndex);
	DefineSimpleUnit(IndexInBetween);
	DefineSimpleUnit(DetectIndex);
	DefineSimpleUnit(Shaper);
	DefineSimpleUnit(FSinOsc);
	DefineSimpleUnit(PSinGrain);
	DefineSimpleUnit(SinOsc);
	DefineSimpleUnit(SinOscFB);
	DefineSimpleUnit(VOsc);
	DefineSimpleUnit(VOsc3);
	DefineSimpleUnit(Osc);
	DefineSimpleUnit(OscN);
	DefineSimpleUnit(COsc);
	DefineSimpleUnit(Formant);
	DefineSimpleUnit(Blip);
	DefineSimpleUnit(Saw);
	DefineSimpleUnit(Pulse);
	DefineDtorUnit(Klang);
	DefineDtorUnit(Klank);

	DefineBufGen("cheby", ChebyFill);
	DefineBufGen("sine1", SineFill1);
	DefineBufGen("sine2", SineFill2);
	DefineBufGen("sine3", SineFill3);
	DefineBufGen("normalize", NormalizeBuf);
	DefineBufGen("wnormalize", NormalizeWaveBuf);
	DefineBufGen("copy", CopyBuf);
	DefineBufGen("cantorFill", CantorFill);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
