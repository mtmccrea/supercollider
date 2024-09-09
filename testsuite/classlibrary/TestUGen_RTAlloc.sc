TestUGen_RTAlloc : UnitTest {

	classvar server;
	classvar asyncTestCounter = 0;

	*initClass {
		passVerbosity = UnitTest.brief;
	}

	*setUpClass {
		server = Server(this.name);
		server.options.sampleRate = 48000;
		server.options.blockSize = 64;
		server.options.memSize = 2 ** 13; // scsynth default
		// - tests fail with memSize < 256 (GVerb fails)
		// - testing with memSize >= 2 ** 20 would require multiple allocations
		//   (allocation size itself can overflow size_t for some UGens)
		server.bootSync;
	}

	*tearDownClass {
		server.quit.remove;
	}

	// Async helpers

	awaitSynthOutput { |synthFn, runtimeSamples=1, timeout=1|
		var cond = CondVar(), synthOutput = nil;
		var runtimeSeconds = runtimeSamples / server.sampleRate;
		// allow extra timeout for longer tests (e.g. testFFT with big memSize)
		if (runtimeSeconds > timeout) {
			timeout = runtimeSeconds + timeout
		};
		synthFn.loadToFloatArray(runtimeSeconds, server) { |out|
			synthOutput = out;
			cond.signalOne();
		};
		cond.waitFor(timeout) { synthOutput.notNil };
		^synthOutput;
	}

	assertAllocPass { |finishCond, name, func, runtimeSamples|
		fork {
			var out = this.awaitSynthOutput(func, runtimeSamples);
			if (out.isNil) {
				this.assert(false, "% allocPass test should complete".format(name))
			} {
				this.assert(out.sum != 0,
					"% allocPass test should produce non-zero output".format(name)
				);
			};

			asyncTestCounter = asyncTestCounter + 1;
			finishCond.signalOne;
		};
	}

	assertAllocFail { |name, func, runtimeSamples|
		var out = this.awaitSynthOutput(func, runtimeSamples);
		if (out.isNil) {
			this.assert(false, "% allocFail test should complete".format(name))
		} {
			this.assertEquals(out.sum, 0,
				"% allocFail test should produce zero output".format(name)
			);
		};
	}

	// mem calc helpers

	memSizeKb { ^server.options.memSize }
	memSizeFloats { ^this.memSizeKb * (1024 / 4) }
	memSizeSeconds { ^this.memSizeFloats / server.sampleRate }

	// Tests

	test_allocFail_clearUnit {
		var output = this.awaitSynthOutput {
			DelayN.ar(DC.ar(1), this.memSizeSeconds * 2)
		};
		this.assertEquals(output[0], 0.0,
			"a failed RTAlloc should clear unit outputs")
	}

	test_allocFail_continueProcessing_sameNode {
		var output = this.awaitSynthOutput {
			DelayN.ar(DC.ar(1), this.memSizeSeconds * 2); DC.ar(1)
		};
		this.assertEquals(output[0], 1.0,
			"a failed RTAlloc should not block other units in the same node")
	}

	test_allocFail_continueProcessing_otherNodes {
		var output = nil;
		var dummyOutBus = Bus.audio(server, 1);
		var memFill = { DelayN.ar(Silent.ar, this.memSizeSeconds * 2) }.play(server, dummyOutBus);
		server.sync;
		output = this.awaitSynthOutput { DC.ar(1) };
		memFill.free; dummyOutBus.free;
		this.assertEquals(output[0], 1.0,
			"a failed RTAlloc should not block other nodes in the processing chain")
	}

	test_allocFail_setDoneFlag {
		var output = this.awaitSynthOutput {
			K2A.ar(Done.kr(DelayN.ar(Silent.ar, this.memSizeSeconds * 2)))
		};
		this.assert(output[0] > 0,
			"a failed RTAlloc should set the Done flag")
	}

	// TEST ALL UGENS
	// for each UGen in /server/plugins which calls RTAlloc
	// - test successful alloc: should output non-zero values
	// - test failed alloc: should only output zeros

	// NOTE: we try to keep runtime for each test to the bare minimum, which varies
	// from UGen to UGen (some need only 1 sample, some need 4 or more)
	// Runtimes for each UGen are kept the same in both allocPass and allocFails tests

	// FFT UGens needs separate testing (see below)
	// Some UGens were left out (see end of file)

	test_allUGens_allocPass {
		// allocate as little memory as possible for each UGen
		var allocSeconds = 1 / server.sampleRate; // 1 sample
		var blockSize = server.options.blockSize;
		var convBuf, convBufFrames;
		var condVar = CondVar();
		var numTests = 0, waitResult = false;
		var tests, convBufs;

		asyncTestCounter = 0; // reset global counter

		/* Test batch 1: general memory allocation */
		tests = [
			// ["UGenName", funcToTest, runtimeSamples]
			["AllpassC", { AllpassC.ar(DC.ar(1), allocSeconds, 0) }, 4],
			["AllpassL", { AllpassL.ar(DC.ar(1), allocSeconds, 0) }, 4],
			["AllpassN", { AllpassN.ar(DC.ar(1), allocSeconds, 0) }, 4],
			["CombC", { CombC.ar(DC.ar(1), allocSeconds, 0) }, 4],
			["CombL", { CombL.ar(DC.ar(1), allocSeconds, 0) }, 4],
			["CombN", { CombN.ar(DC.ar(1), allocSeconds, 0) }, 4],
			["DelayC", { DelayC.ar(DC.ar(1), allocSeconds, 0) }, 4], // fail on DelayC
			["DelayL", { DelayL.ar(DC.ar(1), allocSeconds, 0) }, 4],
			["DelayN", { DelayN.ar(DC.ar(1), allocSeconds, 0) }, 4],
			["Pluck", { Pluck.ar(DC.ar(1), 1, allocSeconds, 1e-4) }, 4],
			["RunningSum", { RunningSum.ar(DC.ar(1), 1) }, 1],
			["PitchShift", { PitchShift.ar(DC.ar(1), allocSeconds) }, 4],
			["Pitch", { K2A.ar(Pitch.kr(SinOsc.ar)[0]).round }, 2],
			["Limiter", { Limiter.ar(DC.ar(1), dur: allocSeconds) }, 2],
			["Normalizer", { Normalizer.ar(DC.ar(1), dur: allocSeconds) }, 2],

			["GrainBuf", { GrainBuf.ar(1, 1, 1, LocalBuf(1).set(1), maxGrains: 2)}, 1],
			["GrainFM", { GrainFM.ar(1, 1, 1, maxGrains: 2) }, 1],
			["GrainIn", { GrainIn.ar(1, 1, in:DC.ar(1), maxGrains:2) }, 1],
			["GrainSin", { GrainSin.ar(1, 1, maxGrains: 2) }, 1],

			["Gendy1", { Gendy1.ar(initCPs: 1) }, 2],
			["Gendy2", { Gendy2.ar(initCPs: 1) }, 2],
			["Gendy3", { Gendy3.ar(initCPs: 2) }, 2],

			["IFFT", { IFFT(LocalBuf(blockSize).set(1)) }, 2],

			["SpecPcile", { K2A.ar(SpecPcile.kr(LocalBuf(blockSize))) }, 2],

			["GVerb", { GVerb.ar(DC.ar(1), 1, maxroomsize: 1) }, 2],

			// -- ONLY PASS UGENS: couldn't reliably test allocFail on these:
			["IEnvGen", { IEnvGen.ar(Env([1,1]), 0) }, 1],
			["PanAz", { PanAz.ar(1, DC.ar(1)) }, 1],
			["Klank", { Klank.ar(`[100!2, nil, 1], DC.ar(1)) }, 1],
			["Klang", { Klang.ar(`[100!2, nil, 1]) }, 1],
			["Dshuf", { Demand.ar(DC.ar(1),0, Dshuf([1])) }, 1],
			// FFT UGens: only tested for allocPass (hard to make these fail before LocalBuf or FFT fails)
			["Onsets", { K2A.ar(Onsets.kr(FFT(LocalBuf(blockSize*2), Impulse.ar(1)))) }, blockSize],
			["Loudness", { K2A.ar(Loudness.kr(FFT(LocalBuf(blockSize*2), SinOsc.ar))) }, blockSize],
			["BeatTrack", { K2A.ar(BeatTrack.kr(FFT(LocalBuf(blockSize*2)))) }, blockSize],
			["Convolution", { Convolution.ar(DC.ar(1), DC.ar(1), blockSize) }, 1];
		];
		numTests = numTests + tests.size;
		tests.do{ |arr, i|
			this.assertAllocPass(condVar, *arr);
		};

		/* Test batch 2: special test for Convolution UGens: requires buffer alloc */

		convBufFrames = blockSize/2;
		convBufs = 4.collect{
			convBuf = Buffer.loadCollection(server, 1!convBufFrames);
		};
		server.sync;
		tests = [
			// ["UGenName", funcToTest, runtimeSamples]
			["Convolution2", { Convolution2.ar(DC.ar(1), convBufs[0], 1, blockSize) }, 1],
			// Convolution2L onset output is not predictable, might sum to zero sporadically: give extra time
			["Convolution2L", { Convolution2L.ar(DC.ar(1), convBufs[1], 1, blockSize) }, blockSize/2],
			["Convolution3", { Convolution3.ar(DC.ar(1), convBufs[2], 1, blockSize) }, 1],
			["StereoConvolution2L", { StereoConvolution2L.ar(DC.ar(1), convBufs[3], convBufs[3], 1, blockSize) }, 1]
		];
		numTests = numTests + tests.size;
		tests.do{ |arr|
			this.assertAllocPass(condVar, *arr);
		};

		// wait for all above tests to finish
		waitResult = condVar.waitFor(5, { asyncTestCounter == numTests });
		this.assert(waitResult, "allocPass TIMEOUT: First and second test group");
		convBufs.do(_.free);

		/* Test batch 3: special test for PartConv: requires larger buffer alloc */

		convBufFrames = blockSize * 8;
		convBuf = Buffer.loadCollection(server, 1!convBufFrames);
		server.sync;
		numTests = numTests + 1;
		this.assertAllocPass(condVar, "PartConv", { PartConv.ar(DC.ar(1), convBufFrames / 2, convBuf) }, convBufFrames / 4);
		waitResult = condVar.waitFor(5, { asyncTestCounter == numTests });
		this.assert(waitResult, "allocPass TIMEOUT: PartConv");
		convBuf.free;

		asyncTestCounter = 0; // reset global counter
	}

	test_allUGens_allocFail {

		var allocSeconds = this.memSizeSeconds;
		var allocSamples = nextPowerOfTwo(this.memSizeFloats);
		var blockSize = server.options.blockSize;

		this.assertAllocFail("AllpassC", { AllpassC.ar(DC.ar(1), allocSeconds) }, 4);
		this.assertAllocFail("AllpassL", { AllpassL.ar(DC.ar(1), allocSeconds) }, 4);
		this.assertAllocFail("AllpassN", { AllpassN.ar(DC.ar(1), allocSeconds) }, 4);
		this.assertAllocFail("CombC", { CombC.ar(DC.ar(1), allocSeconds, 0) }, 4);
		this.assertAllocFail("CombL", { CombL.ar(DC.ar(1), allocSeconds, 0) }, 4);
		this.assertAllocFail("CombN", { CombN.ar(DC.ar(1), allocSeconds, 0) }, 4);
		this.assertAllocFail("DelayC", { DelayC.ar(DC.ar(1), allocSeconds, 0) }, 4);
		this.assertAllocFail("DelayL", { DelayL.ar(DC.ar(1), allocSeconds, 0) }, 4);
		this.assertAllocFail("DelayN", { DelayN.ar(DC.ar(1), allocSeconds, 0) }, 4);
		this.assertAllocFail("Pluck", { Pluck.ar(DC.ar(1), 1, allocSeconds, 0.001) }, 4);
		this.assertAllocFail("RunningSum", { RunningSum.ar(DC.ar(1), allocSamples) }, 1);
		this.assertAllocFail("PitchShift", { PitchShift.ar(SinOsc.ar, allocSeconds) }, 4);
		// alloc big buffer via very small minFreq
		this.assertAllocFail("Pitch", { K2A.ar(Pitch.kr(SinOsc.ar, minFreq: 1 / allocSeconds)[0]).round }, 2);
		this.assertAllocFail("Limiter", { Limiter.ar(DC.ar(1), dur: allocSeconds) }, 2);
		this.assertAllocFail("Normalizer", { Normalizer.ar(DC.ar(1), dur: allocSeconds) }, 2);

		// grain UGens: allocFail with high maxGrains
		// allocated memory = maxGrains * size(GrainStruct)
		// magic numbers are approximated from the sizes of GrainStruct for each UGen (see GrainUGens.cpp)
		this.assertAllocFail("GrainBuf", { GrainBuf.ar(1, 1, 1, LocalBuf(1).set(1), maxGrains: allocSamples/20)}, 1);
		this.assertAllocFail("GrainFM", { GrainFM.ar(1, 1, 1, maxGrains: allocSamples/20) }, 1);
		this.assertAllocFail("GrainIn", { GrainIn.ar(1, 1, in:DC.ar(1), maxGrains: allocSamples/15) }, 1);
		this.assertAllocFail("GrainSin", { GrainSin.ar(1, 1, maxGrains: allocSamples/15) }, 1);

		this.assertAllocFail("Gendy1", { Gendy1.ar(initCPs: allocSamples) }, 2);
		this.assertAllocFail("Gendy2", { Gendy2.ar(initCPs: allocSamples) }, 2);
		this.assertAllocFail("Gendy3", { Gendy3.ar(initCPs: allocSamples) }, 2);

		// these UGens allocate same size as FFT buf
		// allocating memSizeFloats * 0.75 assures alloc fails in UGen and not in LocalBuf
		this.assertAllocFail("IFFT", { IFFT(LocalBuf(this.memSizeFloats * 0.75)) }, 2);
		this.assertAllocFail("SpecPcile", { K2A.ar(SpecPcile.kr(LocalBuf(this.memSizeFloats * 0.75))) }, 2);

		// GVerb needs more tests, as it can fail multiple ways
		// memSizeSeconds * 340 produces an allocation of memSizeFloats floats
		this.assertAllocFail("GVerb", { GVerb.ar(DC.ar(1), 1, maxroomsize: this.memSizeSeconds * 340) }, 2);

		// Convolution UGens don't need a valid buffer when allocation fails
		this.assertAllocFail("Convolution", { Convolution.ar(DC.ar(1), DC.ar(1), allocSamples) }, 1);
		this.assertAllocFail("Convolution2", { Convolution2.ar(DC.ar(1), -1, 1, allocSamples) }, 1);
		this.assertAllocFail("Convolution2L", { Convolution2L.ar(DC.ar(1), -1, 1, allocSamples) }, blockSize/2);
		this.assertAllocFail("Convolution3", { Convolution3.ar(DC.ar(1), -1, 1, allocSamples) }, 1);
		this.assertAllocFail("StereoConvolution2L", {
			var k = LocalBuf(blockSize).set(1); StereoConvolution2L.ar(DC.ar(1), k, k, 1, allocSamples)
		}, 1);
		this.assertAllocFail("PartConv", { PartConv.ar(DC.ar(1), allocSamples, -1) }, server.options.blockSize * 2);
	}

	// special test for LocalBuf alloc
	test_LocalBuf {
		var out = this.awaitSynthOutput({ K2A.ar(LocalBuf(1)) }, 1);
		if (out.isNil) {
			this.assert(false, "successful alloc test should complete");
		} {
			var bufnum = out[0];
			this.assert(bufnum > -1, "should allocate a valid LocalBuf (bufnum: %)".format(bufnum));
		};

		out = this.awaitSynthOutput({ K2A.ar(LocalBuf(this.memSizeFloats * 10)) }, 1);
		if (out.isNil) {
			this.assert(false, "failed alloc test should complete");
		} {
			this.assertEquals(out[0], -1, "should return -1 on alloc fail");
		};
	}

	// FFT UGens: test RTAlloc success and fail
	// FFT normally outputs a stream of -1. When FFT has processed a frame, which typically takes more than a blockSize,
	// on that block it will output the bufnum instead of -1.
	// (It's a way for PV_UGens to know that a new spectral frame is ready and at the same time get the bufnum)
	// In these tests we are scanning the output of FFT for valid bufnums
	test_FFT {
		var blockSize = server.options.blockSize;
		var fftSize = blockSize;
		var out = this.awaitSynthOutput({ FFT(LocalBuf(fftSize), DC.ar(1), hop:1) }, fftSize);
		if (out.isNil) {
			this.assert(false, "successful alloc test should complete (fftSize: %)".format(fftSize));
		} {
			// FFT returns a stream of -1, except when a spectral frame is ready,
			this.assert(out.any(_ > 0), "should return a valid bufnum (bufnum: %)".format(out));
		};
		// FFT should return always -1 when alloc fails
		// - alloc half of memSize, so that LocalBuf can succeed, but FFT can fail
		fftSize = nextPowerOfTwo(this.memSizeFloats * 0.5);
		// - use hop = blockSize / fftSize to make FFT output frames every blockSize
		out = this.awaitSynthOutput({ FFT(LocalBuf(fftSize), DC.ar(1), hop: blockSize / fftSize) }, blockSize);
		out.postln;
		if (out.isNil) {
			this.assert(false, "failed alloc test should complete (fftSize: %)".format(fftSize));
		} {
			this.assert(out.every(_ == -1), "should only return -1 on alloc fail");
		}
	}

	// PV_UGens
	// we need to call .drop(1) on UGen output when testing
	// because all these UGens just output their input for the first block, without calling _next
	// (thus returning the valid bufnum even for failed allocations)
	test_PVUGens {
		var blockSize = server.options.blockSize;
		var smallBuf = Buffer.alloc(server, blockSize);
		var bigBuf = Buffer.alloc(server, nextPowerOfTwo(this.memSizeFloats) * 2);

		var testAllocPass = { |ugen, testFunc|
			var out = this.awaitSynthOutput(testFunc, blockSize * 4);
			if (out.isNil) {
				this.assert(false, "% successful alloc test should complete".format(ugen));
			} {
				var noneIsInvalid = out.drop(1).any(_ == (-1)).not;
				this.assert(noneIsInvalid,
					"% successful alloc should return a valid bufnum".format(ugen)
				);
			};
		};

		var testAllocFail = { |ugen, testFunc|
			var out = this.awaitSynthOutput(testFunc, blockSize * 4);
			if (out.isNil) {
				this.assert(false, "% failed alloc test should complete".format(ugen));
			} {
				var allInvalid = out.drop(1).every(_ == (-1));
				this.assert(allInvalid,
					"% should only return -1 on alloc fail".format(ugen)
				);
			};
		};

		server.sync;

		[PV_RandComb, PV_Diffuser, PV_MagFreeze, PV_BinScramble].do { |ugen|
			testAllocPass.value(ugen, { ugen.new(smallBuf) });
			testAllocFail.value(ugen, { ugen.new(bigBuf) });
		};

		testAllocPass.value(PV_RandWipe, {PV_RandWipe(smallBuf, smallBuf)});
		testAllocFail.value(PV_RandWipe, {PV_RandWipe(bigBuf, bigBuf)});
	}

}


/*
// These UGens are left out:

// BeatTrack2 requires a special test, doesn't work with .loadToFloatArray (defName too long)

// Onsets, Loudness, BeatTrack, and KeyTrack can't be easily tested for allocFail: it's hard not to make FFT fail first
// KeyTrack is also unreliable for allocPass

// MaxLocalBufs is delicate, can leak when calling .increment alone
// { m = MaxLocalBufs(); 128.do{ m.increment } }.play

// AudioControl, Demand, LocalIn (memory depends on numOutputs)
// LagControl, OffsetOut, Out (memory depends on numInputs)
// difficult to fail alloc before "exceeded number of interconnect buffers"

// Poll, Dpoll, SendReply
// memory depends on string size
// difficult to allocFail before `/s_new out of real time memory`
// {Duty.kr(0.5, 0, Dpoll(Dseq([1], inf), String.fill(400000, $a)))}.play
*/
