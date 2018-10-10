// without mul and add.
HoaSphCoeff : UGen {

	*ar { arg l, m, theta = 0.0, phi = 0.0;
        ^this.multiNew('audio', l, m, theta, phi)
    }

	*kr { arg l, m, theta = 0.0, phi = 0.0;
		if (theta.rate == 'audio' or: { phi.rate == 'audio' }) {
			warn(
				format(
					"HoaSphCoeff.kr: argument rate exceeds output rate. theta: %, phi: %",
					theta.rate, phi.rate
				)
			)
		};

        ^this.multiNew('control', l, m, theta, phi)
    }
}