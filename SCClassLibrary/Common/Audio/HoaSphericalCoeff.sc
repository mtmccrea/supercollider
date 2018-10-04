// without mul and add.
HoaSphericalCoeff : UGen {
    *ar { arg l, m, theta = 0.0, phi = 0.0;
        ^this.multiNew('audio', l, m, theta, phi)
    }
    *kr { arg l, m, theta = 0.0, phi = 0.0;
        ^this.multiNew('control', l, m, theta, phi)
    }
}