class JNIMethods {
    private static final String DLL_PATH = System.getProperty("user.dir")
            + "\\plugins\\CISSCO\\Calculate_Magnetic_Moment_3D_Native.dll";

    static {
        Calculate_Magnetic_Moment_3D.logger = new LogManager();
        Calculate_Magnetic_Moment_3D.logger.addInfo("Attempting to load .dll");

        try {
            Calculate_Magnetic_Moment_3D.logger.addVariable("DLL_PATH", DLL_PATH);
            System.load(DLL_PATH);

            Calculate_Magnetic_Moment_3D.logger.addInfo(".dll loaded successfully");
        } catch (Throwable exc) {
            Calculate_Magnetic_Moment_3D.logger.addInfo(".dll loading error message:", exc.toString());
        }
    }

    native void setmVariables(int jm_SubPixels, double jm_R0, double jm_RCenter,
            double jm_CenterX2, double jm_CenterY2, double jm_CenterZ2, double jphaseValue);

    native void setMagMoment(double nm);

    native void setBackPhase(double jBackPhase);

    native void setRealImagNumbers(float[][][] jreal, float[][][] jimag);

    native void setXYZ(double jx, double jy, double jz);

    native void setPhaseXYMatrix(float[][] jphase);

    native void setPhaseXZMatrix(float[][] jphase);

    native void setMagXYMatrix(float[][] jmag);

    native void setMagXZMatrix(float[][] jmag);

    native void setSmallBox(int jxi, int jyi, int jzi, int jxsize, int jysize, int jzsize);

    native void setCenterL(double jx, double jy, double jz);

    native void setCenterM(double jx, double jy, double jz);

    native void setCenterS(double jx, double jy, double jz);

    native void setmR123(double jmr1, double jmr2, double jmr3);

    native void setR123PhaseCalc(double j1, double j2, double j3);

    native void setMagMomentVariables(double jSNR, double je12, double je23, double jB0, double jRChi, double jTE);

    native void setRi(double jri);

    native void setStep6Variables(double jtef, double jtel, double jb0, double jrx);

    native void setSimulatedMatrices(float[][][] jreal, int jsize);

    native void removeBackgroundPhase(double jbPhase);

    native void generateSubpixelArray();

    native String estimateSubpixelCenter();

    native String calculateRealSum();

    native String calculateImagSum();

    native void estBkgAndSpinDensity();

    native void calcSusceptibility();

    native void interpolateVoxelsSIM(int jsize);

    native double SumCircleElementsReal3DSIMMED(int jradii, int jx, int jy, int jz);

    native double equation10(double jp, double jphi_i, double jphi_j);

    native double calculateUncertainty(double je12, double je23);

    native double getSpinDensity();

    native double getBkg();

    native double getSubX();

    native double getSubY();

    native double getSubZ();

    native double getSubXOther();

    native double getSubYOther();

    native double getSubZOther();

    native String calculateMagneticMoment();

    native double getMR1Calc();

    native double getMR2Calc();

    native double getMR3Calc();

    native double getSNR();

    native double getE12();

    native double getE23();

    native double getB0();

    native double getTE();

    native double getRChi();

    native double getRho();

    native double getChi();

    native double getA();

    native double getUncertainty();

    native double getP();

    native double getP0();

    native double getMagMoment();

    native double getResX();

    native double getResY();

    native double getResZ();

    native double getRealSum();

    native double getImagSum();
}
