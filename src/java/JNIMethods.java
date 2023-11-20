class JNIMethods {

    private static JNIMethods instance;

    private JNIMethods() {
    }

    public static synchronized JNIMethods getInstance() {
        if (instance == null)
            instance = new JNIMethods();
        return instance;
    }

    private static final String DLL_PATH = System.getProperty("user.dir")
            + "\\plugins\\CISSCO\\Calculate_Magnetic_Moment_3D_Native.dll";
    private static final String SO_PATH = System.getProperty("user.dir")
            + "/ImageJ/plugins/CISSCO/Calculate_Magnetic_Moment_3D_Native.so";

    static {
        Calculate_Magnetic_Moment_3D.logger = new LogManager();

        String os = System.getProperty("os.name").toLowerCase();

        if (os.contains("win")) {

            Calculate_Magnetic_Moment_3D.logger.addInfo("Attempting to load .dll");

            try {
                Calculate_Magnetic_Moment_3D.logger.addVariable("DLL_PATH", DLL_PATH);
                System.load(DLL_PATH);

                Calculate_Magnetic_Moment_3D.logger.addInfo(".dll loaded successfully");
            } catch (Throwable exc) {
                Calculate_Magnetic_Moment_3D.logger.addInfo(".dll loading error message:", exc.toString());
            }

        } else if (os.contains("nix") || os.contains("nux") || os.contains("mac")) {

            Calculate_Magnetic_Moment_3D.logger.addInfo("Attempting to load .so");

            try {
                Calculate_Magnetic_Moment_3D.logger.addVariable("SO_PATH", SO_PATH);
                System.load(SO_PATH);

                Calculate_Magnetic_Moment_3D.logger.addInfo(".so loaded successfully");
            } catch (Throwable exc) {
                Calculate_Magnetic_Moment_3D.logger.addInfo(".so loading error message:", exc.toString());
            }

        } else {
            Calculate_Magnetic_Moment_3D.logger.addInfo("Unknown OS");
        }
    }

    native void passGenSubpixelValues(double jm_R0, int jm_SubPixels, double jm_RCenter, double jm_BackPhase);

    native void passCalcSubCenterValues(int jsmallX, int jsmallDX, int jsmallY, int jsmallDY, int jsmallZ, int jsmallDZ,
            double jm_RCenter, double jm_R0,
            double jclx, double jcly, double jclz,
            double jcmx, double jcmy, double jcmz,
            double jcsx, double jcsy, double jcsz,
            double jcx2, double jcy2, double jcz2);

    native void passSpinDensValues(double jm_cx, double jm_cy, double jm_cz, double jm_r1, double jm_r2, double jm_r3,
            double jm_bkgphase, double jm_magmom);

    native void passMagMomValues(double jm_r1, double jm_r2, double jm_r3,
            double jm_csx, double jm_csy, double jm_csz,
            double jm_R0, double jm_bkg,
            double jm_RChi, double jm_B0, double jm_TEfirst, double jm_snr,
            double jm_in, double jm_mid, double jm_out,
            double jm_e12, double jm_e23);

    native void passSumValues(double jm_ri, double jm_cx, double jm_cy, double jm_cz);

    native void setmVariables(int jm_SubPixels, double jm_R0, double jm_RCenter,
            double jm_CenterX2, double jm_CenterY2, double jm_CenterZ2, double jphaseValue);

    native void setMagMomentVariables(double jSNR, double je12, double je23, double jB0, double jRChi, double jTE);

    native void setMagMoment(double nm);

    native void setBackPhase(double jBackPhase);

    native void setRealImagNumbers(float[][][] jreal, float[][][] jimag);

    native void setPhaseXYMatrix(float[][] jphase);

    native void setPhaseXZMatrix(float[][] jphase);

    native void setMagXYMatrix(float[][] jmag);

    native void setMagXZMatrix(float[][] jmag);

    native void setCenterL(double jx, double jy, double jz);

    native void setCenterM(double jx, double jy, double jz);

    native void setCenterS(double jx, double jy, double jz);

    native void setStep6Variables(double jtef, double jtel, double jb0, double jrx);

    native void setSimulatedMatrices(float[][][] jreal, float[][][] jimag, int jsize);

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

    native String calculateMagneticMoment();

    native double getMR1Calc();

    native double getMR2Calc();

    native double getMR3Calc();

    native double getChi();

    native double getA();

    native double getUncertainty();

    native double getMagMoment();

    native double getRealSum();

    native double getImagSum();
}
