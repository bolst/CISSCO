/*
 * @author Nicholas Bolton (bolton21@uwindsor.ca)
 * @version 3.11, Nov. 14, 2022
 */

 // imports
 import java.awt.Color;
 import java.awt.event.ActionEvent;
 import java.awt.event.ActionListener;
 import java.util.ArrayList;
 import java.lang.Math;
 import javax.swing.JButton;
 import javax.swing.JCheckBox;
 import javax.swing.JFileChooser;
 import javax.swing.JFrame;
 import javax.swing.JLabel;
 import javax.swing.JOptionPane;
 import javax.swing.JTextField;
 import javax.swing.SwingConstants;
 import javax.swing.filechooser.FileSystemView;
 import net.miginfocom.swing.MigLayout;
  
 // ImageJ tool imports
 import ij.ImagePlus;
 import ij.gui.*;
 import ij.plugin.PlugIn;
 import ij.plugin.frame.RoiManager;
 import ij.process.*;
 import ij.WindowManager;
 import ij.io.Opener;
 
 public class Calculate_Magnetic_Moment_3D implements PlugIn, ActionListener {
 
  private static ImageItem item;
   private static LogManager logger;
   private static ROIS roiImgMag, roiImgMagXZ, roiImgPhase, roiImgPhaseXZ, roiImgV1SE, roiImgV1SEXZ;
   private static JFrame frame;
   private static JCheckBox chkbx_showrc;
   public static JLabel lbl_stepone, lbl_steptwo, lbl_stepthree, lbl_stepfour, lbl_stepfive, lbl_stepsix, lbl_stepseven,
      lbl_r3phaseAct, lbl_r2phaseAct, lbl_r1phaseAct, lbl_r3phaseCalc, lbl_r2phaseCalc,
       lbl_r1phaseCalc, lbl_r3AphaseUnit, lbl_r2AphaseUnit, lbl_r1AphaseUnit, lbl_r3phaseUnit,
       lbl_r2phaseUnit, lbl_r1phaseUnit, lbl_estBkgPhaseVal, lbl_V0Val, lbl_d_V0Val, lbl_magMom,
       lbl_rho0, lbl_rho0val, lbl_dchi, lbl_dchiVal, lbl_a,
       lbl_aVal, lbl_MPcnt, lbl_Ri, lbl_ImRi, lbl_ReRi, lbl_echoDChi, lbl_aSE, lbl_err, lbl_errVal,
       lbl_gridSizeBase, lbl_gridSize, lbl_rho0SEVal, lbl_M, lbl_rcx, lbl_rcy, lbl_rcz, lbl_rczCorrection, lbl_eqPhase,
       lbl_eqPhaseUnit, lbl_rc, lbl_rcUnit, lbl_spx, lbl_spy, lbl_spz, lbl_spzCorrection, lbl_calculated, lbl_actual, lbl_r1,
       lbl_r1unit, lbl_r1phase, lbl_r2, lbl_r2phase, lbl_r3, lbl_r2unit, lbl_r3unit, lbl_r3phase, lbl_estBkgPhase, lbl_estBkgPhaseUnit,
       lbl_magMomUnit, lbl_snr, lbl_eps12, lbl_eps23, lbl_ReRiVal, lbl_ImRiVal, lbl_TEFirst, lbl_TEFirstUnit, lbl_RChi, lbl_RChiUnit, lbl_sigSE,
       lbl_spinCenter, lbl_innerBrack1, lbl_v1se, lbl_innerBrack2, lbl_v2se, lbl_outerBrack1, lbl_comma11, lbl_comma21, lbl_comma31,
       lbl_v2seYVal1, lbl_comma12, lbl_comma22, lbl_comma32, lbl_innerBrack3, lbl_outerBrack2, lbl_V0, lbl_rho0SE, lbl_V0Unit, lbl_echoDChiVal,
       lbl_aSEVal, lbl_comma23, lbl_comma33, lbl_comma24, lbl_comma34;
   public static JTextField txt_eqPhaseRC, txt_rcx, txt_rcy, txt_rcz, txt_rc, txt_r1,
       txt_r2, txt_r3, txt_spx, txt_spy, txt_spz, txt_snrVal, txt_eps12val, txt_eps23val, B0Text,
       txt_RChiVal, txt_magMomVal, txt_Ri, txt_M, txt_TEFirstVal, TE_lastText, txt_spinCenterXVal, txt_spinCenterYVal, txt_spinCenterZVal,
       txt_v1seXVal1, txt_v1seYVal1, txt_v1seZVal1, txt_v1seXVal2, txt_v1seYVal2, txt_v1seZVal2, txt_v2seXVal1,
       txt_v2seYVal1, txt_v2seZVal1, txt_v2seXVal2, txt_v2seYVal2, txt_v2seZVal2, txt_sigSEVal;
   private static JButton btn_loadImages, btn_estCR, btn_genSubpix, btn_estSubC,
       btn_verifyRadii, btn_removeBkg, btn_estBkgDens, btn_loadTE, btn_unk,
       btn_loadspinecho, btn_estRadSpinEcho, btn_redraw, btn_plotX, btn_plotY, btn_plotZ,
       btn_calcMagMom, btn_loadSimImg, btn_sumRi;
   public static ImagePlus subpixelMagImage, subpixelMagImageXZ, subpixelPhaseImage, subpixelPhaseImageXZ, V1SE_XYImage,
       V1SE_XZImage;
   private static RoiManager rois;
   private static String subCenterErrorMessage,
       subMagTitle, subMagXZTitle, subPhaseTitle, subPhaseXZTitle, V1XY_Title, V1XZ_Title, s1MagWindowTitle,
       s1PhaseWindowTitle, s5MagWindowTitle, s5PhaseWindowTitle, s6MagWindowTitle, s6PhaseWindowTitle, s7WindowTitle;
   private static int drawnRectangle_initialX, drawnRectangle_initialY, drawnRectangle_initialZ, innerBox_initialX,
       innerBox_initialY, innerBox_initialZ, drawnRectangle_sizeX, drawnRectangle_sizeY, drawnRectangle_sizeZ,
       innerBox_sizeX, innerBox_sizeY, innerBox_sizeZ, nC, P1C, P2C, percentOfMagnitudeToNeglect, sizeOfSubpixelImage,
       centerX_subpixelCoordinates, centerY_subpixelCoordinates, centerZ_subpixelCoordinates, m_xCenter_subpixel,
       m_yCenter_subpixel, m_zCenter_subpixel, V1SE_x1, V1SE_x2, V1SE_y1,
       V1SE_y2, V1SE_z1, V1SE_z2, V2SE_x1, V2SE_x2, V2SE_y1, V2SE_y2, V2SE_z1, V2SE_z2, spinEchoImage_XCoord,
       spinEchoImage_YCoord, spinEchoImage_ZCoord;
   private static double phaseValue, averageOfMagnitudeImageRectangleCorners, RCenter, Center_L_x, Center_L_y,
       Center_L_z, Center_M_x, Center_M_y, Center_M_z, center_sx, center_sy, center_sz, m_R0, m_R1, m_R2, m_R3, m_Ri,
       estimatedBackgroundPhase, spinDensity, centerX_pixelCoordinates, centerY_pixelCoordinates,
       centerZ_pixelCoordinates, m_xCenter, m_yCenter, m_zCenter, estimatedPValue, R1PhaseCalc, R2PhaseCalc, R3PhaseCalc,
       R1PhaseActual, R2PhaseActual, R3PhaseActual, m_SNR, m_e12, m_e23, m_B0, m_RChi, m_Si, m_Si2, m_TEFirst, m_TELast,
       m_a, m_da, m_dV0, m_snrStandardDeviation, m_p, m_dp;
   private static double[] xPhaseValues_Positive, yPhaseValues_Positive, zPhaseValues_Positive, xPhaseValues_Negative,
       yPhaseValues_Negative, zPhaseValues_Negative, neglectedPVP, neglectedPVN, PVP1, PVN1, PVP2, PVN2;
   private static float[][] subpixelMagMatrix, subpixelMagMatrixXZ, subpixelPhaseMatrix, subpixelPhaseMatrixXZ;
   private static double[][][] innerBox_containsValuesBelowThreshold;
   private static float[][][] croppedMagnitudeValues3D, croppedPhaseValues3D, croppedRealNumbers3D,
       croppedImaginaryNumbers3D;
 
   public static boolean isNearEdge;
   public static boolean subpixelIsGenerated = false;
   public static boolean estimatedXIsFound = false;
   public static boolean estimatedYIsFound = false;
   public static boolean estimatedZIsFound = false;
   public static boolean rCenterIsFound = false;
   public static boolean estimatedBGIsFound = false;
   public static boolean subpixelCenterIsFound = false;
   public static boolean R1IsFound = false;
   public static boolean R2IsFound = false;
   public static boolean R3IsFound = false;
   public static boolean estimateCenterRadii_isClicked = false;
   public static boolean estimateSubpixelCenter_isClicked = false;
   public static boolean magMomentIsFound = false;
   public static boolean secondImagesAreLoaded = false;
   public static boolean VSEisFound = false;
   public static boolean spinEchoImageIsLoaded = false;
   public static boolean simulatedImagesAreLoaded = false;
 
   private static final int grid = 10;
   private static final double m_ROuterFrom = 0.2;
   private static final double m_RMiddleFrom = 0.9;
   private static final double m_RInnerFrom = 2.5;
   private static final double GAMMARBAR = 42.58;
   private static final String ACCEPTED_FILE_TYPE = "nii";
   private static final String RHO = "\u03C1";
   private static final String LABEL_0 = "\u2080";
   private static final String EPSILON = "\u03B5";
   private static final String CHI = "\u03C7";
   private static final String _DELTA = "\u0394";
   private static final String DELTA = "\u03B4";
   private static final String ITALICIZED_I = "\uD835\uDC8A";
   private static final String _SIGMA = "\u03C3";
   private static final String LABEL_0SE = "0,SE";
   private static final String LABEL_1SE = "1,SE";
   private static final String LABEL_2SE = "2,SE";
   private static final String PLUS_MINUS = "\u00B1";
   private static final String CUBED = "\u00B3";
 
   private static final String DLL_PATH = System.getProperty("user.dir")
       + "\\plugins\\CISSCO\\Calculate_Magnetic_Moment_3D_Native.dll";
 
   // needed to call JNI
   static {
     logger = new LogManager();
     logger.addInfo("Attempting to load .dll");
 
     try {
       logger.addVariable("DLL_PATH", DLL_PATH);
       System.load(DLL_PATH);
 
       logger.addInfo(".dll loaded successfully");
     } catch (Throwable exc) {
       logger.addInfo(".dll loading error message:", exc.toString());
     }
   }
 
   /*
    * The following functions are defined for calling C++ functions. You can look
    * at the C++ code and find these functions where they are prefixed with
    * "Java_Calculate_1Magnetic_1Moment_13D_."
    */
 
   native void setmVariables(int jm_SubPixels, double jm_R0, double jm_RCenter, double jm_CenterX2,
       double jm_CenterY2, double jm_CenterZ2, double jphaseValue);
 
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
 
   public static void main(String[] args) {
   }
 
   // Function to run CMM3D (plug-in)
   public void run(String arg) {
    initialize();
    frame.setVisible(true);
   }
 
   @Override
   public void actionPerformed(ActionEvent e) {
 
     /*
      * If the step 1 open button is clicked - this simply just prompts the user to
      * open a magnitude and phase image
      */
     if (e.getSource() == btn_loadImages) {
 
       clearVariables();
       try {
         // Initializing file choosing window
         final JFileChooser initialFileChooserWindow = new JFileChooser(
             FileSystemView.getFileSystemView().getHomeDirectory());
 
         // Setting title of choosing the magnitude file window
         initialFileChooserWindow.setDialogTitle("CHOOSE MAG FILE");
 
         // MAG FILE SELECTION
         // Shows the save dialog
         int statusOfFileChooser = initialFileChooserWindow.showSaveDialog(null);
 
         // opens phase and mag file
         Opener fileOpener = new Opener();
 
         // if the user selects the file
         if (statusOfFileChooser != JFileChooser.APPROVE_OPTION) {
           JOptionPane.showMessageDialog(frame, "Error: No magnitude file selected.");
           logger.addError("No mag file");
         } else {
 
           // putting chosen file into magText
           String mag = initialFileChooserWindow.getSelectedFile().getAbsolutePath();
           logger.addVariable("mag", mag);
           // String magFile = mag.substring(mag.lastIndexOf("\\") + 1);
           String magFile = initialFileChooserWindow.getSelectedFile().getName();
           s1MagWindowTitle = WindowManager.makeUniqueName(magFile);
           String magFileFormat = mag.substring(mag.lastIndexOf(".") + 1);
 
           if (magFileFormat.compareTo(ACCEPTED_FILE_TYPE) == 0) {
             fileOpener.open(mag);
             WindowManager.getCurrentImage().setTitle(s1MagWindowTitle);
             // magFileIsLoaded = true;
             // magText.setText(mag);
           } else {
             JOptionPane.showMessageDialog(frame, "Error: image type not .nii");
           }
         }
 
         // PHASE FILE SELECTION
         // Shows the save dialog
         initialFileChooserWindow.setDialogTitle("CHOOSE PHASE FILE");
         // if the user selects the file
         statusOfFileChooser = initialFileChooserWindow.showSaveDialog(null);
 
         if (statusOfFileChooser != JFileChooser.APPROVE_OPTION) {
           JOptionPane.showMessageDialog(frame, "Error: No phase file selected.");
         } else {
 
           // putting chosen file into phaseText
           String phase = initialFileChooserWindow.getSelectedFile().getAbsolutePath();
           logger.addVariable("phase", phase);
           // String phaseFile = phase.substring(phase.lastIndexOf("\\") + 1);
           String phaseFile = initialFileChooserWindow.getSelectedFile().getName();
           s1PhaseWindowTitle = WindowManager.makeUniqueName(phaseFile);
           String phaseFileFormat = phase.substring(phase.lastIndexOf(".") + 1);
 
           if (phaseFileFormat.compareTo(ACCEPTED_FILE_TYPE) == 0) {
             fileOpener.open(phase);
             WindowManager.getCurrentImage().setTitle(s1PhaseWindowTitle);
             // phaseFileIsLoaded = true;
             // phaseText.setText(phase);
           } else {
             JOptionPane.showMessageDialog(frame, "Error: image type not .nii");
           }
         }
       } catch (Exception exc) {
         logger.addWarning(exc.toString());
         JOptionPane.showMessageDialog(frame, exc.toString());
       }
     }
 
     /*
      * If the estimate center button is clicked:
      * The program takes the ROI that the user draws.
      * The program finds the minimum voxel value in this box and this is called
      * Center_L.
      * Then, the program takes all the values in the user-drawn ROI that are below
      * the % defined by the user (defaulted to 50). These values are put into a new
      * box, or the inner box. The center of this box is Center_M
      * The program then sums the X Y and Z planes and determines where each plane is
      * minimal. Where the X Y and Z planes are minimal, this is called center_s.
      * This center_s is fed into a function to find "RCenter". This is the estimated
      * radius of the object around the estimated center.
      * The center is estimated along with RCenter by the end of this.
      */
     else if (e.getSource() == btn_estCR) {
 
      // updateVariables();

       try{

        item = new ImageItem(s1MagWindowTitle, 
                             s1PhaseWindowTitle, 
                             Integer.parseInt(txt_M.getText()), 
                             Double.parseDouble(txt_eqPhaseRC.getText()));
        item.calcCenterL();
        item.calcCenterM();
        item.calcCenterS();
        estimatedBackgroundPhase = item.estBkg();
        center_sx = item.centerSX();
        center_sy = item.centerSY();
        center_sz = item.centerSZ();
    
         // Setting center to GUI unless user already has an inputted center point
         if (txt_rcx.getText().isEmpty()) {
           txt_rcx.setText(String.valueOf(center_sx));
         }
         estimatedXIsFound = true;
 
         if (txt_rcy.getText().isEmpty()) {
           txt_rcy.setText(String.valueOf(center_sy));
         }
         estimatedYIsFound = true;
 
         if (txt_rcz.getText().isEmpty()) {
           txt_rcz.setText(String.valueOf(center_sz + 1));
         }
         estimatedZIsFound = true;
 
         // Getting x y and z centers from GUI
         item.setCenterSX(Double.parseDouble(txt_rcx.getText()));
         item.setCenterSY(Double.parseDouble(txt_rcy.getText()));
         item.setCenterSZ(Double.parseDouble(txt_rcz.getText()) - 1);
         m_xCenter = item.centerSX();
         m_yCenter = item.centerSY();
         m_zCenter = item.centerSZ();
 
         // Getting RCenter from estimated xyz
         //RCenter = estimateRCenter((int)(m_xCenter), (int)(m_yCenter), (int)(m_zCenter));
         RCenter = item.estimateRCenter();

         // Updating GUI
         txt_rc.setText(String.valueOf(Math.round(RCenter * 10.0) / 10.0));
         rCenterIsFound = true;

         item.calcR0123();
         m_R0 = item.m_R0();
         m_R1 = item.m_R1();
         m_R2 = item.m_R2();
         m_R3 = item.m_R3();

         setMagMoment(Double.parseDouble(txt_eqPhaseRC.getText()) * Math.pow(RCenter, 3));
 
         // estimateCenterRadii_isClicked = true;

        }catch(Exception exc){
          JOptionPane.showMessageDialog(frame, exc);
         }
       }
 
     /*
      * If the generate subpixel button is clicked:
      * The program takes the values of distance m_R0 from the center and puts it in
      * a box. So the size of the box will be 2 * m_R0 + 1.
      * The box is then interpolated into 1000 subvoxels per voxel. This will be
      * displayed as the subpixel images. One for each orthogonal view (XY and XZ)
      * per image (phase and mag).
      * So since we have the mag and phase matrices, we can calculate the real and
      * imaginary matrices, since real numbers are defined by (mag * cos(phase)) and
      * imaginary numbers are defined by (mag * sin(phase)).
      */
     else if (e.getSource() == btn_genSubpix) {
 
       updateVariables();
 
       // if no RCenter is inputted
       if (txt_rc.getText().isEmpty()) {
         JOptionPane.showMessageDialog(frame, "Error: No RCenter found.");
       }
       // if no XYZ center is inputted
       else if (txt_rcx.getText().isEmpty() || 
                txt_rcy.getText().isEmpty() || 
                txt_rcz.getText().isEmpty()) {
         JOptionPane.showMessageDialog(frame, "Error: No XYZ Center found.");
       }
       // If object is too close to edge of image
       else if (!isNearEdge) {
         JOptionPane.showMessageDialog(frame, "Error: Object too close to edge of image.");
       } else {
 
         // Getting mag and phase images
         ImagePlus magnitudeImage = WindowManager.getImage(s1MagWindowTitle);
         ImagePlus phaseImage = WindowManager.getImage(s1PhaseWindowTitle);
 
         // Size of new images
         sizeOfSubpixelImage = (int) ((2 * m_R0 + 1) * grid);
 
         // Initializing 3D arrays to give to C++
         croppedMagnitudeValues3D = new float[sizeOfSubpixelImage][sizeOfSubpixelImage][sizeOfSubpixelImage];
         croppedPhaseValues3D = new float[sizeOfSubpixelImage][sizeOfSubpixelImage][sizeOfSubpixelImage];
         croppedRealNumbers3D = new float[sizeOfSubpixelImage][sizeOfSubpixelImage][sizeOfSubpixelImage];
         croppedImaginaryNumbers3D = new float[sizeOfSubpixelImage][sizeOfSubpixelImage][sizeOfSubpixelImage];
 
         // Initializing arrays that will be displayed in all 4 images
         subpixelMagMatrix = new float[sizeOfSubpixelImage][sizeOfSubpixelImage];
         subpixelMagMatrixXZ = new float[sizeOfSubpixelImage][sizeOfSubpixelImage];
         subpixelPhaseMatrix = new float[sizeOfSubpixelImage][sizeOfSubpixelImage];
         subpixelPhaseMatrixXZ = new float[sizeOfSubpixelImage][sizeOfSubpixelImage];
 
         // Initial points of subpixel images
         int xi = (int)(m_xCenter) - (int) m_R0;
         int yi = (int)(m_yCenter) - (int) m_R0;
         int zi = (int)(m_zCenter) - (int) m_R0;
         double tempmag;
         double tempphase;
 
         // Getting 3D array data to give to C++ and for images
         for (int k = zi; k <= zi + 2 * (int) m_R0; k++) {
           magnitudeImage.setSlice(k + 1);
           phaseImage.setSlice(k + 1);
           for (int j = yi; j <= yi + 2 * (int) m_R0; j++) {
             for (int i = xi; i <= xi + 2 * (int) m_R0; i++) {
               tempphase = phaseImage.getProcessor().getPixelValue(i, j);
               tempmag = magnitudeImage.getProcessor().getPixelValue(i, j);
 
               croppedMagnitudeValues3D[i - xi][j - yi][k - zi] = (float) tempmag;
               croppedPhaseValues3D[i - xi][j - yi][k - zi] = (float) tempphase;
 
               croppedRealNumbers3D[i - xi][j - yi][k - zi] = (float) (tempmag * Math.cos(tempphase));
               croppedImaginaryNumbers3D[i - xi][j - yi][k - zi] = (float) (tempmag * Math.sin(tempphase));
             }
           }
         }
 
         // Passing necessary variables to C++
         setmVariables(grid, m_R0, RCenter, (int)(m_xCenter), (int)(m_yCenter),
             (int)(m_zCenter), phaseValue);
         setBackPhase(estimatedBackgroundPhase);
         setRealImagNumbers(croppedRealNumbers3D, croppedImaginaryNumbers3D);
 
         // Generate subpixel in C++
         generateSubpixelArray();
 
         int croppedImageSize = 2 * (int) m_R0 + 1;
         int Nfinal = grid;
 
         // Interpolating XY subpixel phase matrix
         for (int m = 0; m < croppedImageSize; m++) {
           for (int k = 0; k < croppedImageSize; k++) {
             for (int i = m * Nfinal; i < (m + 1) * Nfinal; i++) {
               for (int j = k * Nfinal; j < (k + 1) * Nfinal; j++) {
                 subpixelPhaseMatrix[j][i] = croppedPhaseValues3D[k][m][(int) m_R0];
               }
             }
           }
         }
 
         // Interpolating XY subpixel mag matrix
         for (int m = 0; m < croppedImageSize; m++) {
           for (int k = 0; k < croppedImageSize; k++) {
             for (int i = m * Nfinal; i < (m + 1) * Nfinal; i++) {
               for (int j = k * Nfinal; j < (k + 1) * Nfinal; j++) {
                 subpixelMagMatrix[j][i] = croppedMagnitudeValues3D[k][m][(int) m_R0];
               }
             }
           }
         }
 
         // Interpolating XZ subpixel phase matrix
         for (int m = 0; m < croppedImageSize; m++) {
           for (int k = 0; k < croppedImageSize; k++) {
             for (int i = m * Nfinal; i < (m + 1) * Nfinal; i++) {
               for (int j = k * Nfinal; j < (k + 1) * Nfinal; j++) {
                 subpixelPhaseMatrixXZ[j][i] = croppedPhaseValues3D[k][(int) m_R0][m];
               }
             }
           }
         }
 
         // Interpolating XZ subpixel mag matrix
         for (int m = 0; m < croppedImageSize; m++) {
           for (int k = 0; k < croppedImageSize; k++) {
             for (int i = m * Nfinal; i < (m + 1) * Nfinal; i++) {
               for (int j = k * Nfinal; j < (k + 1) * Nfinal; j++) {
                 subpixelMagMatrixXZ[j][i] = croppedMagnitudeValues3D[k][(int) m_R0][m];
               }
             }
           }
         }
 
         // Creating ImageProcessors
         ImageProcessor IP_subpixelMagImage = new FloatProcessor(sizeOfSubpixelImage, sizeOfSubpixelImage);
         ImageProcessor IP_subpixelMagImageXZ = new FloatProcessor(sizeOfSubpixelImage, sizeOfSubpixelImage);
         ImageProcessor IP_subpixelPhaseImage = new FloatProcessor(sizeOfSubpixelImage, sizeOfSubpixelImage);
         ImageProcessor IP_subpixelPhaseImageXZ = new FloatProcessor(sizeOfSubpixelImage, sizeOfSubpixelImage);
 
         // Adding mag and phase data to image processors
         for (int i = 0; i < sizeOfSubpixelImage; i++) {
           for (int j = 0; j < sizeOfSubpixelImage; j++) {
             IP_subpixelMagImage.putPixelValue(i, j, subpixelMagMatrix[i][j]);
             IP_subpixelMagImageXZ.putPixelValue(i, j, subpixelMagMatrixXZ[i][j]);
             IP_subpixelPhaseImage.putPixelValue(i, j, subpixelPhaseMatrix[i][j]);
             IP_subpixelPhaseImageXZ.putPixelValue(i, j, subpixelPhaseMatrixXZ[i][j]);
           }
         }
 
         // Setting mag and phase subpixel image titles
         subMagTitle = "Subpixel Mag Image";
         subMagXZTitle = "Subpixel Mag Image XZ";
         subPhaseTitle = "Subpixel Phase Image";
         subPhaseXZTitle = "Subpixel Phase Image XZ";
 
         // Creating new ImagePlus objects with respective ImageProcessors to be
         // displayed
         subpixelMagImage = new ImagePlus(subMagTitle, IP_subpixelMagImage);
         subpixelMagImageXZ = new ImagePlus(subMagXZTitle, IP_subpixelMagImageXZ);
         subpixelPhaseImage = new ImagePlus(subPhaseTitle, IP_subpixelPhaseImage);
         subpixelPhaseImageXZ = new ImagePlus(subPhaseXZTitle, IP_subpixelPhaseImageXZ);
 
         // Displaying images
         subpixelMagImage.show();
         subpixelMagImageXZ.show();
         subpixelPhaseImage.show();
         subpixelPhaseImageXZ.show();
 
         subpixelIsGenerated = true;
       }
     }
 
     /*
      * If the remove background phase button is clicked:
      * The program takes the subpixel XY and XZ matrices and removes the estimated
      * background phase
      */
     else if (e.getSource() == btn_removeBkg) {
 
       updateVariables();
 
       if (!estimatedBGIsFound) {
         JOptionPane.showMessageDialog(frame, "Error: No background phase found.");
       } else {
         // Removing BG phase in C++
         removeBackgroundPhase(estimatedBackgroundPhase);
 
         // Removing BG phase in Java
         removeBGPhase(subpixelPhaseMatrix);
         removeBGPhase(subpixelPhaseMatrixXZ);
         // removeBGPhase(estimatedBackgroundPhase);
         JOptionPane.showMessageDialog(frame, "Removed Background Phase!");
       }
     }
 
     /*
      * If the estimate subpixel center button is clicked:
      * The program feeds the subpixel matrices into C++ and C++ handles the
      * calculations using Amoeba and other functions from Numerical Recipes
      */
     else if (e.getSource() == btn_estSubC) {
 
       updateVariables();
 
       // condition for program to continue, must have generated subpixel and estimated
       // a center and RCenter
       boolean condition = estimatedXIsFound && estimatedYIsFound && estimatedZIsFound && rCenterIsFound
           && estimatedBGIsFound && subpixelIsGenerated;
 
       if (condition) {
 
         // JOptionPane.showMessageDialog(frame, "Finding subpixel center...");
 
         // Passing necessary data to C++
         setXYZ(Double.parseDouble(txt_rcx.getText()), Double.parseDouble(txt_rcy.getText()),
             Double.parseDouble(txt_rcz.getText()));
         setPhaseXYMatrix(subpixelPhaseMatrix);
         setSmallBox(innerBox_initialX, innerBox_initialY, innerBox_initialZ, innerBox_sizeX, innerBox_sizeY,
             innerBox_sizeZ);
         setCenterL(Center_L_x, Center_L_y, Center_L_z);
         setCenterM(Center_M_x, Center_M_y, Center_M_z);
         setCenterS(center_sx, center_sy, center_sz);
 
         // Calculating subpixel center, if there are no errors then the returned string
         // will be empty
         subCenterErrorMessage = estimateSubpixelCenter();
 
         if (subCenterErrorMessage.compareTo("") == 0) {
 
           // Getting estimated subpixel centers from C++ in terms of pixels
           centerX_pixelCoordinates = getSubX();
           centerY_pixelCoordinates = getSubY();
           centerZ_pixelCoordinates = getSubZ();
 
           // Converting pixel values to subpixel values
           centerX_subpixelCoordinates = (int) pixelToSubpixel(centerX_pixelCoordinates, 0);
           centerY_subpixelCoordinates = (int) pixelToSubpixel(centerY_pixelCoordinates, 1);
           centerZ_subpixelCoordinates = (int) pixelToSubpixel(centerZ_pixelCoordinates, 2);
 
           // Setting pixel values to GUI, because the coordinates are relative to the
           // uploaded images
           txt_spx.setText(String.valueOf(Math.round(centerX_pixelCoordinates * 100.0) / 100.0));
           txt_spy.setText(String.valueOf(Math.round(centerY_pixelCoordinates * 100.0) / 100.0));
           txt_spz.setText(String.valueOf(Math.round(centerZ_pixelCoordinates * 100.0) / 100.0 + 1.0));
           subpixelCenterIsFound = true;
 
           // ---------- begin to put ROIs on images
 
           /*
            * The following is used for displaying the centers and RCenter on the images.
            * The program uses another class I built called ROIS that basically handles
            * ImageJ's RoiManager in a better way. ImageJ makes it very annoying to display
            * multiple ROI's on an image, so this way the code is much easier to read and
            * is not repetitively long. You can read the ROIS Javadoc, as it is already
            * written out, and may have to be compiled to view.
            */
 
           // Creating a new ROIS for the XY mag image
           roiImgMag = new ROIS("MXY");
           // Adding the center as a point to the list
           roiImgMag.addPointROI("Center", centerX_subpixelCoordinates, centerY_subpixelCoordinates);
           if (chkbx_showrc.isSelected()) {
             // If the box is selected, adds RCenter circle ROI to the list
             roiImgMag.addCircleROI("RCenter", centerX_subpixelCoordinates, centerY_subpixelCoordinates, RCenter * 10.0);
           }
           // Displays the ROIS on the image
           roiImgMag.displayROIS();
 
           // The next 3 blocks of code follow the same logic as described above
 
           roiImgMagXZ = new ROIS("MXZ");
           roiImgMagXZ.addPointROI("Center", centerX_subpixelCoordinates, centerZ_subpixelCoordinates);
           if (chkbx_showrc.isSelected()) {
             roiImgMagXZ.addCircleROI("RCenter", centerX_subpixelCoordinates, centerZ_subpixelCoordinates,
                 RCenter * 10.0);
           }
           roiImgMagXZ.displayROIS();
 
           roiImgPhase = new ROIS("PXY");
           roiImgPhase.addPointROI("Center", centerX_subpixelCoordinates, centerY_subpixelCoordinates);
           if (chkbx_showrc.isSelected()) {
             roiImgPhase.addCircleROI("RCenter", centerX_subpixelCoordinates, centerY_subpixelCoordinates,
                 RCenter * 10.0);
           }
           roiImgPhase.displayROIS();
 
           roiImgPhaseXZ = new ROIS("PXZ");
           roiImgPhaseXZ.addPointROI("Center", centerX_subpixelCoordinates, centerZ_subpixelCoordinates);
           if (chkbx_showrc.isSelected()) {
             roiImgPhaseXZ.addCircleROI("RCenter", centerX_subpixelCoordinates, centerZ_subpixelCoordinates,
                 RCenter * 10.0);
           }
           roiImgPhaseXZ.displayROIS();
 
           // ---------- end to put ROIS on images
 
           // Update center_s to be the estimated subpixel center
           center_sx = centerX_pixelCoordinates;
           center_sy = centerY_pixelCoordinates;
           center_sz = centerZ_pixelCoordinates;
 
           estimateSubpixelCenter_isClicked = true;
         } else {
           // Display the error message in a message box if needed
           JOptionPane.showMessageDialog(frame, subCenterErrorMessage);
         }
       } else {
         // Display error message if insufficient data is in program
         JOptionPane.showMessageDialog(frame, "Error: Subpixel images not generated and/or step 2 was not completed");
       }
     }
 
     /*
      * If the redraw center button is clicked:
      * This part of the program is all GUI. This basically just checks the subpixel
      * center and RCenter values in the GUI and updates the subpixel image ROIs to
      * correspond
      */
     else if (e.getSource() == btn_redraw) {
 
       updateVariables();
 
       boolean condition = subpixelIsGenerated && subpixelCenterIsFound;
 
       if (condition) {
 
         // This is basically a repetition from the btn_estSubC code
 
         roiImgMag.clear();
         roiImgMag.addPointROI("Center", centerX_subpixelCoordinates, centerY_subpixelCoordinates);
         if (chkbx_showrc.isSelected()) {
           roiImgMag.addCircleROI("RCenter", centerX_subpixelCoordinates, centerY_subpixelCoordinates, RCenter * 10.0);
         }
         roiImgMag.displayROIS();
 
         roiImgMagXZ.clear();
         roiImgMagXZ.addPointROI("Center", centerX_subpixelCoordinates, centerZ_subpixelCoordinates);
         if (chkbx_showrc.isSelected()) {
           roiImgMagXZ.addCircleROI("RCenter", centerX_subpixelCoordinates, centerZ_subpixelCoordinates, RCenter * 10.0);
         }
         roiImgMagXZ.displayROIS();
 
         roiImgPhase.clear();
         roiImgPhase.addPointROI("Center", centerX_subpixelCoordinates, centerY_subpixelCoordinates);
         if (chkbx_showrc.isSelected()) {
           roiImgPhase.addCircleROI("RCenter", centerX_subpixelCoordinates, centerY_subpixelCoordinates, RCenter * 10.0);
         }
         roiImgPhase.displayROIS();
 
         roiImgPhaseXZ.clear();
         roiImgPhaseXZ.addPointROI("Center", centerX_subpixelCoordinates, centerZ_subpixelCoordinates);
         if (chkbx_showrc.isSelected()) {
           roiImgPhaseXZ.addCircleROI("RCenter", centerX_subpixelCoordinates, centerZ_subpixelCoordinates,
               RCenter * 10.0);
         }
         roiImgPhaseXZ.displayROIS();
       } else {
         JOptionPane.showMessageDialog(frame, "Error: Subpixel center not found");
       }
     }
 
     /*
      * If the plot X button is clicked:
      * This displays a graph of phase values along the x-axis through the subpixel
      * center
      */
     else if (e.getSource() == btn_plotX) {
 
       updateVariables();
 
       boolean condition = subpixelCenterIsFound && subpixelIsGenerated;
 
       if (condition) {
 
         // Initializing plot object for graph
         Plot xPlot = new Plot("X-Profile", "Location", "Phase");
 
         // Initializing ArrayLists for graph values on each axis
         ArrayList<Double> intensity = new ArrayList<Double>();
         ArrayList<Double> location = new ArrayList<Double>();
 
         // Setting colour components of graph
         xPlot.setColor(Color.RED, Color.BLACK);
 
         try {
 
           // Adding phase values to both ArrayLists
           for (int i = (int)(m_xCenter) - (int) m_R0,
               c = 0; i <= (int)(m_xCenter) + (int) m_R0; i++, c++) {
             intensity
                 .add((double) subpixelPhaseMatrix[(int) pixelToSubpixel(i, 0)][(int)m_yCenter_subpixel]);
             location.add((double) (i - (int)(m_xCenter)));
             if (i != (int)(m_xCenter) - (int) m_R0)
               xPlot.drawLine(location.get(c - 1), intensity.get(c - 1), location.get(c), intensity.get(c));
           }
         } catch (Exception exc) {
           JOptionPane.showMessageDialog(frame, exc.toString());
         }
 
         try {
           // Adding points to graph and displaying
           xPlot.addPoints(location, intensity, ij.gui.Plot.DIAMOND);
           xPlot.show();
         } catch (Exception exc) {
           JOptionPane.showMessageDialog(frame, exc.toString());
         }
       } else {
         JOptionPane.showMessageDialog(frame, "Error: Subpixel images not generated");
       }
     }
 
     /*
      * If the plot Y button is clicked:
      * This displays a graph of phase values along the y-axis through the subpixel
      * center
      */
     else if (e.getSource() == btn_plotY) {
 
       updateVariables();
 
       boolean condition = subpixelCenterIsFound && subpixelIsGenerated;
 
       if (condition) {
 
         // Initializing plot object for graph
         Plot yPlot = new Plot("Y-Profile", "Location", "Phase");
 
         // Initializing ArrayLists for graph values on each axis
         ArrayList<Double> intensity = new ArrayList<Double>();
         ArrayList<Double> location = new ArrayList<Double>();
 
         // Setting colour components of graph
         yPlot.setColor(Color.BLUE, Color.BLACK);
 
         try {
 
           // Adding phase values to both ArrayLists
           for (int i = (int)(m_yCenter) - (int) m_R0,
               c = 0; i <= (int)(m_yCenter) + m_R0; i++, c++) {
             intensity
                 .add((double) subpixelPhaseMatrix[(int)m_xCenter_subpixel][(int) pixelToSubpixel(i, 1)]);
             location.add((double) (i - (int)(m_yCenter)));
             if (i != (int)(m_yCenter) - (int) m_R0)
               yPlot.drawLine(location.get(c - 1), intensity.get(c - 1), location.get(c), intensity.get(c));
           }
         } catch (Exception exc) {
           JOptionPane.showMessageDialog(frame, exc.toString());
         }
 
         try {
           // Adding points to graph and displaying
           yPlot.addPoints(location, intensity, ij.gui.Plot.DIAMOND);
           yPlot.show();
         } catch (Exception exc) {
           JOptionPane.showMessageDialog(frame, exc.toString());
         }
       } else {
         JOptionPane.showMessageDialog(frame, "Error: Subpixel images not generated");
       }
     }
 
     /*
      * If the plot Z button is clicked:
      * This displays a graph of phase values along the z-axis through the subpixel
      * center
      */
     else if (e.getSource() == btn_plotZ) {
 
       updateVariables();
 
       boolean condition = subpixelCenterIsFound && subpixelIsGenerated;
 
       if (condition) {
 
         // Initializing plot object for graph
         Plot zPlot = new Plot("Z-Profile", "Location", "Phase");
 
         // Initializing ArrayLists for graph values on each axis
         ArrayList<Double> intensity = new ArrayList<Double>();
         ArrayList<Double> location = new ArrayList<Double>();
 
         // Setting colour components of graph
         zPlot.setColor(Color.GREEN, Color.BLACK);
 
         try {
 
           // Adding phase values to both ArrayLists
           for (int i = (int)(m_zCenter) - (int) m_R0,
               c = 0; i <= (int)(m_zCenter) + m_R0; i++, c++) {
             intensity
                 .add((double) subpixelPhaseMatrixXZ[(int)m_xCenter_subpixel][(int) pixelToSubpixel(i, 2)]);
             location.add((double) (i - (int)(m_zCenter)));
             if (i != (int)(m_zCenter) - (int) m_R0) {
               zPlot.drawLine(location.get(c - 1), intensity.get(c - 1), location.get(c), intensity.get(c));
             }
           }
         } catch (Exception exc) {
           JOptionPane.showMessageDialog(frame, exc.toString());
         }
 
         try {
           // Adding points to graph and displaying
           zPlot.addPoints(location, intensity, ij.gui.Plot.DIAMOND);
           zPlot.show();
         } catch (Exception exc) {
           JOptionPane.showMessageDialog(frame, exc.toString());
         }
       } else {
         JOptionPane.showMessageDialog(frame, "Error: Subpixel images not generated");
       }
     }
 
     /*
      * If the verify radii button is clicked:
      * The three radii (m_R1, m_R2, m_R3) are displayed on the images. The program
      * also averages 4 points for each radius, the points are where the radii
      * intercept the equitorial axis
      */
     else if (e.getSource() == btn_verifyRadii) {
 
       updateVariables();
 
       boolean condition = subpixelIsGenerated && subpixelCenterIsFound;
 
       if (condition) {
 
         // Clearing ROI list
         roiImgMag.clear();
         // Adding center ROI to list
         roiImgMag.addPointROI("Center", centerX_subpixelCoordinates, centerY_subpixelCoordinates);
         // Adding R1 to list
         roiImgMag.addCircleROI("MR1", centerX_subpixelCoordinates, centerY_subpixelCoordinates, m_R1 * 10.0);
         // Adding R2 to list
         roiImgMag.addCircleROI("MR2", centerX_subpixelCoordinates, centerY_subpixelCoordinates, m_R2 * 10.0);
         // Adding R3 to list
         roiImgMag.addCircleROI("MR3", centerX_subpixelCoordinates, centerY_subpixelCoordinates, m_R3 * 10.0);
         // Displaying list
         roiImgMag.displayROIS();
 
         // The same logic as above is followed for the next three blocks
 
         roiImgMagXZ.clear();
         roiImgMagXZ.addPointROI("Center", centerX_subpixelCoordinates, centerZ_subpixelCoordinates);
         roiImgMagXZ.addCircleROI("MR1", centerX_subpixelCoordinates, centerZ_subpixelCoordinates, m_R1 * 10.0);
         roiImgMagXZ.addCircleROI("MR2", centerX_subpixelCoordinates, centerZ_subpixelCoordinates, m_R2 * 10.0);
         roiImgMagXZ.addCircleROI("MR3", centerX_subpixelCoordinates, centerZ_subpixelCoordinates, m_R3 * 10.0);
         roiImgMagXZ.displayROIS();
 
         roiImgPhase.clear();
         roiImgPhase.addPointROI("Center", centerX_subpixelCoordinates, centerY_subpixelCoordinates);
         roiImgPhase.addCircleROI("MR1", centerX_subpixelCoordinates, centerY_subpixelCoordinates, m_R1 * 10.0);
         roiImgPhase.addCircleROI("MR2", centerX_subpixelCoordinates, centerY_subpixelCoordinates, m_R2 * 10.0);
         roiImgPhase.addCircleROI("MR3", centerX_subpixelCoordinates, centerY_subpixelCoordinates, m_R3 * 10.0);
         roiImgPhase.displayROIS();
 
         roiImgPhaseXZ.clear();
         roiImgPhaseXZ.addPointROI("Center", centerX_subpixelCoordinates, centerZ_subpixelCoordinates);
         roiImgPhaseXZ.addCircleROI("MR1", centerX_subpixelCoordinates, centerZ_subpixelCoordinates, m_R1 * 10.0);
         roiImgPhaseXZ.addCircleROI("MR2", centerX_subpixelCoordinates, centerZ_subpixelCoordinates, m_R2 * 10.0);
         roiImgPhaseXZ.addCircleROI("MR3", centerX_subpixelCoordinates, centerZ_subpixelCoordinates, m_R3 * 10.0);
         roiImgPhaseXZ.displayROIS();
 
         float radialPhaseMean_mr1 = 0;
         float radialPhaseMean_mr2 = 0;
         float radialPhaseMean_mr3 = 0;
 
         // Summing up all phase values where the radii and equitorial axis intercept
         switch (item.neglectedAxis()) {
           case "x":
             radialPhaseMean_mr1 += subpixelPhaseImage.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerY_subpixelCoordinates + (int) m_R1 * 10);
             logger.addInfo("Got coordinate (" + String.valueOf(centerX_subpixelCoordinates) + ","
                 + String.valueOf(centerY_subpixelCoordinates + (int) m_R1 * 10) + ")");
             radialPhaseMean_mr1 += subpixelPhaseImage.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerY_subpixelCoordinates - (int) m_R1 * 10);
             logger.addInfo("Got coordinate (" + String.valueOf(centerX_subpixelCoordinates) + ","
                 + String.valueOf(centerY_subpixelCoordinates - (int) m_R1 * 10) + ")");
             radialPhaseMean_mr1 += subpixelPhaseImageXZ.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerZ_subpixelCoordinates + (int) m_R1 * 10);
             logger.addInfo("Got coordinate (" + String.valueOf(centerX_subpixelCoordinates) + ","
                 + String.valueOf(centerZ_subpixelCoordinates + (int) m_R1 * 10) + ")");
             radialPhaseMean_mr1 += subpixelPhaseImageXZ.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerZ_subpixelCoordinates - (int) m_R1 * 10);
             logger.addInfo("Got coordinate (" + String.valueOf(centerX_subpixelCoordinates) + ","
                 + String.valueOf(centerZ_subpixelCoordinates - (int) m_R1 * 10) + ")");
 
             radialPhaseMean_mr2 += subpixelPhaseImage.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerY_subpixelCoordinates + (int) m_R2 * 10);
             logger.addInfo("Got coordinate (" + String.valueOf(centerX_subpixelCoordinates) + ","
                 + String.valueOf(centerY_subpixelCoordinates + (int) m_R2 * 10) + ")");
             radialPhaseMean_mr2 += subpixelPhaseImage.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerY_subpixelCoordinates - (int) m_R2 * 10);
             logger.addInfo("Got coordinate (" + String.valueOf(centerX_subpixelCoordinates) + ","
                 + String.valueOf(centerY_subpixelCoordinates - (int) m_R2 * 10) + ")");
             radialPhaseMean_mr2 += subpixelPhaseImageXZ.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerZ_subpixelCoordinates + (int) m_R2 * 10);
             logger.addInfo("Got coordinate (" + String.valueOf(centerX_subpixelCoordinates) + ","
                 + String.valueOf(centerZ_subpixelCoordinates + (int) m_R2 * 10) + ")");
             radialPhaseMean_mr2 += subpixelPhaseImageXZ.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerZ_subpixelCoordinates - (int) m_R2 * 10);
             logger.addInfo("Got coordinate (" + String.valueOf(centerX_subpixelCoordinates) + ","
                 + String.valueOf(centerZ_subpixelCoordinates - (int) m_R2 * 10) + ")");
 
             radialPhaseMean_mr3 += subpixelPhaseImage.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerY_subpixelCoordinates + (int) m_R3 * 10);
             logger.addInfo("Got coordinate (" + String.valueOf(centerX_subpixelCoordinates) + ","
                 + String.valueOf(centerY_subpixelCoordinates + (int) m_R3 * 10) + ")");
             radialPhaseMean_mr3 += subpixelPhaseImage.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerY_subpixelCoordinates - (int) m_R3 * 10);
             logger.addInfo("Got coordinate (" + String.valueOf(centerX_subpixelCoordinates) + ","
                 + String.valueOf(centerY_subpixelCoordinates - (int) m_R3 * 10) + ")");
             radialPhaseMean_mr3 += subpixelPhaseImageXZ.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerZ_subpixelCoordinates + (int) m_R3 * 10);
             logger.addInfo("Got coordinate (" + String.valueOf(centerX_subpixelCoordinates) + ","
                 + String.valueOf(centerZ_subpixelCoordinates + (int) m_R3 * 10) + ")");
             radialPhaseMean_mr3 += subpixelPhaseImageXZ.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerZ_subpixelCoordinates - (int) m_R3 * 10);
             logger.addInfo("Got coordinate (" + String.valueOf(centerX_subpixelCoordinates) + ","
                 + String.valueOf(centerZ_subpixelCoordinates - (int) m_R3 * 10) + ")");
             break;
 
           case "y":
             radialPhaseMean_mr1 += subpixelPhaseImageXZ.getProcessor()
                 .getPixelValue(centerX_subpixelCoordinates + (int) m_R1 * 10, centerZ_subpixelCoordinates);
             radialPhaseMean_mr1 += subpixelPhaseImageXZ.getProcessor()
                 .getPixelValue(centerX_subpixelCoordinates - (int) m_R1 * 10, centerZ_subpixelCoordinates);
             radialPhaseMean_mr1 += subpixelPhaseImageXZ.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerZ_subpixelCoordinates + (int) m_R1 * 10);
             radialPhaseMean_mr1 += subpixelPhaseImageXZ.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerZ_subpixelCoordinates - (int) m_R1 * 10);
 
             radialPhaseMean_mr2 += subpixelPhaseImageXZ.getProcessor()
                 .getPixelValue(centerX_subpixelCoordinates + (int) m_R2 * 10, centerZ_subpixelCoordinates);
             radialPhaseMean_mr2 += subpixelPhaseImageXZ.getProcessor()
                 .getPixelValue(centerX_subpixelCoordinates - (int) m_R2 * 10, centerZ_subpixelCoordinates);
             radialPhaseMean_mr2 += subpixelPhaseImageXZ.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerZ_subpixelCoordinates + (int) m_R2 * 10);
             radialPhaseMean_mr2 += subpixelPhaseImageXZ.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerZ_subpixelCoordinates - (int) m_R2 * 10);
 
             radialPhaseMean_mr3 += subpixelPhaseImageXZ.getProcessor()
                 .getPixelValue(centerX_subpixelCoordinates + (int) m_R3 * 10, centerZ_subpixelCoordinates);
             radialPhaseMean_mr3 += subpixelPhaseImageXZ.getProcessor()
                 .getPixelValue(centerX_subpixelCoordinates - (int) m_R3 * 10, centerZ_subpixelCoordinates);
             radialPhaseMean_mr3 += subpixelPhaseImageXZ.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerZ_subpixelCoordinates + (int) m_R3 * 10);
             radialPhaseMean_mr3 += subpixelPhaseImageXZ.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerZ_subpixelCoordinates - (int) m_R3 * 10);
             break;
 
           case "z":
             radialPhaseMean_mr1 += subpixelPhaseImage.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerY_subpixelCoordinates + (int) m_R1 * 10);
             radialPhaseMean_mr1 += subpixelPhaseImage.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerY_subpixelCoordinates - (int) m_R1 * 10);
             radialPhaseMean_mr1 += subpixelPhaseImage.getProcessor()
                 .getPixelValue(centerX_subpixelCoordinates + (int) m_R1 * 10, centerY_subpixelCoordinates);
             radialPhaseMean_mr1 += subpixelPhaseImage.getProcessor()
                 .getPixelValue(centerX_subpixelCoordinates - (int) m_R1 * 10, centerY_subpixelCoordinates);
 
             radialPhaseMean_mr2 += subpixelPhaseImage.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerY_subpixelCoordinates + (int) m_R2 * 10);
             radialPhaseMean_mr2 += subpixelPhaseImage.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerY_subpixelCoordinates - (int) m_R2 * 10);
             radialPhaseMean_mr2 += subpixelPhaseImage.getProcessor()
                 .getPixelValue(centerX_subpixelCoordinates + (int) m_R2 * 10, centerY_subpixelCoordinates);
             radialPhaseMean_mr2 += subpixelPhaseImage.getProcessor()
                 .getPixelValue(centerX_subpixelCoordinates - (int) m_R2 * 10, centerY_subpixelCoordinates);
 
             radialPhaseMean_mr3 += subpixelPhaseImage.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerY_subpixelCoordinates + (int) m_R3 * 10);
             radialPhaseMean_mr3 += subpixelPhaseImage.getProcessor().getPixelValue(centerX_subpixelCoordinates,
                 centerY_subpixelCoordinates - (int) m_R3 * 10);
             radialPhaseMean_mr3 += subpixelPhaseImage.getProcessor()
                 .getPixelValue(centerX_subpixelCoordinates + (int) m_R3 * 10, centerY_subpixelCoordinates);
             radialPhaseMean_mr3 += subpixelPhaseImage.getProcessor()
                 .getPixelValue(centerX_subpixelCoordinates - (int) m_R3 * 10, centerY_subpixelCoordinates);
             break;
 
           default:
             JOptionPane.showMessageDialog(frame, "Error in Verify Radii: No MRI field axis found");
             break;
         }
 
         // Dividing each phase sum by 4 for average
         radialPhaseMean_mr1 /= 4.0;
         radialPhaseMean_mr2 /= 4.0;
         radialPhaseMean_mr3 /= 4.0;
         // Removing background phase off each phase value
         radialPhaseMean_mr1 -= estimatedBackgroundPhase;
         radialPhaseMean_mr2 -= estimatedBackgroundPhase;
         radialPhaseMean_mr3 -= estimatedBackgroundPhase;
 
         R1PhaseActual = (double) radialPhaseMean_mr1;
         R2PhaseActual = (double) radialPhaseMean_mr2;
         R3PhaseActual = (double) radialPhaseMean_mr3;
 
         lbl_r1phaseAct.setText(String.valueOf(Math.round(radialPhaseMean_mr1 * 100.0) / 100.0));
         lbl_r2phaseAct.setText(String.valueOf(Math.round(radialPhaseMean_mr2 * 100.0) / 100.0));
         lbl_r3phaseAct.setText(String.valueOf(Math.round(radialPhaseMean_mr3 * 100.0) / 100.0));
       } else {
         JOptionPane.showMessageDialog(frame, "Error: Insufficient data to verify radii");
       }
     }
 
     /*
      * If the estimate bkg and spin density button is clicked:
      * The program calls C++ to calculate the background phase and spin density
      * based off of the three radii and their phase values
      */
     else if (e.getSource() == btn_estBkgDens) {
 
       updateVariables();
 
       // Condition for code to run - all radii must be found and the subpixel images
       // and center must be found
       boolean condition = R1IsFound && R2IsFound && R3IsFound && subpixelCenterIsFound && subpixelIsGenerated;
 
       if (condition) {
 
         // Calling C++ to calculate background phase and spin density
         estBkgAndSpinDensity();
 
         // Getting calculated background phase from C++
         estimatedBackgroundPhase = Math.abs(getBkg());
 
         // Getting calculated spin density from C++
         spinDensity = getSpinDensity();
 
         logger.addVariable("estimatedBackgroundPhase", estimatedBackgroundPhase);
         logger.addVariable("spinDensity", spinDensity);
 
         // Setting background phase and spin density to GUI
         lbl_estBkgPhaseVal.setText(String.valueOf(Math.round(estimatedBackgroundPhase * 100.0) / 100.0));
         lbl_rho0val.setText(String.valueOf(Math.round(spinDensity * 100.0) / 100.0));
       } else {
         JOptionPane.showMessageDialog(frame, "Error: Insufficient data to calculate spin density and background phase");
       }
     }
 
     /*
      * If the load simulated images button is clicked:
      * The program prompts the user to open a new set of simulated images. These are
      * used for finding the uncertainty of the magnetic moment and for calculating
      * eij
      */
     else if (e.getSource() == btn_loadSimImg) {
       // Begin loading simulated images
 
       if (WindowManager.getImage(s5MagWindowTitle) != null) {
         WindowManager.getImage(s5MagWindowTitle).close();
       }
 
       if (WindowManager.getImage(s5PhaseWindowTitle) != null) {
         WindowManager.getImage(s5PhaseWindowTitle).close();
       }
 
       updateVariables();
 
       // ---------- Begin to open simulated images, similar to step 1
 
       final JFileChooser simulatedImageChooserWindow = new JFileChooser(
           FileSystemView.getFileSystemView().getHomeDirectory());
       simulatedImageChooserWindow.setDialogTitle("CHOOSE SIMULATED MAG FILE");
       int statusOfFileChooser = simulatedImageChooserWindow.showSaveDialog(null);
 
       Opener fileOpener = new Opener();
 
       if (statusOfFileChooser != JFileChooser.APPROVE_OPTION) {
         JOptionPane.showMessageDialog(frame, "Error: no file selected");
       } else {
         String simMagImageName = simulatedImageChooserWindow.getSelectedFile().getAbsolutePath();
         // String simMagImageFile =
         // simMagImageName.substring(simMagImageName.lastIndexOf("\\") + 1);
         String simMagImageFile = simulatedImageChooserWindow.getSelectedFile().getName();
         s5MagWindowTitle = WindowManager.makeUniqueName(simMagImageFile);
         String simMagImageFileFormat = simMagImageName.substring(simMagImageName.lastIndexOf(".") + 1);
 
         if (simMagImageFileFormat.compareTo(ACCEPTED_FILE_TYPE) == 0) {
           fileOpener.open(simMagImageName);
           WindowManager.getCurrentImage().setTitle(s5MagWindowTitle);
         } else {
           JOptionPane.showMessageDialog(frame, "Error: image type not .nii");
         }
       }
 
       simulatedImageChooserWindow.setDialogTitle("CHOOSE SIMULATED PHASE FILE");
       statusOfFileChooser = simulatedImageChooserWindow.showSaveDialog(null);
 
       if (statusOfFileChooser != JFileChooser.APPROVE_OPTION) {
         JOptionPane.showMessageDialog(frame, "Error: no file selected");
       } else {
         String simPhaseImageName = simulatedImageChooserWindow.getSelectedFile().getAbsolutePath();
         // String simPhaseImageFile =
         // simPhaseImageName.substring(simPhaseImageName.lastIndexOf("\\") + 1);
         String simPhaseImageFile = simulatedImageChooserWindow.getSelectedFile().getName();
         s5PhaseWindowTitle = WindowManager.makeUniqueName(simPhaseImageFile);
         String simPhaseImageFileFormat = simPhaseImageName.substring(simPhaseImageName.lastIndexOf(".") + 1);
 
         if (simPhaseImageFileFormat.compareTo(ACCEPTED_FILE_TYPE) == 0) {
           fileOpener.open(simPhaseImageName);
           WindowManager.getCurrentImage().setTitle(s5PhaseWindowTitle);
         } else {
           JOptionPane.showMessageDialog(frame, "Error: image type not .nii");
         }
       }
 
       // ---------- End to open simulated images
 
       simulatedImagesAreLoaded = (WindowManager.getImage(s5MagWindowTitle) != null)
           && (WindowManager.getImage(s5PhaseWindowTitle) != null);
 
       // End loading simulated images
 
       if (simulatedImagesAreLoaded) {
 
         ImagePlus simulatedMagnitudeImage = WindowManager.getImage(s5MagWindowTitle);
         ImagePlus simulatedPhaseImage = WindowManager.getImage(s5PhaseWindowTitle);
 
         if (simulatedMagnitudeImage.getNSlices() != simulatedPhaseImage.getNSlices()) {
           JOptionPane.showMessageDialog(frame, "Error: different # of slices in images");
         } else if (simulatedMagnitudeImage.getHeight() != simulatedPhaseImage.getHeight()) {
           JOptionPane.showMessageDialog(frame, "Error: different height in images");
         } else if (simulatedMagnitudeImage.getWidth() != simulatedPhaseImage.getWidth()) {
           JOptionPane.showMessageDialog(frame, "Error: different width in images");
         } else if (!R1IsFound) {
           JOptionPane.showMessageDialog(frame, "Error: no m_R1 found");
         } else {
 
           // Size of the matrix that will be fetched from the simulated images
           int sizeOfSimmedMatrices = 2 * ((int) Math.ceil(m_R1) + 1);
 
           // Initializing matrix for simulated images
           float[][][] simulatedRealNumbers = new float[sizeOfSimmedMatrices + 1][sizeOfSimmedMatrices
               + 1][sizeOfSimmedMatrices + 1];
 
           // Getting center of image, since in simulated images the center of the image is
           // the center of the object
           int simmedCenterX = simulatedMagnitudeImage.getWidth() / 2;
           int simmedCenterY = simulatedMagnitudeImage.getHeight() / 2;
           int simmedCenterZ = simulatedMagnitudeImage.getNSlices() / 2;
 
           // Getting bounds for the upcoming for loop: a box around the centers found
           // above
           int k_i = simmedCenterZ - sizeOfSimmedMatrices / 2;
           int k_f = simmedCenterZ + sizeOfSimmedMatrices / 2;
           int j_i = simmedCenterY - sizeOfSimmedMatrices / 2;
           int j_f = simmedCenterY + sizeOfSimmedMatrices / 2;
           int i_i = simmedCenterX - sizeOfSimmedMatrices / 2;
           int i_f = simmedCenterX + sizeOfSimmedMatrices / 2;
 
           // Variables to hold each pixel's mag and phase value, just so code is read
           // easier
           double tempmag;
           double tempphase;
 
           // Loop to get real numbers from simulated images using (mag * cos(phase))
           for (int k = k_i, nk = 0; k <= k_f; k++, nk++) {
             simulatedMagnitudeImage.setSlice(k + 1);
             simulatedPhaseImage.setSlice(k + 1);
             for (int j = j_i, nj = 0; j <= j_f; j++, nj++) {
               for (int i = i_i, ni = 0; i <= i_f; i++, ni++) {
                 tempmag = simulatedMagnitudeImage.getProcessor().getPixelValue(i, j);
                 tempphase = simulatedPhaseImage.getProcessor().getPixelValue(i, j);
                 simulatedRealNumbers[ni][nj][nk] = (float) (tempmag * Math.cos(tempphase));
               }
             }
           }
 
           logger.addInfo("Finished loops");
           logger.addVariable("simmedCenter", simmedCenterX);
           logger.addVariable("m_R1", m_R1);
           logger.addVariable("m_R2", m_R2);
           logger.addVariable("m_R3", m_R3);
 
           // Giving matrix to C++
           setSimulatedMatrices(simulatedRealNumbers, sizeOfSimmedMatrices + 1);
           logger.addInfo("Set matrices to C++");
           // Interpolating matrix in C++
           interpolateVoxelsSIM(sizeOfSimmedMatrices * 10);
           logger.addInfo("Interpolation complete");
 
           // Getting current rho_0 value from GUI
           double rho_0 = Double.parseDouble(lbl_rho0val.getText());
           // double rho_0 = 10;
 
           // Getting center of simulated matrix
           int subCenterSIM = sizeOfSimmedMatrices * 10 / 2;
 
           // Summing real values in each radius in C++
           double S1 = SumCircleElementsReal3DSIMMED((int) (m_R1 * 10.0), subCenterSIM, subCenterSIM, subCenterSIM);
           double S2 = SumCircleElementsReal3DSIMMED((int) (m_R2 * 10.0), subCenterSIM, subCenterSIM, subCenterSIM);
           double S3 = SumCircleElementsReal3DSIMMED((int) (m_R3 * 10.0), subCenterSIM, subCenterSIM, subCenterSIM);
 
           // Getting phi (phase) values from GUI
           double phi1 = Double.parseDouble(lbl_r1phaseCalc.getText());
           double phi2 = Double.parseDouble(lbl_r2phaseCalc.getText());
           double phi3 = Double.parseDouble(lbl_r3phaseCalc.getText());
 
           // Calculating fij using equation 9 (real)
           double Re_f12_REAL = (9.0 * Math.sqrt(3)) / (4.0 * Math.PI * rho_0) * (S1 - S2);
           double Re_f23_REAL = (9.0 * Math.sqrt(3)) / (4.0 * Math.PI * rho_0) * (S2 - S3);
 
           // Calculating fij using equation 10 in C++ (theoretical)
           double Re_f12_THEORETICAL = equation10(m_p, phi1, phi2);
           double Re_f23_THEORETICAL = equation10(m_p, phi2, phi3);
 
           logger.addVariable("S1", S1);
           logger.addVariable("S2", S2);
           logger.addVariable("S3", S3);
           logger.addVariable("Re_f12_REAL", Re_f12_REAL);
           logger.addVariable("Re_f23_REAL", Re_f23_REAL);
           logger.addVariable("Re_f12_THEORETICAL", Re_f12_THEORETICAL);
           logger.addVariable("Re_f23_THEORETICAL", Re_f23_THEORETICAL);
 
           // Calculating eij by doing [(real - theoretical) / theoretical]
           double e12 = Math.abs(Re_f12_REAL - Re_f12_THEORETICAL) / Re_f12_THEORETICAL;
           double e23 = Math.abs(Re_f23_REAL - Re_f23_THEORETICAL) / Re_f23_THEORETICAL;
 
           // Setting eij to GUI
           txt_eps12val.setText(String.valueOf(Math.round(e12 * 100.0) / 100.0));
           txt_eps23val.setText(String.valueOf(Math.round(e23 * 100.0) / 100.0));
 
           // Calculating uncertainty in C++
           double uncertainty = calculateUncertainty(e12, e23);
 
           logger.addVariable("uncertainty", uncertainty);
 
           // Setting uncertainty to GUI
           lbl_errVal.setText(String.valueOf(Math.round(uncertainty * 100.0) / 100.0));
 
         }
       }
     }
 
     /*
      * If the calculate magnetic moment button is clicked:
      * C++ calculates magnetic moment
      */
     else if (e.getSource() == btn_calcMagMom) {
 
       updateVariables();
 
       // Condition for code to run - radii, RCenter, subpixel center and images must
       // be found
       boolean condition = R1IsFound && R2IsFound && R3IsFound && subpixelCenterIsFound && subpixelIsGenerated
           && rCenterIsFound;
 
       if (condition) {
         String errorMessage_Mag;
 
         // Calculating magnetic moment in C++, if function returns empty string then no
         // error message
         errorMessage_Mag = calculateMagneticMoment();
 
         if (errorMessage_Mag.compareTo("") == 0) {
           // Getting various values from C++ as a result of the function to calculate the
           // magnetic moment being ran
           lbl_r1phaseCalc.setText(String.valueOf(Math.round(getMR1Calc() * 100.0) / 100.0));
           lbl_r2phaseCalc.setText(String.valueOf(Math.round(getMR2Calc() * 100.0) / 100.0));
           lbl_r3phaseCalc.setText(String.valueOf(Math.round(getMR3Calc() * 100.0) / 100.0));
           txt_magMomVal.setText(String.valueOf(Math.round(getMagMoment() * 100.0) / 100.0));
           if (getUncertainty() == -1.0) {
             lbl_errVal.setText("");
             JOptionPane.showMessageDialog(frame,
                 "Error: Cannot calculate error\nMake sure SNR, " + EPSILON + "12 and " + EPSILON + "23 are set.");
           } else {
             lbl_errVal.setText(String.valueOf(Math.round(getUncertainty() * 100.0) / 100.0));
           }
           lbl_dchiVal.setText(String.valueOf(Math.round(getChi() * 100.0) / 100.0) + " ppm");
           lbl_aVal.setText(String.valueOf(Math.round(getA() * 100.0) / 100.0) + " pixels");
           lbl_rho0val.setText(String.valueOf(Math.round(getSpinDensity() * 100.0) / 100.0));
           magMomentIsFound = true;
         } else {
           JOptionPane.showMessageDialog(frame, errorMessage_Mag);
         }
       } else {
         JOptionPane.showMessageDialog(frame, "Error: insufficient data to calculate magnetic moment");
       }
     }
 
     /*
      * If the sum button is clicked:
      * C++ calculates the sum of the real and imag numbers within the user-inputted
      * radius on the subpixel images
      */
     else if (e.getSource() == btn_sumRi) {
       updateVariables();
 
       String Imag_errmsg;
       String Real_errmsg;
 
       Imag_errmsg = calculateImagSum();
       Real_errmsg = calculateRealSum();
 
       if (Imag_errmsg.compareTo("") == 0) {
         m_Si = getImagSum();
         lbl_ImRi.setText("S" + ITALICIZED_I + "= " + String.valueOf(Math.round(m_Si * 100.0) / 100.0));
       } else {
         JOptionPane.showMessageDialog(frame, Imag_errmsg);
       }

       if (Real_errmsg.compareTo("") == 0) {
        m_Si2 = getRealSum();
        lbl_ReRi.setText("S" + ITALICIZED_I + "= " + String.valueOf(Math.round(m_Si2 * 100.0) / 100.0));
      } else {
        JOptionPane.showMessageDialog(frame, Real_errmsg);
      }
     }

 
     /*
      * If the load first TE images button is clicked:
      * This idea is basically to get the same object that was uploaded in step 1 but
      * with a smaller echo time. Since we have the equation
      * (TE_First / TE_Last = p_first / p_last), we can use this to solve for unknown
      * variables such as a, g, or delta chi
      */
     /*
      * TODO: Not finished. Currently the code only gets the images but does not do
      * anything with them. We need to consider the cases where echo times are in the
      * originally uploaded image (and ImageJ can recognize these as "hyperstacks")
      * and when they need to be uploaded separately.
      */
     else if (e.getSource() == btn_loadTE) {
 
       // If images are already open, close them
       if (WindowManager.getImage(s6MagWindowTitle) != null) {
         WindowManager.getImage(s6MagWindowTitle).close();
       }
 
       if (WindowManager.getImage(s6PhaseWindowTitle) != null) {
         WindowManager.getImage(s6PhaseWindowTitle).close();
       }
 
       ImagePlus magnitudeImage = WindowManager.getImage(s1MagWindowTitle);
       ImagePlus phaseImage = WindowManager.getImage(s1PhaseWindowTitle);
 
       // ---------- Begin to get first TE images
 
       if (magnitudeImage.isHyperStack() && phaseImage.isHyperStack()) {
         secondImagesAreLoaded = true;
 
         // TODO: handle if second echo is within image that is uploaded in step 1
       } else {
         JOptionPane.showMessageDialog(frame,
             "Error: Image does not contain hyperstacks.\nSelect another image for second echo time.");
 
         try {
 
           final JFileChooser secondEchoChooserWindow = new JFileChooser(
               FileSystemView.getFileSystemView().getHomeDirectory());
           secondEchoChooserWindow.setDialogTitle("CHOOSE MAG FILE");
           int statusOfFileChooser = secondEchoChooserWindow.showSaveDialog(null);
 
           Opener fileOpener = new Opener();
 
           if (statusOfFileChooser != JFileChooser.APPROVE_OPTION) {
             JOptionPane.showMessageDialog(frame, "Error: No magnitude file selected");
             secondImagesAreLoaded = false;
           } else {
             String mag_secondEcho = secondEchoChooserWindow.getSelectedFile().getAbsolutePath();
             // String mag_secondEchoFile =
             // mag_secondEcho.substring(mag_secondEcho.lastIndexOf("\\") + 1);
             String mag_secondEchoFile = secondEchoChooserWindow.getSelectedFile().getName();
             s6MagWindowTitle = WindowManager.makeUniqueName(mag_secondEchoFile);
             String mag_secondEchoFileFormat = mag_secondEcho.substring(mag_secondEcho.lastIndexOf(".") + 1);
 
             if (mag_secondEchoFileFormat.compareTo(ACCEPTED_FILE_TYPE) == 0) {
               fileOpener.open(mag_secondEcho);
               WindowManager.getCurrentImage().setTitle(s6MagWindowTitle);
               // mag_secondEchoIsLoaded = true;
             } else {
               JOptionPane.showMessageDialog(frame, "Error: image type not .nii");
             }
             secondImagesAreLoaded = true;
           }
 
           secondEchoChooserWindow.setDialogTitle("CHOOSE PHASE FILE");
           statusOfFileChooser = secondEchoChooserWindow.showSaveDialog(null);
 
           if (statusOfFileChooser != JFileChooser.APPROVE_OPTION) {
             JOptionPane.showMessageDialog(frame, "Error: No phase file selected");
             secondImagesAreLoaded = false;
           } else {
             String phase_secondEcho = secondEchoChooserWindow.getSelectedFile().getAbsolutePath();
             // String phase_secondEchoFile =
             // phase_secondEcho.substring(phase_secondEcho.lastIndexOf("\\") + 1);
             String phase_secondEchoFile = secondEchoChooserWindow.getSelectedFile().getName();
             s6PhaseWindowTitle = WindowManager.makeUniqueName(phase_secondEchoFile);
             String phase_secondEchoFileFormat = phase_secondEcho.substring(phase_secondEcho.lastIndexOf(".") + 1);
 
             if (phase_secondEchoFileFormat.compareTo(ACCEPTED_FILE_TYPE) == 0) {
               fileOpener.open(phase_secondEcho);
               WindowManager.getCurrentImage().setTitle(s6PhaseWindowTitle);
               // mag_secondEchoIsLoaded = true;
             } else {
               JOptionPane.showMessageDialog(frame, "Error: image type not .nii");
             }
             secondImagesAreLoaded = true;
           }
         } catch (Exception exc) {
           JOptionPane.showMessageDialog(frame, exc.toString());
         }
       }
 
       // ---------- end to get first TE images
     }
 
     /*
      * If the step 6 (not named yet) button is clicked:
      * The program calculates the susceptibility in C++. It does this by using the
      * scaled image's TE values to estimate a magnetic moment in order to estimate a
      * and dChi
      */
     /*
      * TODO: Not finished. This kind of goes with the load step 6 image button.
      * There should be some link between the two, where there is currently not. Also
      * code was not fully tested and was not confirmed to be working properly.
      */
     else if (e.getSource() == btn_unk) {
       if (secondImagesAreLoaded) {
         updateVariables();
 
         // Calculating susceptibility in C++
         calcSusceptibility();
 
         // Getting values from C++ and setting to GUI
         lbl_aVal.setText(String.valueOf(Math.round(getA() * 100.0) / 100.0) + " pixels");
         lbl_dchiVal.setText(String.valueOf(Math.round(getChi() * 100.0) / 100.0) + " ppm");
       }
     }
 
     /*
      * If the load spin echo image button is selected:
      * The program gets the user to open a spin echo image (see section 2.6 of
      * smoment.pdf for more detail)
      */
     else if (e.getSource() == btn_loadspinecho) {
 
       if (WindowManager.getImage(s7WindowTitle) != null) {
         WindowManager.getImage(s7WindowTitle).close();
       }
 
       updateVariables();
 
       // Begin opening spin echo image
       final JFileChooser spinEchoChooserWindow = new JFileChooser(
           FileSystemView.getFileSystemView().getHomeDirectory());
       spinEchoChooserWindow.setDialogTitle("CHOOSE MAG FILE");
       int statusOfFileChooser = spinEchoChooserWindow.showSaveDialog(null);
 
       Opener fileOpener = new Opener();
 
       if (statusOfFileChooser != JFileChooser.APPROVE_OPTION) {
         JOptionPane.showMessageDialog(frame, "Error: No file selected");
         spinEchoImageIsLoaded = false;
       } else {
         String spinEchoImage = spinEchoChooserWindow.getSelectedFile().getAbsolutePath();
         // String spinEchoImageFile =
         // spinEchoImage.substring(spinEchoImage.lastIndexOf("\\") + 1);
         String spinEchoImageFile = spinEchoChooserWindow.getSelectedFile().getName();
         s7WindowTitle = WindowManager.makeUniqueName(spinEchoImageFile);
         String spinEchoImageFileFormat = spinEchoImage.substring(spinEchoImage.lastIndexOf(".") + 1);
 
         if (spinEchoImageFileFormat.compareTo(ACCEPTED_FILE_TYPE) == 0) {
           fileOpener.open(spinEchoImage);
           WindowManager.getCurrentImage().setTitle(s7WindowTitle);
           spinEchoImageIsLoaded = true;
         } else {
           JOptionPane.showMessageDialog(frame, "Error: image type not .nii");
         }
       }
       // End opening spin echo image
 
     }
 
     /*
      * If the estimate object radius from spin echo button is selected:
      * Program gets user-inputted box and uses equations in section 2.6 to solve for
      * radius.
      * User defines a box V1, and another box V2, and a center. The program creates
      * new images based off these boxes so that they can be viewed. These new images
      * are used for the equations in section 2.6.
      */
     /*
      * TODO: Not finished. Code was not tested in depth and is not confirmed to be
      * working properly.
      */
     else if (e.getSource() == btn_estRadSpinEcho) {
 
       if (spinEchoImageIsLoaded) {
         updateVariables();
 
         // If images are already open, close them
         if (WindowManager.getImage(V1XY_Title) != null) {
           V1SE_XYImage.close();
         }
 
         if (WindowManager.getImage(V1XZ_Title) != null) {
           V1SE_XZImage.close();
         }
 
         // Getting spin echo image
         ImagePlus spinEchoImg = WindowManager.getImage(s7WindowTitle);
 
         // Getting volume of boxes
         int V1SE_size = (V1SE_x2 - V1SE_x1) * (V1SE_y2 - V1SE_y1) * (V1SE_z2 - V1SE_z1);
         int V2SE_size = (V2SE_x2 - V2SE_x1) * (V2SE_y2 - V2SE_y1) * (V2SE_z2 - V2SE_z1);
 
         // Getting center from GUI
         int VSE_centerX = spinEchoImage_XCoord;
         int VSE_centerY = spinEchoImage_YCoord;
         int VSE_centerZ = spinEchoImage_ZCoord;
 
         // V1 must be bigger than V2 and center must be inside both
         if (V2SE_size > V1SE_size) {
           JOptionPane.showMessageDialog(frame, "Error: V2 cannot be bigger than V1");
         } else if (VSE_centerX < V1SE_x1 || VSE_centerX > V1SE_x2) {
           JOptionPane.showMessageDialog(frame, "Error: center x coordinate must be on image");
         } else if (VSE_centerY < V1SE_y1 || VSE_centerY > V1SE_y2) {
           JOptionPane.showMessageDialog(frame, "Error: center y coordinate must be on image");
         } else if (VSE_centerZ < V1SE_z1 || VSE_centerZ > V1SE_z2) {
           JOptionPane.showMessageDialog(frame, "Error: center z coordinate must be on image");
         } else {
 
           // Initializing ImageProcessor objects to add data values to
           ImageProcessor IP_V1SE_XY = new FloatProcessor(V1SE_x2 - V1SE_x1, V1SE_y2 - V1SE_y1);
           ImageProcessor IP_V1SE_XZ = new FloatProcessor(V1SE_x2 - V1SE_x1, V1SE_z2 - V1SE_z1);
 
           // Adding data to new V1 image XY. Essentially just 'cropping' the V1 image
           spinEchoImg.setSlice(VSE_centerZ + 1);
           for (int i = V1SE_x1; i < V1SE_x2; i++) {
             for (int j = V1SE_y1; j < V1SE_y2; j++) {
               IP_V1SE_XY.putPixelValue(i - V1SE_x1, j - V1SE_y1, spinEchoImg.getProcessor().getPixelValue(i, j));
             }
           }
 
           // Same thing but for XZ
           for (int i = V1SE_x1; i < V1SE_x2; i++) {
             for (int j = V1SE_z1; j < V1SE_z2; j++) {
               spinEchoImg.setSlice(j + 1);
               IP_V1SE_XZ.putPixelValue(i - V1SE_x1, j - V1SE_z1,
                   spinEchoImg.getProcessor().getPixelValue(i, VSE_centerY));
             }
           }
 
           // Setting image titles
           V1XY_Title = "V1SE XY (V2SE displayed as yellow box)";
           V1XZ_Title = "V1SE XZ (V2SE displayed as yellow box)";
 
           // Creating new ImagePlus objects to display
           V1SE_XYImage = new ImagePlus(V1XY_Title, IP_V1SE_XY);
           V1SE_XZImage = new ImagePlus(V1XZ_Title, IP_V1SE_XZ);
 
           // Displaying
           V1SE_XYImage.show();
           V1SE_XZImage.show();
 
           // Adding V2 box to V1 XY images
           roiImgV1SE = new ROIS("MAG2XY");
           roiImgV1SE.addPointROI("Center", Math.abs(VSE_centerX - V1SE_x1), Math.abs(VSE_centerY - V1SE_y1));
           roiImgV1SE.addRectangle("V1SE", Math.abs(V1SE_x1 - V2SE_x1), Math.abs(V1SE_y1 - V2SE_y1), V2SE_x2 - V2SE_x1,
               V2SE_y2 - V2SE_y1);
           roiImgV1SE.displayROIS();
 
           // Same thing but for XZ
           roiImgV1SEXZ = new ROIS("MAG2XZ");
           roiImgV1SEXZ.addPointROI("Center", Math.abs(VSE_centerX - V1SE_x1), Math.abs(VSE_centerZ - V1SE_z1));
           roiImgV1SEXZ.addRectangle("V2SE", Math.abs(V1SE_x1 - V2SE_x1), Math.abs(V1SE_z1 - V2SE_z1), V2SE_x2 - V2SE_x1,
               V2SE_z2 - V2SE_z1);
           roiImgV1SEXZ.displayROIS();
 
           VSEisFound = true;
 
           // Variables for the sum and volume of the boxes
           double S1_SE = 0.0;
           double S2_SE = 0.0;
           double V1_SE = 0.0;
           double V2_SE = 0.0;
 
           // Getting S1 sum and volume
           for (int k = V1SE_z1; k < V1SE_z2; k++) {
             spinEchoImg.setSlice(k + 1);
             for (int j = V1SE_y1; j < V1SE_y2; j++) {
               for (int i = V1SE_x1; i < V1SE_x2; i++, V1_SE++) {
                 S1_SE += spinEchoImg.getProcessor().getPixelValue(i, j);
               }
             }
           }
 
           // Getting S2 sum and volume
           for (int k = V2SE_z1; k < V2SE_z2; k++) {
             spinEchoImg.setSlice(k + 1);
             for (int j = V2SE_y1; j < V2SE_y2; j++) {
               for (int i = V2SE_x1; i < V2SE_x2; i++, V2_SE++) {
                 S2_SE += spinEchoImg.getProcessor().getPixelValue(i, j);
               }
             }
           }
 
           logger.addVariable("S1_se", S1_SE);
           logger.addVariable("S2_se", S2_SE);
           logger.addVariable("V1_SE", V1_SE);
           logger.addVariable("V2_SE", V2_SE);
 
           // Equation 18
           double V0 = (S1_SE * V2_SE - S2_SE * V1_SE) / (S1_SE - S2_SE);
 
           // Adding to GUI
           lbl_V0Val.setText(String.valueOf(Math.round(V0 * 100.0) / 100.0));
           logger.addVariable("V0", V0);
 
           // Equation 16
           m_a = Math.pow((V0 * 3.0) / (4.0 * Math.PI), 1.0 / 3.0);
 
           // Equation 17
           double rho_SE0 = S2_SE / (V2_SE - V0);
           // Setting to GUI
           lbl_rho0SEVal.setText(String.valueOf(Math.round(rho_SE0 * 100.0) / 100.0));
           logger.addVariable("rho_SE0", rho_SE0);
 
           double deltaV = 1;
 
           // Defining signal to noise ratio
           double SNR_SE = rho_SE0 / m_snrStandardDeviation;
 
           // Equation 19
           m_dV0 = (Math.sqrt(deltaV) / SNR_SE) * Math.sqrt(V2_SE + Math.pow(V2_SE - V0, 2) / (V1_SE - V2_SE));
 
           // Adding v0 and dv0 to GUI
           lbl_V0Val.setText(String.valueOf(Math.round(V0 * 100.0) / 100.0) + " " + PLUS_MINUS + " "
               + String.valueOf(Math.round(m_dV0 * 100.0) / 100.0));
           logger.addVariable("m_dV0", m_dV0);
 
           // Error for a - was derived with Norman
           m_da = (m_a * m_dV0) / (3 * V0);
           // Adding to GUI
           lbl_aSE.setText(String.valueOf(Math.round(m_a * 100.0) / 100.0) + " " + PLUS_MINUS + " "
               + String.valueOf(Math.round(m_da * 100.0) / 100.0) + " pixels");
           logger.addVariable("m_a", m_a);
           logger.addVariable("m_da", m_da);
 
           // p = ga^3 can be rewritten to get dChi
           double dChi = (2.0 * m_p) / (GAMMARBAR * m_B0 * m_TELast * V0);
           // And its error
           double d_dChi = dChi * Math.sqrt(Math.pow(m_dp / m_p, 2) + Math.pow(m_dV0 / V0, 2));
 
           // Setting to GUI
           lbl_echoDChi.setText(String.valueOf(Math.round(dChi * 100.0) / 100.0) + " " + PLUS_MINUS + " "
               + String.valueOf(Math.round(d_dChi * 100.0) / 100.0) + " ppm");
         }
       } else {
         JOptionPane.showMessageDialog(frame, "Error: no spin echo image found");
       }
     }
   }
 

   /*
    * Function to remove background phase
    *
    * @param phaseVals 3D array of phase values
    */
   public void removeBGPhase(double[][][] phaseVals) {
 
     // Background phase must have a value to continue
     if (estimatedBGIsFound) {
 
       // Removing background phase from all data points in matrix
       for (int k = drawnRectangle_initialZ; k < drawnRectangle_initialZ + drawnRectangle_sizeZ; k++) {
         for (int i = drawnRectangle_initialX; i < drawnRectangle_initialX + drawnRectangle_sizeX; i++) {
           for (int j = drawnRectangle_initialY; j < drawnRectangle_initialY + drawnRectangle_sizeY; j++) {
             phaseVals[i][j][k] = Math.abs(phaseVals[i][j][k] - estimatedBackgroundPhase);
           }
         }
       }
     } else {
       JOptionPane.showMessageDialog(frame, "Error: No background phase found.");
     }
 
     return;
   }
 
   /*
    * Function to remove background phase
    *
    * @param phaseVals 2D array of phase values
    */
   public void removeBGPhase(float[][] phaseVals) {
 
     if (estimatedBGIsFound) {
       for (int i = 0; i < sizeOfSubpixelImage; i++) {
         for (int j = 0; j < sizeOfSubpixelImage; j++) {
           if (phaseVals[i][j] > 0) {
             phaseVals[i][j] = phaseVals[i][j] - (float) estimatedBackgroundPhase;
           } else if (phaseVals[i][j] < 0) {
             phaseVals[i][j] = phaseVals[i][j] + (float) estimatedBackgroundPhase;
           }
         }
       }
     } else {
       JOptionPane.showMessageDialog(frame, "Error: No background phase found.");
     }
 
     return;
   }
 
   /*
    * Function to convert voxel coordinate to subvoxel image coordinates
    *
    * @param coordinate the voxel position
    * 
    * @param axisFlag flag for what axis coordinate is in(0 for X, 1 for Y, 2 for
    * Z)
    * 
    * @return the subpixel coordinate
    */
   public static double pixelToSubpixel(double coordinate, int axisFlag) {
     // Calculated by taking total size of cropped image and dividing it by 2 - this
     // will give the center which is also the estimated center of the image in step
     // 2
     // The + 4.9 moves it from the middle of the subpixel to the maximum edge
     // because this is how ImageJ handles the voxel locations
     double subCenter = (2 * m_R0 + 1) * (10.0 / 2.0);
     double subpixelCoordinate = 0.0;
 
     // if x axis
     if (axisFlag == 0) {
       subpixelCoordinate = subCenter + (coordinate - (double) (int)(m_xCenter)) * 10.0;
     }
     if (axisFlag == 1) {
       subpixelCoordinate = subCenter + (coordinate - (double) (int)(m_yCenter)) * 10.0;
     }
     if (axisFlag == 2) {
       subpixelCoordinate = subCenter + (coordinate - (double) (int)(m_zCenter)) * 10.0;
     }
 
     return subpixelCoordinate;
   }
 
   /*
    * Function to update global variables in Java
    * Updates the variables to whatever values are in the respective textboxes
    * 
    * This function was written so that every time a button is clicked the
    * variables that go with each textbox can be updated - there are also some
    * functions in here that format the textboxes so that they have to be
    * doubles/ints.
    * This is so that there isn't an extra 20-40 lines of code in each action
    * handler.
    * Most global variables have the prefix "m_", similar to how it was written in
    * C++.
    */
   public void updateVariables() {
 
     try {
 
       if (WindowManager.getImage(s7WindowTitle) == null) {
         spinEchoImageIsLoaded = false;
       }
 
       if (WindowManager.getImage(s6MagWindowTitle) == null || WindowManager.getImage(s6PhaseWindowTitle) == null) {
         secondImagesAreLoaded = false;
       }
 
       if (WindowManager.getImage(s5MagWindowTitle) == null || WindowManager.getImage(s5PhaseWindowTitle) == null) {
         simulatedImagesAreLoaded = false;
       }
 
       if (!(WindowManager.getImage(subMagTitle) == null) && !(WindowManager.getImage(subMagXZTitle) == null)
           && !(WindowManager.getImage(subPhaseTitle) == null) && !(WindowManager.getImage(subPhaseXZTitle) == null)) {
         subpixelIsGenerated = true;
       } else {
         subpixelIsGenerated = false;
       }
 
       if (txt_sigSEVal.getText().isEmpty()) {
         txt_sigSEVal.setText("1.0");
         m_snrStandardDeviation = 1.0;
       } else {
         formatTextToDouble(txt_sigSEVal);
         m_snrStandardDeviation = Double.parseDouble(txt_sigSEVal.getText());
       }
 
       if (txt_snrVal.getText().isEmpty()) {
         txt_snrVal.setText("1.0");
         m_SNR = 1.0;
       } else {
         formatTextToDouble(txt_snrVal);
         m_SNR = Double.parseDouble(txt_snrVal.getText());
       }
 
       if (txt_eps12val.getText().isEmpty()) {
         txt_eps12val.setText("0.0");
         m_e12 = 0.0;
       } else {
         formatTextToDouble(txt_eps12val);
         m_e12 = Double.parseDouble(txt_eps12val.getText());
       }
 
       if (txt_eps23val.getText().isEmpty()) {
         txt_eps23val.setText("0.0");
         m_e23 = 0.0;
       } else {
         formatTextToDouble(txt_eps23val);
         m_e23 = Double.parseDouble(txt_eps23val.getText());
       }
 
       if (B0Text.getText().isEmpty()) {
         B0Text.setText("0.0");
         m_B0 = 0.0;
       } else {
         formatTextToDouble(B0Text);
         m_B0 = Double.parseDouble(B0Text.getText());
       }
       /*
        * if (XverText.getText().isEmpty()) {
        * XverText.setText("1.0");
        * m_X = 1.0;
        * } else {
        * m_X = Double.parseDouble(XverText.getText());
        * }
        *
        * if (YverText.getText().isEmpty()) {
        * YverText.setText("1.0");
        * m_Y = 1.0;
        * } else {
        * m_Y = Double.parseDouble(YverText.getText());
        * }
        *
        * if (ZverText.getText().isEmpty()) {
        * ZverText.setText("1.0");
        * m_Z = 1.0;
        * } else {
        * m_Z = Double.parseDouble(ZverText.getText());
        * }
        */
 
       if (txt_TEFirstVal.getText().isEmpty()) {
         txt_TEFirstVal.setText("0.0");
         m_TEFirst = 0.0;
       } else {
         formatTextToDouble(txt_TEFirstVal);
         m_TEFirst = Double.parseDouble(txt_TEFirstVal.getText()) / 1000.0;
       }
 
       if (TE_lastText.getText().isEmpty()) {
         TE_lastText.setText("0.0");
         m_TELast = 0.0;
       } else {
         formatTextToDouble(TE_lastText);
         m_TELast = Double.parseDouble(TE_lastText.getText()) / 1000.0;
       }
 
       if (txt_RChiVal.getText().isEmpty()) {
         txt_RChiVal.setText("0.0");
         m_RChi = 0.0;
       } else {
         formatTextToDouble(txt_RChiVal);
         m_RChi = Double.parseDouble(txt_RChiVal.getText());
       }
 
       if (txt_magMomVal.getText().isEmpty()) {
         txt_magMomVal.setText("0.0");
         m_p = 0.0;
       } else {
         formatTextToDouble(txt_magMomVal);
         m_p = Double.parseDouble(txt_magMomVal.getText());
       }
 
       if (lbl_errVal.getText().isEmpty()) {
         m_dp = 0.0;
       } else {
         m_dp = Double.parseDouble(lbl_errVal.getText());
       }
 
       if (txt_eqPhaseRC.getText().isEmpty()) {
         txt_eqPhaseRC.setText("1.0");
         phaseValue = 1.0;
       } else {
         formatTextToDouble(txt_eqPhaseRC);
         phaseValue = Double.parseDouble(txt_eqPhaseRC.getText());
       }
 
       if (txt_Ri.getText().isEmpty()) {
         txt_Ri.setText("0");
       } else {
         formatTextToDouble(txt_Ri);
         m_Ri = Double.parseDouble(txt_Ri.getText());
       }
 
       if (txt_rcx.getText().isEmpty()) {
         estimatedXIsFound = false;
       } else {
         estimatedXIsFound = true;
         formatTextToDouble(txt_rcx);
         m_xCenter = (int) Double.parseDouble(txt_rcx.getText());
         m_xCenter_subpixel = (int) pixelToSubpixel((int)(m_xCenter), 0);
       }
 
       if (txt_rcy.getText().isEmpty()) {
         estimatedYIsFound = false;
       } else {
         estimatedYIsFound = true;
         formatTextToDouble(txt_rcy);
         m_yCenter = (int) Double.parseDouble(txt_rcy.getText());
         m_yCenter_subpixel = (int) pixelToSubpixel((int)(m_yCenter), 1);
       }
 
       if (txt_rcz.getText().isEmpty()) {
         estimatedZIsFound = false;
       } else {
         estimatedZIsFound = true;
         formatTextToDouble(txt_rcz);
         m_zCenter = (int) Double.parseDouble(txt_rcz.getText()) - 1;
         m_zCenter_subpixel = (int) pixelToSubpixel((int)(m_zCenter), 2);
       }
 
       if (txt_M.getText().isEmpty()) {
         txt_M.setText(String.valueOf(50));
         percentOfMagnitudeToNeglect = 50;
       } else {
         formatTextToInteger(txt_M);
         percentOfMagnitudeToNeglect = Integer.parseInt(txt_M.getText());
       }
 
       if (txt_rc.getText().isEmpty()) {
         rCenterIsFound = false;
       } else {
         rCenterIsFound = true;
         formatTextToDouble(txt_rc);
         RCenter = Double.parseDouble(txt_rc.getText());
         setmVariables(grid, m_R0, RCenter, (double) (int)(m_xCenter), (double) (int)(m_yCenter),
             (double) (int)(m_zCenter), phaseValue);
       }
 
       if (lbl_estBkgPhaseVal.getText().isEmpty()) {
         estimatedBGIsFound = false;
       } else {
         estimatedBGIsFound = true;
         logger.addVariable("estimatedBackgroundPhase", estimatedBackgroundPhase);
         setBackPhase(estimatedBackgroundPhase);
       }
 
       if (txt_spx.getText().isEmpty() || txt_spy.getText().isEmpty() || txt_spz.getText().isEmpty()) {
         subpixelCenterIsFound = false;
       } else {
         subpixelCenterIsFound = true;
 
         formatTextToDouble(txt_spx);
         formatTextToDouble(txt_spy);
         formatTextToDouble(txt_spz);
 
         centerX_pixelCoordinates = Double.parseDouble(txt_spx.getText());
         logger.addVariable("centerX_pixelCoordinates", centerX_pixelCoordinates);
         centerY_pixelCoordinates = Double.parseDouble(txt_spy.getText());
         logger.addVariable("centerY_pixelCoordinates", centerY_pixelCoordinates);
         centerZ_pixelCoordinates = Double.parseDouble(txt_spz.getText()) - 1.0;
         logger.addVariable("centerZ_pixelCoordinates", centerZ_pixelCoordinates);
         centerX_subpixelCoordinates = (int) pixelToSubpixel(centerX_pixelCoordinates, 0);
         logger.addVariable("centerX_subpixelCoordinates", centerX_subpixelCoordinates);
         centerY_subpixelCoordinates = (int) pixelToSubpixel(centerY_pixelCoordinates, 1);
         logger.addVariable("centerY_subpixelCoordinates", centerY_subpixelCoordinates);
         centerZ_subpixelCoordinates = (int) pixelToSubpixel(centerZ_pixelCoordinates, 2);
         logger.addVariable("centerZ_subpixelCoordinates", centerZ_subpixelCoordinates);
       }
 
       if (txt_r1.getText().isEmpty()) {
         R1IsFound = false;
       } else {
         R1IsFound = true;
         formatTextToDouble(txt_r1);
         m_R1 = Double.parseDouble(txt_r1.getText());
       }
 
       if (txt_r2.getText().isEmpty()) {
         R2IsFound = false;
       } else {
         R2IsFound = true;
         formatTextToDouble(txt_r2);
         m_R2 = Double.parseDouble(txt_r2.getText());
       }
 
       if (txt_r3.getText().isEmpty()) {
         R3IsFound = false;
       } else {
         R3IsFound = true;
         formatTextToDouble(txt_r3);
         m_R3 = Double.parseDouble(txt_r3.getText());
       }
 
       if (txt_spinCenterXVal.getText().isEmpty() || txt_spinCenterXVal.getText().compareTo("0") == 0) {
         if (estimatedXIsFound) {
           txt_spinCenterXVal.setText(txt_rcx.getText());
         } else {
           txt_spinCenterXVal.setText("0");
         }
         formatTextToInteger(txt_spinCenterXVal);
         spinEchoImage_XCoord = Integer.parseInt(txt_spinCenterXVal.getText());
       } else {
         formatTextToInteger(txt_spinCenterXVal);
         spinEchoImage_XCoord = Integer.parseInt(txt_spinCenterXVal.getText());
       }
 
       if (txt_spinCenterYVal.getText().isEmpty() || txt_spinCenterYVal.getText().compareTo("0") == 0) {
         if (estimatedYIsFound) {
           txt_spinCenterYVal.setText(txt_rcy.getText());
         } else {
           txt_spinCenterYVal.setText("0");
         }
         formatTextToInteger(txt_spinCenterYVal);
         spinEchoImage_YCoord = Integer.parseInt(txt_spinCenterYVal.getText());
       } else {
         formatTextToInteger(txt_spinCenterYVal);
         spinEchoImage_YCoord = Integer.parseInt(txt_spinCenterYVal.getText());
       }
 
       if (txt_spinCenterZVal.getText().isEmpty() || txt_spinCenterZVal.getText().compareTo("1") == 0) {
         if (estimatedZIsFound) {
           txt_spinCenterZVal.setText(txt_rcz.getText());
         } else {
           txt_spinCenterZVal.setText("1");
         }
         formatTextToInteger(txt_spinCenterZVal);
         spinEchoImage_ZCoord = Integer.parseInt(txt_spinCenterZVal.getText()) - 1;
       } else {
         formatTextToInteger(txt_spinCenterZVal);
         spinEchoImage_ZCoord = Integer.parseInt(txt_spinCenterZVal.getText()) - 1;
       }
 
       if (txt_v1seXVal1.getText().isEmpty()) {
         txt_v1seXVal1.setText("0");
         V1SE_x1 = 0;
       } else {
         formatTextToInteger(txt_v1seXVal1);
         V1SE_x1 = Integer.parseInt(txt_v1seXVal1.getText());
       }
 
       if (txt_v1seXVal2.getText().isEmpty()) {
         txt_v1seXVal2.setText("0");
         V1SE_x2 = 0;
       } else {
         formatTextToInteger(txt_v1seXVal2);
         V1SE_x2 = Integer.parseInt(txt_v1seXVal2.getText());
       }
 
       if (txt_v1seYVal1.getText().isEmpty()) {
         txt_v1seYVal1.setText("0");
         V1SE_y1 = 0;
       } else {
         formatTextToInteger(txt_v1seYVal1);
         V1SE_y1 = Integer.parseInt(txt_v1seYVal1.getText());
       }
 
       if (txt_v1seYVal2.getText().isEmpty()) {
         txt_v1seYVal2.setText("0");
         V1SE_y2 = 0;
       } else {
         formatTextToInteger(txt_v1seYVal2);
         V1SE_y2 = Integer.parseInt(txt_v1seYVal2.getText());
       }
 
       if (txt_v1seZVal1.getText().isEmpty()) {
         txt_v1seZVal1.setText("1");
         V1SE_z1 = 0;
       } else {
         formatTextToInteger(txt_v1seZVal1);
         V1SE_z1 = Integer.parseInt(txt_v1seZVal1.getText()) - 1;
       }
 
       if (txt_v1seZVal2.getText().isEmpty()) {
         txt_v1seZVal2.setText("1");
         V1SE_z2 = 0;
       } else {
         formatTextToInteger(txt_v1seZVal2);
         V1SE_z2 = Integer.parseInt(txt_v1seZVal2.getText()) - 1;
       }
 
       if (txt_v2seXVal1.getText().isEmpty()) {
         txt_v2seXVal1.setText("0");
         V2SE_x1 = 0;
       } else {
         formatTextToInteger(txt_v2seXVal1);
         V2SE_x1 = Integer.parseInt(txt_v2seXVal1.getText());
       }
 
       if (txt_v2seXVal2.getText().isEmpty()) {
         txt_v2seXVal2.setText("0");
         V2SE_x2 = 0;
       } else {
         formatTextToInteger(txt_v2seXVal2);
         V2SE_x2 = Integer.parseInt(txt_v2seXVal2.getText());
       }
 
       if (txt_v2seYVal1.getText().isEmpty()) {
         txt_v2seYVal1.setText("0");
         V2SE_y1 = 0;
       } else {
         formatTextToInteger(txt_v2seYVal1);
         V2SE_y1 = Integer.parseInt(txt_v2seYVal1.getText());
       }
 
       if (txt_v2seYVal2.getText().isEmpty()) {
         txt_v2seYVal2.setText("0");
         V2SE_y2 = 0;
       } else {
         formatTextToInteger(txt_v2seYVal2);
         V2SE_y2 = Integer.parseInt(txt_v2seYVal2.getText());
       }
 
       if (txt_v2seZVal1.getText().isEmpty()) {
         txt_v2seZVal1.setText("1");
         V2SE_z1 = 0;
       } else {
         formatTextToInteger(txt_v2seZVal1);
         V2SE_z1 = Integer.parseInt(txt_v2seZVal1.getText()) - 1;
       }
 
       if (txt_v2seZVal2.getText().isEmpty()) {
         txt_v2seZVal2.setText("1");
         V2SE_z2 = 0;
       } else {
         formatTextToInteger(txt_v2seZVal2);
         V2SE_z2 = Integer.parseInt(txt_v2seZVal2.getText()) - 1;
       }
 
       if (estimateCenterRadii_isClicked) {
         setXYZ(m_xCenter, m_yCenter, m_zCenter);
         setmR123(m_R1, m_R2, m_R3);
         estimatedPValue = phaseValue * Math.pow(RCenter, 3);
         R1PhaseCalc = estimatedPValue / Math.pow(m_R1, 3);
         R2PhaseCalc = estimatedPValue / Math.pow(m_R2, 3);
         R3PhaseCalc = estimatedPValue / Math.pow(m_R3, 3);
 
         setR123PhaseCalc(R1PhaseCalc, R2PhaseCalc, R3PhaseCalc);
         setSmallBox(innerBox_initialX, innerBox_initialY, innerBox_initialZ, innerBox_sizeX, innerBox_sizeY,
             innerBox_sizeZ);
         setCenterL(Center_L_x, Center_L_y, Center_L_z);
         setCenterM(Center_M_x, Center_M_y, Center_M_z);
         setCenterS(center_sx, center_sy, center_sz);
       }
 
       setMagMomentVariables(m_SNR, m_e12, m_e23, m_B0, m_RChi, m_TELast);
       setRi(m_Ri);
 
       if (subpixelIsGenerated) {
         setPhaseXYMatrix(subpixelPhaseMatrix);
         setPhaseXZMatrix(subpixelPhaseMatrixXZ);
         setMagXYMatrix(subpixelMagMatrix);
         setMagXZMatrix(subpixelMagMatrixXZ);
         setRealImagNumbers(croppedRealNumbers3D, croppedImaginaryNumbers3D);
       }
 
       if (estimateSubpixelCenter_isClicked) {
         setXYZ(centerX_pixelCoordinates, centerY_pixelCoordinates, centerZ_pixelCoordinates);
       }
 
       if (!(txt_TEFirstVal.getText().isEmpty()) && !(TE_lastText.getText().isEmpty()) && !(B0Text.getText().isEmpty())
           && !(txt_RChiVal.getText().isEmpty())) {
         setStep6Variables(m_TEFirst, m_TELast, m_B0, m_RChi);
       }
 
       logger.addInfo("Updated variables");
     } catch (Exception exc) {
       JOptionPane.showMessageDialog(frame, "Error in updateVariables():\n" + exc.toString());
     }
   }
 
   /*
    * Function to reset CMM3D
    * Used when loading new images
    * Clears the GUI and local variables
    */
   public void clearVariables() {
     try {
 
       txt_spinCenterXVal.setText("0");
       txt_spinCenterYVal.setText("0");
       txt_spinCenterZVal.setText("1");
 
       txt_v1seXVal1.setText("0");
       txt_v1seXVal2.setText("0");
       txt_v1seYVal1.setText("0");
       txt_v1seYVal2.setText("0");
       txt_v1seZVal1.setText("1");
       txt_v1seZVal2.setText("1");
 
       txt_v2seXVal1.setText("0");
       txt_v2seXVal2.setText("0");
       txt_v2seYVal1.setText("0");
       txt_v2seYVal2.setText("0");
       txt_v2seZVal1.setText("1");
       txt_v2seZVal2.setText("1");
 
       txt_sigSEVal.setText("1.0");
 
       lbl_V0Val.setText("");
 
       lbl_aSE.setText("");
 
       lbl_echoDChi.setText("");
 
       // lbl_d_V0Val.setText("");
 
       lbl_rho0SEVal.setText("");
 
       txt_snrVal.setText("1.0");
       m_SNR = 1.0;
 
       txt_eps12val.setText("0.0");
       m_e12 = 0.0;
 
       txt_eps23val.setText("0.0");
       m_e23 = 0.0;
 
       B0Text.setText("0.0");
       m_B0 = 0.0;
 
       /*
        * XverText.setText("1.0");
        * m_X = 1.0;
        * YverText.setText("1.0");
        *
        * m_Y = 1.0;
        * ZverText.setText("1.0");
        *
        * m_Z = 1.0;
        */
 
       txt_sigSEVal.setText("1.0");
       m_snrStandardDeviation = 1.0;
 
       txt_TEFirstVal.setText("0.0");
       m_TEFirst = 0.0;
 
       TE_lastText.setText("0.0");
       m_TELast = 0.0;
 
       txt_RChiVal.setText("0.0");
       m_RChi = 0.0;
 
       txt_Ri.setText("0");
       m_Ri = 0;
 
       txt_magMomVal.setText("");
 
       lbl_rho0val.setText("");
 
       lbl_errVal.setText("");
 
       lbl_aVal.setText("");
 
       lbl_dchiVal.setText("");
 
       lbl_ImRi.setText("Si=");
 
       lbl_ReRi.setText("Si=");
 
       m_Si = Double.NaN;
 
       m_Si2 = Double.NaN;
 
       txt_M.setText("50");
       percentOfMagnitudeToNeglect = 50;
 
       txt_eqPhaseRC.setText("1.0");
       phaseValue = 1.0;
 
       txt_rc.setText("");
       RCenter = Double.NaN;
 
       lbl_estBkgPhaseVal.setText("");
       estimatedBackgroundPhase = Double.NaN;
 
       txt_spx.setText("");
       centerX_pixelCoordinates = Double.NaN;
       centerX_subpixelCoordinates = 0;
       txt_spy.setText("");
       centerY_pixelCoordinates = Double.NaN;
       centerY_subpixelCoordinates = 0;
       txt_spz.setText("");
       centerZ_pixelCoordinates = Double.NaN;
       centerZ_subpixelCoordinates = 0;
 
       txt_r1.setText("");
       m_R1 = Double.NaN;
       txt_r2.setText("");
       m_R2 = Double.NaN;
       txt_r3.setText("");
       m_R3 = Double.NaN;
 
       lbl_r3phaseAct.setText("");
       lbl_r3phaseCalc.setText("");
       lbl_r2phaseAct.setText("");
       lbl_r2phaseCalc.setText("");
       lbl_r1phaseAct.setText("");
       lbl_r1phaseCalc.setText("");
 
       txt_rcx.setText("");
       m_xCenter = 0;
       m_xCenter_subpixel = 0;
       txt_rcy.setText("");
       m_yCenter = 0;
       m_yCenter_subpixel = 0;
       txt_rcz.setText("");
       m_zCenter = 0;
       m_zCenter_subpixel = 0;
 
       estimatedPValue = Double.NaN;
       R1PhaseCalc = Double.NaN;
       R2PhaseCalc = Double.NaN;
       R3PhaseCalc = Double.NaN;
 
       chkbx_showrc.setSelected(false);
 
       // magText.setText("");
 
       // phaseText.setText("");
 
       if (subpixelIsGenerated) {
         if (WindowManager.getImage(subMagTitle) != null)
           WindowManager.getImage(subMagTitle).close();
         if (WindowManager.getImage(subMagXZTitle) != null)
           WindowManager.getImage(subMagXZTitle).close();
         if (WindowManager.getImage(subPhaseTitle) != null)
           WindowManager.getImage(subPhaseTitle).close();
         if (WindowManager.getImage(subPhaseXZTitle) != null)
           WindowManager.getImage(subPhaseXZTitle).close();
       }
 
       if (WindowManager.getImage(V1XY_Title) != null) {
         WindowManager.getImage(V1XY_Title).close();
       }
       if (WindowManager.getImage(V1XZ_Title) != null) {
         WindowManager.getImage(V1XZ_Title).close();
       }
 
       WindowManager.closeAllWindows();
 
       subpixelIsGenerated = false;
       estimatedXIsFound = false;
       estimatedYIsFound = false;
       estimatedZIsFound = false;
       rCenterIsFound = false;
       estimatedBGIsFound = false;
       subpixelCenterIsFound = false;
       R1IsFound = false;
       R2IsFound = false;
       R3IsFound = false;
       estimateCenterRadii_isClicked = false;
       estimateSubpixelCenter_isClicked = false;
 
       logger.addInfo("Cleared variables");
     } catch (Exception exc) {
       JOptionPane.showMessageDialog(frame, "Error in clearVariables():\n" + exc.toString());
     }
   }
 
   /*
    * Function to make sure a textbox has to be a double
    * 
    * @param text textbox to convert
    */
   public void formatTextToDouble(JTextField text) {
     if (!(text.getText().contains("."))) {
       text.setText(text.getText() + ".0");
     }
   }
 
   /*
    * Function to make sure a textbox has to be an integer
    * 
    * @param text textbox to convert
    */
   public void formatTextToInteger(JTextField text) {
     if (text.getText().contains(".")) {
       text.setText(text.getText().substring(0, text.getText().indexOf(".")));
     }
   }
 
   /*
    * Built using Eclipse and MigLayout. Much easier to manage then before.
    * If this is to be edited I advise against changing code directly.
    * Copy and paste this function into Eclipse (or any other IDE) and
    * make sure you have some GUI plug-in configured so that you can view it
    * without compiling + running every time. Also make sure to have MigLayout
    * configured.
    */
   private static void initialize() {
    try{
      frame = new JFrame("Calculate Magnetic Moment 3D");
      frame.setAlwaysOnTop(false);
      frame.setBounds(100, 100, 800, 700);
      frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE); 
		  frame.getContentPane().setLayout(new MigLayout("", "[][][grow][45px:45.00px:45px][45px:45px:45px,grow][45px:45px:45px,grow][48.00,grow][]", "[][][][][][][][][][][][][][][][][][][][][][][]"));
  }
    catch(Exception e){
      logger.addError(e.toString());
    }

		lbl_stepone = new JLabel("1.");
		frame.getContentPane().add(lbl_stepone, "cell 0 0");
		
		btn_loadImages = new JButton("Load Magnitude and Phase Images");
		frame.getContentPane().add(btn_loadImages, "cell 1 0,growx");
		btn_loadImages.addActionListener(new Calculate_Magnetic_Moment_3D());
    
		lbl_steptwo = new JLabel("2.");
		frame.getContentPane().add(lbl_steptwo, "cell 0 1");
		
		btn_estCR = new JButton("Estimate Center/Radii");
		frame.getContentPane().add(btn_estCR, "cell 1 1,growx");
		btn_estCR.addActionListener(new Calculate_Magnetic_Moment_3D());
    
		lbl_M = new JLabel("|M|:");
		frame.getContentPane().add(lbl_M, "flowx,cell 2 1,alignx left");
		
		lbl_rcx = new JLabel("x=");
		frame.getContentPane().add(lbl_rcx, "flowx,cell 3 1");
		
		lbl_rcy = new JLabel("y=");
		frame.getContentPane().add(lbl_rcy, "flowx,cell 4 1,alignx left");
		
		lbl_rcz = new JLabel("z=");
		frame.getContentPane().add(lbl_rcz, "flowx,cell 5 1,alignx left");
		
		lbl_rczCorrection = new JLabel("-1");
		frame.getContentPane().add(lbl_rczCorrection, "cell 6 1,alignx left");
		
		lbl_eqPhase = new JLabel("Equatorial Phase at RCenter=");
		frame.getContentPane().add(lbl_eqPhase, "flowx,cell 1 2");
		
		lbl_eqPhaseUnit = new JLabel("radian(s)");
		frame.getContentPane().add(lbl_eqPhaseUnit, "cell 2 2");
		
		lbl_rc = new JLabel("RCenter=");
		frame.getContentPane().add(lbl_rc, "cell 3 2,alignx trailing");
		
		txt_rc = new JTextField();
		txt_rc.setColumns(3);
		frame.getContentPane().add(txt_rc, "cell 4 2,growx");
		
		lbl_rcUnit = new JLabel("pixels");
		frame.getContentPane().add(lbl_rcUnit, "cell 5 2");
		
		lbl_stepthree = new JLabel("3.");
		frame.getContentPane().add(lbl_stepthree, "cell 0 3");
		
		btn_genSubpix = new JButton("Generate Subpixel Grid/Data");
		frame.getContentPane().add(btn_genSubpix, "cell 1 3,growx");
		btn_genSubpix.addActionListener(new Calculate_Magnetic_Moment_3D());
    
		btn_removeBkg = new JButton("Remove Bkg");
		frame.getContentPane().add(btn_removeBkg, "cell 2 3");
		btn_removeBkg.addActionListener(new Calculate_Magnetic_Moment_3D());
    
		chkbx_showrc = new JCheckBox("Show RCenter");
		chkbx_showrc.setVerticalAlignment(SwingConstants.TOP);
		frame.getContentPane().add(chkbx_showrc, "cell 3 3");
		
		lbl_gridSize = new JLabel("Grid Size:");
		frame.getContentPane().add(lbl_gridSize, "flowx,cell 7 3");
		
		lbl_stepfour = new JLabel("4.");
		frame.getContentPane().add(lbl_stepfour, "cell 0 4");
		
		btn_estSubC = new JButton("Estimate Subpixel Center");
		frame.getContentPane().add(btn_estSubC, "cell 1 4,growx");
		btn_estSubC.addActionListener(new Calculate_Magnetic_Moment_3D());
    
		txt_rcx = new JTextField();
		frame.getContentPane().add(txt_rcx, "cell 3 1");
		txt_rcx.setColumns(4);
		
		txt_rcy = new JTextField();
		frame.getContentPane().add(txt_rcy, "cell 4 1");
		txt_rcy.setColumns(4);
		
		txt_rcz = new JTextField();
		frame.getContentPane().add(txt_rcz, "cell 5 1");
		txt_rcz.setColumns(4);
		
		txt_M = new JTextField();
		txt_M.setText("50");
		frame.getContentPane().add(txt_M, "cell 2 1,alignx left");
		txt_M.setColumns(2);
		
		lbl_MPcnt = new JLabel("%");
		frame.getContentPane().add(lbl_MPcnt, "cell 2 1");
		
		txt_eqPhaseRC = new JTextField();
		txt_eqPhaseRC.setText("1.0");
		frame.getContentPane().add(txt_eqPhaseRC, "cell 1 2");
		txt_eqPhaseRC.setColumns(5);
		
		lbl_spx = new JLabel("x=");
		frame.getContentPane().add(lbl_spx, "flowx,cell 3 4");
		
		txt_spx = new JTextField();
		lbl_spx.setLabelFor(txt_spx);
		frame.getContentPane().add(txt_spx, "cell 3 4");
		txt_spx.setColumns(4);
		
		lbl_spy = new JLabel("y=");
		frame.getContentPane().add(lbl_spy, "flowx,cell 4 4");
		
		txt_spy = new JTextField();
		lbl_spy.setLabelFor(txt_spy);
		txt_spy.setColumns(4);
		frame.getContentPane().add(txt_spy, "cell 4 4");
		
		lbl_spz = new JLabel("z=");
		frame.getContentPane().add(lbl_spz, "flowx,cell 5 4");
		
		txt_spz = new JTextField();
		lbl_spz.setLabelFor(txt_spz);
		txt_spz.setColumns(4);
		frame.getContentPane().add(txt_spz, "cell 5 4");
		
		lbl_spzCorrection = new JLabel("-1");
		lbl_spzCorrection.setLabelFor(txt_spz);
		frame.getContentPane().add(lbl_spzCorrection, "cell 6 4");
		
		btn_redraw = new JButton("Redraw Center");
		frame.getContentPane().add(btn_redraw, "flowx,cell 1 5");
		btn_redraw.addActionListener(new Calculate_Magnetic_Moment_3D());
    
		btn_verifyRadii = new JButton("Verify Radii");
		frame.getContentPane().add(btn_verifyRadii, "cell 1 5");
		btn_verifyRadii.addActionListener(new Calculate_Magnetic_Moment_3D());
    
		lbl_calculated = new JLabel("Calculated");
		frame.getContentPane().add(lbl_calculated, "cell 3 5");
		
		lbl_actual = new JLabel("Actual");
		frame.getContentPane().add(lbl_actual, "cell 5 5");
		
		lbl_r1 = new JLabel("R1=");
		frame.getContentPane().add(lbl_r1, "flowx,cell 1 6");
		
		txt_r1 = new JTextField();
		frame.getContentPane().add(txt_r1, "cell 1 6");
		txt_r1.setColumns(5);
		
		lbl_r1unit = new JLabel("pixels");
		frame.getContentPane().add(lbl_r1unit, "cell 1 6");
		
		lbl_r1phase = new JLabel("R1 Corresponding Phase:");
		frame.getContentPane().add(lbl_r1phase, "cell 2 6");
		
		lbl_r1phaseCalc = new JLabel("");
		frame.getContentPane().add(lbl_r1phaseCalc, "cell 3 6");
		
		lbl_r1phaseUnit = new JLabel("radians");
		frame.getContentPane().add(lbl_r1phaseUnit, "cell 4 6");
		
		lbl_r1phaseAct = new JLabel("");
		frame.getContentPane().add(lbl_r1phaseAct, "cell 5 6");
		
		lbl_r1AphaseUnit = new JLabel("radians");
		frame.getContentPane().add(lbl_r1AphaseUnit, "cell 6 6");
		
		btn_plotX = new JButton("Plot X Phase Profiles");
		frame.getContentPane().add(btn_plotX, "cell 7 6");
		btn_plotX.addActionListener(new Calculate_Magnetic_Moment_3D());
    
		lbl_r2 = new JLabel("R2=");
		frame.getContentPane().add(lbl_r2, "flowx,cell 1 7");
		
		lbl_r2phase = new JLabel("R2 Corresponding Phase:");
		frame.getContentPane().add(lbl_r2phase, "cell 2 7");
		
		lbl_r2phaseCalc = new JLabel("");
		frame.getContentPane().add(lbl_r2phaseCalc, "cell 3 7");
		
		lbl_r2phaseUnit = new JLabel("radians");
		frame.getContentPane().add(lbl_r2phaseUnit, "cell 4 7");
		
		lbl_r2phaseAct = new JLabel("");
		frame.getContentPane().add(lbl_r2phaseAct, "cell 5 7");
		
		lbl_r2AphaseUnit = new JLabel("radians");
		frame.getContentPane().add(lbl_r2AphaseUnit, "cell 6 7");
		
		btn_plotY = new JButton("Plot Y Phase Profiles");
		frame.getContentPane().add(btn_plotY, "cell 7 7");
		btn_plotY.addActionListener(new Calculate_Magnetic_Moment_3D());
    
		lbl_r3 = new JLabel("R3=");
		frame.getContentPane().add(lbl_r3, "flowx,cell 1 8");
		
		txt_r2 = new JTextField();
		txt_r2.setColumns(5);
		frame.getContentPane().add(txt_r2, "cell 1 7");
		
		txt_r3 = new JTextField();
		txt_r3.setColumns(5);
		frame.getContentPane().add(txt_r3, "cell 1 8");
		
		lbl_r2unit = new JLabel("pixels");
		frame.getContentPane().add(lbl_r2unit, "cell 1 7");
		
		lbl_r3unit = new JLabel("pixels");
		frame.getContentPane().add(lbl_r3unit, "cell 1 8");
		
		lbl_r3phase = new JLabel("R3 Corresponding Phase:");
		frame.getContentPane().add(lbl_r3phase, "cell 2 8");
		
		lbl_r3phaseCalc = new JLabel("");
		frame.getContentPane().add(lbl_r3phaseCalc, "cell 3 8");
		
		lbl_r3phaseUnit = new JLabel("radians");
		frame.getContentPane().add(lbl_r3phaseUnit, "cell 4 8");
		
		lbl_r3phaseAct = new JLabel("");
		frame.getContentPane().add(lbl_r3phaseAct, "cell 5 8");
		
		lbl_r3AphaseUnit = new JLabel("radians");
		frame.getContentPane().add(lbl_r3AphaseUnit, "cell 6 8");
		
		lbl_gridSizeBase = new JLabel("10" + CUBED);
		frame.getContentPane().add(lbl_gridSizeBase, "cell 7 3");
		
		btn_plotZ = new JButton("Plot Z Phase Profiles");
		frame.getContentPane().add(btn_plotZ, "cell 7 8");
		btn_plotZ.addActionListener(new Calculate_Magnetic_Moment_3D());
    
		btn_estBkgDens = new JButton("Estimate Bkg & " + RHO + "0");
		frame.getContentPane().add(btn_estBkgDens, "cell 1 9,growx");
		btn_estBkgDens.addActionListener(new Calculate_Magnetic_Moment_3D());
    
		lbl_rho0 = new JLabel(RHO + "0 =");
		frame.getContentPane().add(lbl_rho0, "flowx,cell 2 9");
		
		lbl_rho0val = new JLabel("");
		frame.getContentPane().add(lbl_rho0val, "cell 2 9");
		
		lbl_estBkgPhase = new JLabel("Estimated Background Phase =");
		frame.getContentPane().add(lbl_estBkgPhase, "cell 3 9");
		
		lbl_estBkgPhaseVal = new JLabel("");
		frame.getContentPane().add(lbl_estBkgPhaseVal, "cell 6 9");
		
		lbl_estBkgPhaseUnit = new JLabel("radians");
		frame.getContentPane().add(lbl_estBkgPhaseUnit, "cell 7 9");
		
		lbl_stepfive = new JLabel("5.");
		frame.getContentPane().add(lbl_stepfive, "cell 0 10");
		
		btn_calcMagMom = new JButton("Calculate Magnetic Moment");
		frame.getContentPane().add(btn_calcMagMom, "cell 1 10,growx");
		btn_calcMagMom.addActionListener(new Calculate_Magnetic_Moment_3D());
    
		lbl_magMom = new JLabel("|p|=");
		frame.getContentPane().add(lbl_magMom, "flowx,cell 2 10");
		
		btn_loadSimImg = new JButton("Load Simulated Images");
		frame.getContentPane().add(btn_loadSimImg, "cell 1 11,growx");
		btn_loadSimImg.addActionListener(new Calculate_Magnetic_Moment_3D());
    
		txt_magMomVal = new JTextField();
		frame.getContentPane().add(txt_magMomVal, "cell 2 10");
		txt_magMomVal.setColumns(4);
		
		lbl_magMomUnit = new JLabel("radians*pixel" + CUBED);
		frame.getContentPane().add(lbl_magMomUnit, "cell 2 10");
		
		lbl_snr = new JLabel("SNR=");
		frame.getContentPane().add(lbl_snr, "flowx,cell 1 12");
		
		txt_snrVal = new JTextField();
		txt_snrVal.setText("1.0");
		frame.getContentPane().add(txt_snrVal, "cell 1 12");
		txt_snrVal.setColumns(3);
		
		lbl_eps12 = new JLabel(EPSILON + "12 =");
		frame.getContentPane().add(lbl_eps12, "flowx,cell 1 13");
		
		txt_eps12val = new JTextField();
		txt_eps12val.setText("0.0");
		txt_eps12val.setColumns(3);
		frame.getContentPane().add(txt_eps12val, "cell 1 13");
		
		lbl_ReRi = new JLabel("Real(S" + ITALICIZED_I + ") =");
		frame.getContentPane().add(lbl_ReRi, "flowx,cell 2 13");
		
		lbl_eps23 = new JLabel(EPSILON + "23 =");
		frame.getContentPane().add(lbl_eps23, "flowx,cell 1 14");
		
		txt_eps23val = new JTextField();
		txt_eps23val.setText("0.0");
		txt_eps23val.setColumns(3);
		frame.getContentPane().add(txt_eps23val, "cell 1 14");
		
		lbl_err = new JLabel("      " + DELTA + RHO + "/" + RHO + " =");
		lbl_err.setHorizontalAlignment(SwingConstants.RIGHT);
		frame.getContentPane().add(lbl_err, "cell 1 13");
		
		lbl_errVal = new JLabel("");
		frame.getContentPane().add(lbl_errVal, "cell 1 13");
		
		lbl_Ri = new JLabel("R" + ITALICIZED_I + " =");
		frame.getContentPane().add(lbl_Ri, "flowx,cell 2 12");
		
		txt_Ri = new JTextField();
		txt_Ri.setColumns(4);
		frame.getContentPane().add(txt_Ri, "cell 2 12");
		
		lbl_ImRi = new JLabel("Imag(S" + ITALICIZED_I + ") =");
		frame.getContentPane().add(lbl_ImRi, "flowx,cell 2 14");
		
		lbl_ReRiVal = new JLabel("");
		frame.getContentPane().add(lbl_ReRiVal, "cell 2 13");
		
		lbl_ImRiVal = new JLabel("");
		frame.getContentPane().add(lbl_ImRiVal, "cell 2 14");
		
		btn_sumRi = new JButton("Sum");
		frame.getContentPane().add(btn_sumRi, "cell 2 12");
		btn_sumRi.addActionListener(new Calculate_Magnetic_Moment_3D());
    
		lbl_stepsix = new JLabel("6.");
		frame.getContentPane().add(lbl_stepsix, "cell 0 15");
		
		btn_loadTE = new JButton("Load First TE Images");
		frame.getContentPane().add(btn_loadTE, "cell 1 15,growx");
		btn_loadTE.addActionListener(new Calculate_Magnetic_Moment_3D());
    
		lbl_TEFirst = new JLabel("TE_first =");
		frame.getContentPane().add(lbl_TEFirst, "flowx,cell 2 15");
		
		lbl_dchi = new JLabel(_DELTA + CHI + " =");
		frame.getContentPane().add(lbl_dchi, "flowx,cell 4 15");
		
		btn_unk = new JButton("TODO: Name");
		frame.getContentPane().add(btn_unk, "cell 1 16,growx");
		btn_unk.addActionListener(new Calculate_Magnetic_Moment_3D());
    
		txt_TEFirstVal = new JTextField();
		frame.getContentPane().add(txt_TEFirstVal, "cell 2 15");
		txt_TEFirstVal.setColumns(4);
		
		lbl_TEFirstUnit = new JLabel("ms");
		frame.getContentPane().add(lbl_TEFirstUnit, "cell 2 15");
		
		lbl_RChi = new JLabel("      R" + _DELTA + CHI + " =");
		frame.getContentPane().add(lbl_RChi, "flowx,cell 2 16");
		
		txt_RChiVal = new JTextField();
		txt_RChiVal.setText("0.0");
		frame.getContentPane().add(txt_RChiVal, "cell 2 16");
		txt_RChiVal.setColumns(4);
		
		lbl_RChiUnit = new JLabel("pixels");
		frame.getContentPane().add(lbl_RChiUnit, "cell 2 16");
		
		lbl_dchiVal = new JLabel("");
		frame.getContentPane().add(lbl_dchiVal, "cell 4 15");
		
		lbl_a = new JLabel("  a =");
		frame.getContentPane().add(lbl_a, "flowx,cell 4 16");
		
		lbl_aVal = new JLabel("");
		frame.getContentPane().add(lbl_aVal, "cell 4 16");
		
		lbl_stepseven = new JLabel("7.");
		frame.getContentPane().add(lbl_stepseven, "cell 0 17");
		
		btn_loadspinecho = new JButton("Load Spin Echo Images");
		frame.getContentPane().add(btn_loadspinecho, "cell 1 17,growx");
		btn_loadspinecho.addActionListener(new Calculate_Magnetic_Moment_3D());
    
		lbl_sigSE = new JLabel("        " + _SIGMA + "SE =");
		frame.getContentPane().add(lbl_sigSE, "flowx,cell 2 17,alignx left");
		
		txt_sigSEVal = new JTextField();
		frame.getContentPane().add(txt_sigSEVal, "cell 2 17");
		txt_sigSEVal.setColumns(4);
		
		lbl_spinCenter = new JLabel("Second Image Center:  (");
		frame.getContentPane().add(lbl_spinCenter, "flowx,cell 1 18,alignx trailing");
		
		txt_v1seXVal1 = new JTextField();
		txt_v1seXVal1.setColumns(3);
		frame.getContentPane().add(txt_v1seXVal1, "flowx,cell 2 18,alignx center");
		
		lbl_innerBrack1 = new JLabel("-1 )");
		frame.getContentPane().add(lbl_innerBrack1, "flowx,cell 3 18");
		
		lbl_v1se = new JLabel("V1,SE Box Coordinates:  (");
		frame.getContentPane().add(lbl_v1se, "flowx,cell 1 19,alignx trailing");
		
		txt_v1seXVal1 = new JTextField();
		txt_v1seXVal1.setColumns(3);
		frame.getContentPane().add(txt_v1seXVal1, "flowx,cell 2 19,alignx center");
		
		lbl_innerBrack2 = new JLabel("-1 )        (");
		frame.getContentPane().add(lbl_innerBrack2, "cell 3 19,alignx left");
		
		txt_v1seXVal2 = new JTextField();
		frame.getContentPane().add(txt_v1seXVal2, "flowx,cell 4 19,alignx center");
		txt_v1seXVal2.setColumns(3);
		
		txt_v1seYVal2 = new JTextField();
		frame.getContentPane().add(txt_v1seYVal2, "flowx,cell 5 19,alignx center");
		txt_v1seYVal2.setColumns(3);
		
		txt_v1seZVal2 = new JTextField();
		txt_v1seZVal2.setColumns(3);
		frame.getContentPane().add(txt_v1seZVal2, "cell 6 19,alignx center");
		
		lbl_outerBrack1 = new JLabel("-1 )");
		frame.getContentPane().add(lbl_outerBrack1, "cell 7 19");
		
		lbl_v2se = new JLabel("V2,SE Box Coordinates:  (");
		frame.getContentPane().add(lbl_v2se, "flowx,cell 1 20,alignx trailing");
		
		txt_v2seXVal1 = new JTextField();
		txt_v2seXVal1.setColumns(3);
		frame.getContentPane().add(txt_v2seXVal1, "flowx,cell 2 20,alignx center");
		
		lbl_comma11 = new JLabel(",");
		frame.getContentPane().add(lbl_comma11, "cell 2 18");
		
		lbl_comma21 = new JLabel(",");
		frame.getContentPane().add(lbl_comma21, "cell 2 19");
		
		lbl_comma31 = new JLabel(",");
		frame.getContentPane().add(lbl_comma31, "cell 2 20");
		
		txt_v1seYVal1 = new JTextField();
		txt_v1seYVal1.setColumns(3);
		frame.getContentPane().add(txt_v1seYVal1, "cell 2 18");
		
		txt_v1seYVal1 = new JTextField();
		txt_v1seYVal1.setColumns(3);
		frame.getContentPane().add(txt_v1seYVal1, "cell 2 19,alignx center");
		
		txt_v2seYVal1 = new JTextField();
		txt_v2seYVal1.setColumns(3);
		frame.getContentPane().add(txt_v2seYVal1, "cell 2 20");
		
		lbl_comma12 = new JLabel(",");
		frame.getContentPane().add(lbl_comma12, "cell 2 18");
		
		lbl_comma22 = new JLabel(",");
		frame.getContentPane().add(lbl_comma22, "cell 2 19");
		
		lbl_comma32 = new JLabel(",");
		frame.getContentPane().add(lbl_comma32, "cell 2 20");
		
		txt_v1seZVal1 = new JTextField();
		frame.getContentPane().add(txt_v1seZVal1, "cell 2 18");
		txt_v1seZVal1.setColumns(3);
		
		txt_v1seZVal1 = new JTextField();
		txt_v1seZVal1.setColumns(3);
		frame.getContentPane().add(txt_v1seZVal1, "cell 2 19");
		
		txt_v2seZVal1 = new JTextField();
		txt_v2seZVal1.setColumns(3);
		frame.getContentPane().add(txt_v2seZVal1, "cell 2 20");
		
		lbl_innerBrack3 = new JLabel("-1 )        (");
		frame.getContentPane().add(lbl_innerBrack3, "cell 3 20,alignx left");
		
		txt_v2seXVal2 = new JTextField();
		txt_v2seXVal2.setColumns(3);
		frame.getContentPane().add(txt_v2seXVal2, "flowx,cell 4 20,alignx center");
		
		txt_v2seYVal2 = new JTextField();
		txt_v2seYVal2.setColumns(3);
		frame.getContentPane().add(txt_v2seYVal2, "flowx,cell 5 20,alignx center");
		
		txt_v2seZVal2 = new JTextField();
		txt_v2seZVal2.setColumns(3);
		frame.getContentPane().add(txt_v2seZVal2, "cell 6 20,alignx center");
		
		lbl_outerBrack2 = new JLabel("-1 )");
		frame.getContentPane().add(lbl_outerBrack2, "cell 7 20");
		
		btn_estRadSpinEcho = new JButton("Estimate Object Radius From Spin Echo");
		frame.getContentPane().add(btn_estRadSpinEcho, "cell 1 21");
		btn_estRadSpinEcho.addActionListener(new Calculate_Magnetic_Moment_3D());
    
		lbl_V0 = new JLabel("V0 =");
		frame.getContentPane().add(lbl_V0, "flowx,cell 2 21");
		
		lbl_rho0SE = new JLabel(RHO +"0,SE =");
		frame.getContentPane().add(lbl_rho0SE, "flowx,cell 3 21");
		
		lbl_echoDChi = new JLabel(_DELTA + CHI + " =");
		frame.getContentPane().add(lbl_echoDChi, "flowx,cell 2 22");
		
		lbl_V0Val = new JLabel("     ");
		frame.getContentPane().add(lbl_V0Val, "cell 2 21");
		
		lbl_V0Unit = new JLabel("pixels" + CUBED);
		frame.getContentPane().add(lbl_V0Unit, "cell 2 21");
		
		lbl_echoDChiVal = new JLabel("");
		frame.getContentPane().add(lbl_echoDChiVal, "cell 2 22");
		
		lbl_rho0SEVal = new JLabel("");
		frame.getContentPane().add(lbl_rho0SEVal, "cell 3 21");
		
		lbl_aSE = new JLabel("a =");
		frame.getContentPane().add(lbl_aSE, "flowx,cell 3 22");
		
		lbl_aSEVal = new JLabel("");
		frame.getContentPane().add(lbl_aSEVal, "cell 3 22");
		
		lbl_comma23 = new JLabel(",");
		frame.getContentPane().add(lbl_comma23, "cell 4 19");
		
		lbl_comma33 = new JLabel(",");
		frame.getContentPane().add(lbl_comma33, "cell 4 20");
		
		lbl_comma24 = new JLabel(",");
		frame.getContentPane().add(lbl_comma24, "cell 5 19");
		
		lbl_comma34 = new JLabel(",");
		frame.getContentPane().add(lbl_comma34, "cell 5 20");
	}
 
 }