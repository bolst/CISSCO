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

  JNIMethods jni = new JNIMethods();

  public static JFrame frame;
  public static JCheckBox chkbx_showrc;
  public static JLabel lbl_stepone, lbl_steptwo, lbl_stepthree, lbl_stepfour, lbl_stepfive, lbl_stepsix,
      lbl_stepseven, lbl_r3phaseAct, lbl_r2phaseAct, lbl_r1phaseAct, lbl_r3phaseCalc, lbl_r2phaseCalc,
      lbl_r1phaseCalc, lbl_r3AphaseUnit, lbl_r2AphaseUnit, lbl_r1AphaseUnit, lbl_r3phaseUnit, lbl_r2phaseUnit,
      lbl_r1phaseUnit, lbl_estBkgPhaseVal, lbl_V0Val, lbl_d_V0Val, lbl_magMom, lbl_rho0, lbl_rho0val, lbl_dchi,
      lbl_dchiVal, lbl_a, lbl_aVal, lbl_MPcnt, lbl_Ri, lbl_ImRi, lbl_ReRi, lbl_echoDChi, lbl_aSE, lbl_err,
      lbl_errVal, lbl_gridSizeBase, lbl_gridSize, lbl_rho0SEVal, lbl_M, lbl_rcx, lbl_rcy, lbl_rcz,
      lbl_rczCorrection, lbl_eqPhase, lbl_eqPhaseUnit, lbl_rc, lbl_rcUnit, lbl_spx, lbl_spy, lbl_spz,
      lbl_spzCorrection, lbl_calculated, lbl_actual, lbl_r1, lbl_r1unit, lbl_r1phase, lbl_r2, lbl_r2phase, lbl_r3,
      lbl_r2unit, lbl_r3unit, lbl_r3phase, lbl_estBkgPhase, lbl_estBkgPhaseUnit, lbl_magMomUnit, lbl_snr,
      lbl_eps12, lbl_eps23, lbl_ReRiVal, lbl_ImRiVal, lbl_TEFirst, lbl_TEFirstUnit, lbl_RChi, lbl_RChiUnit,
      lbl_sigSE, lbl_spinCenter, lbl_innerBrack1, lbl_v1se, lbl_innerBrack2, lbl_v2se, lbl_outerBrack1,
      lbl_comma11, lbl_comma21, lbl_comma31, lbl_v2seYVal1, lbl_comma12, lbl_comma22, lbl_comma32,
      lbl_innerBrack3, lbl_outerBrack2, lbl_V0, lbl_rho0SE, lbl_V0Unit, lbl_echoDChiVal, lbl_aSEVal, lbl_comma23,
      lbl_comma33, lbl_comma24, lbl_comma34;
  public static DefTextField txt_eqPhaseRC, txt_rcx, txt_rcy, txt_rcz, txt_rc, txt_r1,
      txt_r2, txt_r3, txt_spx, txt_spy, txt_spz, txt_snrVal, txt_eps12val, txt_eps23val, txt_B0Val,
      txt_RChiVal, txt_magMomVal, txt_Ri, txt_M, txt_TEFirstVal, txt_TELastVal, txt_spinCenterXVal,
      txt_spinCenterYVal, txt_spinCenterZVal, txt_v1seXVal1, txt_v1seYVal1, txt_v1seZVal1, txt_v1seXVal2,
      txt_v1seYVal2, txt_v1seZVal2, txt_v2seXVal1, txt_v2seYVal1, txt_v2seZVal1, txt_v2seXVal2, txt_v2seYVal2,
      txt_v2seZVal2, txt_sigSEVal;
  public static JButton btn_loadImages, btn_estCR, btn_genSubpix, btn_estSubC, btn_verifyRadii, btn_removeBkg,
      btn_estBkgDens, btn_loadTE, btn_unk, btn_loadspinecho, btn_estRadSpinEcho, btn_redraw, btn_plotX, btn_plotY,
      btn_plotZ, btn_calcMagMom, btn_loadSimImg, btn_sumRi;

  private static ImageItem item;
  public static LogManager logger;
  private static ROIS roiImgMag, roiImgMagXZ, roiImgPhase, roiImgPhaseXZ, roiImgV1SE, roiImgV1SEXZ;
  public static ImagePlus subpixelMagImage, subpixelMagImageXZ, subpixelPhaseImage, subpixelPhaseImageXZ, V1SE_XYImage,
      V1SE_XZImage;
  private static RoiManager rois;
  private static String subCenterErrorMessage,
      subMagTitle, subMagXZTitle, subPhaseTitle, subPhaseXZTitle, V1XY_Title, V1XZ_Title, s1MagWindowTitle,
      s1PhaseWindowTitle, s5MagWindowTitle, s5PhaseWindowTitle, s6MagWindowTitle, s6PhaseWindowTitle, s7WindowTitle;
  private static int drawnRectangle_initialX, drawnRectangle_initialY, drawnRectangle_initialZ, drawnRectangle_sizeX,
      drawnRectangle_sizeY, drawnRectangle_sizeZ;

  private static double m_R0;
  private static float[][] subpixelMagMatrix, subpixelMagMatrixXZ, subpixelPhaseMatrix, subpixelPhaseMatrixXZ;

  public static boolean estimateCenterRadii_isClicked = false;

  private static final int grid = 10;
  private static final double m_ROuterFrom = 0.2;
  private static final double m_RMiddleFrom = 0.9;
  private static final double m_RInnerFrom = 2.5;
  private static final double GAMMARBAR = 42.58;
  private static final String ACCEPTED_FILE_TYPE = "nii";
  private static final String ITALICIZED_I = "\uD835\uDC8A";
  private static final String PLUS_MINUS = "\u00B1";

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

      // clearVariables();
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

      updateVariables();

      try {

        int irand = 0;
        logger.addInfo("estcr", irand++);
        item = new ImageItem(s1MagWindowTitle,
            s1PhaseWindowTitle,
            Integer.parseInt(txt_M.getText()),
            Double.parseDouble(txt_eqPhaseRC.getText()));
        logger.addInfo("estcr", irand++);
        item.calcCenterL();
        logger.addInfo("estcr", irand++);
        item.calcCenterM();
        logger.addInfo("estcr", irand++);
        item.calcCenterS();
        logger.addInfo("estcr", irand++);
        lbl_estBkgPhaseVal.setText(String.valueOf(Math.round(item.estBkg() * 100.0) / 100.0));
        double center_sx = item.centerS().get(0);
        double center_sy = item.centerS().get(1);
        double center_sz = item.centerS().get(2);
        logger.addInfo("csssss", item.centerS());

        // Setting center to GUI unless user already has an inputted center point
        if (txt_rcx.isDefault()) {
          txt_rcx.setText(String.valueOf(center_sx));
        }

        if (txt_rcy.isDefault()) {
          txt_rcy.setText(String.valueOf(center_sy));
        }

        if (txt_rcz.isDefault()) {
          txt_rcz.setText(String.valueOf(center_sz + 1));
        }

        // Getting x y and z centers from GUI
        item.setCenterSX(Double.parseDouble(txt_rcx.getText()));
        item.setCenterSY(Double.parseDouble(txt_rcy.getText()));
        item.setCenterSZ(Double.parseDouble(txt_rcz.getText()) - 1);

        // Getting RCenter from estimated xyz
        // RCenter = estimateRCenter((int)(item.centerS().get(0)),
        // (int)(item.centerS().get(1)),
        // (int)(item.centerS().get(2)));
        double RCenter = item.estimateRCenter();

        // Updating GUI
        txt_rc.setText(String.valueOf(Math.round(RCenter * 10.0) / 10.0));

        item.calcR0123();
        m_R0 = item.m_R0();

        jni.setmVariables(grid, m_R0, RCenter,
            Double.parseDouble(txt_rcx.getText()),
            Double.parseDouble(txt_rcy.getText()),
            Double.parseDouble(txt_rcx.getText()) - 1.0,
            Double.parseDouble(txt_eqPhaseRC.getText()));
        jni.setMagMoment(Double.parseDouble(txt_eqPhaseRC.getText()) * Math.pow(RCenter, 3));

        estimateCenterRadii_isClicked = true;

      } catch (Exception exc) {
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
      else if (!item.isNearEdge) {
        JOptionPane.showMessageDialog(frame, "Error: Object too close to edge of image.");
      } else {

        logger.addVariable("s1magwt", s1MagWindowTitle);
        logger.addVariable("s2phasewt", s1PhaseWindowTitle);
        // Getting mag and phase images
        ImagePlus magnitudeImage = WindowManager.getImage(s1MagWindowTitle);
        ImagePlus phaseImage = WindowManager.getImage(s1PhaseWindowTitle);
        logger.addInfo("Got images");

        // Size of new images
        int size_subpixelimg = (int) ((2 * m_R0 + 1) * grid);
        logger.addVariable("sizeofspi", size_subpixelimg);

        // Initializing 3D arrays to give to C++
        float[][][] croppedMagnitudeValues3D = new float[size_subpixelimg][size_subpixelimg][size_subpixelimg];
        float[][][] croppedPhaseValues3D = new float[size_subpixelimg][size_subpixelimg][size_subpixelimg];
        float[][][] croppedRealNumbers3D = new float[size_subpixelimg][size_subpixelimg][size_subpixelimg];
        float[][][] croppedImaginaryNumbers3D = new float[size_subpixelimg][size_subpixelimg][size_subpixelimg];

        // Initializing arrays that will be displayed in all 4 images
        subpixelMagMatrix = new float[size_subpixelimg][size_subpixelimg];
        subpixelMagMatrixXZ = new float[size_subpixelimg][size_subpixelimg];
        subpixelPhaseMatrix = new float[size_subpixelimg][size_subpixelimg];
        subpixelPhaseMatrixXZ = new float[size_subpixelimg][size_subpixelimg];
        logger.addInfo("made arrays");

        // Initial points of subpixel images
        int xi = (item.centerS().get(0).intValue()) - (int) m_R0;
        int yi = (item.centerS().get(1).intValue()) - (int) m_R0;
        int zi = (item.centerS().get(2).intValue()) - (int) m_R0;
        logger.addVariable("xi", xi);
        logger.addVariable("yi", yi);
        logger.addVariable("zi", zi);
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
        double RCenter = Double.parseDouble(txt_rc.getText());
        jni.setmVariables(grid, m_R0, RCenter,
            Double.parseDouble(txt_rcx.getText()),
            Double.parseDouble(txt_rcy.getText()),
            Double.parseDouble(txt_rcx.getText()) - 1.0,
            Double.parseDouble(txt_eqPhaseRC.getText()));
        jni.setBackPhase(item.bkgPhase);
        jni.setRealImagNumbers(croppedRealNumbers3D, croppedImaginaryNumbers3D);

        // Generate subpixel in C++
        jni.generateSubpixelArray();

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
        ImageProcessor IP_subpixelMagImage = new FloatProcessor(size_subpixelimg, size_subpixelimg);
        ImageProcessor IP_subpixelMagImageXZ = new FloatProcessor(size_subpixelimg, size_subpixelimg);
        ImageProcessor IP_subpixelPhaseImage = new FloatProcessor(size_subpixelimg, size_subpixelimg);
        ImageProcessor IP_subpixelPhaseImageXZ = new FloatProcessor(size_subpixelimg, size_subpixelimg);

        // Adding mag and phase data to image processors
        for (int i = 0; i < size_subpixelimg; i++) {
          for (int j = 0; j < size_subpixelimg; j++) {
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

        // Passing matrices to C++
        jni.setRealImagNumbers(croppedRealNumbers3D, croppedImaginaryNumbers3D);
        jni.setPhaseXYMatrix(subpixelPhaseMatrix);
        jni.setPhaseXZMatrix(subpixelPhaseMatrixXZ);
        jni.setMagXYMatrix(subpixelMagMatrix);
        jni.setMagXZMatrix(subpixelMagMatrixXZ);
      }
    }

    /*
     * If the remove background phase button is clicked:
     * The program takes the subpixel XY and XZ matrices and removes the estimated
     * background phase
     */
    else if (e.getSource() == btn_removeBkg) {

      updateVariables();

      if (!!lbl_estBkgPhaseVal.getText().isEmpty()) {
        JOptionPane.showMessageDialog(frame, "Error: No background phase found.");
      } else

      {
        // Removing BG phase in C++
        jni.removeBackgroundPhase(item.bkgPhase);

        // Removing BG phase in Java
        removeBGPhase(subpixelPhaseMatrix);
        removeBGPhase(subpixelPhaseMatrixXZ);
        // removeBGPhase(item.bkgPhase);
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
      double RCenter = Double.parseDouble(txt_rc.getText());
      jni.setmVariables(grid, m_R0, RCenter,
          Double.parseDouble(txt_rcx.getText()),
          Double.parseDouble(txt_rcy.getText()),
          Double.parseDouble(txt_rcx.getText()) - 1.0,
          Double.parseDouble(txt_eqPhaseRC.getText()));

      // condition for program to continue, must have generated subpixel and estimated
      // a center and RCenter
      boolean condition = !lbl_rcx.getText().isEmpty() && !lbl_rcy.getText().isEmpty() && !lbl_rcz.getText().isEmpty()
          && !lbl_rc.getText().isEmpty()
          && !lbl_estBkgPhaseVal.getText().isEmpty()
          && (WindowManager.getImage(subMagTitle) != null && WindowManager.getImage(subMagXZTitle) != null
              && WindowManager.getImage(subPhaseTitle) != null && WindowManager.getImage(subPhaseXZTitle) != null);

      if (condition) {

        logger.addInfo("rcx", txt_rcx.getText());
        logger.addInfo("rcy", txt_rcy.getText());
        logger.addInfo("rcz", txt_rcz.getText() + "-1");
        logger.addInfo("ibix", item.roi_mag_belowM_xi);
        logger.addInfo("ibiy", item.roi_mag_belowM_yi);
        logger.addInfo("ibiz", item.roi_mag_belowM_zi);
        logger.addInfo("ibsx", item.roi_mag_belowM_Dx);
        logger.addInfo("ibsy", item.roi_mag_belowM_Dy);
        logger.addInfo("ibsz", item.roi_mag_belowM_Dz);
        logger.addInfo("cl", item.centerL());
        logger.addInfo("cm", item.centerM());
        logger.addInfo("cs", item.centerS());

        // Passing necessary data to C++
        jni.setXYZ(Double.parseDouble(txt_rcx.getText()), Double.parseDouble(txt_rcy.getText()),
            Double.parseDouble(txt_rcz.getText()) - 1.0);
        jni.setPhaseXYMatrix(subpixelPhaseMatrix);
        jni.setSmallBox(item.roi_mag_belowM_xi, item.roi_mag_belowM_yi, item.roi_mag_belowM_zi, item.roi_mag_belowM_Dx,
            item.roi_mag_belowM_Dy,
            item.roi_mag_belowM_Dz);
        jni.setCenterL(item.centerL().get(0), item.centerL().get(1), item.centerL().get(2));
        jni.setCenterM(item.centerM().get(0), item.centerM().get(1), item.centerM().get(2));
        jni.setCenterS(item.centerS().get(0), item.centerS().get(1), item.centerS().get(2));

        // Calculating subpixel center, if there are no errors then the returned string
        // will be empty
        subCenterErrorMessage = jni.estimateSubpixelCenter();
        logger.addInfo("hi");

        if (subCenterErrorMessage.compareTo("") == 0) {

          // Getting estimated subpixel centers from C++ in terms of pixels
          double centerX_pixelCoordinates = jni.getSubX();
          double centerY_pixelCoordinates = jni.getSubY();
          double centerZ_pixelCoordinates = jni.getSubZ();

          // Update center_s to be the estimated subpixel center
          item.centerS().set(0, centerX_pixelCoordinates);
          item.centerS().set(1, centerY_pixelCoordinates);
          item.centerS().set(2, centerZ_pixelCoordinates);

          // Setting pixel values to GUI, because the coordinates are relative to the
          // uploaded images
          txt_spx.setText(String.valueOf(Math.round(centerX_pixelCoordinates * 100.0) / 100.0));
          txt_spy.setText(String.valueOf(Math.round(centerY_pixelCoordinates * 100.0) / 100.0));
          txt_spz.setText(String.valueOf(Math.round(centerZ_pixelCoordinates * 100.0) / 100.0 + 1.0));

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
          roiImgMag.addPointROI("Center", (int) pixelToSubpixel(item.centerS().get(0), 0),
              (int) pixelToSubpixel(item.centerS().get(1), 1));
          if (chkbx_showrc.isSelected()) {
            // If the box is selected, adds RCenter circle ROI to the list
            roiImgMag.addCircleROI("RCenter", (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(1), 1), RCenter * 10.0);
          }
          // Displays the ROIS on the image
          roiImgMag.displayROIS();

          // The next 3 blocks of code follow the same logic as described above

          roiImgMagXZ = new ROIS("MXZ");
          roiImgMagXZ.addPointROI("Center", (int) pixelToSubpixel(item.centerS().get(0), 0),
              (int) pixelToSubpixel(item.centerS().get(2), 2));
          if (chkbx_showrc.isSelected()) {
            roiImgMagXZ.addCircleROI("RCenter", (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(2), 2),
                RCenter * 10.0);
          }
          roiImgMagXZ.displayROIS();

          roiImgPhase = new ROIS("PXY");
          roiImgPhase.addPointROI("Center", (int) pixelToSubpixel(item.centerS().get(0), 0),
              (int) pixelToSubpixel(item.centerS().get(1), 1));
          if (chkbx_showrc.isSelected()) {
            roiImgPhase.addCircleROI("RCenter", (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(1), 1),
                RCenter * 10.0);
          }
          roiImgPhase.displayROIS();

          roiImgPhaseXZ = new ROIS("PXZ");
          roiImgPhaseXZ.addPointROI("Center", (int) pixelToSubpixel(item.centerS().get(0), 0),
              (int) pixelToSubpixel(item.centerS().get(2), 2));
          if (chkbx_showrc.isSelected()) {
            roiImgPhaseXZ.addCircleROI("RCenter", (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(2), 2),
                RCenter * 10.0);
          }
          roiImgPhaseXZ.displayROIS();

          // ---------- end to put ROIS on images

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
      double RCenter = Double.parseDouble(txt_rc.getText());
      jni.setmVariables(grid, m_R0, RCenter,
          Double.parseDouble(txt_rcx.getText()),
          Double.parseDouble(txt_rcy.getText()),
          Double.parseDouble(txt_rcx.getText()) - 1.0,
          Double.parseDouble(txt_eqPhaseRC.getText()));

      boolean condition = (WindowManager.getImage(subMagTitle) != null && WindowManager.getImage(subMagXZTitle) != null
          && WindowManager.getImage(subPhaseTitle) != null && WindowManager.getImage(subPhaseXZTitle) != null)
          && !(txt_spx.getText().isEmpty() || txt_spy.getText().isEmpty() || txt_spz.getText().isEmpty());

      if (condition) {

        // This is basically a repetition from the btn_estSubC code

        roiImgMag.clear();
        roiImgMag.addPointROI("Center", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(1), 1));
        if (chkbx_showrc.isSelected()) {
          roiImgMag.addCircleROI("RCenter", (int) pixelToSubpixel(item.centerS().get(0), 0),
              (int) pixelToSubpixel(item.centerS().get(1), 1), RCenter * 10.0);
        }
        roiImgMag.displayROIS();

        roiImgMagXZ.clear();
        roiImgMagXZ.addPointROI("Center", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(2), 2));
        if (chkbx_showrc.isSelected()) {
          roiImgMagXZ.addCircleROI("RCenter", (int) pixelToSubpixel(item.centerS().get(0), 0),
              (int) pixelToSubpixel(item.centerS().get(2), 2), RCenter * 10.0);
        }
        roiImgMagXZ.displayROIS();

        roiImgPhase.clear();
        roiImgPhase.addPointROI("Center", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(1), 1));
        if (chkbx_showrc.isSelected()) {
          roiImgPhase.addCircleROI("RCenter", (int) pixelToSubpixel(item.centerS().get(0), 0),
              (int) pixelToSubpixel(item.centerS().get(1), 1), RCenter * 10.0);
        }
        roiImgPhase.displayROIS();

        roiImgPhaseXZ.clear();
        roiImgPhaseXZ.addPointROI("Center", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(2), 2));
        if (chkbx_showrc.isSelected()) {
          roiImgPhaseXZ.addCircleROI("RCenter", (int) pixelToSubpixel(item.centerS().get(0), 0),
              (int) pixelToSubpixel(item.centerS().get(2), 2),
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

      boolean condition = !(txt_spx.getText().isEmpty() || txt_spy.getText().isEmpty() || txt_spz.getText().isEmpty())
          && (WindowManager.getImage(subMagTitle) != null && WindowManager.getImage(subMagXZTitle) != null
              && WindowManager.getImage(subPhaseTitle) != null && WindowManager.getImage(subPhaseXZTitle) != null);

      if (condition) {

        // Initializing plot object for graph
        Plot xPlot = new Plot("X-Profile", "Location", "Phase");

        // Initializing ArrayLists for graph values on each axis
        ArrayList<Double> intensity = new ArrayList<Double>();
        ArrayList<Double> location = new ArrayList<Double>();

        // Setting colour components of graph
        xPlot.setColor(Color.RED, Color.BLACK);

        try {

          int x_px = item.centerS().get(0).intValue();
          int y_spx = (int) pixelToSubpixel(item.centerS().get(1).intValue(), 1);
          // Adding phase values to both ArrayLists
          for (int i = x_px - (int) m_R0,
              c = 0; i <= x_px + (int) m_R0; i++, c++) {
            intensity
                .add((double) subpixelPhaseMatrix[(int) pixelToSubpixel(i, 0)][y_spx]);
            location.add((double) (i - x_px));
            if (i != x_px - (int) m_R0)
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

      boolean condition = !(txt_spx.getText().isEmpty() || txt_spy.getText().isEmpty() || txt_spz.getText().isEmpty())
          && (WindowManager.getImage(subMagTitle) != null && WindowManager.getImage(subMagXZTitle) != null
              && WindowManager.getImage(subPhaseTitle) != null && WindowManager.getImage(subPhaseXZTitle) != null);

      if (condition) {

        // Initializing plot object for graph
        Plot yPlot = new Plot("Y-Profile", "Location", "Phase");

        // Initializing ArrayLists for graph values on each axis
        ArrayList<Double> intensity = new ArrayList<Double>();
        ArrayList<Double> location = new ArrayList<Double>();

        // Setting colour components of graph
        yPlot.setColor(Color.BLUE, Color.BLACK);

        try {

          int x_spx = (int) pixelToSubpixel(item.centerS().get(0).intValue(), 0);
          int y_px = item.centerS().get(1).intValue();
          // Adding phase values to both ArrayLists
          for (int i = y_px - (int) m_R0,
              c = 0; i <= y_px + m_R0; i++, c++) {
            intensity
                .add((double) subpixelPhaseMatrix[x_spx][(int) pixelToSubpixel(i, 1)]);
            location.add((double) (i - y_px));
            if (i != y_px - (int) m_R0)
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

      boolean condition = !(txt_spx.getText().isEmpty() || txt_spy.getText().isEmpty() || txt_spz.getText().isEmpty())
          && (WindowManager.getImage(subMagTitle) != null && WindowManager.getImage(subMagXZTitle) != null
              && WindowManager.getImage(subPhaseTitle) != null && WindowManager.getImage(subPhaseXZTitle) != null);

      if (condition) {

        // Initializing plot object for graph
        Plot zPlot = new Plot("Z-Profile", "Location", "Phase");

        // Initializing ArrayLists for graph values on each axis
        ArrayList<Double> intensity = new ArrayList<Double>();
        ArrayList<Double> location = new ArrayList<Double>();

        // Setting colour components of graph
        zPlot.setColor(Color.GREEN, Color.BLACK);

        try {

          int x_spx = (int) pixelToSubpixel(item.centerS().get(0).intValue(), 0);
          int z_px = item.centerS().get(2).intValue();

          // Adding phase values to both ArrayLists
          for (int i = z_px - (int) m_R0,
              c = 0; i <= z_px + m_R0; i++, c++) {
            intensity
                .add((double) subpixelPhaseMatrixXZ[x_spx][(int) pixelToSubpixel(i, 2)]);
            location.add((double) (i - z_px));
            if (i != z_px - (int) m_R0) {
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
      double m_R1 = Double.parseDouble(txt_r1.getText());
      double m_R2 = Double.parseDouble(txt_r2.getText());
      double m_R3 = Double.parseDouble(txt_r3.getText());

      boolean condition = (WindowManager.getImage(subMagTitle) != null && WindowManager.getImage(subMagXZTitle) != null
          && WindowManager.getImage(subPhaseTitle) != null && WindowManager.getImage(subPhaseXZTitle) != null)
          && !(txt_spx.getText().isEmpty() || txt_spy.getText().isEmpty() || txt_spz.getText().isEmpty());

      if (condition) {

        // Clearing ROI list
        roiImgMag.clear();
        // Adding center ROI to list
        roiImgMag.addPointROI("Center", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(1), 1));
        // Adding R1 to list
        roiImgMag.addCircleROI("MR1", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(1), 1),
            m_R1 * 10.0);
        // Adding R2 to list
        roiImgMag.addCircleROI("MR2", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(1), 1),
            m_R2 * 10.0);
        // Adding R3 to list
        roiImgMag.addCircleROI("MR3", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(1), 1),
            m_R3 * 10.0);
        // Displaying list
        roiImgMag.displayROIS();

        // The same logic as above is followed for the next three blocks

        roiImgMagXZ.clear();
        roiImgMagXZ.addPointROI("Center", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(2), 2));
        roiImgMagXZ.addCircleROI("MR1", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(2), 2),
            m_R1 * 10.0);
        roiImgMagXZ.addCircleROI("MR2", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(2), 2),
            m_R2 * 10.0);
        roiImgMagXZ.addCircleROI("MR3", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(2), 2),
            m_R3 * 10.0);
        roiImgMagXZ.displayROIS();

        roiImgPhase.clear();
        roiImgPhase.addPointROI("Center", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(1), 1));
        roiImgPhase.addCircleROI("MR1", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(1), 1),
            m_R1 * 10.0);
        roiImgPhase.addCircleROI("MR2", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(1), 1),
            m_R2 * 10.0);
        roiImgPhase.addCircleROI("MR3", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(1), 1),
            m_R3 * 10.0);
        roiImgPhase.displayROIS();

        roiImgPhaseXZ.clear();
        roiImgPhaseXZ.addPointROI("Center", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(2), 2));
        roiImgPhaseXZ.addCircleROI("MR1", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(2), 2),
            m_R1 * 10.0);
        roiImgPhaseXZ.addCircleROI("MR2", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(2), 2),
            m_R2 * 10.0);
        roiImgPhaseXZ.addCircleROI("MR3", (int) pixelToSubpixel(item.centerS().get(0), 0),
            (int) pixelToSubpixel(item.centerS().get(2), 2),
            m_R3 * 10.0);
        roiImgPhaseXZ.displayROIS();

        float R1_phase_actual = 0;
        float R2_phase_actual = 0;
        float R3_phase_actual = 0;

        // Summing up all phase values where the radii and equitorial axis intercept
        switch (item.neglectedAxis()) {
          case "x":
            R1_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(1), 1) + (int) m_R1 * 10);
            logger.addInfo("Got coordinate (" + String.valueOf((int) pixelToSubpixel(item.centerS().get(0), 0)) + ","
                + String.valueOf((int) pixelToSubpixel(item.centerS().get(1), 1) + (int) m_R1 * 10) + ")");
            R1_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(1), 1) - (int) m_R1 * 10);
            logger.addInfo("Got coordinate (" + String.valueOf((int) pixelToSubpixel(item.centerS().get(0), 0)) + ","
                + String.valueOf((int) pixelToSubpixel(item.centerS().get(1), 1) - (int) m_R1 * 10) + ")");
            R1_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(2), 2) + (int) m_R1 * 10);
            logger.addInfo("Got coordinate (" + String.valueOf((int) pixelToSubpixel(item.centerS().get(0), 0)) + ","
                + String.valueOf((int) pixelToSubpixel(item.centerS().get(2), 2) + (int) m_R1 * 10) + ")");
            R1_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(2), 2) - (int) m_R1 * 10);
            logger.addInfo("Got coordinate (" + String.valueOf((int) pixelToSubpixel(item.centerS().get(0), 0)) + ","
                + String.valueOf((int) pixelToSubpixel(item.centerS().get(2), 2) - (int) m_R1 * 10) + ")");

            R2_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(1), 1) + (int) m_R2 * 10);
            logger.addInfo("Got coordinate (" + String.valueOf((int) pixelToSubpixel(item.centerS().get(0), 0)) + ","
                + String.valueOf((int) pixelToSubpixel(item.centerS().get(1), 1) + (int) m_R2 * 10) + ")");
            R2_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(1), 1) - (int) m_R2 * 10);
            logger.addInfo("Got coordinate (" + String.valueOf((int) pixelToSubpixel(item.centerS().get(0), 0)) + ","
                + String.valueOf((int) pixelToSubpixel(item.centerS().get(1), 1) - (int) m_R2 * 10) + ")");
            R2_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(2), 2) + (int) m_R2 * 10);
            logger.addInfo("Got coordinate (" + String.valueOf((int) pixelToSubpixel(item.centerS().get(0), 0)) + ","
                + String.valueOf((int) pixelToSubpixel(item.centerS().get(2), 2) + (int) m_R2 * 10) + ")");
            R2_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(2), 2) - (int) m_R2 * 10);
            logger.addInfo("Got coordinate (" + String.valueOf((int) pixelToSubpixel(item.centerS().get(0), 0)) + ","
                + String.valueOf((int) pixelToSubpixel(item.centerS().get(2), 2) - (int) m_R2 * 10) + ")");

            R3_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(1), 1) + (int) m_R3 * 10);
            logger.addInfo("Got coordinate (" + String.valueOf((int) pixelToSubpixel(item.centerS().get(0), 0)) + ","
                + String.valueOf((int) pixelToSubpixel(item.centerS().get(1), 1) + (int) m_R3 * 10) + ")");
            R3_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(1), 1) - (int) m_R3 * 10);
            logger.addInfo("Got coordinate (" + String.valueOf((int) pixelToSubpixel(item.centerS().get(0), 0)) + ","
                + String.valueOf((int) pixelToSubpixel(item.centerS().get(1), 1) - (int) m_R3 * 10) + ")");
            R3_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(2), 2) + (int) m_R3 * 10);
            logger.addInfo("Got coordinate (" + String.valueOf((int) pixelToSubpixel(item.centerS().get(0), 0)) + ","
                + String.valueOf((int) pixelToSubpixel(item.centerS().get(2), 2) + (int) m_R3 * 10) + ")");
            R3_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(2), 2) - (int) m_R3 * 10);
            logger.addInfo("Got coordinate (" + String.valueOf((int) pixelToSubpixel(item.centerS().get(0), 0)) + ","
                + String.valueOf((int) pixelToSubpixel(item.centerS().get(2), 2) - (int) m_R3 * 10) + ")");
            break;

          case "y":
            R1_phase_actual += subpixelPhaseImageXZ.getProcessor()
                .getPixelValue((int) pixelToSubpixel(item.centerS().get(0), 0) + (int) m_R1 * 10,
                    (int) pixelToSubpixel(item.centerS().get(2), 2));
            R1_phase_actual += subpixelPhaseImageXZ.getProcessor()
                .getPixelValue((int) pixelToSubpixel(item.centerS().get(0), 0) - (int) m_R1 * 10,
                    (int) pixelToSubpixel(item.centerS().get(2), 2));
            R1_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(2), 2) + (int) m_R1 * 10);
            R1_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(2), 2) - (int) m_R1 * 10);

            R2_phase_actual += subpixelPhaseImageXZ.getProcessor()
                .getPixelValue((int) pixelToSubpixel(item.centerS().get(0), 0) + (int) m_R2 * 10,
                    (int) pixelToSubpixel(item.centerS().get(2), 2));
            R2_phase_actual += subpixelPhaseImageXZ.getProcessor()
                .getPixelValue((int) pixelToSubpixel(item.centerS().get(0), 0) - (int) m_R2 * 10,
                    (int) pixelToSubpixel(item.centerS().get(2), 2));
            R2_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(2), 2) + (int) m_R2 * 10);
            R2_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(2), 2) - (int) m_R2 * 10);

            R3_phase_actual += subpixelPhaseImageXZ.getProcessor()
                .getPixelValue((int) pixelToSubpixel(item.centerS().get(0), 0) + (int) m_R3 * 10,
                    (int) pixelToSubpixel(item.centerS().get(2), 2));
            R3_phase_actual += subpixelPhaseImageXZ.getProcessor()
                .getPixelValue((int) pixelToSubpixel(item.centerS().get(0), 0) - (int) m_R3 * 10,
                    (int) pixelToSubpixel(item.centerS().get(2), 2));
            R3_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(2), 2) + (int) m_R3 * 10);
            R3_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(2), 2) - (int) m_R3 * 10);
            break;

          case "z":
            R1_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(1), 1) + (int) m_R1 * 10);
            R1_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(1), 1) - (int) m_R1 * 10);
            R1_phase_actual += subpixelPhaseImage.getProcessor()
                .getPixelValue((int) pixelToSubpixel(item.centerS().get(0), 0) + (int) m_R1 * 10,
                    (int) pixelToSubpixel(item.centerS().get(1), 1));
            R1_phase_actual += subpixelPhaseImage.getProcessor()
                .getPixelValue((int) pixelToSubpixel(item.centerS().get(0), 0) - (int) m_R1 * 10,
                    (int) pixelToSubpixel(item.centerS().get(1), 1));

            R2_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(1), 1) + (int) m_R2 * 10);
            R2_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(1), 1) - (int) m_R2 * 10);
            R2_phase_actual += subpixelPhaseImage.getProcessor()
                .getPixelValue((int) pixelToSubpixel(item.centerS().get(0), 0) + (int) m_R2 * 10,
                    (int) pixelToSubpixel(item.centerS().get(1), 1));
            R2_phase_actual += subpixelPhaseImage.getProcessor()
                .getPixelValue((int) pixelToSubpixel(item.centerS().get(0), 0) - (int) m_R2 * 10,
                    (int) pixelToSubpixel(item.centerS().get(1), 1));

            R3_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(1), 1) + (int) m_R3 * 10);
            R3_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(
                (int) pixelToSubpixel(item.centerS().get(0), 0),
                (int) pixelToSubpixel(item.centerS().get(1), 1) - (int) m_R3 * 10);
            R3_phase_actual += subpixelPhaseImage.getProcessor()
                .getPixelValue((int) pixelToSubpixel(item.centerS().get(0), 0) + (int) m_R3 * 10,
                    (int) pixelToSubpixel(item.centerS().get(1), 1));
            R3_phase_actual += subpixelPhaseImage.getProcessor()
                .getPixelValue((int) pixelToSubpixel(item.centerS().get(0), 0) - (int) m_R3 * 10,
                    (int) pixelToSubpixel(item.centerS().get(1), 1));
            break;

          default:
            JOptionPane.showMessageDialog(frame, "Error in Verify Radii: No MRI field axis found");
            break;
        }

        // Dividing each phase sum by 4 for average
        R1_phase_actual /= 4.0;
        R2_phase_actual /= 4.0;
        R3_phase_actual /= 4.0;
        // Removing background phase off each phase value
        R1_phase_actual -= item.bkgPhase;
        R2_phase_actual -= item.bkgPhase;
        R3_phase_actual -= item.bkgPhase;

        lbl_r1phaseAct.setText(String.valueOf(Math.round(R1_phase_actual * 100.0) / 100.0));
        lbl_r2phaseAct.setText(String.valueOf(Math.round(R2_phase_actual * 100.0) / 100.0));
        lbl_r3phaseAct.setText(String.valueOf(Math.round(R3_phase_actual * 100.0) / 100.0));
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
      boolean condition = !(txt_r1.getText().isEmpty()) && !(txt_r2.getText().isEmpty())
          && !(txt_r3.getText().isEmpty())
          && !(txt_spx.getText().isEmpty() || txt_spy.getText().isEmpty() || txt_spz.getText().isEmpty())
          && (WindowManager.getImage(subMagTitle) != null && WindowManager.getImage(subMagXZTitle) != null
              && WindowManager.getImage(subPhaseTitle) != null && WindowManager.getImage(subPhaseXZTitle) != null);

      if (condition) {

        // Calling C++ to calculate background phase and spin density
        jni.estBkgAndSpinDensity();

        // Getting calculated background phase from C++
        item.bkgPhase = Math.abs(jni.getBkg());

        // Getting calculated spin density from C++
        double spinDensity = jni.getSpinDensity();

        logger.addVariable("estimatedBackgroundPhase", item.bkgPhase);
        logger.addVariable("spinDensity", spinDensity);

        // Setting background phase and spin density to GUI
        lbl_estBkgPhaseVal.setText(String.valueOf(Math.round(item.bkgPhase * 100.0) / 100.0));
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
      double snr = Double.parseDouble(txt_snrVal.getText());
      double e12 = Double.parseDouble(txt_eps12val.getText());
      double e23 = Double.parseDouble(txt_eps23val.getText());
      double B0 = Double.parseDouble(txt_B0Val.getText());
      double R_Chi = Double.parseDouble(txt_RChiVal.getText());
      double TElast = Double.parseDouble(txt_TELastVal.getText());
      jni.setMagMomentVariables(snr, e12, e23, B0, R_Chi, TElast);

      double m_R1 = Double.parseDouble(txt_r1.getText());
      double m_R2 = Double.parseDouble(txt_r2.getText());
      double m_R3 = Double.parseDouble(txt_r3.getText());

      // ---------- Begin to open simulated images, similar to step 1

      final JFileChooser simulatedImageChooserWindow = new JFileChooser(
          FileSystemView.getFileSystemView()
              .getHomeDirectory());
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

      // End loading simulated images

      if ((WindowManager.getImage(s5MagWindowTitle) != null) && (WindowManager.getImage(s5PhaseWindowTitle) != null)) {

        ImagePlus simulatedMagnitudeImage = WindowManager.getImage(s5MagWindowTitle);
        ImagePlus simulatedPhaseImage = WindowManager.getImage(s5PhaseWindowTitle);

        if (simulatedMagnitudeImage.getNSlices() != simulatedPhaseImage.getNSlices()) {
          JOptionPane.showMessageDialog(frame, "Error: different # of slices in images");
        } else if (simulatedMagnitudeImage.getHeight() != simulatedPhaseImage.getHeight()) {
          JOptionPane.showMessageDialog(frame, "Error: different height in images");
        } else if (simulatedMagnitudeImage.getWidth() != simulatedPhaseImage.getWidth()) {
          JOptionPane.showMessageDialog(frame, "Error: different width in images");
        } else if (!!(txt_r1.getText().isEmpty())) {
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
          jni.setSimulatedMatrices(simulatedRealNumbers, sizeOfSimmedMatrices + 1);
          logger.addInfo("Set matrices to C++");
          // Interpolating matrix in C++
          jni.interpolateVoxelsSIM(sizeOfSimmedMatrices * 10);
          logger.addInfo("Interpolation complete");

          // Getting current rho_0 value from GUI
          double rho_0 = Double.parseDouble(lbl_rho0val.getText());
          // double rho_0 = 10;

          // Getting center of simulated matrix
          int subCenterSIM = sizeOfSimmedMatrices * 10 / 2;

          // Summing real values in each radius in C++
          double S1 = jni.SumCircleElementsReal3DSIMMED((int) (m_R1 * 10.0), subCenterSIM, subCenterSIM, subCenterSIM);
          double S2 = jni.SumCircleElementsReal3DSIMMED((int) (m_R2 * 10.0), subCenterSIM, subCenterSIM, subCenterSIM);
          double S3 = jni.SumCircleElementsReal3DSIMMED((int) (m_R3 * 10.0), subCenterSIM, subCenterSIM, subCenterSIM);

          // Getting phi (phase) values from GUI
          double phi1 = Double.parseDouble(lbl_r1phaseCalc.getText());
          double phi2 = Double.parseDouble(lbl_r2phaseCalc.getText());
          double phi3 = Double.parseDouble(lbl_r3phaseCalc.getText());

          // Calculating fij using equation 9 (real)
          double Re_f12_REAL = (9.0 * Math.sqrt(3)) / (4.0 * Math.PI * rho_0) * (S1 - S2);
          double Re_f23_REAL = (9.0 * Math.sqrt(3)) / (4.0 * Math.PI * rho_0) * (S2 - S3);

          // Calculating fij using equation 10 in C++ (theoretical)
          double mag_moment = Double.parseDouble(txt_magMomVal.getText());
          double Re_f12_THEORETICAL = jni.equation10(mag_moment, phi1, phi2);
          double Re_f23_THEORETICAL = jni.equation10(mag_moment, phi2, phi3);

          logger.addVariable("S1", S1);
          logger.addVariable("S2", S2);
          logger.addVariable("S3", S3);
          logger.addVariable("Re_f12_REAL", Re_f12_REAL);
          logger.addVariable("Re_f23_REAL", Re_f23_REAL);
          logger.addVariable("Re_f12_THEORETICAL", Re_f12_THEORETICAL);
          logger.addVariable("Re_f23_THEORETICAL", Re_f23_THEORETICAL);

          // Calculating eij by doing [(real - theoretical) / theoretical]
          double se12 = Math.abs(Re_f12_REAL - Re_f12_THEORETICAL) / Re_f12_THEORETICAL;
          double se23 = Math.abs(Re_f23_REAL - Re_f23_THEORETICAL) / Re_f23_THEORETICAL;

          // Setting eij to GUI
          txt_eps12val.setText(String.valueOf(Math.round(se12 * 100.0) / 100.0));
          txt_eps23val.setText(String.valueOf(Math.round(se23 * 100.0) / 100.0));

          // Calculating uncertainty in C++
          double uncertainty = jni.calculateUncertainty(se12, se23);

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
      boolean condition = !(txt_r1.getText().isEmpty()) && !(txt_r2.getText().isEmpty())
          && !(txt_r3.getText().isEmpty())
          && !(txt_spx.getText().isEmpty() || txt_spy.getText().isEmpty() || txt_spz.getText().isEmpty())
          && (WindowManager.getImage(subMagTitle) != null && WindowManager.getImage(subMagXZTitle) != null
              && WindowManager.getImage(subPhaseTitle) != null && WindowManager.getImage(subPhaseXZTitle) != null)
          && !txt_rc.getText().isEmpty();

      if (condition) {
        String errorMessage_Mag;

        // Calculating magnetic moment in C++, if function returns empty string then no
        // error message
        errorMessage_Mag = jni.calculateMagneticMoment();

        if (errorMessage_Mag.compareTo("") == 0) {
          // Getting various values from C++ as a result of the function to calculate the
          // magnetic moment being ran
          lbl_r1phaseCalc.setText(String.valueOf(Math.round(jni.getMR1Calc() * 100.0) / 100.0));
          lbl_r2phaseCalc.setText(String.valueOf(Math.round(jni.getMR2Calc() * 100.0) / 100.0));
          lbl_r3phaseCalc.setText(String.valueOf(Math.round(jni.getMR3Calc() * 100.0) / 100.0));
          txt_magMomVal.setText(String.valueOf(Math.round(jni.getMagMoment() * 100.0) / 100.0));
          if (jni.getUncertainty() == -1.0) {
            lbl_errVal.setText("");
            JOptionPane.showMessageDialog(frame,
                "<html>Error: Cannot calculate error\nMake sure SNR, &epsilon;12 and &epsilon;23 are set.");
          } else {
            lbl_errVal.setText(String.valueOf(Math.round(jni.getUncertainty() * 100.0) / 100.0));
          }
          lbl_dchiVal.setText(String.valueOf(Math.round(jni.getChi() * 100.0) / 100.0) + " ppm");
          lbl_aVal.setText(String.valueOf(Math.round(jni.getA() * 100.0) / 100.0) + " pixels");
          lbl_rho0val.setText(String.valueOf(Math.round(jni.getSpinDensity() * 100.0) / 100.0));
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
      double m_Ri = Double.parseDouble(txt_Ri.getText());
      jni.setRi(m_Ri);

      String Imag_errmsg;
      String Real_errmsg;

      Imag_errmsg = jni.calculateImagSum();
      Real_errmsg = jni.calculateRealSum();

      if (Imag_errmsg.compareTo("") == 0) {
        double m_Si = jni.getImagSum();
        lbl_ImRi.setText("S" + ITALICIZED_I + "= " + String.valueOf(Math.round(m_Si * 100.0) / 100.0));
      } else {
        JOptionPane.showMessageDialog(frame, Imag_errmsg);
      }

      if (Real_errmsg.compareTo("") == 0) {
        double m_Si2 = jni.getRealSum();
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
          }

          secondEchoChooserWindow.setDialogTitle("CHOOSE PHASE FILE");
          statusOfFileChooser = secondEchoChooserWindow.showSaveDialog(null);

          if (statusOfFileChooser != JFileChooser.APPROVE_OPTION) {
            JOptionPane.showMessageDialog(frame, "Error: No phase file selected");
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
      if (!(WindowManager.getImage(s6MagWindowTitle) == null || WindowManager.getImage(s6PhaseWindowTitle) == null)) {
        updateVariables();

        // Calculating susceptibility in C++
        jni.calcSusceptibility();

        // Getting values from C++ and setting to GUI
        lbl_aVal.setText(String.valueOf(Math.round(jni.getA() * 100.0) / 100.0) + " pixels");
        lbl_dchiVal.setText(String.valueOf(Math.round(jni.getChi() * 100.0) / 100.0) + " ppm");
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

      // Begin opening spin echo image
      final JFileChooser spinEchoChooserWindow = new JFileChooser(
          FileSystemView.getFileSystemView().getHomeDirectory());
      spinEchoChooserWindow.setDialogTitle("CHOOSE MAG FILE");
      int statusOfFileChooser = spinEchoChooserWindow.showSaveDialog(null);

      Opener fileOpener = new Opener();

      if (statusOfFileChooser != JFileChooser.APPROVE_OPTION) {
        JOptionPane.showMessageDialog(frame, "Error: No file selected");
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

      if (WindowManager.getImage(s7WindowTitle) != null) {
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
        int V1SE_x1 = Integer.parseInt(txt_v1seXVal1.getText());
        int V1SE_x2 = Integer.parseInt(txt_v1seXVal2.getText());
        int V1SE_y1 = Integer.parseInt(txt_v1seYVal1.getText());
        int V1SE_y2 = Integer.parseInt(txt_v1seYVal2.getText());
        int V1SE_z1 = Integer.parseInt(txt_v1seZVal1.getText());
        int V1SE_z2 = Integer.parseInt(txt_v1seZVal2.getText());
        int V2SE_x1 = Integer.parseInt(txt_v2seXVal1.getText());
        int V2SE_x2 = Integer.parseInt(txt_v2seXVal2.getText());
        int V2SE_y1 = Integer.parseInt(txt_v2seYVal1.getText());
        int V2SE_y2 = Integer.parseInt(txt_v2seYVal2.getText());
        int V2SE_z1 = Integer.parseInt(txt_v2seZVal1.getText());
        int V2SE_z2 = Integer.parseInt(txt_v2seZVal2.getText());

        int V1SE_size = (V1SE_x2 - V1SE_x1) * (V1SE_y2 - V1SE_y1) * (V1SE_z2 - V1SE_z1);
        int V2SE_size = (V2SE_x2 - V2SE_x1) * (V2SE_y2 - V2SE_y1) * (V2SE_z2 - V2SE_z1);

        // Getting center from GUI
        int VSE_centerX = Integer.parseInt(txt_spinCenterXVal.getText());
        int VSE_centerY = Integer.parseInt(txt_spinCenterYVal.getText());
        int VSE_centerZ = Integer.parseInt(txt_spinCenterZVal.getText()) - 1;

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
          double a = Math.pow((V0 * 3.0) / (4.0 * Math.PI), 1.0 / 3.0);

          // Equation 17
          double rho_SE0 = S2_SE / (V2_SE - V0);
          // Setting to GUI
          lbl_rho0SEVal.setText(String.valueOf(Math.round(rho_SE0 * 100.0) / 100.0));
          logger.addVariable("rho_SE0", rho_SE0);

          double deltaV = 1;

          // Defining signal to noise ratio
          double snrStandardDeviation = Double.parseDouble(txt_sigSEVal.getText());
          double SNR_SE = rho_SE0 / snrStandardDeviation;

          // Equation 19
          double dV0 = (Math.sqrt(deltaV) / SNR_SE) * Math.sqrt(V2_SE + Math.pow(V2_SE - V0, 2) / (V1_SE - V2_SE));

          // Adding v0 and dv0 to GUI
          lbl_V0Val.setText(String.valueOf(Math.round(V0 * 100.0) / 100.0) + " " + PLUS_MINUS + " "
              + String.valueOf(Math.round(dV0 * 100.0) / 100.0));
          logger.addVariable("dV0", dV0);

          // Error for a - was derived with Norman
          double da = (a * dV0) / (3 * V0);
          // Adding to GUI
          lbl_aSE.setText(String.valueOf(Math.round(a * 100.0) / 100.0) + " " + PLUS_MINUS + " "
              + String.valueOf(Math.round(da * 100.0) / 100.0) + " pixels");
          logger.addVariable("a", a);
          logger.addVariable("da", da);

          double B0 = Double.parseDouble(txt_B0Val.getText());
          double TELast = Double.parseDouble(txt_TELastVal.getText()) / 1000.0;
          double mag_moment = Double.parseDouble(txt_magMomVal.getText());
          double d_mag_moment = Double.parseDouble(lbl_errVal.getText());

          // p = ga^3 can be rewritten to get dChi
          double dChi = (2.0 * mag_moment) / (GAMMARBAR * B0 * TELast * V0);
          // And its error
          double d_dChi = dChi * Math.sqrt(Math.pow(d_mag_moment / mag_moment, 2) + Math.pow(dV0 / V0, 2));

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
    if (!lbl_estBkgPhaseVal.getText().isEmpty()) {

      // Removing background phase from all data points in matrix
      for (int k = drawnRectangle_initialZ; k < drawnRectangle_initialZ + drawnRectangle_sizeZ; k++) {
        for (int i = drawnRectangle_initialX; i < drawnRectangle_initialX + drawnRectangle_sizeX; i++) {
          for (int j = drawnRectangle_initialY; j < drawnRectangle_initialY + drawnRectangle_sizeY; j++) {
            phaseVals[i][j][k] = Math.abs(phaseVals[i][j][k] - item.bkgPhase);
          }
        }
      }
    } else

    {
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

    if (!lbl_estBkgPhaseVal.getText().isEmpty()) {
      int subpix_size = (int) ((2 * m_R0 + 1) * grid);
      for (int i = 0; i < subpix_size; i++) {
        for (int j = 0; j < subpix_size; j++) {
          if (phaseVals[i][j] > 0) {
            phaseVals[i][j] = phaseVals[i][j] - (float) item.bkgPhase;
          } else if (phaseVals[i][j] < 0) {
            phaseVals[i][j] = phaseVals[i][j] + (float) item.bkgPhase;
          }
        }
      }
    } else

    {
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
      subpixelCoordinate = subCenter + (coordinate - (double) (item.centerS().get(0).intValue())) * 10.0;
    }
    if (axisFlag == 1) {
      subpixelCoordinate = subCenter + (coordinate - (double) (item.centerS().get(1).intValue())) * 10.0;
    }
    if (axisFlag == 2) {
      subpixelCoordinate = subCenter + (coordinate - (double) (item.centerS().get(2).intValue())) * 10.0;
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

      if (estimateCenterRadii_isClicked) {
        logger.addVariable("mxyz", item.centerS().get(0));
        logger.addVariable("mxyz", item.centerS().get(1));
        logger.addVariable("mxyz", item.centerS().get(2));

        jni.setXYZ(item.centerS().get(0), item.centerS().get(1), item.centerS().get(2));

        jni.setmR123(Double.parseDouble(txt_r1.getText()), Double.parseDouble(txt_r2.getText()),
            Double.parseDouble(txt_r3.getText()));

        double RCenter = Double.parseDouble(txt_rc.getText());
        double phaseValue = Double.parseDouble(txt_eqPhaseRC.getText());
        double estimatedPValue = phaseValue * Math.pow(RCenter, 3);
        double R1PhaseCalc = estimatedPValue / Math.pow(Double.parseDouble(txt_r1.getText()), 3);
        double R2PhaseCalc = estimatedPValue / Math.pow(Double.parseDouble(txt_r2.getText()), 3);
        double R3PhaseCalc = estimatedPValue / Math.pow(Double.parseDouble(txt_r3.getText()), 3);

        jni.setR123PhaseCalc(R1PhaseCalc, R2PhaseCalc, R3PhaseCalc);
        jni.setSmallBox(item.roi_mag_belowM_xi, item.roi_mag_belowM_yi, item.roi_mag_belowM_zi, item.roi_mag_belowM_Dx,
            item.roi_mag_belowM_Dy,
            item.roi_mag_belowM_Dz);
        jni.setCenterL(item.centerL().get(0), item.centerL().get(1), item.centerL().get(2));
        jni.setCenterM(item.centerM().get(0), item.centerM().get(1), item.centerM().get(2));

      }

      double snr = Double.parseDouble(txt_snrVal.getText());
      double e12 = Double.parseDouble(txt_eps12val.getText());
      double e23 = Double.parseDouble(txt_eps23val.getText());
      double B0 = Double.parseDouble(txt_B0Val.getText());
      double R_Chi = Double.parseDouble(txt_RChiVal.getText());
      double TElast = Double.parseDouble(txt_TELastVal.getText());
      jni.setMagMomentVariables(snr, e12, e23, B0, R_Chi, TElast);

      B0 = Double.parseDouble(txt_B0Val.getText());
      double RChi = Double.parseDouble(txt_RChiVal.getText());
      double TEFirst = Double.parseDouble(txt_TEFirstVal.getText()) / 1000.0;
      double TELast = Double.parseDouble(txt_TELastVal.getText()) / 1000.0;
      jni.setStep6Variables(TEFirst, TELast, B0, RChi);

      logger.addInfo("Updated variables");
    } catch (

    Exception exc) {
      JOptionPane.showMessageDialog(frame, "Error in updateVariables():\n" + exc.toString());
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
    frame = new JFrame("Calculate Magnetic Moment 3D");
    frame.setAlwaysOnTop(false);
    frame.setBounds(100, 100, 900, 700);
    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    frame.getContentPane().setLayout(
        new MigLayout("",
            "[][][grow][][45px:45.00px:45px][45px:45px:45px,grow][45px:45px:45px,grow][48.00,grow][]",
            "[][][][][][][][][][][][][][][][][][][][][][][][]"));

    lbl_stepone = new JLabel("1.");
    frame.getContentPane().add(lbl_stepone, "cell 0 0");

    btn_loadImages = new JButton("Load Magnitude and Phase Images");
    frame.getContentPane().add(btn_loadImages, "cell 1 0,growx");
    // btn_loadImages.addActionListener(new Calculate_Magnetic_Moment_3D());

    lbl_steptwo = new JLabel("2.");
    frame.getContentPane().add(lbl_steptwo, "cell 0 1");

    btn_estCR = new JButton("Estimate Center/Radii");
    frame.getContentPane().add(btn_estCR, "cell 1 1,growx");

    lbl_M = new JLabel("|M%|:");
    frame.getContentPane().add(lbl_M, "flowx,cell 4 1,alignx right");

    txt_M = new DefTextField("50");
    // txt_M.setText("50");
    frame.getContentPane().add(txt_M, "cell 5 1,alignx left");
    txt_M.setColumns(2);

    lbl_eqPhase = new JLabel("Equatorial Phase at RCenter=");
    frame.getContentPane().add(lbl_eqPhase, "flowx,cell 1 2");

    lbl_eqPhaseUnit = new JLabel("radian(s)");
    frame.getContentPane().add(lbl_eqPhaseUnit, "cell 2 2");

    lbl_rc = new JLabel("RCenter=");
    frame.getContentPane().add(lbl_rc, "cell 4 2,alignx trailing");

    txt_rc = new DefTextField();
    txt_rc.setColumns(3);
    frame.getContentPane().add(txt_rc, "cell 5 2,growx");

    lbl_rcUnit = new JLabel("pixels");
    frame.getContentPane().add(lbl_rcUnit, "cell 6 2");

    lbl_stepthree = new JLabel("3.");
    frame.getContentPane().add(lbl_stepthree, "cell 0 3");

    btn_genSubpix = new JButton("Generate Subpixel Grid/Data");
    frame.getContentPane().add(btn_genSubpix, "cell 1 3,growx");

    btn_removeBkg = new JButton("Remove Bkg");
    frame.getContentPane().add(btn_removeBkg, "cell 2 3");

    chkbx_showrc = new JCheckBox("Show RCenter");
    chkbx_showrc.setVerticalAlignment(SwingConstants.TOP);
    frame.getContentPane().add(chkbx_showrc, "cell 4 3");

    lbl_gridSize = new JLabel("<html>Grid Size: 10<sup>3</sup></html>");
    frame.getContentPane().add(lbl_gridSize, "flowx,cell 8 3");

    lbl_stepfour = new JLabel("4.");
    frame.getContentPane().add(lbl_stepfour, "cell 0 4");

    btn_estSubC = new JButton("Estimate Subpixel Center");
    frame.getContentPane().add(btn_estSubC, "cell 1 4,growx");

    txt_eqPhaseRC = new DefTextField("1.0");
    frame.getContentPane().add(txt_eqPhaseRC, "cell 1 2");
    txt_eqPhaseRC.setColumns(5);

    lbl_spx = new JLabel("x=");
    frame.getContentPane().add(lbl_spx, "flowx,cell 2 4");

    btn_redraw = new JButton("Redraw Center");
    frame.getContentPane().add(btn_redraw, "flowx,cell 1 5");

    btn_verifyRadii = new JButton("Verify Radii");
    frame.getContentPane().add(btn_verifyRadii, "cell 1 5");

    lbl_calculated = new JLabel("Calculated");
    frame.getContentPane().add(lbl_calculated, "cell 4 5");

    lbl_actual = new JLabel("Actual");
    frame.getContentPane().add(lbl_actual, "cell 6 5");

    lbl_r1 = new JLabel("R1=");
    frame.getContentPane().add(lbl_r1, "flowx,cell 1 6");

    txt_r1 = new DefTextField();
    frame.getContentPane().add(txt_r1, "cell 1 6");
    txt_r1.setColumns(5);

    lbl_r1unit = new JLabel("pixels");
    frame.getContentPane().add(lbl_r1unit, "cell 1 6");

    lbl_r1phase = new JLabel("R1 Corresponding Phase:");
    frame.getContentPane().add(lbl_r1phase, "cell 2 6");

    lbl_r1phaseCalc = new JLabel("");
    frame.getContentPane().add(lbl_r1phaseCalc, "cell 4 6");

    lbl_r1phaseUnit = new JLabel("radians");
    frame.getContentPane().add(lbl_r1phaseUnit, "cell 5 6");

    lbl_r1phaseAct = new JLabel("");
    frame.getContentPane().add(lbl_r1phaseAct, "cell 6 6");

    lbl_r1AphaseUnit = new JLabel("radians");
    frame.getContentPane().add(lbl_r1AphaseUnit, "cell 7 6");

    btn_plotX = new JButton("Plot X Phase Profiles");
    frame.getContentPane().add(btn_plotX, "cell 8 6");

    lbl_r2 = new JLabel("R2=");
    frame.getContentPane().add(lbl_r2, "flowx,cell 1 7");

    lbl_r2phase = new JLabel("R2 Corresponding Phase:");
    frame.getContentPane().add(lbl_r2phase, "cell 2 7");

    lbl_r2phaseCalc = new JLabel("");
    frame.getContentPane().add(lbl_r2phaseCalc, "cell 4 7");

    lbl_r2phaseUnit = new JLabel("radians");
    frame.getContentPane().add(lbl_r2phaseUnit, "cell 5 7");

    lbl_r2phaseAct = new JLabel("");
    frame.getContentPane().add(lbl_r2phaseAct, "cell 6 7");

    lbl_r2AphaseUnit = new JLabel("radians");
    frame.getContentPane().add(lbl_r2AphaseUnit, "cell 7 7");

    btn_plotY = new JButton("Plot Y Phase Profiles");
    frame.getContentPane().add(btn_plotY, "cell 8 7");

    lbl_r3 = new JLabel("R3=");
    frame.getContentPane().add(lbl_r3, "flowx,cell 1 8");

    txt_r2 = new DefTextField();
    txt_r2.setColumns(5);
    frame.getContentPane().add(txt_r2, "cell 1 7");

    txt_r3 = new DefTextField();
    txt_r3.setColumns(5);
    frame.getContentPane().add(txt_r3, "cell 1 8");

    lbl_r2unit = new JLabel("pixels");
    frame.getContentPane().add(lbl_r2unit, "cell 1 7");

    lbl_r3unit = new JLabel("pixels");
    frame.getContentPane().add(lbl_r3unit, "cell 1 8");

    lbl_r3phase = new JLabel("R3 Corresponding Phase:");
    frame.getContentPane().add(lbl_r3phase, "cell 2 8");

    lbl_r3phaseCalc = new JLabel("");
    frame.getContentPane().add(lbl_r3phaseCalc, "cell 4 8");

    lbl_r3phaseUnit = new JLabel("radians");
    frame.getContentPane().add(lbl_r3phaseUnit, "cell 5 8");

    lbl_r3phaseAct = new JLabel("");
    frame.getContentPane().add(lbl_r3phaseAct, "cell 6 8");

    lbl_r3AphaseUnit = new JLabel("radians");
    frame.getContentPane().add(lbl_r3AphaseUnit, "cell 7 8");

    btn_plotZ = new JButton("Plot Z Phase Profiles");
    frame.getContentPane().add(btn_plotZ, "cell 8 8");

    btn_estBkgDens = new JButton("<html>Estimate Bkg & &rho;0");
    frame.getContentPane().add(btn_estBkgDens, "cell 1 9,growx");

    lbl_rho0 = new JLabel("<html>&rho;0 =</html>");
    frame.getContentPane().add(lbl_rho0, "flowx,cell 2 9");

    lbl_rho0val = new JLabel("");
    frame.getContentPane().add(lbl_rho0val, "cell 2 9");

    lbl_estBkgPhase = new JLabel("Estimated Background Phase =");
    frame.getContentPane().add(lbl_estBkgPhase, "cell 4 9");

    lbl_estBkgPhaseVal = new JLabel("");
    frame.getContentPane().add(lbl_estBkgPhaseVal, "cell 7 9");

    lbl_estBkgPhaseUnit = new JLabel("radians");
    frame.getContentPane().add(lbl_estBkgPhaseUnit, "cell 8 9");

    lbl_stepfive = new JLabel("5.");
    frame.getContentPane().add(lbl_stepfive, "cell 0 10");

    btn_calcMagMom = new JButton("Calculate Magnetic Moment");
    frame.getContentPane().add(btn_calcMagMom, "cell 1 10,growx");

    lbl_magMom = new JLabel("|p|=");
    frame.getContentPane().add(lbl_magMom, "flowx,cell 2 10");

    btn_loadSimImg = new JButton("Load Simulated Images");
    frame.getContentPane().add(btn_loadSimImg, "cell 1 11,growx");

    txt_magMomVal = new DefTextField();
    frame.getContentPane().add(txt_magMomVal, "cell 2 10");
    txt_magMomVal.setColumns(4);

    lbl_magMomUnit = new JLabel("<html>radians*pixel<sup>3</sup></html>");
    frame.getContentPane().add(lbl_magMomUnit, "cell 2 10");

    lbl_snr = new JLabel("SNR=");
    frame.getContentPane().add(lbl_snr, "flowx,cell 1 12");

    txt_snrVal = new DefTextField("1.0");
    frame.getContentPane().add(txt_snrVal, "cell 1 12");
    txt_snrVal.setColumns(3);

    lbl_eps12 = new JLabel("<html>&epsilon;12 =</html>");
    frame.getContentPane().add(lbl_eps12, "flowx,cell 1 13");

    txt_eps12val = new DefTextField();
    txt_eps12val.setColumns(3);
    frame.getContentPane().add(txt_eps12val, "cell 1 13");

    lbl_ReRi = new JLabel("Real(S" + ITALICIZED_I + ") =");
    frame.getContentPane().add(lbl_ReRi, "flowx,cell 2 13");

    lbl_eps23 = new JLabel("<html>&epsilon;23 =</html>");
    frame.getContentPane().add(lbl_eps23, "flowx,cell 1 14");

    txt_eps23val = new DefTextField();
    txt_eps23val.setColumns(3);
    frame.getContentPane().add(txt_eps23val, "cell 1 14");

    lbl_err = new JLabel("<html>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&delta&rho/&rho =</html>");
    lbl_err.setHorizontalAlignment(SwingConstants.RIGHT);
    frame.getContentPane().add(lbl_err, "cell 1 13");

    lbl_errVal = new JLabel("");
    frame.getContentPane().add(lbl_errVal, "cell 1 13");

    lbl_Ri = new JLabel("R" + ITALICIZED_I + " =");
    frame.getContentPane().add(lbl_Ri, "flowx,cell 2 12");

    txt_Ri = new DefTextField();
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

    JLabel lbl_B0 = new JLabel("B0  =");
    frame.getContentPane().add(lbl_B0, "flowx,cell 1 15");

    txt_B0Val = new DefTextField();
    frame.getContentPane().add(txt_B0Val, "cell 1 15");
    txt_B0Val.setColumns(3);

    JLabel lbl_B0Unit = new JLabel("T");
    frame.getContentPane().add(lbl_B0Unit, "cell 1 15");

    JLabel lbl_TElast = new JLabel("    TE_last =");
    frame.getContentPane().add(lbl_TElast, "cell 1 15,alignx right");

    lbl_stepsix = new JLabel("6.");
    frame.getContentPane().add(lbl_stepsix, "cell 0 16");

    btn_loadTE = new JButton("Load First TE Images");
    frame.getContentPane().add(btn_loadTE, "cell 1 16,growx");

    lbl_TEFirst = new JLabel("TE_first =");
    frame.getContentPane().add(lbl_TEFirst, "flowx,cell 2 16");

    lbl_dchi = new JLabel("<html>&Delta;&Chi =</html>");
    frame.getContentPane().add(lbl_dchi, "flowx,cell 5 16");

    btn_unk = new JButton("TODO: Name");
    frame.getContentPane().add(btn_unk, "cell 1 17,growx");

    txt_TEFirstVal = new DefTextField();
    frame.getContentPane().add(txt_TEFirstVal, "cell 2 16");
    txt_TEFirstVal.setColumns(4);

    lbl_TEFirstUnit = new JLabel("ms");
    frame.getContentPane().add(lbl_TEFirstUnit, "cell 2 16");

    lbl_RChi = new JLabel("<html>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;R<sub>&Delta&Chi</sub> =</html>");
    frame.getContentPane().add(lbl_RChi, "flowx,cell 2 17");

    txt_RChiVal = new DefTextField();
    frame.getContentPane().add(txt_RChiVal, "cell 2 17");
    txt_RChiVal.setColumns(4);

    lbl_RChiUnit = new JLabel("pixels");
    frame.getContentPane().add(lbl_RChiUnit, "cell 2 17");

    lbl_dchiVal = new JLabel("");
    frame.getContentPane().add(lbl_dchiVal, "cell 5 16");

    lbl_a = new JLabel("  a =");
    frame.getContentPane().add(lbl_a, "flowx,cell 5 17");

    lbl_aVal = new JLabel("");
    frame.getContentPane().add(lbl_aVal, "cell 5 17");

    lbl_stepseven = new JLabel("7.");
    frame.getContentPane().add(lbl_stepseven, "cell 0 18");

    btn_loadspinecho = new JButton("Load Spin Echo Images");
    frame.getContentPane().add(btn_loadspinecho, "cell 1 18,growx");

    lbl_sigSE = new JLabel("<html>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&sigma<sub>SE</sub> =</html>");
    frame.getContentPane().add(lbl_sigSE, "flowx,cell 2 18,alignx left");

    txt_sigSEVal = new DefTextField();
    frame.getContentPane().add(txt_sigSEVal, "cell 2 18");
    txt_sigSEVal.setColumns(4);

    lbl_spinCenter = new JLabel("Second Image Center:  (");
    frame.getContentPane().add(lbl_spinCenter, "flowx,cell 1 19,alignx trailing");

    txt_spinCenterXVal = new DefTextField();
    txt_spinCenterXVal.setColumns(4);
    frame.getContentPane().add(txt_spinCenterXVal, "flowx,cell 2 19,alignx left");

    // lbl_v1se = new JLabel("V1,SE Box Coordinates: (");
    lbl_v1se = new JLabel("<html>V<sub>1,SE</sub> Box Coordinates:  (</html>");
    frame.getContentPane().add(lbl_v1se, "flowx,cell 1 20,alignx trailing");

    txt_v1seXVal1 = new DefTextField();
    txt_v1seXVal1.setColumns(4);
    frame.getContentPane().add(txt_v1seXVal1, "flowx,cell 2 20,alignx left");

    JLabel lbl_brack23 = new JLabel("(");
    frame.getContentPane().add(lbl_brack23, "cell 3 20,alignx trailing");

    txt_v1seXVal2 = new DefTextField();
    frame.getContentPane().add(txt_v1seXVal2, "cell 4 20,alignx center");
    txt_v1seXVal2.setColumns(4);

    JLabel lbl_brack24 = new JLabel("-1 )");
    frame.getContentPane().add(lbl_brack24, "cell 7 20");

    lbl_v2se = new JLabel("<html>V<sub>2,SE</sub> Box Coordinates:  (</html>");
    frame.getContentPane().add(lbl_v2se, "flowx,cell 1 21,alignx trailing");

    txt_v2seXVal1 = new DefTextField();
    txt_v2seXVal1.setColumns(4);
    frame.getContentPane().add(txt_v2seXVal1, "flowx,cell 2 21,alignx left");

    lbl_comma11 = new JLabel(",");
    frame.getContentPane().add(lbl_comma11, "cell 2 19");

    lbl_comma21 = new JLabel(",");
    frame.getContentPane().add(lbl_comma21, "cell 2 20");

    lbl_comma31 = new JLabel(",");
    frame.getContentPane().add(lbl_comma31, "cell 2 21");

    txt_spinCenterYVal = new DefTextField();
    txt_spinCenterYVal.setColumns(4);
    frame.getContentPane().add(txt_spinCenterYVal, "cell 2 19");

    txt_v1seYVal1 = new DefTextField();
    txt_v1seYVal1.setColumns(4);
    frame.getContentPane().add(txt_v1seYVal1, "cell 2 20,alignx center");

    txt_v2seYVal1 = new DefTextField();
    txt_v2seYVal1.setColumns(4);
    frame.getContentPane().add(txt_v2seYVal1, "cell 2 21");

    lbl_comma12 = new JLabel(",");
    frame.getContentPane().add(lbl_comma12, "cell 2 19");

    lbl_comma22 = new JLabel(",");
    frame.getContentPane().add(lbl_comma22, "cell 2 20");

    lbl_comma32 = new JLabel(",");
    frame.getContentPane().add(lbl_comma32, "cell 2 21");

    txt_spinCenterZVal = new DefTextField();
    frame.getContentPane().add(txt_spinCenterZVal, "cell 2 19");
    txt_spinCenterZVal.setColumns(4);

    txt_v1seZVal1 = new DefTextField();
    txt_v1seZVal1.setColumns(4);
    frame.getContentPane().add(txt_v1seZVal1, "cell 2 20");

    txt_v2seZVal1 = new DefTextField();
    txt_v2seZVal1.setColumns(4);
    frame.getContentPane().add(txt_v2seZVal1, "cell 2 21");

    JLabel lbl_brack33 = new JLabel("(");
    frame.getContentPane().add(lbl_brack33, "cell 3 21");

    txt_v2seXVal2 = new DefTextField();
    txt_v2seXVal2.setColumns(4);
    frame.getContentPane().add(txt_v2seXVal2, "cell 4 21,alignx center");

    JLabel lbl_brack34 = new JLabel("-1 )");
    frame.getContentPane().add(lbl_brack34, "cell 7 21");

    btn_estRadSpinEcho = new JButton("Estimate Object Radius From Spin Echo");
    frame.getContentPane().add(btn_estRadSpinEcho, "cell 1 22");

    lbl_V0 = new JLabel("<html>V<sub>0</sub> =</html>");
    frame.getContentPane().add(lbl_V0, "flowx,cell 2 22");

    lbl_rho0SE = new JLabel("<html>&rho<sub>0,SE</sub> =</html>");
    frame.getContentPane().add(lbl_rho0SE, "flowx,cell 4 22");

    lbl_echoDChi = new JLabel("<html>&Delta&Chi =</html>");
    frame.getContentPane().add(lbl_echoDChi, "flowx,cell 2 23");

    lbl_V0Val = new JLabel("     ");
    frame.getContentPane().add(lbl_V0Val, "cell 2 22");

    lbl_V0Unit = new JLabel("<html>pixels<sup>3</sup></html>");
    frame.getContentPane().add(lbl_V0Unit, "cell 2 22,aligny top");

    lbl_echoDChiVal = new JLabel("");
    frame.getContentPane().add(lbl_echoDChiVal, "cell 2 23");

    lbl_rho0SEVal = new JLabel("");
    frame.getContentPane().add(lbl_rho0SEVal, "cell 4 22");

    lbl_aSE = new JLabel("a =");
    frame.getContentPane().add(lbl_aSE, "flowx,cell 4 23");

    lbl_aSEVal = new JLabel("");
    frame.getContentPane().add(lbl_aSEVal, "cell 4 23");

    lbl_comma23 = new JLabel(",");
    frame.getContentPane().add(lbl_comma23, "flowx,cell 5 20");

    lbl_comma33 = new JLabel(",");
    frame.getContentPane().add(lbl_comma33, "flowx,cell 5 21");

    lbl_comma24 = new JLabel(",");
    frame.getContentPane().add(lbl_comma24, "flowx,cell 6 20");

    lbl_comma34 = new JLabel(",");
    frame.getContentPane().add(lbl_comma34, "flowx,cell 6 21");
    lbl_spx.setLabelFor(txt_spx);

    txt_spx = new DefTextField();
    frame.getContentPane().add(txt_spx, "cell 2 4");
    txt_spx.setColumns(4);

    lbl_spy = new JLabel("y=");
    frame.getContentPane().add(lbl_spy, "cell 2 4");

    txt_spy = new DefTextField();
    txt_spy.setColumns(4);
    frame.getContentPane().add(txt_spy, "cell 2 4");
    lbl_spy.setLabelFor(txt_spy);

    lbl_spz = new JLabel("z=");
    frame.getContentPane().add(lbl_spz, "cell 2 4");

    txt_spz = new DefTextField();
    txt_spz.setColumns(4);
    frame.getContentPane().add(txt_spz, "cell 2 4");
    lbl_spz.setLabelFor(txt_spz);

    lbl_spzCorrection = new JLabel("-1");
    frame.getContentPane().add(lbl_spzCorrection, "cell 2 4");
    lbl_spzCorrection.setLabelFor(txt_spz);

    lbl_innerBrack1 = new JLabel("-1 )");
    frame.getContentPane().add(lbl_innerBrack1, "cell 2 19");

    JLabel lbl_brack22 = new JLabel("-1 )");
    frame.getContentPane().add(lbl_brack22, "cell 2 20,alignx left");

    JLabel lbl_brack32 = new JLabel("-1 )");
    frame.getContentPane().add(lbl_brack32, "cell 2 21,alignx left");

    txt_v1seYVal2 = new DefTextField();
    frame.getContentPane().add(txt_v1seYVal2, "cell 5 20,alignx center");
    txt_v1seYVal2.setColumns(4);

    txt_v2seYVal2 = new DefTextField();
    txt_v2seYVal2.setColumns(4);
    frame.getContentPane().add(txt_v2seYVal2, "cell 5 21,alignx center");

    txt_v1seZVal2 = new DefTextField();
    txt_v1seZVal2.setColumns(4);
    frame.getContentPane().add(txt_v1seZVal2, "cell 6 20,alignx center");

    txt_v2seZVal2 = new DefTextField();
    txt_v2seZVal2.setColumns(4);
    frame.getContentPane().add(txt_v2seZVal2, "cell 6 21,alignx center");

    lbl_rcx = new JLabel("x=");
    frame.getContentPane().add(lbl_rcx, "flowx,cell 2 1");

    txt_rcx = new DefTextField();
    frame.getContentPane().add(txt_rcx, "cell 2 1");
    txt_rcx.setColumns(4);

    lbl_rcy = new JLabel("y=");
    frame.getContentPane().add(lbl_rcy, "cell 2 1,alignx left");

    txt_rcy = new DefTextField();
    frame.getContentPane().add(txt_rcy, "cell 2 1");
    txt_rcy.setColumns(4);

    lbl_rcz = new JLabel("z=");
    frame.getContentPane().add(lbl_rcz, "cell 2 1,alignx left");

    txt_rcz = new DefTextField();
    frame.getContentPane().add(txt_rcz, "cell 2 1");
    txt_rcz.setColumns(4);

    lbl_rczCorrection = new JLabel("-1");
    frame.getContentPane().add(lbl_rczCorrection, "cell 2 1,alignx left");

    txt_TELastVal = new DefTextField();
    frame.getContentPane().add(txt_TELastVal, "cell 1 15,alignx right");
    txt_TELastVal.setColumns(3);

    JLabel lbl_TELastUnit = new JLabel("ms");
    lbl_TELastUnit.setHorizontalAlignment(SwingConstants.RIGHT);
    frame.getContentPane().add(lbl_TELastUnit, "cell 1 15,alignx right");
    frame.getContentPane().add(lbl_TELastUnit, "cell 1 15,alignx right");

    // ActionListeners
    btn_loadImages.addActionListener(new Calculate_Magnetic_Moment_3D());
    btn_estCR.addActionListener(new Calculate_Magnetic_Moment_3D());
    btn_genSubpix.addActionListener(new Calculate_Magnetic_Moment_3D());
    btn_removeBkg.addActionListener(new Calculate_Magnetic_Moment_3D());
    btn_estSubC.addActionListener(new Calculate_Magnetic_Moment_3D());
    btn_redraw.addActionListener(new Calculate_Magnetic_Moment_3D());
    btn_verifyRadii.addActionListener(new Calculate_Magnetic_Moment_3D());
    btn_plotX.addActionListener(new Calculate_Magnetic_Moment_3D());
    btn_plotY.addActionListener(new Calculate_Magnetic_Moment_3D());
    btn_plotZ.addActionListener(new Calculate_Magnetic_Moment_3D());
    btn_estBkgDens.addActionListener(new Calculate_Magnetic_Moment_3D());
    btn_calcMagMom.addActionListener(new Calculate_Magnetic_Moment_3D());
    btn_loadSimImg.addActionListener(new Calculate_Magnetic_Moment_3D());
    btn_sumRi.addActionListener(new Calculate_Magnetic_Moment_3D());
    btn_loadTE.addActionListener(new Calculate_Magnetic_Moment_3D());
    btn_unk.addActionListener(new Calculate_Magnetic_Moment_3D());
    btn_loadspinecho.addActionListener(new Calculate_Magnetic_Moment_3D());
    btn_estRadSpinEcho.addActionListener(new Calculate_Magnetic_Moment_3D());
  }
}