/*
 * @author Nicholas Bolton (bolton21@uwindsor.ca)
 * @version 3.12, May 5, 2023
 */

import ij.gui.OvalRoi;
import ij.gui.PointRoi;
import ij.gui.Roi;
import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.lang.Math;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.filechooser.FileSystemView;
import java.io.File;
import java.io.IOException;
import java.util.Scanner;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.*;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import ij.process.*;
import ij.WindowManager;
import ij.io.Opener;

public class Calculate_Magnetic_Moment_3D implements PlugIn {

  private static JNIMethods jni = new JNIMethods();
  public static GUI gui;

  private static ImageItem item;
  public static LogManager logger;
  private static ROIS roiImgMag, roiImgMagXZ, roiImgPhase, roiImgPhaseXZ, roiImgV1SE, roiImgV1SEXZ;
  public static ImagePlus subpixelMagImage, subpixelMagImageXZ, subpixelPhaseImage, subpixelPhaseImageXZ, V1SE_XYImage,
      V1SE_XZImage;
  private static String subCenterErrorMessage,
      V1XY_Title, V1XZ_Title, s1MagWindowTitle,
      s1PhaseWindowTitle, s5MagWindowTitle, s5PhaseWindowTitle, s6MagWindowTitle, s6PhaseWindowTitle, s7WindowTitle;

  private static float[][] subpixelPhaseMatrix, subpixelPhaseMatrixXZ;
  private static float[][][] croppedRealNumbers3D, croppedImaginaryNumbers3D;

  private static final int grid = 10;
  private static final double GAMMARBAR = 42.58;
  private static final String ACCEPTED_FILE_TYPE = "nii";
  private static final String PLUS_MINUS = "\u00B1";
  private static final String subMagTitle = "Subpixel Mag Image";
  private static final String subMagXZTitle = "Subpixel Mag Image XZ";
  private static final String subPhaseTitle = "Subpixel Phase Image";
  private static final String subPhaseXZTitle = "Subpixel Phase Image XZ";

  // set this to true if debugging in IntelliJ
  private static final boolean DEBUG = false;

  public static void main(String[] args) {

    if (DEBUG) {
      // set the plugins.dir property to make the plugin appear in the Plugins menu
      Class<?> clazz = Calculate_Magnetic_Moment_3D.class;
      String url = clazz.getResource("/" + clazz.getName().replace('.', '/') +
          ".class").toString();
      String pluginsDir = url.substring("file:".length(), url.length() -
          clazz.getName().length() - ".class".length());
      System.setProperty("plugins.dir", pluginsDir);

      // start ImageJ
      new ImageJ();

      // run the plugin
      IJ.runPlugIn(clazz.getName(), "");
      // IJ.runPlugIn(clazz.getName(),"batch"); // Running the plugin to test batch
      // processing.
    }

  }

  // Function to run CMM3D (plug-in)
  public void run(String arg) {
    gui = new GUI();
  }

  // =====================================================================================
  // "Load Magnitude and Phase Images"
  // =====================================================================================
  public static void load_mag_phase_images() {

    s1MagWindowTitle = loadImagesIJ(s1MagWindowTitle);
    s1PhaseWindowTitle = loadImagesIJ(s1PhaseWindowTitle);

    /*
     * try {
     * s1MagWindowTitle = loadImages(s1MagWindowTitle, "mag");
     * s1PhaseWindowTitle = loadImages(s1PhaseWindowTitle, "phase");
     * } catch (IOException exc) {
     * return;
     * }
     */
  }

  // =====================================================================================
  // "Estimate Center/Radii"
  // =====================================================================================
  public static void est_center_rad() {

    try {

      // if a center was already calculated and a new ROI is drawn, clear previously
      // calculated center
      if (item != null) {
        ImageItem new_item = new ImageItem(s1MagWindowTitle,
            s1PhaseWindowTitle,
            (int) Double.parseDouble(gui.ltf_M.getValue()),
            Double.parseDouble(gui.ltf_eqPhase.getValue()));

        // if the ROI was changed
        if (new_item.roi_xi != item.roi_xi ||
            new_item.roi_yi != item.roi_yi ||
            new_item.roi_zi != item.roi_zi ||
            new_item.roi_dx != item.roi_dx ||
            new_item.roi_dy != item.roi_dy ||
            new_item.roi_dz != item.roi_dz) {
          gui.ltf_rcx.setValue("0.0");
          gui.ltf_rcy.setValue("0.0");
          gui.ltf_rcz.setValue("0.0");

          new_item = null; // don't need this anymore
        }
      }

      // new ImageItem object, passing in the Window titles, along
      // with M% and equitorial phase value
      item = new ImageItem(s1MagWindowTitle,
          s1PhaseWindowTitle,
          (int) Double.parseDouble(gui.ltf_M.getValue()),
          Double.parseDouble(gui.ltf_eqPhase.getValue()));

      // calculating each center
      item.calcCenterL();
      jni.setCenterL(item.centerL().get(0), item.centerL().get(1), item.centerL().get(2));
      item.calcCenterM();
      jni.setCenterM(item.centerM().get(0), item.centerM().get(1), item.centerM().get(2));
      item.calcCenterS();

      // setting GUI bkg phase to estimate bkg phase
      gui.ll_estBkgPhase.setValue(String.valueOf(Math.round(item.estBkg() * 100.0) / 100.0));

      // grabbing center_s from item
      double center_sx = item.centerS().get(0);
      double center_sy = item.centerS().get(1);
      double center_sz = item.centerS().get(2);
      logger.addInfo("cs", item.centerS());

      // Setting center to GUI unless user already has an inputted center point
      if (gui.ltf_rcx.getValueTF().isDefault()) {
        gui.ltf_rcx.setValue(String.valueOf(center_sx));
      }

      if (gui.ltf_rcy.getValueTF().isDefault()) {
        gui.ltf_rcy.setValue(String.valueOf(center_sy));
      }

      if (gui.ltf_rcz.getValueTF().isDefault()) {
        gui.ltf_rcz.setValue(String.valueOf(center_sz + 1));
      }

      // Getting x y and z centers from GUI
      item.setCenterSX(Double.parseDouble(gui.ltf_rcx.getValue()));
      item.setCenterSY(Double.parseDouble(gui.ltf_rcy.getValue()));
      item.setCenterSZ(Double.parseDouble(gui.ltf_rcz.getValue()) - 1);

      // Getting RCenter from estimated xyz
      // RCenter = estimateRCenter((int)(item.centerS().get(0)),
      // (int)(item.centerS().get(1)),
      // (int)(item.centerS().get(2)));

      // removing background phase then estimating rcenter
      item.remove_bkg();
      double RCenter = item.estimateRCenter();

      // Updating GUI
      gui.ltf_rc.setValue(String.valueOf(Math.round(RCenter * 10.0) / 10.0));

      item.calcR0123();
      double m_R0 = item.m_R0();

      jni.setmVariables(grid, m_R0, RCenter,
          Double.parseDouble(gui.ltf_rcx.getValue()),
          Double.parseDouble(gui.ltf_rcy.getValue()),
          Double.parseDouble(gui.ltf_rcz.getValue()) - 1.0,
          Double.parseDouble(gui.ltf_eqPhase.getValue()));
      jni.setMagMoment(Double.parseDouble(gui.ltf_eqPhase.getValue()) *
          Math.pow(RCenter, 3));

      // setting slice to center
      item.setSlice((int) Double.parseDouble(gui.ltf_rcz.getValue()));

    } catch (Exception exc) {
      JOptionPane.showMessageDialog(gui.frame, exc);
    }

  }

  // =====================================================================================
  // "Generate Subpixel Grid/Data"
  // =====================================================================================
  public static void gen_subpix() {
    // updateVariables();
    float[][] subpixelMagMatrix, subpixelMagMatrixXZ;

    // pass step 2 estimated center to item. This is needed for subpixel/pixel
    // coordinate conversions
    double c2x = Double.parseDouble(gui.ltf_rcx.getValue());
    double c2y = Double.parseDouble(gui.ltf_rcy.getValue());
    double c2z = Double.parseDouble(gui.ltf_rcz.getValue()) - 1.0;
    item.subpix_image_center = new Triplet<Double>(c2x, c2y, c2z);
    logger.addVariable("subpix image center", item.subpix_image_center.toString());

    // if no RCenter is inputted
    if (gui.ltf_rc.getValue().isEmpty()) {
      JOptionPane.showMessageDialog(gui.frame, "Error: No RCenter found.");
    }
    // if no XYZ center is inputted
    else if (gui.ltf_rcx.getValue().isEmpty() ||
        gui.ltf_rcy.getValue().isEmpty() ||
        gui.ltf_rcz.getValue().isEmpty()) {
      JOptionPane.showMessageDialog(gui.frame, "Error: No XYZ Center found.");
    }
    // If object is too close to edge of image
    else if (!item.isNearEdge) {
      JOptionPane.showMessageDialog(gui.frame, "Error: Object too close to edge of image.");
    } else {

      logger.addVariable("s1magwt", s1MagWindowTitle);
      logger.addVariable("s2phasewt", s1PhaseWindowTitle);
      // Getting mag and phase images
      ImagePlus magnitudeImage = WindowManager.getImage(s1MagWindowTitle);
      ImagePlus phaseImage = WindowManager.getImage(s1PhaseWindowTitle);
      logger.addInfo("Got images");

      double m_R0 = item.m_R0();

      // Size of new images
      int size_subpixelimg = (int) ((2 * m_R0 + 1) * grid);
      logger.addVariable("sizeofspi", size_subpixelimg);

      // Initializing 3D arrays to give to C++
      float[][][] croppedMagnitudeValues3D = new float[size_subpixelimg][size_subpixelimg][size_subpixelimg];
      float[][][] croppedPhaseValues3D = new float[size_subpixelimg][size_subpixelimg][size_subpixelimg];
      croppedRealNumbers3D = new float[size_subpixelimg][size_subpixelimg][size_subpixelimg];
      croppedImaginaryNumbers3D = new float[size_subpixelimg][size_subpixelimg][size_subpixelimg];

      // Initializing arrays that will be displayed in all 4 images
      subpixelMagMatrix = new float[size_subpixelimg][size_subpixelimg];
      subpixelMagMatrixXZ = new float[size_subpixelimg][size_subpixelimg];
      subpixelPhaseMatrix = new float[size_subpixelimg][size_subpixelimg];
      subpixelPhaseMatrixXZ = new float[size_subpixelimg][size_subpixelimg];
      logger.addInfo("made arrays");

      // Initial points of subpixel images
      int xi = (int) (c2x) - (int) m_R0;
      int yi = (int) (c2y) - (int) m_R0;
      int zi = (int) (c2z) - (int) m_R0;
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
      double RCenter = Double.parseDouble(gui.ltf_rc.getValue());
      jni.passGenSubpixelValues(m_R0, grid, RCenter, item.bkgPhase);
      jni.setRealImagNumbers(croppedRealNumbers3D, croppedImaginaryNumbers3D);

      // Generate subpixel in C++
      jni.generateSubpixelArray();
      logger.addVariable("bkgphase", item.bkgPhase);

      int croppedImageSize = 2 * (int) m_R0 + 1;
      int Nfinal = grid;

      // Interpolating XY subpixel phase matrix
      for (int m = 0; m < croppedImageSize * Nfinal; m++) {
        for (int k = 0; k < croppedImageSize * Nfinal; k++) {
          subpixelPhaseMatrix[k][m] = croppedPhaseValues3D[(int) (k / Nfinal)][(int) (m / Nfinal)][(int) m_R0];
        }
      }

      // Interpolating XY subpixel mag matrix
      for (int m = 0; m < croppedImageSize * Nfinal; m++) {
        for (int k = 0; k < croppedImageSize * Nfinal; k++) {
          subpixelMagMatrix[k][m] = croppedMagnitudeValues3D[(int) (k / Nfinal)][(int) (m / Nfinal)][(int) m_R0];
        }
      }

      // Interpolating XZ subpixel mag matrix
      for (int m = 0; m < croppedImageSize * Nfinal; m++) {
        for (int k = 0; k < croppedImageSize * Nfinal; k++) {
          subpixelMagMatrixXZ[k][m] = croppedMagnitudeValues3D[(int) (k / Nfinal)][(int) m_R0][(int) (m / Nfinal)];
        }
      }

      // Interpolating XZ subpixel phase matrix
      for (int m = 0; m < croppedImageSize * Nfinal; m++) {
        for (int k = 0; k < croppedImageSize * Nfinal; k++) {
          subpixelPhaseMatrixXZ[k][m] = croppedPhaseValues3D[(int) (k / Nfinal)][(int) m_R0][(int) (m / Nfinal)];
        }
      }

      // Creating ImageProcessors
      ImageProcessor IP_subpixelMagImage = new FloatProcessor(size_subpixelimg, size_subpixelimg);
      ImageProcessor IP_subpixelMagImageXZ = new FloatProcessor(size_subpixelimg, size_subpixelimg);
      ImageProcessor IP_subpixelPhaseImage = new FloatProcessor(size_subpixelimg, size_subpixelimg);
      ImageProcessor IP_subpixelPhaseImageXZ = new FloatProcessor(size_subpixelimg, size_subpixelimg);

      // Adding mag and phase data to image processors (with no bkg phase)
      for (int i = 0; i < size_subpixelimg; i++) {
        for (int j = 0; j < size_subpixelimg; j++) {
          IP_subpixelMagImage.putPixelValue(i, j, subpixelMagMatrix[i][j]);
          IP_subpixelMagImageXZ.putPixelValue(i, j, subpixelMagMatrixXZ[i][j]);
          IP_subpixelPhaseImage.putPixelValue(i, j, subpixelPhaseMatrix[i][j] - item.bkgPhase);
          IP_subpixelPhaseImageXZ.putPixelValue(i, j, subpixelPhaseMatrixXZ[i][j] - item.bkgPhase);
        }
      }

      // Creating new ImagePlus objects with respective ImageProcessors to be
      // displayed
      subpixelMagImage = new ImagePlus(subMagTitle, IP_subpixelMagImage);
      subpixelMagImageXZ = new ImagePlus(subMagXZTitle, IP_subpixelMagImageXZ);
      subpixelPhaseImage = new ImagePlus(subPhaseTitle, IP_subpixelPhaseImage);
      subpixelPhaseImageXZ = new ImagePlus(subPhaseXZTitle, IP_subpixelPhaseImageXZ);

      JOptionPane.showMessageDialog(gui.frame, "Background Phase "
          + String.valueOf(Math.round(item.bkgPhase * 100.0) / 100.0)
          + " Has Been Removed");

      // Displaying images
      subpixelMagImage.show();
      subpixelMagImageXZ.show();
      subpixelPhaseImage.show();
      subpixelPhaseImageXZ.show();

      // Passing matrices to C++
      jni.setPhaseXYMatrix(subpixelPhaseMatrix);
      jni.setPhaseXZMatrix(subpixelPhaseMatrixXZ);
      jni.setMagXYMatrix(subpixelMagMatrix);
      jni.setMagXZMatrix(subpixelMagMatrixXZ);

      // resetting slice to center
      magnitudeImage.setSlice(item.centerS().get(2).intValue() + 1);
      phaseImage.setSlice(item.centerS().get(2).intValue() + 1);

      // creating ROIs for later use
      roiImgMag = new ROIS(subMagTitle);
      roiImgMagXZ = new ROIS(subMagXZTitle);
      roiImgPhase = new ROIS(subPhaseTitle);
      roiImgPhaseXZ = new ROIS(subPhaseXZTitle);
    }
  }

  // =====================================================================================
  // "Remove Bkg"
  // =====================================================================================
  public static void remove_bkg() {

    // if subpixel images are not generated
    if (subpixelMagImage == null || subpixelMagImageXZ == null || subpixelPhaseImage == null
        || subpixelPhaseImageXZ == null) {
      JOptionPane.showMessageDialog(gui.frame, "Error: subpixel images not generated");
      return;
    }

    double bkg_phase = 0.0;

    // not sure how this warning would ever present itself, but including as a
    // safeguard just in case
    if (item == null) {
      JOptionPane.showMessageDialog(gui.frame, "Warning: no background phase found. Continuing with a value of 0.0");
    } else {
      bkg_phase = item.bkgPhase;
    }

    // Removing BG phase in C++
    jni.removeBackgroundPhase(bkg_phase);

    JOptionPane.showMessageDialog(gui.frame, "Background Phase "
        + String.valueOf(Math.round(bkg_phase * 100.0) / 100.0)
        + " Has Been Removed");
  }

  // =====================================================================================
  // "Estimate Subpixel Center"
  // =====================================================================================
  public static void est_subpix_ctr() {

    // updateVariables();
    // jni.setmVariables(grid, m_R0, RCenter,
    // Double.parseDouble(gui.ltf_rcx.getValue()),
    // Double.parseDouble(gui.ltf_rcy.getValue()),
    // Double.parseDouble(gui.ltf_rcz.getValue()) - 1.0,
    // Double.parseDouble(gui.ltf_eqPhase.getValue()));

    // condition for program to continue, must have generated subpixel and estimated
    // a center and RCenter
    boolean condition = !gui.ltf_rcx.getValue().isEmpty() && !gui.ltf_rcy.getValue().isEmpty()
        && !gui.ltf_rcz.getValue().isEmpty()
        && !gui.ltf_rc.getValue().isEmpty()
        && !gui.ll_estBkgPhase.getValue().isEmpty()
        && (WindowManager.getImage(subMagTitle) != null && WindowManager.getImage(subMagXZTitle) != null
            && WindowManager.getImage(subPhaseTitle) != null && WindowManager.getImage(subPhaseXZTitle) != null);

    // if the condition to proceed is false, display message and return
    if (!condition) {
      JOptionPane.showMessageDialog(gui.frame, "Error: Subpixel images not generated and/or step 2 was not completed");
      return;
    }

    // passing needed values to C++

    double RCenter = Double.parseDouble(gui.ltf_rc.getValue());
    double m_R0 = item.m_R0();
    double csx = Double.parseDouble(gui.ltf_rcx.getValue());
    double csy = Double.parseDouble(gui.ltf_rcy.getValue());
    double csz = Double.parseDouble(gui.ltf_rcz.getValue()) - 1.0;
    int c2x = (int) csx;
    int c2y = (int) csy;
    int c2z = (int) csz;
    jni.passCalcSubCenterValues(item.roi_mag_belowM_xi, item.roi_mag_belowM_dx, item.roi_mag_belowM_yi,
        item.roi_mag_belowM_dy, item.roi_mag_belowM_zi, item.roi_mag_belowM_dz, RCenter, m_R0,
        item.centerL().get(0), item.centerL().get(1), item.centerL().get(2),
        item.centerM().get(0), item.centerM().get(1), item.centerM().get(2),
        csx, csy, csz,
        (double) c2x, (double) c2y, (double) c2z);

    // Calculating subpixel center, if there are no errors then the returned string
    // will be empty
    subCenterErrorMessage = jni.estimateSubpixelCenter();

    // if C++ sends back an error message
    if (subCenterErrorMessage.compareTo("") != 0) {
      JOptionPane.showMessageDialog(gui.frame, subCenterErrorMessage);
      return;
    }

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
    gui.ltf_spx.setValue(String.valueOf(Math.round(centerX_pixelCoordinates * 100.0) / 100.0));
    gui.ltf_spy.setValue(String.valueOf(Math.round(centerY_pixelCoordinates * 100.0) / 100.0));
    gui.ltf_spz.setValue(String.valueOf(Math.round(centerZ_pixelCoordinates * 100.0) / 100.0 + 1.0));

    // ---------- begin to put ROIs on images

    /*
     * The following is used for displaying the centers and RCenter on the images.
     * The program uses another class I built called ROIS that basically handles
     * ImageJ's RoiManager in a better way. ImageJ makes it very annoying to display
     * multiple ROI's on an image, so this way the code is much easier to read and
     * is not repetitively long. You can read the ROIS Javadoc, as it is already
     * written out, and may have to be compiled to view.
     */

    int sub_x = (int) pixelToSubpixel(centerX_pixelCoordinates, 0);
    int sub_y = (int) pixelToSubpixel(centerY_pixelCoordinates, 1);
    int sub_z = (int) pixelToSubpixel(centerZ_pixelCoordinates, 2);
    logger.addVariable("subxyz (est_subpix)",
        String.valueOf(sub_x) + ' ' + String.valueOf(sub_y) + ' ' + String.valueOf(sub_z));

    // Creating a new ROIS for the XY mag image
    roiImgMag.clear();
    // Adding the center as a point to the list
    roiImgMag.addPointROI(sub_x, sub_y);
    // if (gui.chkbx_showrc.isSelected()) {
    // // If the box is selected, adds RCenter circle ROI to the list
    // roiImgMag.addCircleROI(sub_x, sub_y, RCenter * 10.0);
    // }
    // Displays the ROIS on the image
    roiImgMag.displayROIS();

    roiImgMagXZ.clear();
    roiImgMagXZ.addPointROI(sub_x, sub_z);
    // if (gui.chkbx_showrc.isSelected()) {
    // roiImgMagXZ.addCircleROI(sub_x, sub_z, RCenter * 10.0);
    // }
    roiImgMagXZ.displayROIS();

    roiImgPhase.clear();
    roiImgPhase.addPointROI(sub_x, sub_y);
    // if (gui.chkbx_showrc.isSelected()) {
    // roiImgPhase.addCircleROI(sub_x, sub_y, RCenter * 10.0);
    // }
    roiImgPhase.displayROIS();

    roiImgPhaseXZ.clear();
    roiImgPhaseXZ.addPointROI(sub_x, sub_z);
    // if (gui.chkbx_showrc.isSelected()) {
    // roiImgPhaseXZ.addCircleROI(sub_x, sub_z, RCenter * 10.0);
    // }
    roiImgPhaseXZ.displayROIS();

    // ---------- end to put ROIS on images

  }

  // =====================================================================================
  // "Redraw Center"
  // =====================================================================================
  public static void redraw_center() {
    // updateVariables();

    // we need this to be true in order to have sufficient values in the GUI
    boolean condition = (WindowManager.getImage(subMagTitle) != null && WindowManager.getImage(subMagXZTitle) != null
        && WindowManager.getImage(subPhaseTitle) != null && WindowManager.getImage(subPhaseXZTitle) != null)
        && !(gui.ltf_spx.getValue().isEmpty() || gui.ltf_spy.getValue().isEmpty() || gui.ltf_spz.getValue().isEmpty());

    // if the condition doesn't pass, tell the user then return
    if (!condition) {
      JOptionPane.showMessageDialog(gui.frame, "Error: subpixel images/center not generated/calculated");
      return;
    }

    // get rcenter from GUI
    double RCenter = Double.parseDouble(gui.ltf_rc.getValue());

    // update item center based off current subcenter GUI values
    item.setCenterSX(Double.parseDouble(gui.ltf_spx.getValue()));
    item.setCenterSY(Double.parseDouble(gui.ltf_spy.getValue()));
    item.setCenterSZ(Double.parseDouble(gui.ltf_spz.getValue()) - 1.0);
    logger.addVariable("Center redrawn as", String.valueOf(Double.parseDouble(gui.ltf_spx.getValue()))
        + ',' + String.valueOf(Double.parseDouble(gui.ltf_spy.getValue())) + ','
        + String.valueOf(Double.parseDouble(gui.ltf_spz.getValue()) - 1.0));
    logger.addVariable("Center returned values", String.valueOf(item.centerS().get(0)) + ','
        + String.valueOf(item.centerS().get(1)) + ',' + String.valueOf(item.centerS().get(2)));

    // This is basically a repetition from the btn_estSubC code

    // getting sub center
    int sub_x = (int) pixelToSubpixel(item.centerS().get(0), 0);
    int sub_y = (int) pixelToSubpixel(item.centerS().get(1), 1);
    int sub_z = (int) pixelToSubpixel(item.centerS().get(2), 2);
    logger.addVariable("subxyz (in redraw)",
        String.valueOf(sub_x) + ',' + String.valueOf(sub_y) + ',' + String.valueOf(sub_z));

    // Creating a new ROIS for the XY mag image
    // roiImgMag = new ROIS(subMagTitle);
    // Adding the center as a point to the list
    // roiImgMag.addPointROI("Center", sub_x, sub_y);
    // if (gui.chkbx_showrc.isSelected()) {
    // If the box is selected, adds RCenter circle ROI to the list
    // roiImgMag.addCircleROI("RCenter", sub_x, sub_y, RCenter * 10.0);
    // }
    // Displays the ROIS on the image
    // roiImgMag.displayROIS();

    roiImgMag.clear();
    roiImgMag.addPointROI(sub_x, sub_y);
    // if (gui.chkbx_showrc.isSelected()) {
    // roiImgMag.addCircleROI(sub_x, sub_y, RCenter * 10.0);
    // }
    roiImgMag.displayROIS();

    // The next 3 blocks of code follow the same logic as described above

    roiImgMagXZ.clear();
    roiImgMagXZ.addPointROI(sub_x, sub_z);
    // if (gui.chkbx_showrc.isSelected()) {
    // roiImgMagXZ.addCircleROI(sub_x, sub_z, RCenter * 10.0);
    // }
    roiImgMagXZ.displayROIS();

    roiImgPhase.clear();
    roiImgPhase.addPointROI(sub_x, sub_y);
    // if (gui.chkbx_showrc.isSelected()) {
    // roiImgPhase.addCircleROI(sub_x, sub_y, RCenter * 10.0);
    // }
    roiImgPhase.displayROIS();

    roiImgPhaseXZ.clear();
    roiImgPhaseXZ.addPointROI(sub_x, sub_z);
    // if (gui.chkbx_showrc.isSelected()) {
    // roiImgPhaseXZ.addCircleROI(sub_x, sub_z, RCenter * 10.0);
    // }
    roiImgPhaseXZ.displayROIS();

    // ---------- end to put ROIS on images
  }

  // =====================================================================================
  // "Show RCenter"
  // =====================================================================================
  public static void showRCenter() {

    boolean condition = (WindowManager.getImage(subMagTitle) != null && WindowManager.getImage(subMagXZTitle) != null
        && WindowManager.getImage(subPhaseTitle) != null && WindowManager.getImage(subPhaseXZTitle) != null)
        && !(gui.ltf_rcx.getValue().isEmpty() || gui.ltf_rcy.getValue().isEmpty() || gui.ltf_rcz.getValue().isEmpty()
            && !(gui.ltf_rc.getValue().isEmpty()));

    if (!condition) {
      JOptionPane.showMessageDialog(gui.frame, "Error: subpixel images and/or center not found");
      return;
    }

    double csx = Double.parseDouble(gui.ltf_rcx.getValue());
    double csy = Double.parseDouble(gui.ltf_rcy.getValue());
    double csz = Double.parseDouble(gui.ltf_rcz.getValue()) - 1.0;

    // get rcenter from GUI
    double RCenter = Double.parseDouble(gui.ltf_rc.getValue());

    roiImgMag.clear();
    roiImgMagXZ.clear();
    roiImgPhase.clear();
    roiImgPhaseXZ.clear();

    // convert center to subpixel coordinates
    int sub_x = (int) pixelToSubpixel(csx, 0);
    int sub_y = (int) pixelToSubpixel(csy, 1);
    int sub_z = (int) pixelToSubpixel(csz, 2);

    // add point ROIS to each image
    roiImgMag.addPointROI(sub_x, sub_y);
    roiImgMagXZ.addPointROI(sub_x, sub_z);
    roiImgPhase.addPointROI(sub_x, sub_y);
    roiImgPhaseXZ.addPointROI(sub_x, sub_z);

    // add circle ROIS to images
    roiImgMag.addCircleROI(sub_x, sub_y, RCenter * 10.0);
    roiImgMagXZ.addCircleROI(sub_x, sub_z, RCenter * 10.0);
    roiImgPhase.addCircleROI(sub_x, sub_y, RCenter * 10.0);
    roiImgPhaseXZ.addCircleROI(sub_x, sub_z, RCenter * 10.0);

    // display ROIS
    roiImgMag.displayROIS();
    roiImgMagXZ.displayROIS();
    roiImgPhase.displayROIS();
    roiImgPhaseXZ.displayROIS();
  }

  // =====================================================================================
  // "Verify Radii"
  // =====================================================================================
  public static void verify_radii() {
    // updateVariables();

    boolean condition = (WindowManager.getImage(subMagTitle) != null && WindowManager.getImage(subMagXZTitle) != null
        && WindowManager.getImage(subPhaseTitle) != null && WindowManager.getImage(subPhaseXZTitle) != null)
        && !(gui.ltf_spx.getValue().isEmpty() || gui.ltf_spy.getValue().isEmpty() || gui.ltf_spz.getValue().isEmpty());

    // if condition doesn't pass then tell user and return
    if (!condition) {
      JOptionPane.showMessageDialog(gui.frame, "Error: Insufficient data to verify radii");
      return;
    }

    double m_R1 = Double.parseDouble(gui.ltf_r1.getValue());
    double m_R2 = Double.parseDouble(gui.ltf_r2.getValue());
    double m_R3 = Double.parseDouble(gui.ltf_r3.getValue());

    // update item center based off current subcenter GUI values
    item.setCenterSX(Double.parseDouble(gui.ltf_spx.getValue()));
    item.setCenterSY(Double.parseDouble(gui.ltf_spy.getValue()));
    item.setCenterSZ(Double.parseDouble(gui.ltf_spz.getValue()) - 1.0);

    int sub_x = (int) pixelToSubpixel(item.centerS().get(0), 0);
    int sub_y = (int) pixelToSubpixel(item.centerS().get(1), 1);
    int sub_z = (int) pixelToSubpixel(item.centerS().get(2), 2);
    logger.addVariable("subxyz", String.valueOf(sub_x) + ',' + String.valueOf(sub_y) + ',' + String.valueOf(sub_z));

    roiImgMag.clear();
    roiImgMag.addPointROI(sub_x, sub_y);
    roiImgMag.addCircleROI(sub_x, sub_y, m_R1 * 10.0);
    roiImgMag.addCircleROI(sub_x, sub_y, m_R2 * 10.0);
    roiImgMag.addCircleROI(sub_x, sub_y, m_R3 * 10.0);
    roiImgMag.displayROIS();

    roiImgMagXZ.clear();
    roiImgMagXZ.addPointROI(sub_x, sub_z);
    roiImgMagXZ.addCircleROI(sub_x, sub_z, m_R1 * 10.0);
    roiImgMagXZ.addCircleROI(sub_x, sub_z, m_R2 * 10.0);
    roiImgMagXZ.addCircleROI(sub_x, sub_z, m_R3 * 10.0);
    roiImgMagXZ.displayROIS();

    roiImgPhase.clear();
    roiImgPhase.addPointROI(sub_x, sub_y);
    roiImgPhase.addCircleROI(sub_x, sub_y, m_R1 * 10.0);
    roiImgPhase.addCircleROI(sub_x, sub_y, m_R2 * 10.0);
    roiImgPhase.addCircleROI(sub_x, sub_y, m_R3 * 10.0);
    roiImgPhase.displayROIS();

    roiImgPhaseXZ.clear();
    roiImgPhaseXZ.addPointROI(sub_x, sub_z);
    roiImgPhaseXZ.addCircleROI(sub_x, sub_z, m_R1 * 10.0);
    roiImgPhaseXZ.addCircleROI(sub_x, sub_z, m_R2 * 10.0);
    roiImgPhaseXZ.addCircleROI(sub_x, sub_z, m_R3 * 10.0);
    roiImgPhaseXZ.displayROIS();

    float R1_phase_actual = 0;
    float R2_phase_actual = 0;
    float R3_phase_actual = 0;

    // Summing up all phase values where the radii and equitorial axis intercept
    switch (item.MRIAxis()) {
      case X:

        R1_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(sub_x, sub_y + (int) (m_R1 * 10));
        logger.addVariable("got coordinate", String.valueOf(sub_x) + ',' + String.valueOf(sub_y + (int) (m_R1 * 10)));

        R1_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(sub_x, sub_y - (int) (m_R1 * 10));
        logger.addVariable("got coordinate", String.valueOf(sub_x) + ',' + String.valueOf(sub_y - (int) (m_R1 * 10)));

        R1_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(sub_x, sub_z + (int) (m_R1 * 10));
        logger.addVariable("got coordinate", String.valueOf(sub_x) + ',' + String.valueOf(sub_z + (int) (m_R1 * 10)));

        R1_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(sub_x, sub_z - (int) (m_R1 * 10));
        logger.addVariable("got coordinate", String.valueOf(sub_x) + ',' + String.valueOf(sub_z - (int) (m_R1 * 10)));

        R2_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(sub_x, sub_y + (int) (m_R2 * 10));
        R2_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(sub_x, sub_y - (int) (m_R2 * 10));
        R2_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(sub_x, sub_z + (int) (m_R2 * 10));
        R2_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(sub_x, sub_z - (int) (m_R2 * 10));

        R3_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(sub_x, sub_y + (int) (m_R3 * 10));
        R3_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(sub_x, sub_y - (int) (m_R3 * 10));
        R3_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(sub_x, sub_z + (int) (m_R3 * 10));
        R3_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(sub_x, sub_z - (int) (m_R3 * 10));
        break;

      case Y:
        R1_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(sub_x + (int) (m_R1 * 10), sub_z);
        R1_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(sub_x - (int) (m_R1 * 10), sub_z);
        R1_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(sub_x, sub_z + (int) (m_R1 * 10));
        R1_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(sub_x, sub_z - (int) (m_R1 * 10));

        R2_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(sub_x + (int) (m_R2 * 10), sub_z);
        R2_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(sub_x - (int) (m_R2 * 10), sub_z);
        R2_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(sub_x, sub_z + (int) (m_R2 * 10));
        R2_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(sub_x, sub_z - (int) (m_R2 * 10));

        R3_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(sub_x + (int) (m_R3 * 10), sub_z);
        R3_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(sub_x - (int) (m_R3 * 10), sub_z);
        R3_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(sub_x, sub_z + (int) (m_R3 * 10));
        R3_phase_actual += subpixelPhaseImageXZ.getProcessor().getPixelValue(sub_x, sub_z - (int) (m_R3 * 10));
        break;

      case Z:
        R1_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(sub_x, sub_y + (int) (m_R1 * 10));
        R1_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(sub_x, sub_y - (int) (m_R1 * 10));
        R1_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(sub_x + (int) (m_R1 * 10), sub_y);
        R1_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(sub_x - (int) (m_R1 * 10), sub_y);

        R2_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(sub_x, sub_y + (int) (m_R2 * 10));
        R2_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(sub_x, sub_y - (int) (m_R2 * 10));
        R2_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(sub_x + (int) (m_R2 * 10), sub_y);
        R2_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(sub_x - (int) (m_R2 * 10), sub_y);

        R3_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(sub_x, sub_y + (int) (m_R3 * 10));
        R3_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(sub_x, sub_y - (int) (m_R3 * 10));
        R3_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(sub_x + (int) (m_R3 * 10), sub_y);
        R3_phase_actual += subpixelPhaseImage.getProcessor().getPixelValue(sub_x - (int) (m_R3 * 10), sub_y);
        break;

      default:
        JOptionPane.showMessageDialog(gui.frame, "Error in Verify Radii: No MRI field axis found");
        break;
    }

    // Dividing each phase sum by 4 for average
    R1_phase_actual /= 4.0;
    R2_phase_actual /= 4.0;
    R3_phase_actual /= 4.0;
    // Removing background phase off each phase value
    // R1_phase_actual -= item.bkgPhase;
    // R2_phase_actual -= item.bkgPhase;
    // R3_phase_actual -= item.bkgPhase;

    gui.lbl_r1phaseAct.setText(String.valueOf(Math.round(R1_phase_actual * 100.0) / 100.0));
    gui.lbl_r2phaseAct.setText(String.valueOf(Math.round(R2_phase_actual * 100.0) / 100.0));
    gui.lbl_r3phaseAct.setText(String.valueOf(Math.round(R3_phase_actual * 100.0) / 100.0));
  }

  // =====================================================================================
  // "Estimate Bkg & rho_0"
  // =====================================================================================
  public static void est_bkg_dens() {
    // updateVariables();

    // Condition for code to run - all radii must be found and the subpixel images
    // and center must be found
    boolean condition = !(gui.ltf_r1.getValue().isEmpty()) && !(gui.ltf_r2.getValue().isEmpty())
        && !(gui.ltf_r3.getValue().isEmpty())
        && !(gui.ltf_spx.getValue().isEmpty() || gui.ltf_spy.getValue().isEmpty() || gui.ltf_spz.getValue().isEmpty())
        && !gui.ltf_magMom.getValue().isEmpty()
        && (WindowManager.getImage(subMagTitle) != null && WindowManager.getImage(subMagXZTitle) != null
            && WindowManager.getImage(subPhaseTitle) != null && WindowManager.getImage(subPhaseXZTitle) != null);

    if (!condition) {
      JOptionPane.showMessageDialog(gui.frame,
          "Error: Insufficient data to calculate spin density and background phase");
      return;
    }

    // passing required values to C++
    double mag_mom = Double.parseDouble(gui.ltf_magMom.getValue());
    double r1 = Double.parseDouble(gui.ltf_r1.getValue());
    double r2 = Double.parseDouble(gui.ltf_r2.getValue());
    double r3 = Double.parseDouble(gui.ltf_r3.getValue());
    double cx = Double.parseDouble(gui.ltf_spx.getValue());
    double cy = Double.parseDouble(gui.ltf_spy.getValue());
    double cz = Double.parseDouble(gui.ltf_spz.getValue()) - 1.0;
    logger.addVariable("bkgphase", item.bkgPhase);
    jni.passSpinDensValues(cx, cy, cz, r1, r2, r3, item.bkgPhase, mag_mom);

    // Calling C++ to calculate background phase and spin density
    jni.estBkgAndSpinDensity();

    // Getting calculated background phase from C++
    item.bkgPhase = Math.abs(jni.getBkg());

    // Getting calculated spin density from C++
    double spinDensity = jni.getSpinDensity();

    logger.addVariable("estimatedBackgroundPhase", item.bkgPhase);
    logger.addVariable("spinDensity", spinDensity);

    // Setting background phase and spin density to GUI
    gui.ll_estBkgPhase.setValue(String.valueOf(Math.round(item.bkgPhase * 100.0) / 100.0));
    gui.ll_rho0.setValue(String.valueOf(Math.round(spinDensity * 100.0) / 100.0));
  }

  /*
   * If the plot X button is clicked:
   * This displays a graph of phase values along the x-axis through the subpixel
   * center
   */
  public static void plot_x() {
    // updateVariables();

    boolean condition = !(gui.ltf_spx.getValue().isEmpty() || gui.ltf_spy.getValue().isEmpty()
        || gui.ltf_spz.getValue().isEmpty())
        && (WindowManager.getImage(subMagTitle) != null && WindowManager.getImage(subMagXZTitle) != null
            && WindowManager.getImage(subPhaseTitle) != null && WindowManager.getImage(subPhaseXZTitle) != null);

    if (!condition) {
      JOptionPane.showMessageDialog(gui.frame, "Error: Subpixel images/center not generated/calculated");
      return;
    }

    // Initializing plot object for graph
    Plot xPlot = new Plot("X-Profile", "Location", "Phase");

    // Initializing ArrayLists for graph values on each axis
    ArrayList<Double> intensity = new ArrayList<Double>();
    ArrayList<Double> location = new ArrayList<Double>();

    // Setting colour components of graph
    xPlot.setColor(Color.RED, Color.BLACK);

    try {

      // int x_px = item.centerS().get(0).intValue();
      int x_px = (int) Double.parseDouble(gui.ltf_rcx.getValue());
      // int y_spx = (int) pixelToSubpixel(item.centerS().get(1).intValue(), 1);
      int y_spx = (int) pixelToSubpixel((int) Double.parseDouble(gui.ltf_rcy.getValue()), 1);
      logger.addVariable("x_px, y_spx", String.valueOf(x_px) + ',' + String.valueOf(y_spx));

      // Adding phase values to both ArrayLists
      double m_R0 = item.m_R0();
      logger.addVariable("m_R0", m_R0);

      for (int i = x_px - (int) m_R0,
          c = 0; i <= x_px + (int) m_R0; i++, c++) {
        logger.addVariable("uhhhhh",
            String.valueOf((int) pixelToSubpixel(i, 0)) + ',' + y_spx + ',' + String.valueOf(i - x_px));

        intensity.add((double) subpixelPhaseMatrix[(int) pixelToSubpixel(i, 0)][y_spx]);
        location.add((double) (i - x_px));

        if (i != x_px - (int) m_R0)
          xPlot.drawLine(location.get(c - 1), intensity.get(c - 1), location.get(c), intensity.get(c));

      }

    } catch (Exception exc) {
      JOptionPane.showMessageDialog(gui.frame, exc.toString());
      return;
    }

    try {
      // Adding points to graph and displaying
      xPlot.addPoints(location, intensity, ij.gui.Plot.DIAMOND);
      xPlot.show();
    } catch (Exception exc) {
      JOptionPane.showMessageDialog(gui.frame, exc.toString());
      return;
    }
  }

  /*
   * If the plot Y button is clicked:
   * This displays a graph of phase values along the y-axis through the subpixel
   * center
   */
  public static void plot_y() {
    // updateVariables();

    boolean condition = !(gui.ltf_spx.getValue().isEmpty() || gui.ltf_spy.getValue().isEmpty()
        || gui.ltf_spz.getValue().isEmpty())
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
        double m_R0 = item.m_R0();
        for (int i = y_px - (int) m_R0,
            c = 0; i <= y_px + m_R0; i++, c++) {
          intensity
              .add((double) subpixelPhaseMatrix[x_spx][(int) pixelToSubpixel(i, 1)]);
          location.add((double) (i - y_px));
          if (i != y_px - (int) m_R0)
            yPlot.drawLine(location.get(c - 1), intensity.get(c - 1), location.get(c), intensity.get(c));
        }
      } catch (Exception exc) {
        JOptionPane.showMessageDialog(gui.frame, exc.toString());
      }

      try {
        // Adding points to graph and displaying
        yPlot.addPoints(location, intensity, ij.gui.Plot.DIAMOND);
        yPlot.show();
      } catch (Exception exc) {
        JOptionPane.showMessageDialog(gui.frame, exc.toString());
      }
    } else {
      JOptionPane.showMessageDialog(gui.frame, "Error: Subpixel images not generated");
    }
  }

  /*
   * If the plot Z button is clicked:
   * This displays a graph of phase values along the z-axis through the subpixel
   * center
   */
  public static void plot_z() {
    // updateVariables();

    boolean condition = !(gui.ltf_spx.getValue().isEmpty() || gui.ltf_spy.getValue().isEmpty()
        || gui.ltf_spz.getValue().isEmpty())
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
        double m_R0 = item.m_R0();
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
        JOptionPane.showMessageDialog(gui.frame, exc.toString());
      }

      try {
        // Adding points to graph and displaying
        zPlot.addPoints(location, intensity, ij.gui.Plot.DIAMOND);
        zPlot.show();
      } catch (Exception exc) {
        JOptionPane.showMessageDialog(gui.frame, exc.toString());
      }
    } else {
      JOptionPane.showMessageDialog(gui.frame, "Error: Subpixel images not generated");
    }
  }

  // =====================================================================================
  // "Calculate Magnetic Moment"
  // =====================================================================================
  public static void calc_mag_moment() {
    // updateVariables();

    // Condition for code to run - radii, RCenter, subpixel center and images must
    // be found
    boolean condition = !(gui.ltf_r1.getValue().isEmpty()) && !(gui.ltf_r2.getValue().isEmpty())
        && !(gui.ltf_r3.getValue().isEmpty())
        && !(gui.ltf_spx.getValue().isEmpty() || gui.ltf_spy.getValue().isEmpty() || gui.ltf_spz.getValue().isEmpty())
        && (WindowManager.getImage(subMagTitle) != null && WindowManager.getImage(subMagXZTitle) != null
            && WindowManager.getImage(subPhaseTitle) != null && WindowManager.getImage(subPhaseXZTitle) != null)
        && !gui.ltf_rc.getValue().isEmpty();

    if (!condition) {
      JOptionPane.showMessageDialog(gui.frame, "Error: insufficient data to calculate magnetic moment");
      return;
    }

    String errorMessage_Mag;

    double r1 = Double.parseDouble(gui.ltf_r1.getValue());
    double r2 = Double.parseDouble(gui.ltf_r2.getValue());
    double r3 = Double.parseDouble(gui.ltf_r3.getValue());
    double csx = Double.parseDouble(gui.ltf_spx.getValue());
    double csy = Double.parseDouble(gui.ltf_spy.getValue());
    double csz = Double.parseDouble(gui.ltf_spz.getValue()) - 1.0;
    double mr0 = item.m_R0();
    double bkg_phase = item.bkgPhase;
    double rchi = Double.parseDouble(gui.ltf_RChi.getValue());
    double b0 = Double.parseDouble(gui.ltf_B0.getValue());
    double tefirst = Double.parseDouble(gui.ltf_TEFirst.getValue());
    double snr = Double.parseDouble(gui.ltf_snr.getValue());
    double outer_phase = Double.parseDouble(gui.lbl_r1phaseCalc.getText());
    double middle_phase = Double.parseDouble(gui.lbl_r2phaseCalc.getText());
    double inner_phase = Double.parseDouble(gui.lbl_r3phaseCalc.getText());
    double e12 = Double.parseDouble(gui.ltf_eps12.getValue());
    double e23 = Double.parseDouble(gui.ltf_eps23.getValue());
    jni.passMagMomValues(r1, r2, r3, csx, csy, csz, mr0, bkg_phase, rchi, b0, tefirst, snr, inner_phase, middle_phase,
        outer_phase, e12, e23);

    // Calculating magnetic moment in C++, if function returns empty string then no
    // error message
    errorMessage_Mag = jni.calculateMagneticMoment();

    if (errorMessage_Mag.compareTo("") != 0) {
      JOptionPane.showMessageDialog(gui.frame, errorMessage_Mag);
      return;
    }

    // Getting various values from C++ as a result of the function to calculate the
    // magnetic moment being ran
    gui.lbl_r1phaseCalc.setText(String.valueOf(Math.round(jni.getMR1Calc() * 100.0) / 100.0));
    gui.lbl_r2phaseCalc.setText(String.valueOf(Math.round(jni.getMR2Calc() * 100.0) / 100.0));
    gui.lbl_r3phaseCalc.setText(String.valueOf(Math.round(jni.getMR3Calc() * 100.0) / 100.0));
    gui.ltf_magMom.setValue(String.valueOf(Math.round(jni.getMagMoment() * 100.0) / 100.0));
    if (jni.getUncertainty() == -1.0) {
      gui.ll_momenterror.setValue("");
      JOptionPane.showMessageDialog(gui.frame,
          "<html>Error: Cannot calculate error, make sure SNR, &eta;12 and &eta;23 are set.</html>");
    } else {
      gui.ll_momenterror.setValue(String.valueOf(Math.round(jni.getUncertainty() * 100.0) / 100.0));
    }
    // gui.ll_dChi.setValue(String.valueOf(Math.round(jni.getChi() * 100.0) /
    // 100.0));
    // gui.ll_a.setValue(String.valueOf(Math.round(jni.getA() * 100.0) / 100.0));
    // gui.ll_rho0.setValue(String.valueOf(Math.round(jni.getSpinDensity() * 100.0)
    // / 100.0));
  }

  /*
   * If the load simulated images button is clicked:
   * The program prompts the user to open a new set of simulated images. These are
   * used for finding the uncertainty of the magnetic moment and for calculating
   * eij
   */
  public static void load_sim_images() {

    double snr = Double.parseDouble(gui.ltf_snr.getValue());
    double e12 = Double.parseDouble(gui.ltf_eps12.getValue());
    double e23 = Double.parseDouble(gui.ltf_eps23.getValue());
    double B0 = Double.parseDouble(gui.ltf_B0.getValue());
    double R_Chi = Double.parseDouble(gui.ltf_RChi.getValue());
    double TElast = Double.parseDouble(gui.ltf_TELast.getValue());
    jni.setMagMomentVariables(snr, e12, e23, B0, R_Chi, TElast);

    double m_R1 = Double.parseDouble(gui.ltf_r1.getValue());
    double m_R2 = Double.parseDouble(gui.ltf_r2.getValue());
    double m_R3 = Double.parseDouble(gui.ltf_r3.getValue());

    s5MagWindowTitle = loadImagesIJ(s5MagWindowTitle);
    s5PhaseWindowTitle = loadImagesIJ(s5PhaseWindowTitle);
    /*
     * try {
     * s5MagWindowTitle = loadImages(s5MagWindowTitle, "mag");
     * s5PhaseWindowTitle = loadImages(s5PhaseWindowTitle, "phase");
     * } catch (IOException exc) {
     * return;
     * }
     */

    // End loading simulated images

    if ((WindowManager.getImage(s5MagWindowTitle) != null) && (WindowManager.getImage(s5PhaseWindowTitle) != null)) {

      ImagePlus simulatedMagnitudeImage = WindowManager.getImage(s5MagWindowTitle);
      ImagePlus simulatedPhaseImage = WindowManager.getImage(s5PhaseWindowTitle);

      if (simulatedMagnitudeImage.getNSlices() != simulatedPhaseImage.getNSlices()) {
        JOptionPane.showMessageDialog(gui.frame, "Error: different # of slices in images");
      } else if (simulatedMagnitudeImage.getHeight() != simulatedPhaseImage.getHeight()) {
        JOptionPane.showMessageDialog(gui.frame, "Error: different height in images");
      } else if (simulatedMagnitudeImage.getWidth() != simulatedPhaseImage.getWidth()) {
        JOptionPane.showMessageDialog(gui.frame, "Error: different width in images");
      } else if (gui.ltf_r1.getValue().isEmpty()) {
        JOptionPane.showMessageDialog(gui.frame, "Error: no m_R1 found");
      } else {

        // Size of the matrix that will be fetched from the simulated images
        int sizeOfSimmedMatrices = 2 * ((int) Math.ceil(m_R1) + 1);

        // Initializing matrix for simulated images
        float[][][] simulatedRealNumbers = new float[sizeOfSimmedMatrices + 1][sizeOfSimmedMatrices
            + 1][sizeOfSimmedMatrices + 1];
        float[][][] simulatedImagNumbers = new float[sizeOfSimmedMatrices + 1][sizeOfSimmedMatrices
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
              simulatedImagNumbers[ni][nj][nk] = (float) (tempmag * Math.sin(tempphase));
            }
          }
        }

        logger.addInfo("Finished loops");
        logger.addVariable("simmedCenter", simmedCenterX);
        logger.addVariable("m_R1", m_R1);
        logger.addVariable("m_R2", m_R2);
        logger.addVariable("m_R3", m_R3);

        // Giving matrix to C++
        jni.setSimulatedMatrices(simulatedRealNumbers, simulatedImagNumbers, sizeOfSimmedMatrices + 1);
        logger.addInfo("Set matrices to C++");
        // Interpolating matrix in C++
        jni.interpolateVoxelsSIM(sizeOfSimmedMatrices * 10);
        logger.addInfo("Interpolation complete");

        // Getting current rho_0 value from GUI
        double rho_0 = Double.parseDouble(gui.ll_rho0.getValue());
        // double rho_0 = 10;

        // Getting center of simulated matrix
        int subCenterSIM = sizeOfSimmedMatrices * 10 / 2;

        // Summing real values in each radius in C++
        double S1 = jni.SumCircleElementsReal3DSIMMED((int) (m_R1 * 10.0), subCenterSIM, subCenterSIM, subCenterSIM);
        double S2 = jni.SumCircleElementsReal3DSIMMED((int) (m_R2 * 10.0), subCenterSIM, subCenterSIM, subCenterSIM);
        double S3 = jni.SumCircleElementsReal3DSIMMED((int) (m_R3 * 10.0), subCenterSIM, subCenterSIM, subCenterSIM);

        // Getting phi (phase) values from GUI
        double phi1 = Double.parseDouble(gui.lbl_r1phaseCalc.getText());
        double phi2 = Double.parseDouble(gui.lbl_r2phaseCalc.getText());
        double phi3 = Double.parseDouble(gui.lbl_r3phaseCalc.getText());

        // Calculating fij using equation 9 (real)
        double Re_f12_REAL = (9.0 * Math.sqrt(3)) / (4.0 * Math.PI * rho_0) * (S1 - S2);
        double Re_f23_REAL = (9.0 * Math.sqrt(3)) / (4.0 * Math.PI * rho_0) * (S2 - S3);

        // Calculating fij using equation 10 in C++ (theoretical)
        double mag_moment = Double.parseDouble(gui.ltf_magMom.getValue());
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
        gui.ltf_eps12.setValue(String.valueOf(Math.round(se12 * 100.0) / 100.0));
        gui.ltf_eps23.setValue(String.valueOf(Math.round(se23 * 100.0) / 100.0));

        // Calculating uncertainty in C++
        double uncertainty = jni.calculateUncertainty(se12, se23);

        logger.addVariable("uncertainty", uncertainty);

        // Setting uncertainty to GUI
        gui.ll_momenterror.setValue(String.valueOf(Math.round(uncertainty * 100.0) / 100.0));

      }
    }
  }

  /*
   * If the sum button is clicked:
   * C++ calculates the sum of the real and imag numbers within the user-inputted
   * radius on the subpixel images
   */
  public static void sum_ri() {
    // updateVariables();

    double m_Ri = Double.parseDouble(gui.ltf_Ri.getValue());
    double csx = Double.parseDouble(gui.ltf_spx.getValue());
    double csy = Double.parseDouble(gui.ltf_spy.getValue());
    double csz = Double.parseDouble(gui.ltf_spz.getValue()) - 1.0;
    jni.passSumValues(m_Ri, csx, csy, csz);

    String imag_msg = jni.calculateImagSum();
    String real_msg = jni.calculateRealSum();

    if (imag_msg.compareTo("") == 0) {
      gui.ll_ImRi.setValue(String.valueOf(Math.round(jni.getImagSum() * 100.0) / 100.0));
    } else {
      JOptionPane.showMessageDialog(gui.frame, imag_msg);
    }

    if (real_msg.compareTo("") == 0) {
      gui.ll_ReRi.setValue(String.valueOf(Math.round(jni.getRealSum() * 100.0) / 100.0));
    } else {
      JOptionPane.showMessageDialog(gui.frame, real_msg);
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
  public static void load_first_TEImgs() {

    ImagePlus magnitudeImage = WindowManager.getImage(s1MagWindowTitle);
    ImagePlus phaseImage = WindowManager.getImage(s1PhaseWindowTitle);

    // ---------- Begin to get first TE images

    if (magnitudeImage.isHyperStack() && phaseImage.isHyperStack()) {

      // TODO: handle if second echo is within image that is uploaded in step 1
    } else {
      JOptionPane.showMessageDialog(gui.frame,
          "Error: Image does not contain hyperstacks.\nSelect another image for second echo time.");

      s6MagWindowTitle = loadImagesIJ(s6MagWindowTitle);
      s6PhaseWindowTitle = loadImagesIJ(s6PhaseWindowTitle);
      /*
       * try {
       * s6MagWindowTitle = loadImages(s6MagWindowTitle, "mag");
       * s6PhaseWindowTitle = loadImages(s6PhaseWindowTitle, "phase");
       * } catch (IOException exc) {
       * return;
       * }
       */
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
  public static void noname() {
    if (!(WindowManager.getImage(s6MagWindowTitle) == null || WindowManager.getImage(s6PhaseWindowTitle) == null)) {
      // TODO: handle no updateVariables() call
      // updateVariables();

      // Calculating susceptibility in C++
      jni.calcSusceptibility();

      // Getting values from C++ and setting to GUI
      gui.ll_a.setValue(String.valueOf(Math.round(jni.getA() * 100.0) / 100.0));
      gui.ll_dChi.setValue(String.valueOf(Math.round(jni.getChi() * 100.0) / 100.0));
    }
  }

  /*
   * If the load spin echo image button is selected:
   * The program gets the user to open a spin echo image (see section 2.6 of
   * smoment.pdf for more detail)
   */
  public static void load_spin_echoImgs() {

    s7WindowTitle = loadImagesIJ(s7WindowTitle);
    /*
     * try {
     * s7WindowTitle = loadImages(s7WindowTitle, "mag");
     * } catch (IOException exc) {
     * return;
     * }
     */

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
  public static void est_radius_spin_echo() {

    if (WindowManager.getImage(s7WindowTitle) == null) {
      JOptionPane.showMessageDialog(gui.frame, "Error: no spin echo image found");
      return;
    }

    // TODO: handle no updateVariables() call
    // updateVariables();

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
    int V1SE_x1 = (int) Double.parseDouble(gui.ltf_v1seX1.getValue());
    int V1SE_x2 = (int) Double.parseDouble(gui.ltf_v1seX2.getValue());
    int V1SE_y1 = (int) Double.parseDouble(gui.ltf_v1seY1.getValue());
    int V1SE_y2 = (int) Double.parseDouble(gui.ltf_v1seY2.getValue());
    int V1SE_z1 = (int) Double.parseDouble(gui.ltf_v1seZ1.getValue()) - 1;
    int V1SE_z2 = (int) Double.parseDouble(gui.ltf_v1seZ2.getValue()) - 1;
    int V2SE_x1 = (int) Double.parseDouble(gui.ltf_v2seX1.getValue());
    int V2SE_x2 = (int) Double.parseDouble(gui.ltf_v2seX2.getValue());
    int V2SE_y1 = (int) Double.parseDouble(gui.ltf_v2seY1.getValue());
    int V2SE_y2 = (int) Double.parseDouble(gui.ltf_v2seY2.getValue());
    int V2SE_z1 = (int) Double.parseDouble(gui.ltf_v2seZ1.getValue()) - 1;
    int V2SE_z2 = (int) Double.parseDouble(gui.ltf_v2seZ2.getValue()) - 1;

    int V1SE_size = (V1SE_x2 - V1SE_x1) * (V1SE_y2 - V1SE_y1) * (V1SE_z2 - V1SE_z1);
    int V2SE_size = (V2SE_x2 - V2SE_x1) * (V2SE_y2 - V2SE_y1) * (V2SE_z2 - V2SE_z1);

    // Getting center from GUI
    int VSE_centerX = (int) Double.parseDouble(gui.ltf_secondImgX.getValue());
    int VSE_centerY = (int) Double.parseDouble(gui.ltf_secondImgY.getValue());
    int VSE_centerZ = (int) Double.parseDouble(gui.ltf_secondImgZ.getValue()) - 1;

    // V1 must be bigger than V2 and center must be inside both
    if (V2SE_size > V1SE_size) {
      JOptionPane.showMessageDialog(gui.frame, "Error: V2 cannot be bigger than V1");
    } else if (VSE_centerX < V1SE_x1 || VSE_centerX > V1SE_x2) {
      JOptionPane.showMessageDialog(gui.frame, "Error: center x coordinate must be on image");
    } else if (VSE_centerY < V1SE_y1 || VSE_centerY > V1SE_y2) {
      JOptionPane.showMessageDialog(gui.frame, "Error: center y coordinate must be on image");
    } else if (VSE_centerZ < V1SE_z1 || VSE_centerZ > V1SE_z2) {
      JOptionPane.showMessageDialog(gui.frame, "Error: center z coordinate must be on image");
    } else {

      // Initializing ImageProcessor objects to add data values to
      ImageProcessor IP_V1SE_XY = new FloatProcessor(V1SE_x2 - V1SE_x1, V1SE_y2 - V1SE_y1);
      ImageProcessor IP_V1SE_XZ = new FloatProcessor(V1SE_x2 - V1SE_x1, V1SE_z2 - V1SE_z1);

      // Adding data to new V1 image XY. Essentially just 'cropping' the V1 image
      spinEchoImg.setSlice(VSE_centerZ + 1);
      for (int i = V1SE_x1; i <= V1SE_x2; i++) {
        for (int j = V1SE_y1; j <= V1SE_y2; j++) {
          IP_V1SE_XY.putPixelValue(i - V1SE_x1, j - V1SE_y1, spinEchoImg.getProcessor().getPixelValue(i, j));
        }
      }

      // Same thing but for XZ
      for (int i = V1SE_x1; i <= V1SE_x2; i++) {
        for (int j = V1SE_z1; j <= V1SE_z2; j++) {
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
      roiImgV1SE = new ROIS(V1XY_Title);
      roiImgV1SE.addPointROI(Math.abs(VSE_centerX - V1SE_x1), Math.abs(VSE_centerY - V1SE_y1));
      roiImgV1SE.addRectangle(Math.abs(V1SE_x1 - V2SE_x1), Math.abs(V1SE_y1 - V2SE_y1), V2SE_x2 - V2SE_x1,
          V2SE_y2 - V2SE_y1);
      roiImgV1SE.displayROIS();

      // Same thing but for XZ
      roiImgV1SEXZ = new ROIS(V1XZ_Title);
      roiImgV1SEXZ.addPointROI(Math.abs(VSE_centerX - V1SE_x1), Math.abs(VSE_centerZ - V1SE_z1));
      roiImgV1SEXZ.addRectangle(Math.abs(V1SE_x1 - V2SE_x1), Math.abs(V1SE_z1 - V2SE_z1), V2SE_x2 - V2SE_x1,
          V2SE_z2 - V2SE_z1);
      roiImgV1SEXZ.displayROIS();

      // Variables for the sum and volume of the boxes
      double S1_SE = 0.0;
      double S2_SE = 0.0;
      double V1_SE = 0.0;
      double V2_SE = 0.0;

      logger.addVariable("v1se ranges", String.valueOf(V1SE_x1) + " " + String.valueOf(V1SE_x2));
      logger.addVariable("v1se ranges", String.valueOf(V1SE_y1) + " " + String.valueOf(V1SE_y2));
      logger.addVariable("v1se ranges", String.valueOf(V1SE_z1) + " " + String.valueOf(V1SE_z2));

      // Getting S1 sum and volume
      for (int k = V1SE_z1; k <= V1SE_z2; k++) {
        spinEchoImg.setSlice(k + 1);
        for (int j = V1SE_y1; j <= V1SE_y2; j++) {
          for (int i = V1SE_x1; i <= V1SE_x2; i++, V1_SE++) {
            S1_SE += spinEchoImg.getProcessor().getPixelValue(i, j);
          }
        }
      }

      logger.addVariable("v2se ranges", String.valueOf(V2SE_x1) + " " + String.valueOf(V2SE_x2));
      logger.addVariable("v2se ranges", String.valueOf(V2SE_y1) + " " + String.valueOf(V2SE_y2));
      logger.addVariable("v2se ranges", String.valueOf(V2SE_z1) + " " + String.valueOf(V2SE_z2));

      // Getting S2 sum and volume
      for (int k = V2SE_z1; k <= V2SE_z2; k++) {
        spinEchoImg.setSlice(k + 1);
        for (int j = V2SE_y1; j <= V2SE_y2; j++) {
          for (int i = V2SE_x1; i <= V2SE_x2; i++, V2_SE++) {
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
      gui.ll_V0.setValue(String.valueOf(Math.round(V0 * 100.0) / 100.0));
      logger.addVariable("V0", V0);

      // Equation 16
      double a = Math.pow((V0 * 3.0) / (4.0 * Math.PI), 1.0 / 3.0);

      // Equation 17
      double rho_SE0 = S2_SE / (V2_SE - V0);
      // Setting to GUI
      gui.ll_rho0SE.setValue(String.valueOf(Math.round(rho_SE0 * 100.0) / 100.0));
      logger.addVariable("rho_SE0", rho_SE0);

      double deltaV = 1.0;

      // Defining signal to noise ratio
      double snrStandardDeviation = Double.parseDouble(gui.ltf_sigSE.getValue());
      double SNR_SE = rho_SE0 / snrStandardDeviation;

      // Equation 19
      double dV0 = (Math.sqrt(deltaV) / SNR_SE) * Math.sqrt(V2_SE + Math.pow(V2_SE - V0, 2) / (V1_SE - V2_SE));

      // Adding v0 and dv0 to GUI
      gui.ll_V0.setValue(String.valueOf(Math.round(V0 * 100.0) / 100.0) + " " + PLUS_MINUS + " "
          + String.valueOf(Math.round(dV0 * 100.0) / 100.0));
      logger.addVariable("dV0", dV0);

      // Error for a - was derived with Dr. Cheng
      double da = (a * dV0) / (3 * V0);
      // Adding to GUI
      gui.ll_aSE.setValue(String.valueOf(Math.round(a * 100.0) / 100.0) + " " + PLUS_MINUS + " "
          + String.valueOf(Math.round(da * 100.0) / 100.0));
      logger.addVariable("a", a);
      logger.addVariable("da", da);

      double B0 = Double.parseDouble(gui.ltf_B0.getValue());
      double TELast = Double.parseDouble(gui.ltf_TELast.getValue()) / 1000.0;
      double mag_moment = Double.parseDouble(gui.ltf_magMom.getValue());
      double d_mag_moment = Double.parseDouble(gui.ll_momenterror.getValue());

      // p = ga^3 can be rewritten to get dChi
      double dChi = (2.0 * mag_moment) / (GAMMARBAR * B0 * TELast * V0);
      // And its error
      double d_dChi = dChi * Math.sqrt(Math.pow(d_mag_moment / mag_moment, 2) + Math.pow(dV0 / V0, 2));

      // Setting to GUI
      gui.ll_dChiSE.setValue(String.valueOf(Math.round(dChi * 100.0) / 100.0) + " " + PLUS_MINUS + " "
          + String.valueOf(Math.round(d_dChi * 100.0) / 100.0));
    }
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
  private static double pixelToSubpixel(double coordinate, int axisFlag) {
    // Calculated by taking total size of cropped image and dividing it by 2 - this
    // will give the center which is also the estimated center of the image in step
    // 2
    // The + 4.9 moves it from the middle of the subpixel to the maximum edge
    // because this is how ImageJ handles the voxel locations
    double m_R0 = item.m_R0();
    double subCenter = (2 * m_R0 + 1) * (10.0 / 2.0);
    double subpixelCoordinate = 0.0;

    subpixelCoordinate = subCenter + (coordinate - (double) (item.subpix_image_center.get(axisFlag))) * 10.0;

    return subpixelCoordinate;
  }

  /**
   * Function to load an image
   * 
   * @param title title to be displayed on ImageJ image window
   * @param type  type of file to be opened - will be displayed on the chooser
   *              window
   * @return the unique title of the image
   * @throws IOException if file is not a valid type
   */
  private static String loadImages(String title, String type) throws IOException {

    String retval;

    if (WindowManager.getImage(title) != null)
      WindowManager.getImage(title).close();

    JFileChooser chooser;
    File setPathFile = new File("pth.txt");
    if (setPathFile.exists()) {
      try (Scanner scnr = new Scanner(setPathFile)) {
        String set_path = scnr.nextLine();
        chooser = new JFileChooser(set_path);
      } catch (Exception exc) {
        logger.addWarning(exc.toString());
        chooser = new JFileChooser(FileSystemView.getFileSystemView().getHomeDirectory());
      }
    } else {
      logger.addVariable("home", FileSystemView.getFileSystemView().getHomeDirectory());
      chooser = new JFileChooser(FileSystemView.getFileSystemView().getHomeDirectory());
    }
    chooser.setDialogTitle("Open " + type);
    int status = chooser.showSaveDialog(null);

    Opener fOpener = new Opener();

    // if user does not choose a file
    if (status != JFileChooser.APPROVE_OPTION)
      return "";

    String img_path = chooser.getSelectedFile().getAbsolutePath();
    String filename = chooser.getSelectedFile().getName();
    retval = WindowManager.makeUniqueName(filename);
    String format = img_path.substring(img_path.lastIndexOf(".") + 1);

    if (format.compareTo(ACCEPTED_FILE_TYPE) != 0)
      throw new IOException("File type not accepted");

    fOpener.open(img_path);
    WindowManager.getCurrentImage().setTitle(retval);

    return retval;
  }

  private static String loadImagesIJ(String title) {
    if (WindowManager.getImage(title) != null)
      WindowManager.getImage(title).close();

    IJ.open();
    title = WindowManager.getCurrentImage().getTitle();
    title = WindowManager.makeUniqueName(title);
    WindowManager.getCurrentImage().setTitle(title);

    return title;
  }
}