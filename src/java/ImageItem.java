import java.awt.Rectangle;
import javax.swing.JOptionPane;
import ij.ImagePlus;
import ij.gui.*;
import ij.WindowManager;

public class ImageItem {

  private ImagePlus mag_img, phase_img;
  private Vec3<Double> center_l;
  private Vec3<Double> center_m;
  private Vec3<Double> center_s;
  public Vec3<Double> subImageCenter;
  public double bkgPhase, estimatedBkgPhase;
  private double roi_mag_belowM_sumX, roi_mag_belowM_sumY, roi_mag_belowM_sumZ;
  public int roi_xi, roi_yi, roi_zi, roi_dx, roi_dy, roi_dz;
  public int roi_mag_belowM_xi, roi_mag_belowM_yi, roi_mag_belowM_zi,
      roi_mag_belowM_dx, roi_mag_belowM_dy, roi_mag_belowM_dz;
  public int M;
  private double RCenter;
  private double phaseValue;
  private Axis mriAxis = Axis.UNKNOWN;
  private double m_R0;
  private double estMagMoment;
  public double R1PhaseCalc, R2PhaseCalc, R3PhaseCalc;
  private double[][][] phase_nobkg;

  public int echoImageIndex;

  private final int grid = 10;
  private final double m_ROuterFrom = 0.2;
  private final double m_RMiddleFrom = 0.9;
  private final double m_RInnerFrom = 2.5;

  public boolean isNearEdge = false;

  public boolean imagesOpen(String mag_title, String phase_title) {
    return WindowManager.getImage(mag_title) != null && WindowManager.getImage(phase_title) != null;
  }

  public ImageItem(String magTitle, String phaseTitle, int M_pct, double pV) {
    M = M_pct;
    phaseValue = pV;

    // If mag and phase images are not open
    if (!imagesOpen(magTitle, phaseTitle)) {
      throw new IllegalStateException("No magnitude or phase file open.");
    }

    // Getting image instances
    try {
      mag_img = WindowManager.getImage(magTitle);
      phase_img = WindowManager.getImage(phaseTitle);

      echoImageIndex = (mag_img.getRoi() != null) ? mag_img.getFrame() : phase_img.getFrame();
      Calculate_Magnetic_Moment_3D.logger.addVariable("echo image index", echoImageIndex);

      if (mag_img.getRoi() == null)
        mag_img.setT(echoImageIndex);
      else
        phase_img.setT(echoImageIndex);

    } catch (Exception exc) {
      JOptionPane.showMessageDialog(Calculate_Magnetic_Moment_3D.gui.frame, exc);
    }

    // If no ROI is set
    if (mag_img.getRoi() == null && phase_img.getRoi() == null) {
      throw new IllegalStateException("No ROI found");
    }

    // getting ROI
    Roi user_roi = (mag_img.getRoi() != null) ? mag_img.getRoi() : phase_img.getRoi();

    // Putting ROI into rectangle
    Rectangle user_roi_rectangle = user_roi.getBounds();

    // Getting ROI parameters
    roi_xi = (int) user_roi.getXBase();
    roi_yi = (int) user_roi.getYBase();
    roi_zi = (mag_img.getRoi() != null) ? mag_img.getSlice() - 1 : phase_img.getSlice() - 1;
    // square box dimensions are the larger of roi_dx and roi_dy
    roi_dx = roi_dy = roi_dz = Math.max((int) user_roi_rectangle.getWidth(),
        (int) user_roi_rectangle.getHeight());

    roi_zi = roi_zi - roi_dz / 2;

    Calculate_Magnetic_Moment_3D.logger.addVariable("roi_d", new Vec3<>(roi_dx, roi_dy, roi_dz).toString());

    // get current slice from whatever image ROI is on

    Calculate_Magnetic_Moment_3D.logger.addVariable("roi before",
        new Vec3<>(roi_xi, roi_yi, roi_zi).toString() + ',' + new Vec3<>(roi_dx, roi_dy, roi_dz));
    // making sure roi parameters fit onto image
    roi_xi = goodROI(mag_img.getWidth(), roi_xi, roi_dx).get(0);
    roi_dx = goodROI(mag_img.getWidth(), roi_xi, roi_dx).get(1);
    roi_yi = goodROI(mag_img.getHeight(), roi_yi, roi_dy).get(0);
    roi_dy = goodROI(mag_img.getHeight(), roi_yi, roi_dy).get(1);
    // saving old initial z position because of how goodROI(int,int,int) is
    // written...
    // having a negative initial z (which is only possible in the slice direction)
    // has a direct correlation with how roi_dz is calculated
    int old_roi_zi = roi_zi;
    roi_zi = goodROI(mag_img.getNSlices(), roi_zi, roi_dz).get(0);
    roi_dz = goodROI(mag_img.getNSlices(), old_roi_zi, roi_dz).get(1);

    // Setting new ROI parameters to rectangle object
    user_roi_rectangle.setBounds(roi_xi, roi_yi, roi_dx,
        roi_dy);

    // Setting new ROI/rectangle to image
    if (mag_img.getRoi() != null)
      mag_img.setRoi(user_roi_rectangle);
    else
      phase_img.setRoi(user_roi_rectangle);

    Calculate_Magnetic_Moment_3D.logger.addVariable("roi after",
        new Vec3<>(roi_xi, roi_yi, roi_zi).toString() + ',' + new Vec3<>(roi_dx, roi_dy, roi_dz));

    // ----------------------------------------------------------------------------
    // Finding average mag intensity of box corners (2b of step2.md)
    // ----------------------------------------------------------------------------

    // setting slice to starting ROI point in the z
    int currFrame = mag_img.getFrame();
    Calculate_Magnetic_Moment_3D.logger.addVariable("curr frame", currFrame);

    Vec3<Integer> corner = new Vec3<Integer>(roi_xi, roi_yi, roi_zi);
    Vec3<Integer> roiSizeNoEdge = new Vec3<Integer>(roi_xi + roi_dx - 1, roi_yi + roi_dy - 1, roi_zi + roi_dz - 1);
    double avgCorners = ImageMethods.sumBoxCorners(mag_img, corner, roiSizeNoEdge, echoImageIndex) / 8.0;

    Calculate_Magnetic_Moment_3D.logger.addVariable("avgCorners", avgCorners);

    // ----------------------------------------------------------------------------
    // Finding region that contains values below M% (2d of step2.md)
    // ----------------------------------------------------------------------------

    // 3d array mapping the inner box that will be made
    double[][][] roi_mag_belowM = new double[roi_dx][roi_dy][roi_dz];
    // flooding roi_mag_belowM with -1's so valid points can
    // be easily distingushed
    for (int i = 0; i < roi_dx; i++) {
      for (int j = 0; j < roi_dy; j++) {
        for (int k = 0; k < roi_dz; k++) {
          roi_mag_belowM[i][j][k] = -1.0;
        }
      }
    }

    // a new ROI will be defined with all the values in the user's ROI that are less
    // than this value:
    double newbox_max = avgCorners * (double) M / 100.0;
    Calculate_Magnetic_Moment_3D.logger.addVariable("M% val", newbox_max);

    // ----------------------------------------------------------------------------
    // Finding corners of new ROI (1/2 of 2e in step2.md)
    // ----------------------------------------------------------------------------
    roi_mag_belowM_xi = roi_xi + roi_dx;
    roi_mag_belowM_yi = roi_yi + roi_dy;
    roi_mag_belowM_zi = roi_zi + roi_dz;
    roi_mag_belowM_dx = 0;
    roi_mag_belowM_dy = 0;
    roi_mag_belowM_dz = 0;

    // getting all pixel values that are below value
    for (int dx = 0; dx < roi_dx; dx++) {
      for (int dy = 0; dy < roi_dy; dy++) {
        for (int dz = 0; dz < roi_dz; dz++) {
          double mag_intensity = ImageMethods.getVoxelValue(mag_img, roi_xi + dx, roi_yi + dy, roi_zi + dz,
              echoImageIndex);

          if (mag_intensity < newbox_max) // if below value
          {
            // store value in array
            roi_mag_belowM[dx][dy][dz] = mag_intensity;

            // getting smallest coordinate
            if (roi_xi + dx < roi_mag_belowM_xi) {
              roi_mag_belowM_xi = roi_xi + dx;
            }
            if (roi_yi + dy < roi_mag_belowM_yi) {
              roi_mag_belowM_yi = roi_yi + dy;
            }
            if (roi_zi + dz < roi_mag_belowM_zi) {
              roi_mag_belowM_zi = roi_zi + dz;
            }

            // getting largest delta
            if (dx + roi_xi - roi_mag_belowM_xi > roi_mag_belowM_dx) {
              roi_mag_belowM_dx = dx + roi_xi - roi_mag_belowM_xi;
            }
            if (dy + roi_yi - roi_mag_belowM_yi > roi_mag_belowM_dy) {
              roi_mag_belowM_dy = dy + roi_yi - roi_mag_belowM_yi;
            }
            if (dz + roi_zi - roi_mag_belowM_zi > roi_mag_belowM_dz) {
              roi_mag_belowM_dz = dz + roi_zi - roi_mag_belowM_zi;
            }
          }
        }
      }
    }

    Calculate_Magnetic_Moment_3D.logger.addVariable("roi xi yi zi",
        new Vec3<Integer>(roi_mag_belowM_xi, roi_mag_belowM_yi, roi_mag_belowM_zi).toString());
    Calculate_Magnetic_Moment_3D.logger.addVariable("roi dx dy dz",
        new Vec3<Integer>(roi_mag_belowM_dx, roi_mag_belowM_dy, roi_mag_belowM_dz).toString());

    // Remainder of step2 implementation is in other functions in this class

  }

  // =====================================================================================
  // Function to estimate background phase by averaging corners of user-drawn ROI
  // This function also initializes the array in which these removed-bkg phase
  // values are stored in
  // =====================================================================================
  public double estimateBkgPhase() {

    // Initialize 3d array for phase values w/ removed background phase
    phase_nobkg = new double[roi_dx + 1][roi_dy + 1][roi_dz + 1];
    // Populate array with corresponding values
    for (int iz = 0; iz <= roi_dz; iz++) {
      for (int iy = 0; iy <= roi_dy; iy++) {
        for (int ix = 0; ix <= roi_dx; ix++) {
          phase_nobkg[ix][iy][iz] = ImageMethods.getVoxelValue(phase_img, roi_xi + ix, roi_yi + iy, roi_zi + iz,
              echoImageIndex);
        }
      }
    }

    Vec3<Integer> corner1 = new Vec3<Integer>(roi_xi, roi_yi, roi_zi);
    Vec3<Integer> corner2 = new Vec3<Integer>(roi_xi + roi_dx - 1, roi_yi + roi_dy - 1, roi_zi + roi_dz - 1);
    bkgPhase = ImageMethods.sumBoxCorners(phase_img, corner1, corner2, echoImageIndex) / 8.0;

    estimatedBkgPhase = bkgPhase;
    Calculate_Magnetic_Moment_3D.logger.addVariable("estimated bkg phase", estimatedBkgPhase);

    return bkgPhase;

  }

  // =====================================================================================
  // Function to remove current background phase from phase_nobkg array
  // =====================================================================================
  public void removeBkgPhase() {
    for (int iz = 0; iz <= roi_dz; iz++) {
      for (int iy = 0; iy <= roi_dy; iy++) {
        for (int ix = 0; ix <= roi_dx; ix++) {
          phase_nobkg[ix][iy][iz] -= bkgPhase;
        }
      }
    }
  }

  // =====================================================================================
  // Function to access phase values with removed background
  // =====================================================================================
  public double phaseNoBkg(int i, int j, int k) {
    int x = i - roi_xi;
    int y = j - roi_yi;
    int z = k - roi_zi;

    if (x < 0 || x > roi_dx)
      return Double.NaN;
    if (y < 0 || y > roi_dy)
      return Double.NaN;
    if (z < 0 || z > roi_dz)
      return Double.NaN;

    return phase_nobkg[x][y][z];
  }

  public double estimateRCenter() {
    int csx = center_s.get(0).intValue();
    int csy = center_s.get(1).intValue();
    int csz = center_s.get(2).intValue();
    Calculate_Magnetic_Moment_3D.logger.addVariable("csx", csx);
    Calculate_Magnetic_Moment_3D.logger.addVariable("csy", csy);
    Calculate_Magnetic_Moment_3D.logger.addVariable("csz", csz);

    // Phase values for each direction
    int xp_size = roi_xi + roi_dx - csx + 1;
    int yp_size = roi_yi + roi_dy - csy + 1;
    int zp_size = roi_zi + roi_dz - csz + 1;
    int xn_size = csx - roi_xi + 1;
    int yn_size = csy - roi_yi + 1;
    int zn_size = csz - roi_zi + 1;
    double[] phaseVals_xPos = new double[xp_size];
    double[] phaseVals_yPos = new double[yp_size];
    double[] phaseVals_zPos = new double[zp_size];
    double[] phaseVals_xNeg = new double[xn_size];
    double[] phaseVals_yNeg = new double[yn_size];
    double[] phaseVals_zNeg = new double[zn_size];

    // putting phase values into array, following respectively for other 5 loops
    for (int i = 0; i < xp_size; i++) {
      // phaseVals_xPos[i] = Math.abs((double)
      // getImageProcessor(phase_img).getPixelValue(i + csx,
      // csy));
      phaseVals_xPos[i] = Math.abs(phaseNoBkg(csx + i, csy, csz));
    }
    for (int i = 0; i < xn_size; i++) {
      // phaseVals_xNeg[i] = Math.abs((double)
      // getImageProcessor(phase_img).getPixelValue(csx - i,
      // csy));
      phaseVals_xNeg[i] = Math.abs(phaseNoBkg(csx - i, csy, csz));
    }
    for (int j = 0; j < yp_size; j++) {
      // phaseVals_yPos[j] = Math.abs((double)
      // getImageProcessor(phase_img).getPixelValue(csx, j +
      // csy));
      phaseVals_yPos[j] = Math.abs(phaseNoBkg(csx, csy + j, csz));
    }
    for (int j = 0; j < yn_size; j++) {
      // phaseVals_yNeg[j] = Math.abs((double)
      // getImageProcessor(phase_img).getPixelValue(csx, csy
      // - j));
      phaseVals_yNeg[j] = Math.abs(phaseNoBkg(csx, csy - j, csz));
    }
    for (int k = 0; k < zp_size; k++) {
      // phase_img.setSlice(k + csz + 1);
      // phaseVals_zPos[k] = Math.abs((double)
      // getImageProcessor(phase_img).getPixelValue(csx,
      // csy));
      phaseVals_zPos[k] = Math.abs(phaseNoBkg(csx, csy, csz + k));
    }
    for (int k = 0; k < zn_size; k++) {
      // phase_img.setSlice(csz - k + 1);
      // phaseVals_zNeg[k] = Math.abs((double)
      // getImageProcessor(phase_img).getPixelValue(csx,
      // csy));
      phaseVals_zNeg[k] = Math.abs(phaseNoBkg(csx, csy, csz - k));
    }

    // negating all the values that are below the threshold (because phase values
    // can be random within the object)
    Vec3<Integer> roiCorner1 = new Vec3<Integer>(roi_xi, roi_yi, roi_zi);
    Vec3<Integer> roiCorner2 = new Vec3<Integer>(roi_xi + roi_dx - 1, roi_yi + roi_dy - 1, roi_zi + roi_dz - 1);
    ImageMethods.removeValuesBelow(phaseVals_xPos, mag_img, center_s, roiCorner1, roiCorner2, echoImageIndex, M, Axis.X,
        true);
    ImageMethods.removeValuesBelow(phaseVals_xNeg, mag_img, center_s, roiCorner1, roiCorner2, echoImageIndex, M, Axis.X,
        false);
    ImageMethods.removeValuesBelow(phaseVals_yPos, mag_img, center_s, roiCorner1, roiCorner2, echoImageIndex, M, Axis.Y,
        true);
    ImageMethods.removeValuesBelow(phaseVals_yNeg, mag_img, center_s, roiCorner1, roiCorner2, echoImageIndex, M, Axis.Y,
        false);
    ImageMethods.removeValuesBelow(phaseVals_zPos, mag_img, center_s, roiCorner1, roiCorner2, echoImageIndex, M, Axis.Z,
        true);
    ImageMethods.removeValuesBelow(phaseVals_zNeg, mag_img, center_s, roiCorner1, roiCorner2, echoImageIndex, M, Axis.Z,
        false);
    Calculate_Magnetic_Moment_3D.logger.addInfo("negated");
    Calculate_Magnetic_Moment_3D.logger.addVariable("roiCorner1 in M% removal", roiCorner1.toString());
    Calculate_Magnetic_Moment_3D.logger.addVariable("roiCorner2 in M% removal", roiCorner2.toString());

    // Putting equitorial plane and MRI field direction into new variables (since we
    // technically don't know which direction is going to be the MRI field
    // direction)
    // neglected will be MRI and PV will be equitorials
    double[] phaseVals_pos;
    double[] phaseVals_neg;
    int center_mri_axis;
    // String neglectedAxis;
    mriAxis = ImageMethods.calculateMRIAxis(phaseVals_xPos, phaseVals_xNeg, phaseVals_yPos, phaseVals_yNeg,
        phaseVals_zPos,
        phaseVals_zNeg);
    Calculate_Magnetic_Moment_3D.gui.ldd_MriAxis.setValue(mriAxis);
    switch (mriAxis) {
      case X:
        phaseVals_pos = phaseVals_xPos;
        phaseVals_neg = phaseVals_xNeg;
        center_mri_axis = csx;
        break;

      case Y:
        phaseVals_pos = phaseVals_yPos;
        phaseVals_neg = phaseVals_yNeg;
        center_mri_axis = csy;
        break;

      case Z:
        phaseVals_pos = phaseVals_zPos;
        phaseVals_neg = phaseVals_zNeg;
        center_mri_axis = csz;
        break;

      default:
        // TODO: throw warning
        phaseVals_pos = phaseVals_zPos;
        phaseVals_neg = phaseVals_zNeg;
        center_mri_axis = csz;
        break;
    }

    Calculate_Magnetic_Moment_3D.logger.addVariable("mriAxis", mriAxis);

    // calculate RCenter and new center
    double[] RCenterAndNewCenter = ImageMethods.calculateRCenter(center_mri_axis, mriAxis, phaseVals_pos, phaseVals_neg,
        phaseValue);
    RCenter = RCenterAndNewCenter[0];
    double newCenter = RCenterAndNewCenter[1];
    double newCenter1f = Math.round(newCenter * 10.0) / 10.0; // new center rounded to one decimal place
    Calculate_Magnetic_Moment_3D.logger.addVariable("newCenter rounded", newCenter1f);

    // Set the center_s MRI axis value to the average of the found radii + the
    // found center point of the MRI field axis
    switch (mriAxis) {
      case X:
        center_s.set(0, newCenter1f);
        Calculate_Magnetic_Moment_3D.gui.ltf_rcx.setValue(String.valueOf(newCenter1f));
        break;
      case Y:
        center_s.set(1, newCenter1f);
        Calculate_Magnetic_Moment_3D.gui.ltf_rcy.setValue(String.valueOf(newCenter1f));
        break;
      case Z:
        center_s.set(2, newCenter1f);
        Calculate_Magnetic_Moment_3D.gui.ltf_rcz.setValue(String.valueOf(newCenter1f));
        break;
      default:
        Calculate_Magnetic_Moment_3D.logger.addInfo("Unable to determine MRI axis");
    }

    return RCenter;
  }

  public void calcR0123() {
    double factor = Math.pow(RCenter, 3.0) * phaseValue;
    double Practical_Phase_Limit = 0.02;
    double exp = 1.0 / 3.0;
    double suggR1, suggR2, suggR3;
    double base;

    base = factor / m_RInnerFrom;
    suggR3 = Math.pow(base, exp);
    suggR3 = suggR3 * (double) grid;
    suggR3 = Math.ceil(suggR3) / (double) grid;
    double m_R3 = suggR3;

    base = factor / m_RMiddleFrom;
    suggR2 = Math.pow(base, exp);
    suggR2 = suggR2 * (double) grid;
    suggR2 = Math.ceil(suggR2) / (double) grid;
    double m_R2 = suggR2;

    base = factor / m_ROuterFrom;
    suggR1 = Math.pow(base, exp);
    suggR1 = suggR1 * (double) grid;
    suggR1 = Math.ceil(suggR1) / (double) grid;
    double m_R1 = suggR1;

    Calculate_Magnetic_Moment_3D.gui.ltf_r1.setValue(String.valueOf(m_R1));
    Calculate_Magnetic_Moment_3D.gui.ltf_r2.setValue(String.valueOf(m_R2));
    Calculate_Magnetic_Moment_3D.gui.ltf_r3.setValue(String.valueOf(m_R3));

    /*
     * Finding m_R0. It is determined by the calculations below.
     *
     * m_R0 is an estimation of the p-value.
     *
     * m_R0 is also going to be the distance from the center that will be generated
     * in the images. So 2 * m_R0 + 1 is the size of the new images
     *
     * If the distance between the center and the edge of the image is too close
     * (less than 4), send an error message to the GUI
     *
     * If the distance is between 4-8, then that distance will be m_R0
     * Other than that, m_R0 must be between 8 and 21
     */
    m_R0 = Math.pow((factor / Practical_Phase_Limit), exp);

    if (m_R0 > 21.0)
      m_R0 = 21.0;

    if (m_R0 < 8.0)
      m_R0 = 8.0;

    int xDistanceFromEdge1 = mag_img.getWidth() - 1 - center_s.get(0).intValue();
    int xDistanceFromEdge2 = center_s.get(0).intValue();
    int yDistanceFromEdge1 = mag_img.getHeight() - 1 - center_s.get(1).intValue();
    int yDistanceFromEdge2 = center_s.get(1).intValue();
    int zDistanceFromEdge1 = mag_img.getNSlices() - 1 - center_s.get(2).intValue();
    int zDistanceFromEdge2 = center_s.get(2).intValue();

    isNearEdge = true;

    if ((xDistanceFromEdge1 >= 4) && (xDistanceFromEdge1 < 8)) {
      if (xDistanceFromEdge1 < m_R0) {
        m_R0 = xDistanceFromEdge1;
      }
    } else if (xDistanceFromEdge1 < 4) {
      isNearEdge = false;
    }

    if ((xDistanceFromEdge2 >= 4) && (xDistanceFromEdge2 < 8)) {
      if (xDistanceFromEdge2 < m_R0) {
        m_R0 = xDistanceFromEdge2;
      }
    } else if (xDistanceFromEdge2 < 4) {
      isNearEdge = false;
    }

    if ((yDistanceFromEdge1 >= 4) && (yDistanceFromEdge1 < 8)) {
      if (yDistanceFromEdge1 < m_R0) {
        m_R0 = yDistanceFromEdge1;
      }
    } else if (yDistanceFromEdge1 < 4) {
      isNearEdge = false;
    }

    if ((yDistanceFromEdge2 >= 4) && (yDistanceFromEdge2 < 8)) {
      if (yDistanceFromEdge2 < m_R0) {
        m_R0 = yDistanceFromEdge2;
      }
    } else if (yDistanceFromEdge2 < 4) {
      isNearEdge = false;
    }

    if ((zDistanceFromEdge1 >= 4) && (zDistanceFromEdge1 < 8)) {
      if (zDistanceFromEdge1 < m_R0) {
        m_R0 = zDistanceFromEdge1;
      }
    } else if (zDistanceFromEdge1 < 4) {
      isNearEdge = false;
    }

    if ((zDistanceFromEdge2 >= 4) && (zDistanceFromEdge2 < 8)) {
      if (zDistanceFromEdge2 < m_R0) {
        m_R0 = zDistanceFromEdge2;
      }
    } else if (zDistanceFromEdge2 < 4) {
      isNearEdge = false;
    }

    Calculate_Magnetic_Moment_3D.logger.addVariable("xdistedge1", xDistanceFromEdge1);
    Calculate_Magnetic_Moment_3D.logger.addVariable("xdistedge2", xDistanceFromEdge2);
    Calculate_Magnetic_Moment_3D.logger.addVariable("ydistedge1", yDistanceFromEdge1);
    Calculate_Magnetic_Moment_3D.logger.addVariable("ydistedge2", yDistanceFromEdge2);
    Calculate_Magnetic_Moment_3D.logger.addVariable("zdistedge1", zDistanceFromEdge1);
    Calculate_Magnetic_Moment_3D.logger.addVariable("zdistedge2", zDistanceFromEdge2);

    m_R0 = Math.ceil(m_R0);

    Calculate_Magnetic_Moment_3D.logger.addVariable("m_R0", m_R0);

    estMagMoment = phaseValue * Math.pow(RCenter, 3);

    R1PhaseCalc = estMagMoment / Math.pow(m_R1, 3);
    R2PhaseCalc = estMagMoment / Math.pow(m_R2, 3);
    R3PhaseCalc = estMagMoment / Math.pow(m_R3, 3);
    Calculate_Magnetic_Moment_3D.gui.lbl_r1phaseCalc.setText(String.valueOf(Math.round(R1PhaseCalc * 100.0) / 100.0));
    Calculate_Magnetic_Moment_3D.gui.lbl_r2phaseCalc.setText(String.valueOf(Math.round(R2PhaseCalc * 100.0) / 100.0));
    Calculate_Magnetic_Moment_3D.gui.lbl_r3phaseCalc.setText(String.valueOf(Math.round(R3PhaseCalc * 100.0) / 100.0));

    // setMagMoment(phaseValue * Math.pow(RCenter, 3));
  }

  // =====================================================================================
  // Function to calculate center_l, i.e. the coordinate of the minimum mag value
  // =====================================================================================
  public void calcCenterL() {

    float min = Float.POSITIVE_INFINITY;
    center_l = new Vec3<Double>(0.0, 0.0, 0.0);

    for (int i = roi_xi; i < roi_xi + roi_dx; i++) {
      for (int j = roi_yi; j < roi_yi + roi_dy; j++) {
        for (int k = roi_zi; k < roi_zi + roi_dz; k++) {
          if (ImageMethods.getVoxelValue(mag_img, i, j, k, echoImageIndex) < min) {
            min = ImageMethods.getVoxelValue(mag_img, i, j, k, echoImageIndex);
            center_l.set(0, (double) i);
            center_l.set(1, (double) j);
            center_l.set(2, (double) k);
          }
        }
      }
    }

    Calculate_Magnetic_Moment_3D.logger.addVariable("center_l", center_l);
  }

  // =====================================================================================
  // Function to calculate center_m, i.e. center of small ROI
  // =====================================================================================
  public void calcCenterM() {

    // ----------------------------------------------------------------------------
    // Finding center_m (i.e., center of M% region) (2/2 of 2e of step2.md)
    // ----------------------------------------------------------------------------

    center_m = new Vec3<Double>(0.0, 0.0, 0.0);
    int x1 = roi_mag_belowM_xi;
    int y1 = roi_mag_belowM_yi;
    int z1 = roi_mag_belowM_zi;
    int x2 = roi_mag_belowM_xi + roi_mag_belowM_dx;
    int y2 = roi_mag_belowM_yi + roi_mag_belowM_dy;
    int z2 = roi_mag_belowM_zi + roi_mag_belowM_dz;
    Calculate_Magnetic_Moment_3D.logger.addVariable("CM xyz1",
        String.valueOf(x1) + ',' + String.valueOf(y1) + ',' + String.valueOf(z1));
    Calculate_Magnetic_Moment_3D.logger.addVariable("CM xyz2",
        String.valueOf(x2) + ',' + String.valueOf(y2) + ',' + String.valueOf(z2));

    /*
     * if ((x2 - x1) % 2 == 0) {
     * center_m.set(0, (double) (x1 + (x2 - x1) / 2));
     * } else {
     * center_m.set(0, (double) (x1 + (x2 - x1) / 2 + 0.5));
     * }
     *
     * if ((y2 - y1) % 2 == 0) {
     * center_m.set(1, (double) (y1 + (y2 - y1) / 2));
     * } else {
     * center_m.set(1, (double) (y1 + (y2 - y1) / 2 + 0.5));
     * }
     *
     * if ((z2 - z1) % 2 == 0) {
     * center_m.set(2, (double) (z1 + (z2 - z1) / 2));
     * } else {
     * center_m.set(2, (double) (z1 + (z2 - z1) / 2 + 0.5));
     * }
     */

    // So instead of the above calculation, we can simply do the following and get
    // the same result. No matter if the difference is even or odd, the center of
    // the range can be calculated by adding half the difference to the initial
    // coordinate
    // Also note that (and assume x1 and x2 are casted as doubles):
    // x1 + (x2-x1)/2.0 = 0.5*(x1 + x2)
    center_m.set(0, 0.5 * (double) (x1 + x2));
    center_m.set(1, 0.5 * (double) (y1 + y2));
    center_m.set(2, 0.5 * (double) (z1 + z2));

    Calculate_Magnetic_Moment_3D.logger.addVariable("center_m", center_m);
  }

  // =====================================================================================
  // Function to calculate center_s, i.e. coordinate of minimum sum of each plane
  // =====================================================================================
  public void calcCenterS() {

    center_s = new Vec3<Double>(-1.0, -1.0, -1.0);

    // ------ Begin summing planes in innerbox
    // Initializing arrays for sums of planes in small box
    double[] innerBox_sumOfXPlane = new double[roi_mag_belowM_dx + 1];
    double[] innerBox_sumOfYPlane = new double[roi_mag_belowM_dy + 1];
    double[] innerBox_sumOfZPlane = new double[roi_mag_belowM_dz + 1];

    // Summing z plane and putting it into array
    for (int k = 0; k <= roi_mag_belowM_dz; k++) {
      innerBox_sumOfZPlane[k] = 0;
      for (int i = 0; i <= roi_mag_belowM_dx; i++) {
        for (int j = 0; j <= roi_mag_belowM_dy; j++) {
          innerBox_sumOfZPlane[k] += ImageMethods.getVoxelValue(mag_img, i + roi_mag_belowM_xi, j + roi_mag_belowM_yi,
              k + roi_mag_belowM_zi, echoImageIndex);
        }
      }
    }

    // Summing x plane and putting it into array
    for (int i = 0; i <= roi_mag_belowM_dx; i++) {
      innerBox_sumOfXPlane[i] = 0;
      for (int j = 0; j <= roi_mag_belowM_dy; j++) {
        for (int k = 0; k <= roi_mag_belowM_dz; k++) {
          innerBox_sumOfXPlane[i] += ImageMethods.getVoxelValue(mag_img, i + roi_mag_belowM_xi, j + roi_mag_belowM_yi,
              k + roi_mag_belowM_zi, echoImageIndex);
        }
      }
    }

    // Summing y plane and putting it into array
    for (int j = 0; j <= roi_mag_belowM_dy; j++) {
      innerBox_sumOfYPlane[j] = 0;
      for (int k = 0; k <= roi_mag_belowM_dz; k++) {
        for (int i = 0; i <= roi_mag_belowM_dx; i++) {
          innerBox_sumOfYPlane[j] += ImageMethods.getVoxelValue(mag_img, i + roi_mag_belowM_xi, j + roi_mag_belowM_yi,
              k + roi_mag_belowM_zi, echoImageIndex);
        }
      }
    }

    // for finding the minimum summed plane
    roi_mag_belowM_sumX = innerBox_sumOfXPlane[0];
    roi_mag_belowM_sumY = innerBox_sumOfYPlane[0];
    roi_mag_belowM_sumZ = innerBox_sumOfZPlane[0];

    // ------ End summing planes in innerbox

    // ---------- begin to find Center_S

    // Finding minimum x plane and putting it into center_s.get(0)
    for (int i = 0; i <= roi_mag_belowM_dx; i++) {
      if (innerBox_sumOfXPlane[i] < roi_mag_belowM_sumX) {
        roi_mag_belowM_sumX = innerBox_sumOfXPlane[i];
        center_s.set(0, (double) (i + roi_mag_belowM_xi));
      }
    }

    // Finding minimum y plane and putting it into center_s.get(1)
    for (int j = 0; j <= roi_mag_belowM_dy; j++) {
      if (innerBox_sumOfYPlane[j] < roi_mag_belowM_sumY) {
        roi_mag_belowM_sumY = innerBox_sumOfYPlane[j];
        center_s.set(1, (double) (j + roi_mag_belowM_yi));
      }
    }

    // Finding minimum z plane and putting it into center_s.get(2)
    for (int k = 0; k <= roi_mag_belowM_dz; k++) {
      if (innerBox_sumOfZPlane[k] < roi_mag_belowM_sumZ) {
        roi_mag_belowM_sumZ = innerBox_sumOfZPlane[k];
        center_s.set(2, (double) (k + roi_mag_belowM_zi));
      }
    }
    // ---------- end to find Center_S
    Calculate_Magnetic_Moment_3D.logger.addVariable("Center S", center_s.toString());
    Calculate_Magnetic_Moment_3D.logger
        .addVariable("Mag @ Center S", ImageMethods.getVoxelValue(mag_img,
            center_s.get(0).intValue(),
            center_s.get(1).intValue(),
            center_s.get(2).intValue(),
            echoImageIndex));
    Calculate_Magnetic_Moment_3D.logger
        .addVariable("Phase @ Center S", ImageMethods.getVoxelValue(phase_img,
            center_s.get(0).intValue(),
            center_s.get(1).intValue(),
            center_s.get(2).intValue(),
            echoImageIndex));

    // if a component can't be found throw error
    if (center_s.get(0) == -1.0 || center_s.get(1) == -1.0 || center_s.get(2) == -1.0)
      throw new ArithmeticException("Unable to determine center_s completely");

  }

  // =====================================================================================
  // Getter for Center_L
  // =====================================================================================
  public Vec3<Double> centerL() {
    return center_l;
  }

  // =====================================================================================
  // Getter for Center_M
  // =====================================================================================
  public Vec3<Double> centerM() {
    return center_m;
  }

  // =====================================================================================
  // Getter for Center_S
  // =====================================================================================
  public Vec3<Double> centerS() {
    return center_s;
  }

  // =====================================================================================
  // Getter for MRI axis
  // =====================================================================================
  public Axis MRIAxis() {
    mriAxis = Calculate_Magnetic_Moment_3D.gui.ldd_MriAxis.getValue();
    Calculate_Magnetic_Moment_3D.logger.addVariable("mriAxis", mriAxis);
    return mriAxis;
  }

  // =====================================================================================
  // Getter for RCenter
  // =====================================================================================
  public double RCenter() {
    return RCenter;
  }

  // =====================================================================================
  // Getter for R0
  // =====================================================================================
  public Double m_R0() {
    return m_R0;
  }

  // =====================================================================================
  // Setter for Center_Sx
  // =====================================================================================
  public void setCenterSX(double x) {
    center_s.set(0, x);
  }

  // =====================================================================================
  // Setter for Center_Sy
  // =====================================================================================
  public void setCenterSY(double y) {
    center_s.set(1, y);
  }

  // =====================================================================================
  // Setter for Center_Sz
  // =====================================================================================
  public void setCenterSZ(double z) {
    center_s.set(2, z);
  }

  // =====================================================================================
  // Setter for RCenter
  // =====================================================================================
  public void setRCenter(double rc) {
    RCenter = rc;
  }

  /**
   * Ensures the roi size in the nth direction does not extend past the image
   * size in the nth direction
   *
   * @param img_size size of image in nth direction
   * @param point    initial coordinate of ROI in the nth direction
   * @param delta    size of ROI in nth direction
   * @return a Vec3 with contents being the fitted point, the fitted delta,
   *         then null
   */
  private Vec3<Integer> goodROI(int img_size, int point, int delta) {

    // if point + delta extends past image size, restrict delta to only go to the
    // image edge
    if (point + delta > img_size) {
      delta = img_size - point;
      // delta must be odd but our point must account for the change
      point = (delta % 2 == 0) ? (point - 1) : point;
    }

    // if point is negative (this is mainly for if a user picks an ROI near the 0th
    // slice)
    if (point < 0) {
      delta = delta + point;
      point = 0;
    }

    delta = (delta % 2 == 0) ? delta + 1 : delta;

    return new Vec3<Integer>(point, delta, null);
  }

  /**
   * Function to set slice of both mag and phase images
   *
   * @param i slice #
   */
  public void setSlice(int i) {
    mag_img.setSlice(i);
    phase_img.setSlice(i);
  }

  /**
   * Function to set Z of both mag and phase images
   *
   * @param z z to set
   */
  public void setZ(int z) {
    mag_img.setZ(z);
    phase_img.setZ(z);
  }
}