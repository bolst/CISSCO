import java.awt.Rectangle;
import javax.swing.JOptionPane;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.gui.*;
import ij.WindowManager;
import java.util.Arrays;

public class ImageItem {

  private ImagePlus mag_img, phase_img;
  private Triplet<Double> center_l;
  private Triplet<Double> center_m;
  private Triplet<Double> center_s;
  public Triplet<Double> subpix_image_center;
  public double bkgPhase, estimatedBkgPhase;
  private double roi_mag_belowM_sumX, roi_mag_belowM_sumY, roi_mag_belowM_sumZ;
  public int roi_xi, roi_yi, roi_zi, roi_dx, roi_dy, roi_dz;
  public int roi_mag_belowM_xi, roi_mag_belowM_yi, roi_mag_belowM_zi,
      roi_mag_belowM_dx, roi_mag_belowM_dy, roi_mag_belowM_dz;
  public int M;
  private double RCenter;
  private double phaseValue;
  private Axis MRI_axis = Axis.UNKNOWN;
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

  public static enum Axis {
    X, Y, Z, UNKNOWN;

    @Override
    public String toString() {
      switch (this) {
        case X:
          return "X";
        case Y:
          return "Y";
        case Z:
          return "Z";
        case UNKNOWN:
        default:
          return "?";
      }
    }
  }

  public boolean imagesOpen(String mag_title, String phase_title) {
    return WindowManager.getImage(mag_title) != null && WindowManager.getImage(phase_title) != null;
  }

  public Roi getROI(String mag_title, String phase_title) {
    if (!imagesOpen(mag_title, phase_title)) {
      return null;
    }

    ImagePlus mag_image;
    ImagePlus phase_image;
    try {
      mag_image = WindowManager.getImage(mag_title);
      phase_image = WindowManager.getImage(phase_title);
    } catch (Exception exc) {
      return null;
    }

    if (mag_image.getRoi() == null && phase_image.getRoi() == null) {
      return null;
    }

    return (mag_image.getRoi() != null) ? mag_image.getRoi() : phase_image.getRoi();
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

    Calculate_Magnetic_Moment_3D.logger.addVariable("roi_d", new Triplet<>(roi_dx, roi_dy, roi_dz).toString());

    // get current slice from whatever image ROI is on

    Calculate_Magnetic_Moment_3D.logger.addVariable("roi before",
        new Triplet<>(roi_xi, roi_yi, roi_zi).toString() + ',' + new Triplet<>(roi_dx, roi_dy, roi_dz));
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
        new Triplet<>(roi_xi, roi_yi, roi_zi).toString() + ',' + new Triplet<>(roi_dx, roi_dy, roi_dz));

    // ----------------------------------------------------------------------------
    // Finding average mag intensity of box corners (2b of step2.md)
    // ----------------------------------------------------------------------------

    // setting slice to starting ROI point in the z
    int currFrame = mag_img.getFrame();
    Calculate_Magnetic_Moment_3D.logger.addVariable("curr frame", currFrame);

    // summing all corners of ROI box in this slice
    float avgCorners = ImageMethods.getVoxelValue(mag_img, roi_xi, roi_yi, roi_zi)
        + ImageMethods.getVoxelValue(mag_img, roi_xi + (roi_dx - 1), roi_yi, roi_zi)
        + ImageMethods.getVoxelValue(mag_img, roi_xi, roi_yi + (roi_dy - 1), roi_zi)
        + ImageMethods.getVoxelValue(mag_img, roi_xi + (roi_dx - 1),
            roi_xi + (roi_dy - 1), roi_zi);

    // adding these corners of the ROI box in this slice
    avgCorners += ImageMethods.getVoxelValue(mag_img, roi_xi, roi_yi, roi_zi + roi_dz - 1)
        + ImageMethods.getVoxelValue(mag_img, roi_xi + (roi_dx - 1), roi_yi, roi_zi + roi_dz - 1)
        + ImageMethods.getVoxelValue(mag_img, roi_xi, roi_yi + (roi_dy - 1), roi_zi + roi_dz - 1)
        + ImageMethods.getVoxelValue(mag_img, roi_xi + (roi_dx - 1),
            roi_xi + (roi_dy - 1), roi_zi + roi_dz - 1);

    // dividing by 8 for average
    avgCorners /= 8.0;
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
    for (int dz = 0; dz < roi_dz; dz++) {
      for (int dx = 0; dx < roi_dx; dx++) {
        for (int dy = 0; dy < roi_dy; dy++) {
          double mag_intensity = ImageMethods.getVoxelValue(mag_img, roi_xi + dx, roi_yi + dy, roi_zi + dz);

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
        new Triplet<Integer>(roi_mag_belowM_xi, roi_mag_belowM_yi, roi_mag_belowM_zi).toString());
    Calculate_Magnetic_Moment_3D.logger.addVariable("roi dx dy dz",
        new Triplet<Integer>(roi_mag_belowM_dx, roi_mag_belowM_dy, roi_mag_belowM_dz).toString());

    // Remainder of step2 implementation is in other functions in this class

  }

  /*
   * Function to find MRI field direction. Does this by comparing the values in
   * each direction from the estimated center.
   * The two that follow similar trends are going to be on the equitorial plane,
   * while the other will be the MRI field direction.
   * It will iterate along each direction and store the direction with the largest
   * difference between the other two. The most frequent direction
   * will be the MRI field direction
   *
   * @param accuracy The amount of times the code will iterate along each
   * direction
   *
   * @return A string value ("x", "y", or "z") of the MRI field direction
   */

  private Axis calculateMRIAxis(double[] XP,
      double[] XN,
      double[] YP,
      double[] YN,
      double[] ZP,
      double[] ZN) {
    Calculate_Magnetic_Moment_3D.logger.addVariable("XP", Arrays.toString(XP));
    Calculate_Magnetic_Moment_3D.logger.addVariable("XN", Arrays.toString(XN));
    Calculate_Magnetic_Moment_3D.logger.addVariable("YP", Arrays.toString(YP));
    Calculate_Magnetic_Moment_3D.logger.addVariable("YN", Arrays.toString(YN));
    Calculate_Magnetic_Moment_3D.logger.addVariable("ZP", Arrays.toString(ZP));
    Calculate_Magnetic_Moment_3D.logger.addVariable("ZN", Arrays.toString(ZN));

    int max_iter = Arrays.stream(new int[] { XP.length, XN.length, YP.length, YN.length, ZP.length, ZN.length }).min()
        .getAsInt();

    double X = 0.0;
    double Y = 0.0;
    double Z = 0.0;
    for (int i = 0; i < max_iter; i++) {
      if (XP[i] == 0.0 ||
          XN[i] == 0.0 ||
          YP[i] == 0.0 ||
          YN[i] == 0.0 ||
          ZP[i] == 0.0 ||
          ZN[i] == 0.0) {
        continue;
      }

      X += Math.abs(XP[i]) + Math.abs(XN[i]);
      Y += Math.abs(YP[i]) + Math.abs(YN[i]);
      Z += Math.abs(ZP[i]) + Math.abs(ZN[i]);
    }

    double ratioXY = Math.abs(X / Y - 1);
    double ratioXZ = Math.abs(X / Z - 1);
    double ratioYZ = Math.abs(Y / Z - 1);
    double minRatio = Math.min(Math.min(ratioXY, ratioXZ), ratioYZ);

    Axis retval = (minRatio == ratioXY) ? Axis.Z : (minRatio == ratioXZ ? Axis.Y : Axis.X);
    Calculate_Magnetic_Moment_3D.logger.addVariable("mri axis", retval);

    return retval;
  }

  /*
   * Function to void all values below the param threshold.
   *
   * @param values array of values in one direction
   *
   * @param threshold the % to negate, default is 50
   *
   * @param r0 the center coordinate of the direction inputted
   *
   * @param axis The direction
   *
   * @param direction Flag to specify if function should iterate in the + or -
   * direction
   */
  private void removeValuesBelow(double[] values,
      int threshold,
      int r0,
      Axis axis,
      boolean direction) {

    // getting average of ROI
    double avgOfCorners = ImageMethods.getVoxelValue(mag_img, roi_xi, roi_yi, roi_zi)
        + ImageMethods.getVoxelValue(mag_img, roi_xi + roi_dx, roi_yi, roi_zi)
        + ImageMethods.getVoxelValue(mag_img, roi_xi, roi_yi + roi_dy, roi_zi)
        + ImageMethods.getVoxelValue(mag_img, roi_xi + roi_dx, roi_yi + roi_dy, roi_zi);

    avgOfCorners += ImageMethods.getVoxelValue(mag_img, roi_xi, roi_yi, roi_zi + roi_dz)
        + ImageMethods.getVoxelValue(mag_img, roi_xi + roi_dx, roi_yi, roi_zi + roi_dz)
        + ImageMethods.getVoxelValue(mag_img, roi_xi, roi_yi + roi_dy, roi_zi + roi_dz)
        + ImageMethods.getVoxelValue(mag_img, roi_xi + roi_dx, roi_yi + roi_dy, roi_zi + roi_dz);

    avgOfCorners /= 8.0;

    // all values below this variable will be negated
    double maxMagValue = avgOfCorners * (double) threshold / 100.0;

    int intCsx = center_s.get(0).intValue();
    int intCsy = center_s.get(1).intValue();
    int intCsz = center_s.get(2).intValue();

    // setting all values in array below maxMagValue to 0
    switch (axis) {
      case X:
        if (direction) {
          for (int i = 0; i < values.length; i++) {
            if (Math.abs(
                ImageMethods.getVoxelValue(mag_img, i + r0, intCsy, intCsz)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        } else {
          for (int i = 0; i < values.length; i++) {
            if (Math.abs(ImageMethods.getVoxelValue(mag_img, r0 - i, intCsy, intCsz)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        }
      case Y:
        if (direction) {
          for (int i = 0; i < values.length; i++) {
            if (Math.abs(ImageMethods.getVoxelValue(mag_img, intCsx, i + r0, intCsz)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        } else {
          for (int i = 0; i < values.length; i++) {
            if (Math.abs(ImageMethods.getVoxelValue(mag_img, intCsx, r0 - i, intCsz)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        }
      case Z:
        if (direction) {
          for (int i = 0; i < values.length; i++) {
            if (Math.abs(ImageMethods.getVoxelValue(mag_img, intCsx, intCsy, i + r0)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        } else {
          for (int i = 0; i < values.length; i++) {
            if (Math.abs(ImageMethods.getVoxelValue(mag_img, intCsx, intCsy, r0 - i)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        }
      default:
        break;
    }
  }

  /*
   * Helper function used to estimate radius on MRI field axis
   * Uses the theory that on this axis, phase values follow a 2/r^3 trend
   *
   * @return averaged radius in both directions
   */
  private double calculateRC(int center_mri_axis,
      Axis mri_dir,
      double[] phaseVals_pos,
      double[] phaseVals_neg) {
    double r_pos = 1.6;
    double r_neg = 1.6;

    /*
     * The logic below works by finding a voxel phase value that is between the
     * user-inputted phase value.
     * So by default it is 1, so the program finds the voxel that is between two
     * voxels that have greater phase than 1 and less than 1.
     * The program then interpolates the distance between these two voxels to
     * estimate how far from the center a phase value of 1 is
     */

    for (int i = 0; i < phaseVals_pos.length - 1; i++) {
      // this condition is if the current voxel is nested between two voxels greater
      // and less than the user-defined phase value
      if ((phaseVals_pos[i] > phaseValue) && (phaseVals_pos[i + 1] < phaseValue)) {
        // this condition is if the interpolation value is greater than 1.6, because the
        // object radius has to be greater than 1.6
        if (interpolation(phaseValue, phaseVals_pos[i], phaseVals_pos[i + 1], i, i + 1) > 1.6) {
          r_pos = interpolation(phaseValue, phaseVals_pos[i], phaseVals_pos[i + 1], i, i + 1);
        }
        break;
      }
    }

    // Same thing as above just the opposite direction
    for (int i = 0; i < phaseVals_neg.length - 1; i++) {
      if ((phaseVals_neg[i] > phaseValue) && (phaseVals_neg[i + 1] < phaseValue)) {
        if (interpolation(phaseValue, phaseVals_neg[i], phaseVals_neg[i + 1], i, i + 1) > 1.6) {
          r_neg = interpolation(phaseValue, phaseVals_neg[i], phaseVals_neg[i + 1], i, i + 1);
        }
        break;
      }
    }

    // Set the center_s MRI axis value to the average of the found radii + the
    // found center point of the MRI field axis
    if (mri_dir == Axis.X) {
      center_s.set(0, ((center_mri_axis + r_pos) + (center_mri_axis - r_neg)) / 2.0);
      center_s.set(0, Math.round(center_s.get(0) * 10.0) / 10.0);
      Calculate_Magnetic_Moment_3D.gui.ltf_rcx.setValue(String.valueOf(Math.round(center_s.get(0) * 10.0) / 10.0));
    } else if (mri_dir == Axis.Y) {
      center_s.set(1, ((center_mri_axis + r_pos) + (center_mri_axis - r_neg)) / 2.0);
      center_s.set(1, Math.round(center_s.get(1) * 10.0) / 10.0);
      Calculate_Magnetic_Moment_3D.gui.ltf_rcy.setValue(String.valueOf(Math.round(center_s.get(1) * 10.0) / 10.0));
    } else if (mri_dir == Axis.Z) {
      center_s.set(2, ((center_mri_axis + r_pos) + (center_mri_axis - r_neg)) / 2.0);
      center_s.set(2, Math.round(center_s.get(2) * 10.0) / 10.0);
      Calculate_Magnetic_Moment_3D.gui.ltf_rcz.setValue(String.valueOf(Math.round(center_s.get(2) * 10.0) / 10.0));
    } else {
      Calculate_Magnetic_Moment_3D.logger.addInfo("Unable to determine MRI axis");
    }

    // returning equation for RCenter
    RCenter = ((r_pos + r_neg) / 2.0) / Math.cbrt(2);
    return RCenter;
  }

  /*
   * Standard method used for interpolation. Used to find distance to RCenter
   * Phase
   *
   * @param x Target point for desired y value - used as the RCenter Phase
   *
   * @param x1 First phase value that is greater than and adjacent to RCenter
   * Phase
   *
   * @param x2 Second phase value that is less than and adjacent to RCenter Phase
   *
   * @param y1 The first phase value's distance from RCenter
   *
   * @param y2 The second phase value's distance from RCenter
   *
   * @return The interpolated distance to what would be RCenter Phase
   */
  private double interpolation(double x, double x1, double x2, double y1, double y2) {
    return y1 + ((x - x1) / (x2 - x1)) * (y2 - y1);
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
          phase_nobkg[ix][iy][iz] = ImageMethods.getVoxelValue(phase_img, roi_xi + ix, roi_yi + iy, roi_zi);
        }
      }
    }

    // ---------- Begin to find estimated background phase
    bkgPhase = ImageMethods.getVoxelValue(phase_img, roi_xi, roi_yi, roi_zi)
        + ImageMethods.getVoxelValue(phase_img, roi_xi + roi_dx, roi_yi, roi_zi)
        + ImageMethods.getVoxelValue(phase_img, roi_xi, roi_yi + roi_dy, roi_zi)
        + ImageMethods.getVoxelValue(phase_img, roi_xi + roi_dx, roi_yi + roi_dy, roi_zi);

    bkgPhase += ImageMethods.getVoxelValue(phase_img, roi_xi, roi_yi, roi_zi + roi_dz)
        + ImageMethods.getVoxelValue(phase_img, roi_xi + roi_dx, roi_yi, roi_zi + roi_dz)
        + ImageMethods.getVoxelValue(phase_img, roi_xi, roi_yi + roi_dy, roi_zi + roi_dz)
        + ImageMethods.getVoxelValue(phase_img, roi_xi + roi_dx, roi_yi + roi_dy, roi_zi + roi_dz);

    bkgPhase /= 8.0;

    estimatedBkgPhase = bkgPhase;

    return bkgPhase;

    // ---------- end to find background phase
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
  public double phase_noBkg(int i, int j, int k) {
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
    int xn_size = csx - roi_xi + 1;
    int yp_size = roi_yi + roi_dy - csy + 1;
    int yn_size = csy - roi_yi + 1;
    int zp_size = roi_zi + roi_dz - csz + 1;
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
      phaseVals_xPos[i] = Math.abs(phase_noBkg(csx + i, csy, csz));
    }
    for (int i = 0; i < xn_size; i++) {
      // phaseVals_xNeg[i] = Math.abs((double)
      // getImageProcessor(phase_img).getPixelValue(csx - i,
      // csy));
      phaseVals_xNeg[i] = Math.abs(phase_noBkg(csx - i, csy, csz));
    }
    for (int j = 0; j < yp_size; j++) {
      // phaseVals_yPos[j] = Math.abs((double)
      // getImageProcessor(phase_img).getPixelValue(csx, j +
      // csy));
      phaseVals_yPos[j] = Math.abs(phase_noBkg(csx, csy + j, csz));
    }
    for (int j = 0; j < yn_size; j++) {
      // phaseVals_yNeg[j] = Math.abs((double)
      // getImageProcessor(phase_img).getPixelValue(csx, csy
      // - j));
      phaseVals_yNeg[j] = Math.abs(phase_noBkg(csx, csy - j, csz));
    }
    for (int k = 0; k < zp_size; k++) {
      // phase_img.setSlice(k + csz + 1);
      // phaseVals_zPos[k] = Math.abs((double)
      // getImageProcessor(phase_img).getPixelValue(csx,
      // csy));
      phaseVals_zPos[k] = Math.abs(phase_noBkg(csx, csy, csz + k));
    }
    for (int k = 0; k < zn_size; k++) {
      // phase_img.setSlice(csz - k + 1);
      // phaseVals_zNeg[k] = Math.abs((double)
      // getImageProcessor(phase_img).getPixelValue(csx,
      // csy));
      phaseVals_zNeg[k] = Math.abs(phase_noBkg(csx, csy, csz - k));
    }

    // Getting threshold (should be 50 unless user changes it)
    // M = M_pct;

    // negating all the values that are below the threshold (because phase values
    // can be random within the object)
    removeValuesBelow(phaseVals_xPos, M, csx, Axis.X, true);
    removeValuesBelow(phaseVals_xNeg, M, csx, Axis.X, false);
    removeValuesBelow(phaseVals_yPos, M, csy, Axis.Y, true);
    removeValuesBelow(phaseVals_yNeg, M, csy, Axis.Y, false);
    removeValuesBelow(phaseVals_zPos, M, csz, Axis.Z, true);
    removeValuesBelow(phaseVals_zNeg, M, csz, Axis.Z, false);
    Calculate_Magnetic_Moment_3D.logger.addInfo("negated");

    // Putting equitorial plane and MRI field direction into new variables (since we
    // technically don't know which direction is going to be the MRI field
    // direction)
    // neglected will be MRI and PV will be equitorials
    double[] phaseVals_pos;
    double[] phaseVals_neg;
    int center_mri_axis;
    // String neglectedAxis;
    MRI_axis = calculateMRIAxis(phaseVals_xPos, phaseVals_xNeg, phaseVals_yPos, phaseVals_yNeg, phaseVals_zPos,
        phaseVals_zNeg);
    Calculate_Magnetic_Moment_3D.gui.ldd_MriAxis.setValue(MRI_axis);
    switch (MRI_axis) {
      case X:
        phaseVals_pos = phaseVals_xPos;
        phaseVals_neg = phaseVals_xNeg;
        // neglectedAxis = "x";
        center_mri_axis = csx;
        // PVP1 = phaseVals_yPos;
        // PVN1 = phaseVals_yNeg;
        // axis1 = "y";
        // P1C = csy;
        // PVP2 = phaseVals_zPos;
        // PVN2 = phaseVals_zNeg;
        // axis2 = "z";
        // P2C = csz;
        break;

      case Y:
        phaseVals_pos = phaseVals_yPos;
        phaseVals_neg = phaseVals_yNeg;
        // neglectedAxis = "y";
        center_mri_axis = csy;
        // PVP1 = phaseVals_xPos;
        // PVN1 = phaseVals_xNeg;
        // P1C = csx;
        // axis1 = "x";
        // PVP2 = phaseVals_zPos;
        // PVN2 = phaseVals_zNeg;
        // P2C = csz;
        // axis2 = "z";
        break;

      default:
        phaseVals_pos = phaseVals_zPos;
        phaseVals_neg = phaseVals_zNeg;
        // neglectedAxis = "z";
        center_mri_axis = csz;
        // PVP1 = phaseVals_xPos;
        // PVN1 = phaseVals_xNeg;
        // P1C = csx;
        // axis1 = "x";
        // PVP2 = phaseVals_yPos;
        // PVN2 = phaseVals_yNeg;
        // P2C = csy;
        // axis2 = "y";
        break;
    }

    Calculate_Magnetic_Moment_3D.logger.addVariable("MRI_axis", MRI_axis);
    // returning the estimated radius along the MRI field direction
    return calculateRC(center_mri_axis, MRI_axis, phaseVals_pos, phaseVals_neg);
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
    center_l = new Triplet<Double>(0.0, 0.0, 0.0);

    for (int i = roi_xi; i < roi_xi + roi_dx; i++) {
      for (int j = roi_yi; j < roi_yi + roi_dy; j++) {
        for (int k = roi_zi; k < roi_zi + roi_dz; k++) {
          if (ImageMethods.getVoxelValue(mag_img, i, j, k) < min) {
            min = ImageMethods.getVoxelValue(mag_img, i, j, k);
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

    center_m = new Triplet<Double>(0.0, 0.0, 0.0);
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

    center_s = new Triplet<Double>(-1.0, -1.0, -1.0);

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
              k + roi_mag_belowM_zi);
        }
      }
    }

    // Summing x plane and putting it into array
    for (int i = 0; i <= roi_mag_belowM_dx; i++) {
      innerBox_sumOfXPlane[i] = 0;
      for (int j = 0; j <= roi_mag_belowM_dy; j++) {
        for (int k = 0; k <= roi_mag_belowM_dz; k++) {
          innerBox_sumOfXPlane[i] += ImageMethods.getVoxelValue(mag_img, i + roi_mag_belowM_xi, j + roi_mag_belowM_yi,
              k + roi_mag_belowM_zi);
        }
      }
    }

    // Summing y plane and putting it into array
    for (int j = 0; j <= roi_mag_belowM_dy; j++) {
      innerBox_sumOfYPlane[j] = 0;
      for (int k = 0; k <= roi_mag_belowM_dz; k++) {
        for (int i = 0; i <= roi_mag_belowM_dx; i++) {
          innerBox_sumOfYPlane[j] += ImageMethods.getVoxelValue(mag_img, i + roi_mag_belowM_xi, j + roi_mag_belowM_yi,
              k + roi_mag_belowM_zi);
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

    // if a component can't be found throw error
    if (center_s.get(0) == -1.0 || center_s.get(1) == -1.0 || center_s.get(2) == -1.0)
      throw new ArithmeticException("Unable to determine center_s completely");

  }

  // =====================================================================================
  // Getter for Center_L
  // =====================================================================================
  public Triplet<Double> centerL() {
    return center_l;
  }

  // =====================================================================================
  // Getter for Center_M
  // =====================================================================================
  public Triplet<Double> centerM() {
    return center_m;
  }

  // =====================================================================================
  // Getter for Center_S
  // =====================================================================================
  public Triplet<Double> centerS() {
    return center_s;
  }

  // =====================================================================================
  // Getter for MRI axis
  // =====================================================================================
  public Axis MRIAxis() {
    MRI_axis = Calculate_Magnetic_Moment_3D.gui.ldd_MriAxis.getValue();
    Calculate_Magnetic_Moment_3D.logger.addVariable("MRI_Axis", MRI_axis);
    return MRI_axis;
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
   * @return a Triplet with contents being the fitted point, the fitted delta,
   *         then null
   */
  private Triplet<Integer> goodROI(int img_size, int point, int delta) {

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

    return new Triplet<Integer>(point, delta, null);
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