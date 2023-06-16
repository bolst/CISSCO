import java.awt.Rectangle;
import javax.swing.JOptionPane;

import ij.ImagePlus;
import ij.gui.*;
import ij.WindowManager;

public class ImageItem {

  private ImagePlus mag_img, phase_img;
  private Triplet<Double> center_l;
  private Triplet<Double> center_m;
  private Triplet<Double> center_s;
  public Triplet<Double> subpix_image_center;
  public double bkgPhase;
  private double roi_mag_belowM_sumX, roi_mag_belowM_sumY, roi_mag_belowM_sumZ;
  public int roi_xi, roi_yi, roi_zi, roi_dx, roi_dy, roi_dz;
  public int roi_mag_belowM_xi, roi_mag_belowM_yi, roi_mag_belowM_zi,
      roi_mag_belowM_Dx, roi_mag_belowM_Dy, roi_mag_belowM_Dz;
  public int M;
  private double RCenter;
  private double phaseValue;
  private Axis MRI_axis;
  private double m_R0;
  private double estMagMoment;
  private double R1PhaseCalc, R2PhaseCalc, R3PhaseCalc;

  private final int grid = 10;
  private final double m_ROuterFrom = 0.2;
  private final double m_RMiddleFrom = 0.9;
  private final double m_RInnerFrom = 2.5;

  public boolean isNearEdge = false;

  public enum Axis {
    X, Y, Z;

    @Override
    public String toString() {
      if (this == Axis.X) {
        return "X";
      } else if (this == Axis.Y) {
        return "Y";
      }
      return "Z";
    }
  }

  public ImageItem(String magTitle, String phaseTitle, int M_pct, double pV) {
    M = M_pct;
    phaseValue = pV;

    // If mag and phase images are not open
    if (WindowManager.getImage(magTitle) == null || WindowManager.getImage(phaseTitle) == null) {
      throw new IllegalStateException("No magnitude or phase file open.");
    }

    // Getting image instances
    try {
      mag_img = WindowManager.getImage(magTitle);
      phase_img = WindowManager.getImage(phaseTitle);
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
    // square box dimensions are the larger of roi_dx and roi_dy
    roi_dx = roi_dy = roi_dz = Math.max((int) user_roi_rectangle.getWidth(),
        (int) user_roi_rectangle.getHeight());

    Calculate_Magnetic_Moment_3D.logger.addVariable("roi_dz", roi_dz);

    // roi side lengths must be odd
    roi_dx = roi_dy = roi_dz = (roi_dx % 2 == 0) ? (roi_dx + 1) : roi_dx;

    // Setting new ROI to rectangle
    user_roi_rectangle.setBounds(roi_xi, roi_yi, roi_dz,
        roi_dz);

    // Setting new ROI to image
    if (mag_img.getRoi() != null)
      mag_img.setRoi(user_roi_rectangle);
    else
      phase_img.setRoi(user_roi_rectangle);

    Calculate_Magnetic_Moment_3D.logger.addVariable("roi_dz", roi_dz);

    // get current and # of slices from whatever image ROI is on
    int Sn = (mag_img.getRoi() != null) ? mag_img.getNSlices() : phase_img.getNSlices();
    int Si = (mag_img.getRoi() != null) ? mag_img.getCurrentSlice() : phase_img.getCurrentSlice();

    // If current slice is too close to begin/end of img - need to fit dz (slice
    // range) to be within image slice range

    // if ROI extends beyond the 1st slice, then fit dz to be twice
    // the distance from the 1st slice to the current slice (and + 1 for center)
    if (Si - roi_dz / 2 < 1) {
      roi_dz = 2 * Si - 1;
      roi_zi = 0;
    }
    // if ROI extends beyond last slice, then fit dz to be twice the
    // distance from the last slice to the current slice (and +1 for center)
    else if (Si + roi_dz / 2 > Sn) {
      roi_dz = 2 * (Sn - Si) + 1;
      roi_zi = Sn - roi_dz;
    }
    // if slice is ok
    else {
      roi_zi = Si - roi_dz / 2 - 1;
    }

    if (Si - roi_dz / 2 < 1) // if slice is too close to first slice
    {
      roi_dz = 2 * Si - 1;
      roi_zi = 0;
    } else if (Si + roi_dz / 2 > Sn) // if slice is too close to last slice
    {
      roi_dz = 2 * (Sn - Si) + 1;
      roi_zi = Sn - roi_dz;
    } else // if slice is ok
    {
      roi_zi = Si - roi_dz / 2 - 1;
    }

    Calculate_Magnetic_Moment_3D.logger.addVariable("roi_xi", roi_xi);
    Calculate_Magnetic_Moment_3D.logger.addVariable("roi_yi", roi_yi);
    Calculate_Magnetic_Moment_3D.logger.addVariable("roi_zi", roi_zi);
    Calculate_Magnetic_Moment_3D.logger.addVariable("roi_dx", roi_dx);
    Calculate_Magnetic_Moment_3D.logger.addVariable("roi_dy", roi_dy);
    Calculate_Magnetic_Moment_3D.logger.addVariable("roi_dz", roi_dz);

    // ------- Finding innerbox that contains values below threshold

    // setting slice to starting ROI point in the z
    mag_img.setSlice(roi_zi + 1);

    // summing all corners of ROI box in this slice
    float avgCorners = mag_img.getProcessor().getPixelValue(roi_xi,
        roi_yi) +
        mag_img.getProcessor().getPixelValue(roi_xi + roi_dx,
            roi_yi)
        + mag_img.getProcessor().getPixelValue(roi_xi,
            roi_yi + roi_dy)
        + mag_img.getProcessor().getPixelValue(roi_xi + roi_dx,
            roi_xi + roi_dy);

    // setting slice to ending ROI point in the z
    mag_img.setSlice(roi_zi + roi_dz + 1);

    // adding these corners of the ROI box in this slice
    avgCorners += mag_img.getProcessor().getPixelValue(roi_xi, roi_yi)
        + mag_img.getProcessor().getPixelValue(roi_xi + roi_dx, roi_yi)
        + mag_img.getProcessor().getPixelValue(roi_xi, roi_yi + roi_dy)
        + mag_img.getProcessor().getPixelValue(roi_xi + roi_dx, roi_xi + roi_dy);

    // dividing by 8 for average
    avgCorners /= 8.0;

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

    // a new ROI is defined with all the values in the user's ROI that are less than
    // this
    double newbox_max = avgCorners * (double) M / 100.0;

    // getting all pixel values that are below value
    for (int dz = 0; dz < roi_dz; dz++) {
      mag_img.setSlice(roi_zi + dz + 1);
      for (int dx = 0; dx < roi_dx; dx++) {
        for (int dy = 0; dy < roi_dy; dy++) {
          if (mag_img.getProcessor().getPixelValue(roi_xi + dx, roi_yi + dy) < newbox_max) // if below value
          {
            // place value into box
            roi_mag_belowM[dx][dy][dz] = mag_img.getProcessor().getPixelValue(roi_xi + dx, roi_yi + dy);
          }
        }
      }
    }

    roi_mag_belowM_xi = roi_xi + roi_dx;
    roi_mag_belowM_yi = roi_yi + roi_dy;
    roi_mag_belowM_zi = roi_zi + roi_dz;
    roi_mag_belowM_Dx = 0;
    roi_mag_belowM_Dy = 0;
    roi_mag_belowM_Dz = 0;

    // putting original corner of small box into variables
    // stores the smallest indices that have values not equal to -1.0
    for (int dz = 0; dz < roi_dz; dz++) {
      for (int dx = 0; dx < roi_dx; dx++) {
        for (int dy = 0; dy < roi_dy; dy++) {
          if (roi_mag_belowM[dx][dy][dz] != -1.0) {
            if (roi_xi + dx < roi_mag_belowM_xi) {
              roi_mag_belowM_xi = roi_xi + dx;
            }
            if (roi_yi + dy < roi_mag_belowM_yi) {
              roi_mag_belowM_yi = roi_yi + dy;
            }
            if (roi_zi + dz < roi_mag_belowM_zi) {
              roi_mag_belowM_zi = roi_zi + dz;
            }
          }
        }
      }
    }

    // putting size of small box into variables
    // stores the smallest size that contains all values below newbox_max
    for (int dz = 0; dz < roi_dz; dz++) {
      for (int dx = 0; dx < roi_dx; dx++) {
        for (int dy = 0; dy < roi_dy; dy++) {
          if (roi_mag_belowM[dx][dy][dz] != -1.0) {
            if (roi_mag_belowM_Dx < dx + roi_xi - roi_mag_belowM_xi) {
              roi_mag_belowM_Dx = dx + roi_xi - roi_mag_belowM_xi;
            }
            if (roi_mag_belowM_Dy < dy + roi_yi - roi_mag_belowM_yi) {
              roi_mag_belowM_Dy = dy + roi_yi - roi_mag_belowM_yi;
            }
            if (roi_mag_belowM_Dz < dz + roi_zi - roi_mag_belowM_zi) {
              roi_mag_belowM_Dz = dz + roi_zi - roi_mag_belowM_zi;
            }
          }
        }
      }
    }

    // ------ End finding innerbox
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
  private Axis calculateMRIAxis(int accuracy,
      double[] XP,
      double[] XN,
      double[] YP,
      double[] YN,
      double[] ZP,
      double[] ZN) {

    // finding radius where no values are 0.0 (recall values were set to 0.0 if
    // below M%)
    int dr = 0;
    while (dr < XP.length &&
        XP[dr] == 0.0 ||
        XN[dr] == 0.0 ||
        YP[dr] == 0.0 ||
        YN[dr] == 0.0 ||
        ZP[dr] == 0.0 ||
        ZN[dr] == 0.0)
      dr++;

    // Calculating ratio of each direction
    // The non-MRI field directions should have a ratio of ~1
    // The field that is responsible of the ratio not near 1 is the MRI field
    // direction

    // requested accuracy cannot surpass the array lengths
    accuracy = Math.min(accuracy, XP.length - dr - 1);
    double avgX = 0.0;
    double avgY = 0.0;
    double avgZ = 0.0;

    for (int i = dr; i < dr + accuracy; i++) {
      avgX += (Math.abs(XP[i]) + Math.abs(XN[i])) / accuracy;
      avgY += (Math.abs(YP[i]) + Math.abs(YN[i])) / accuracy;
      avgZ += (Math.abs(ZP[i]) + Math.abs(ZN[i])) / accuracy;
    }

    double ratioXY = Math.abs(avgX / avgY - 1);
    double ratioXZ = Math.abs(avgX / avgZ - 1);
    double ratioYZ = Math.abs(avgY / avgZ - 1);
    double minRatio = Math.min(Math.min(ratioXY, ratioXZ), ratioYZ);

    Axis retval = minRatio == ratioXY ? Axis.Z : (minRatio == ratioXZ ? Axis.Y : Axis.X);
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

    // Setting slice to initial z point
    mag_img.setSlice(roi_zi + 1);
    // getting average of ROI
    double avgOfCorners = (double) (mag_img.getProcessor().getPixelValue(roi_xi,
        roi_yi) +
        mag_img.getProcessor().getPixelValue(roi_xi + roi_dx,
            roi_yi)
        +
        mag_img.getProcessor().getPixelValue(roi_xi,
            roi_yi + roi_dy)
        +
        mag_img.getProcessor().getPixelValue(roi_xi + roi_dx,
            roi_yi + roi_dy));

    mag_img.setSlice(roi_zi + roi_dz + 1);
    avgOfCorners += (double) (mag_img.getProcessor().getPixelValue(roi_xi,
        roi_yi) +
        mag_img.getProcessor().getPixelValue(roi_xi + roi_dx,
            roi_yi)
        +
        mag_img.getProcessor().getPixelValue(roi_xi,
            roi_yi + roi_dy)
        +
        mag_img.getProcessor().getPixelValue(roi_xi + roi_dx,
            roi_yi + roi_dy));
    avgOfCorners /= 8.0;

    // all values below this variable will be negated
    double maxMagValue = avgOfCorners * (double) threshold / 100.0;

    mag_img.setSlice(center_s.get(2).intValue() + 1);

    // setting all values in array below maxMagValue to 0
    switch (axis) {
      case X:
        if (direction) {
          for (int i = 0; i < grid + 1; i++) {
            if (Math.abs(
                mag_img.getProcessor().getPixelValue(i + r0, center_s.get(1).intValue())) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        } else {
          for (int i = 0; i < grid + 1; i++) {
            if (Math.abs(mag_img.getProcessor().getPixelValue(r0 - i, center_s.get(1).intValue())) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        }
      case Y:
        if (direction) {
          for (int i = 0; i < grid + 1; i++) {
            if (Math.abs(mag_img.getProcessor().getPixelValue(center_s.get(0).intValue(), i + r0)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        } else {
          for (int i = 0; i < grid + 1; i++) {
            if (Math.abs(mag_img.getProcessor().getPixelValue(center_s.get(0).intValue(), r0 - i)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        }
      case Z:
        if (direction) {
          for (int i = 0; i < grid + 1; i++) {
            mag_img.setSlice(i + r0 + 1);
            if (Math.abs(mag_img.getProcessor().getPixelValue(center_s.get(0).intValue(),
                center_s.get(1).intValue())) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        } else {
          for (int i = 0; i < grid + 1; i++) {
            mag_img.setSlice(r0 - i + 1);
            if (Math.abs(mag_img.getProcessor().getPixelValue(center_s.get(0).intValue(),
                center_s.get(1).intValue())) < maxMagValue) {
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
    double r1 = 1.6;
    double r2 = 1.6;

    /*
     * The logic below works by finding a voxel phase value that is between the
     * user-inputted phase value.
     * So by default it is 1, so the program finds the voxel that is between two
     * voxels that have greater phase than 1 and less than 1.
     * The program then interpolates the distance between these two voxels to
     * estimate how far from the center a phase value of 1 is
     */

    for (int i = 0; i < grid; i++) {
      // this condition is if the current voxel is nested between two voxels greater
      // and less than the user-defined phase value
      if ((phaseVals_pos[i] > phaseValue) && (phaseVals_pos[i + 1] < phaseValue)) {
        // this condition is if the interpolation value is greater than 1.6, because the
        // object radius has to be greater than 1.6
        if (interpolation(phaseValue, phaseVals_pos[i], phaseVals_pos[i + 1], i, i + 1) > 1.6) {
          r1 = interpolation(phaseValue, phaseVals_pos[i], phaseVals_pos[i + 1], i, i + 1);
        }
        break;
      }
    }

    // Same thing as above just the opposite direction
    for (int i = 0; i < grid; i++) {
      if ((phaseVals_neg[i] > phaseValue) && (phaseVals_neg[i + 1] < phaseValue)) {
        if (interpolation(phaseValue, phaseVals_neg[i], phaseVals_neg[i + 1], i, i + 1) > 1.6) {
          r2 = interpolation(phaseValue, phaseVals_neg[i], phaseVals_neg[i + 1], i, i + 1);
        }
        break;
      }
    }

    // Set the center_s MRI axis value to the average of the found radii + the
    // found center point of the MRI field axis
    if (mri_dir == Axis.X) {
      center_s.set(0, ((center_mri_axis + r1) + (center_mri_axis - r2)) / 2.0);
      center_s.set(0, Math.round(center_s.get(0) * 10.0) / 10.0);
      Calculate_Magnetic_Moment_3D.gui.ltf_rcx.setValue(String.valueOf(Math.round(center_s.get(0) * 10.0) / 10.0));
    }
    if (mri_dir == Axis.Y) {
      center_s.set(1, ((center_mri_axis + r1) + (center_mri_axis - r2)) / 2.0);
      center_s.set(1, Math.round(center_s.get(1) * 10.0) / 10.0);
      Calculate_Magnetic_Moment_3D.gui.ltf_rcy.setValue(String.valueOf(Math.round(center_s.get(1) * 10.0) / 10.0));
    }
    if (mri_dir == Axis.Z) {
      center_s.set(2, ((center_mri_axis + r1) + (center_mri_axis - r2)) / 2.0);
      center_s.set(2, Math.round(center_s.get(2) * 10.0) / 10.0);
      Calculate_Magnetic_Moment_3D.gui.ltf_rcz.setValue(String.valueOf(Math.round(center_s.get(2) * 10.0) / 10.0));
    }

    // returning equation for RCenter
    RCenter = ((r1 + r2) / 2.0) / Math.cbrt(2);
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

  public double estimateRCenter() {
    int csx = center_s.get(0).intValue();
    int csy = center_s.get(1).intValue();
    int csz = center_s.get(2).intValue();

    // Phase values for each direction
    double[] phaseVals_xPos = new double[grid + 1];
    double[] phaseVals_yPos = new double[grid + 1];
    double[] phaseVals_zPos = new double[grid + 1];
    double[] phaseVals_xNeg = new double[grid + 1];
    double[] phaseVals_yNeg = new double[grid + 1];
    double[] phaseVals_zNeg = new double[grid + 1];
    Calculate_Magnetic_Moment_3D.logger.addVariable("pvp array size", grid + 1);

    // Setting slice to initial z
    phase_img.setSlice(csz + 1);

    Calculate_Magnetic_Moment_3D.logger.addVariable("csx", csx);
    Calculate_Magnetic_Moment_3D.logger.addVariable("csy", csy);
    Calculate_Magnetic_Moment_3D.logger.addVariable("csz", csz);
    // putting phase values into array, following respectively for other 5 loops
    for (int i = 0; i < grid + 1; i++) {
      phaseVals_xPos[i] = Math.abs((double) phase_img.getProcessor().getPixelValue(i + csx, csy));
    }
    for (int i = 0; i < grid + 1; i++) {
      phaseVals_xNeg[i] = Math.abs((double) phase_img.getProcessor().getPixelValue(csx - i, csy));
    }
    for (int j = 0; j < grid + 1; j++) {
      phaseVals_yPos[j] = Math.abs((double) phase_img.getProcessor().getPixelValue(csx, j + csy));
    }
    for (int j = 0; j < grid + 1; j++) {
      phaseVals_yNeg[j] = Math.abs((double) phase_img.getProcessor().getPixelValue(csx, csy - j));
    }
    for (int k = 0; k < grid + 1; k++) {
      phase_img.setSlice(k + csz + 1);
      phaseVals_zPos[k] = Math.abs((double) phase_img.getProcessor().getPixelValue(csx, csy));
    }
    for (int k = 0; k < grid + 1; k++) {
      phase_img.setSlice(csz - k + 1);
      phaseVals_zNeg[k] = Math.abs((double) phase_img.getProcessor().getPixelValue(csx, csy));
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
    switch (calculateMRIAxis(grid - 2,
        phaseVals_xPos,
        phaseVals_xNeg,
        phaseVals_yPos,
        phaseVals_yNeg,
        phaseVals_zPos,
        phaseVals_zNeg)) {
      case X:
        phaseVals_pos = phaseVals_xPos;
        phaseVals_neg = phaseVals_xNeg;
        // neglectedAxis = "x";
        MRI_axis = Axis.X;
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
        MRI_axis = Axis.Y;
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
        MRI_axis = Axis.Z;
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

  public double estBkg() {
    // ---------- Begin to find estimated background phase
    phase_img.setSlice(roi_zi + 1);
    bkgPhase = Math
        .abs(phase_img.getProcessor().getPixelValue(roi_xi, roi_yi)) +
        Math.abs(phase_img.getProcessor().getPixelValue(roi_xi + roi_dx,
            roi_yi))
        +
        Math.abs(phase_img.getProcessor().getPixelValue(roi_xi,
            roi_yi + roi_dy))
        +
        Math.abs(phase_img.getProcessor().getPixelValue(roi_xi + roi_dx,
            roi_yi + roi_dy));

    phase_img.setSlice(roi_zi + roi_dz + 1);
    bkgPhase += Math
        .abs(phase_img.getProcessor().getPixelValue(roi_xi, roi_yi)) +
        Math.abs(phase_img.getProcessor().getPixelValue(roi_xi + roi_dx,
            roi_yi))
        +
        Math.abs(phase_img.getProcessor().getPixelValue(roi_xi,
            roi_yi + roi_dy))
        +
        Math.abs(phase_img.getProcessor().getPixelValue(roi_xi + roi_dx,
            roi_yi + roi_dy));

    bkgPhase /= 8.0;

    return bkgPhase;

    // ---------- end to find background phase
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
    m_R0 = m_R0 * (double) grid;
    m_R0 = Math.ceil(m_R0) / (double) grid;

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

    // logger.addVariable("m_R0", m_R0);

    estMagMoment = phaseValue * Math.pow(RCenter, 3);

    R1PhaseCalc = estMagMoment / Math.pow(m_R1, 3);
    R2PhaseCalc = estMagMoment / Math.pow(m_R2, 3);
    R3PhaseCalc = estMagMoment / Math.pow(m_R3, 3);
    Calculate_Magnetic_Moment_3D.gui.lbl_r1phaseCalc.setText(String.valueOf(Math.round(R1PhaseCalc * 100.0) / 100.0));
    Calculate_Magnetic_Moment_3D.gui.lbl_r2phaseCalc.setText(String.valueOf(Math.round(R2PhaseCalc * 100.0) / 100.0));
    Calculate_Magnetic_Moment_3D.gui.lbl_r3phaseCalc.setText(String.valueOf(Math.round(R3PhaseCalc * 100.0) / 100.0));

    // setMagMoment(phaseValue * Math.pow(RCenter, 3));

    Calculate_Magnetic_Moment_3D.estimateCenterRadii_isClicked = true;
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
          mag_img.setSlice(k + 1);
          if (mag_img.getProcessor().getPixelValue(i, j) < min) {
            min = mag_img.getProcessor().getPixelValue(i, j);
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
    center_m = new Triplet<Double>(0.0, 0.0, 0.0);

    // Making Center_M the small box to start
    // center_m.set(0, (double) (roi_mag_belowM_xi + roi_mag_belowM_xi +
    // roi_mag_belowM_Dx));
    // center_m.set(1, (double) (roi_mag_belowM_yi + roi_mag_belowM_yi +
    // roi_mag_belowM_Dy));
    // center_m.set(2, (double) (roi_mag_belowM_zi + roi_mag_belowM_zi +
    // roi_mag_belowM_Dz));

    // Center of the box is Center_M

    // I am realizing this chunk of code is literally doing the same thing in
    // both if else statements? Refactored underneath this commented block
    /*
     * if ((int) center_mx % 2 == 0.0) {
     * center_mx /= 2.0;
     * center_mx += 0.5;
     * } else {
     * center_mx = (center_mx + 1.0) / 2.0;
     * }
     * if ((int) center_my % 2.0 == 0.0) {
     * center_my /= 2.0;
     * center_my += 0.5;
     * } else {
     * center_my = (center_mz + 1.0) / 2.0;
     * }
     * if ((int) center_mz % 2.0 == 0.0) {
     * center_mz /= 2.0;
     * center_mz += 0.5;
     * } else {
     * center_mz = (center_mz + 1.0) / 2.0;
     * }
     */

    // center_m.set(0, center_m.get(0) / 2.0 + 0.5);
    // center_m.set(1, center_m.get(1) / 2.0 + 0.5);
    // center_m.set(2, center_m.get(2) / 2.0 + 0.5);

    // center_m is center of small ROI
    center_m.set(0, (double) roi_mag_belowM_xi + 0.5 * (double) (roi_mag_belowM_Dx + 1));
    center_m.set(1, (double) roi_mag_belowM_yi + 0.5 * (double) (roi_mag_belowM_Dy + 1));
    center_m.set(2, (double) roi_mag_belowM_zi + 0.5 * (double) (roi_mag_belowM_Dz + 1));

    Calculate_Magnetic_Moment_3D.logger.addVariable("center_m", center_m);
  }

  // =====================================================================================
  // Function to calculate center_s, i.e. coordinate of minimum sum of each plane
  // =====================================================================================
  public void calcCenterS() {

    center_s = new Triplet<Double>(0.0, 0.0, 0.0);

    // ------ Begin summing planes in innerbox
    // Initializing arrays for sums of planes in small box
    double[] innerBox_sumOfXPlane = new double[roi_mag_belowM_Dx + 1];
    double[] innerBox_sumOfYPlane = new double[roi_mag_belowM_Dy + 1];
    double[] innerBox_sumOfZPlane = new double[roi_mag_belowM_Dz + 1];

    Calculate_Magnetic_Moment_3D.logger.addVariable("rmb xi", roi_mag_belowM_xi);
    Calculate_Magnetic_Moment_3D.logger.addVariable("rmb yi", roi_mag_belowM_yi);
    Calculate_Magnetic_Moment_3D.logger.addVariable("rmb zi", roi_mag_belowM_zi);
    Calculate_Magnetic_Moment_3D.logger.addVariable("rmb dx", roi_mag_belowM_Dx);
    Calculate_Magnetic_Moment_3D.logger.addVariable("rmb dy", roi_mag_belowM_Dy);
    Calculate_Magnetic_Moment_3D.logger.addVariable("rmb dx", roi_mag_belowM_Dz);

    // Summing z plane and putting it into array
    for (int k = 0; k <= roi_mag_belowM_Dz; k++) {
      mag_img.setSlice(k + roi_mag_belowM_zi + 1);
      innerBox_sumOfZPlane[k] = 0;
      for (int i = 0; i <= roi_mag_belowM_Dx; i++) {
        for (int j = 0; j <= roi_mag_belowM_Dy; j++) {
          innerBox_sumOfZPlane[k] += mag_img.getProcessor().getPixelValue(i + roi_mag_belowM_xi, j + roi_mag_belowM_yi);
        }
      }
    }

    // Summing x plane and putting it into array
    for (int i = 0; i <= roi_mag_belowM_Dx; i++) {
      innerBox_sumOfXPlane[i] = 0;
      for (int j = 0; j <= roi_mag_belowM_Dy; j++) {
        for (int k = 0; k <= roi_mag_belowM_Dz; k++) {
          mag_img.setSlice(k + roi_mag_belowM_zi + 1);
          innerBox_sumOfXPlane[i] += mag_img.getProcessor().getPixelValue(i + roi_mag_belowM_xi, j + roi_mag_belowM_yi);
        }
      }
    }

    // Summing y plane and putting it into array
    for (int j = 0; j <= roi_mag_belowM_Dy; j++) {
      innerBox_sumOfYPlane[j] = 0;
      for (int k = 0; k <= roi_mag_belowM_Dz; k++) {
        mag_img.setSlice(k + roi_mag_belowM_zi + 1);
        for (int i = 0; i <= roi_mag_belowM_Dx; i++) {
          innerBox_sumOfYPlane[j] += mag_img.getProcessor().getPixelValue(i + roi_mag_belowM_xi, j + roi_mag_belowM_yi);
        }
      }
    }

    // for finding the minimum summed plane
    roi_mag_belowM_sumX = innerBox_sumOfXPlane[0];
    roi_mag_belowM_sumY = innerBox_sumOfYPlane[0];
    roi_mag_belowM_sumZ = innerBox_sumOfZPlane[0];

    // ------ End summing planes in innerbox

    // ---------- begin to find Center_S

    center_s.set(0, (double) roi_mag_belowM_xi);
    center_s.set(1, (double) roi_mag_belowM_yi);
    center_s.set(2, (double) roi_mag_belowM_zi);

    // Finding minimum x plane and putting it into center_s.get(0)
    for (int i = 0; i <= roi_mag_belowM_Dx; i++) {
      if (innerBox_sumOfXPlane[i] < roi_mag_belowM_sumX) {
        roi_mag_belowM_sumX = innerBox_sumOfXPlane[i];
        center_s.set(0, (double) (i + roi_mag_belowM_xi));
      }
    }

    // Finding minimum y plane and putting it into center_s.get(1)
    for (int j = 0; j <= roi_mag_belowM_Dy; j++) {
      if (innerBox_sumOfYPlane[j] < roi_mag_belowM_sumY) {
        roi_mag_belowM_sumY = innerBox_sumOfYPlane[j];
        center_s.set(1, (double) (j + roi_mag_belowM_yi));
      }
    }

    // Finding minimum z plane and putting it into center_s.get(2)
    for (int k = 0; k <= roi_mag_belowM_Dz; k++) {
      if (innerBox_sumOfZPlane[k] < roi_mag_belowM_sumZ) {
        roi_mag_belowM_sumZ = innerBox_sumOfZPlane[k];
        center_s.set(2, (double) (k + roi_mag_belowM_zi));
      }
    }
    // ---------- end to find Center_S
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

}
