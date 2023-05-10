import java.awt.Rectangle;
import java.util.ArrayList;
import javax.swing.JOptionPane;

import ij.ImagePlus;
import ij.gui.*;
import ij.WindowManager;

public class ImageItem {

  private ImagePlus mag_img, phase_img;
  private Triplet<Double> center_l;
  private Triplet<Double> center_m;
  private Triplet<Double> center_s;
  public double bkgPhase;
  private double roi_mag_belowM_sumX, roi_mag_belowM_sumY, roi_mag_belowM_sumZ;
  public int roi_xi, roi_yi, roi_zi, roi_dx, roi_dy, roi_dz;
  public int roi_mag_belowM_xi, roi_mag_belowM_yi, roi_mag_belowM_zi,
      roi_mag_belowM_Dx, roi_mag_belowM_Dy, roi_mag_belowM_Dz;
  public int M;
  private double RCenter;
  private double phaseValue;
  private String neglectedAxis;
  private double m_R0;
  private double estMagMoment;
  private double R1PhaseCalc, R2PhaseCalc, R3PhaseCalc;

  private final int grid = 10;
  private final double m_ROuterFrom = 0.2;
  private final double m_RMiddleFrom = 0.9;
  private final double m_RInnerFrom = 2.5;

  private final int MIN_ROI_DZ = 5;

  public boolean isNearEdge = false;

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

    // If current slice is too close to begin/end of img - need to fit dz (slice
    // range) to be within image slice range
    if (mag_img.getRoi() != null) // if roi is on mag image
    {
      // if dz range is ok
      if (mag_img.getCurrentSlice() - (roi_dz + 1) / 2 < 0) {
        roi_zi = 0;
        roi_dz = Math.max(2 * mag_img.getCurrentSlice(), MIN_ROI_DZ);
      } else if (mag_img.getCurrentSlice() + (roi_dz + 1) / 2 > mag_img.getNSlices()) {
        roi_dz = Math.max(mag_img.getNSlices() - mag_img.getCurrentSlice(), MIN_ROI_DZ);
        roi_zi = mag_img.getNSlices() - roi_dz - 1;
      }
      // if dz spans past slice range
      else
        roi_zi = mag_img.getCurrentSlice() - (roi_dz + 1) / 2 - 1;

    } else // if roi is on phase image
    {
      if (phase_img.getCurrentSlice() - (roi_dz + 1) / 2 < 0) // if dz range is ok
      {
        roi_zi = 0;
        roi_dz = Math.max(2 * phase_img.getCurrentSlice(), MIN_ROI_DZ);
      } else if (phase_img.getCurrentSlice() + (roi_dz + 1) / 2 > phase_img.getNSlices()) {
        roi_dz = Math.max(phase_img.getNSlices() - phase_img.getCurrentSlice(), MIN_ROI_DZ);
        roi_zi = phase_img.getNSlices() - roi_dz - 1;
      }
      // if dz spans past slice range
      else
        roi_zi = phase_img.getCurrentSlice() - (roi_dz + 1) / 2 - 1;
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

    // getting text of M% (should be 50 unless user changes it)
    /*
     * The program will use this to check what values are below this % of the
     * average of the ROI corners calculated above.
     * These values will be put into a new box for Center_M
     */
    // int M = Integer.parseInt(txt_M.getText());

    // getting all pixel values that are below threshold %
    for (int k = 0; k < roi_dz; k++) {
      mag_img.setSlice(k + roi_zi + 1);
      for (int i = 0; i < roi_dx; i++) {
        for (int j = 0; j < roi_dy; j++) {
          if (mag_img.getProcessor().getPixelValue(i + roi_xi,
              j + roi_yi) < (avgCorners * M / 100.0)) {
            roi_mag_belowM[i][j][k] = mag_img.getProcessor().getPixelValue(i + roi_xi, j + roi_yi);
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
    /*
     * Program takes the smallest index that has a value that isn't -1.0. This is
     * because the box will completely surround where
     * these values are. It will be a box that fits perfectly around the values
     * below M%
     */
    for (int k = 0; k < roi_dz; k++) {
      for (int i = 0; i < roi_dx; i++) {
        for (int j = 0; j < roi_dy; j++) {
          if (roi_mag_belowM[i][j][k] != -1.0) {
            if (i + roi_xi < roi_mag_belowM_xi) {
              roi_mag_belowM_xi = i + roi_zi;
            }
            if (j + roi_yi < roi_mag_belowM_yi) {
              roi_mag_belowM_yi = j + roi_yi;
            }
            if (k + roi_zi < roi_mag_belowM_zi) {
              roi_mag_belowM_zi = k + roi_zi;
            }
          }
        }
      }
    }

    // putting opposite corner into variables
    /*
     * This has the same logic as the above loop, but now it is finding the maximum
     * index that contains values below M% instead of minimum
     */
    for (int k = 0; k < roi_dz; k++) {
      for (int i = 0; i < roi_dx; i++) {
        for (int j = 0; j < roi_dy; j++) {
          if (roi_mag_belowM[i][j][k] != -1.0) {
            if (roi_mag_belowM_Dx < i + roi_xi - roi_mag_belowM_xi) {
              roi_mag_belowM_Dx = i + roi_xi - roi_mag_belowM_xi;
            }
            if (roi_mag_belowM_Dy < j + roi_yi - roi_mag_belowM_yi) {
              roi_mag_belowM_Dy = j + roi_yi - roi_mag_belowM_yi;
            }
            if (roi_mag_belowM_Dz < k + roi_zi - roi_mag_belowM_zi) {
              roi_mag_belowM_Dz = k + roi_zi - roi_mag_belowM_zi;
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
  private String findFaultyAxis(int accuracy,
      double[] XP,
      double[] XN,
      double[] YP,
      double[] YN,
      double[] ZP,
      double[] ZN) {

    int ix = center_s.get(0).intValue();
    int iy = center_s.get(1).intValue();
    int iz = center_s.get(2).intValue();

    double deltaXY = 0.0;
    double deltaXZ = 0.0;
    double deltaYZ = 0.0;
    int accuracyCount = 0;
    double smallestDelta = 0.0;

    ArrayList<String> faultyAxis = new ArrayList<String>();

    while (accuracyCount < accuracy) {
      if ((XP[ix] == 0.0) || (YP[iy] == 0.0)
          || (ZP[iz] == 0.0)) {
        ix++;
        iy++;
        iz++;
        continue;
      }

      deltaXY = Math.abs(XP[ix] - YP[iy]);
      deltaXZ = Math.abs(XP[ix] - ZP[iz]);
      deltaYZ = Math.abs(YP[iy] - ZP[iz]);

      // ternary operator to get minimum value of deltas
      // may seem confusing, but you could always test this line using values of your
      // own to see if the logic works
      smallestDelta = (deltaXY < deltaXZ) ? (deltaXY < deltaYZ ? deltaXY : deltaYZ)
          : (deltaXZ < deltaYZ ? deltaXZ : deltaYZ);

      // the axis with the largest difference is added to the faulty axis list
      if (smallestDelta == deltaXY) {
        faultyAxis.add("z");
      }

      if (smallestDelta == deltaXZ) {
        faultyAxis.add("y");
      }

      if (smallestDelta == deltaYZ) {
        faultyAxis.add("x");
      }

      ix++;
      iy++;
      iz++;
      accuracyCount++;
    }

    accuracyCount = 0;
    /*
     * while (accuracyCount < accuracy) {
     *
     * if ((xPhaseValues_Negative[ix_Neg] == 0.0) || (yPhaseValues_Negative[iy_Neg]
     * == 0.0)
     * || (zPhaseValues_Negative[iz_Neg] == 0.0)) {
     * ix_Neg++;
     * iy_Neg++;
     * iz_Neg++;
     * continue;
     * }
     *
     * deltaXY_Neg = Math.abs(xPhaseValues_Negative[ix_Neg] -
     * yPhaseValues_Negative[iy_Neg]);
     * deltaXZ_Neg = Math.abs(xPhaseValues_Negative[ix_Neg] -
     * zPhaseValues_Negative[iz_Neg]);
     * deltaYZ_Neg = Math.abs(yPhaseValues_Negative[iy_Neg] -
     * zPhaseValues_Negative[iz_Neg]);
     *
     * if ((deltaXY_Neg < deltaXZ_Neg) && (deltaXY_Neg < deltaYZ_Neg)) {
     * faultyNegAxis.add("z");
     *
     * if ((deltaXZ_Neg < deltaXY_Neg) && (deltaXZ_Neg < deltaYZ_Neg)) {
     * faultyNegAxis.add("y");
     * }
     * if ((deltaYZ_Neg < deltaXY_Neg) && (deltaYZ_Neg < deltaXZ_Neg)) {
     * faultyNegAxis.add("x");
     * }
     *
     * ix_Neg++;
     * iy_Neg++;
     * iz_Neg++;
     * accuracyCount++;
     * }
     */
    // axis0Neg = faultyNegAxis.get(0);

    int xCount = 0;
    int yCount = 0;
    int zCount = 0;

    for (int i = 0; i < accuracy; i++) {
      if (faultyAxis.get(i) == "x") {
        xCount++;
      }
      if (faultyAxis.get(i) == "y") {
        yCount++;
      }
      if (faultyAxis.get(i) == "z") {
        zCount++;
      }
    }

    if ((zCount > yCount) && (zCount > xCount)) {
      return "z";
    }

    if ((yCount > xCount) && (yCount > zCount)) {
      return "y";
    }

    return "x";
    // if (axis0 != axis0Neg) {
    // return "n";
    // }
  }

  /*
   * Function to void all values below the param threshold.
   *
   * @param values array of values in one direction
   *
   * @param threshold the % to negate, default is 50
   *
   * @param startingPoint the center coordinate of the direction inputted
   *
   * @param axis The direction
   *
   * @param direction Flag to specify if function should iterate in the + or -
   * direction
   */
  private void negateValues(double[] values,
      int threshold,
      int startingPoint,
      String axis,
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
      case "x":
        if (direction == true) {
          for (int i = startingPoint; i < startingPoint + grid + 1; i++) {
            if (Math.abs(mag_img.getProcessor().getPixelValue(i, center_s.get(1).intValue())) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        }
        if (direction == false) {
          for (int i = startingPoint; i > startingPoint - grid - 1; i--) {
            if (Math.abs(mag_img.getProcessor().getPixelValue(i, center_s.get(1).intValue())) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        }
        break;
      case "y":
        if (direction == true) {
          for (int i = startingPoint; i < startingPoint + grid + 1; i++) {
            if (Math.abs(mag_img.getProcessor().getPixelValue(center_s.get(0).intValue(), i)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        }
        if (direction == false) {
          for (int i = startingPoint; i > startingPoint - grid - 1; i--) {
            if (Math.abs(mag_img.getProcessor().getPixelValue(center_s.get(0).intValue(), i)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        }
        break;
      case "z":
        if (direction == true) {
          for (int i = startingPoint; i < startingPoint + grid + 1; i++) {
            mag_img.setSlice(i + 1);
            if (Math.abs(mag_img.getProcessor().getPixelValue(center_s.get(0).intValue(),
                center_s.get(1).intValue())) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        }
        if (direction == false) {
          for (int i = startingPoint; i > startingPoint - grid - 1; i--) {
            mag_img.setSlice(i + 1);
            if (Math.abs(mag_img.getProcessor().getPixelValue(center_s.get(0).intValue(),
                center_s.get(1).intValue())) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        }
        break;
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
  private double findRadiiNeglected(int nC,
      String neglectedAxis,
      double[] neglectedPVP,
      double[] neglectedPVN) {
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

    for (int i = nC; i < nC + grid; i++) {
      // this condition is if the current voxel is nested between two voxels greater
      // and less than the user-defined phase value
      if ((neglectedPVP[i] > phaseValue) && (neglectedPVP[i + 1] < phaseValue)) {
        // this condition is if the interpolation value is greater than 1.6, because the
        // object radius has to be greater than 1.6
        if (interpolation(phaseValue, neglectedPVP[i], neglectedPVP[i + 1], i - nC, i - nC + 1) > 1.6) {
          r1 = interpolation(phaseValue, neglectedPVP[i], neglectedPVP[i + 1], i - nC, i - nC + 1);
        }
        break;
      }
    }

    // Same thing as above just the opposite direction
    for (int i = nC; i > nC - grid; i--) {
      if ((neglectedPVN[i] > phaseValue) && (neglectedPVN[i - 1] < phaseValue)) {
        if (interpolation(phaseValue, neglectedPVN[i], neglectedPVN[i - 1], nC - i, nC - i + 1) > 1.6) {
          r2 = interpolation(phaseValue, neglectedPVN[i], neglectedPVN[i - 1], nC - i, nC - i + 1);
        }
        break;
      }
    }

    // Set the center_s MRI axis value to the average of the found radii + the
    // found center point of the MRI field axis
    if (neglectedAxis.compareTo("x") == 0) {
      center_s.set(0, ((nC + r1) + (nC - r2)) / 2.0);
      center_s.set(0, Math.round(center_s.get(0) * 10.0) / 10.0);
      Calculate_Magnetic_Moment_3D.gui.ltf_rcx.setValue(String.valueOf(Math.round(center_s.get(0) * 10.0) / 10.0));
    }
    if (neglectedAxis.compareTo("y") == 0) {
      center_s.set(1, ((nC + r1) + (nC - r2)) / 2.0);
      center_s.set(1, Math.round(center_s.get(1) * 10.0) / 10.0);
      Calculate_Magnetic_Moment_3D.gui.ltf_rcy.setValue(String.valueOf(Math.round(center_s.get(1) * 10.0) / 10.0));
    }
    if (neglectedAxis.compareTo("z") == 0) {
      center_s.set(2, ((nC + r1) + (nC - r2)) / 2.0);
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
    // Phase values for each direction
    int csx = center_s.get(0).intValue();
    int csy = center_s.get(1).intValue();
    int csz = center_s.get(2).intValue();
    double[] xPhaseValues_Positive = new double[2 * (csx + grid + 1)];
    double[] yPhaseValues_Positive = new double[2 * (csy + grid + 1)];
    double[] zPhaseValues_Positive = new double[2 * (csz + grid + 1)];
    double[] xPhaseValues_Negative = new double[2 * (csx + grid + 1)];
    double[] yPhaseValues_Negative = new double[2 * (csy + grid + 1)];
    double[] zPhaseValues_Negative = new double[2 * (csz + grid + 1)];
    Calculate_Magnetic_Moment_3D.logger.addVariable("pvp array size", 2 * (csx + grid + 1));

    // Setting slice to initial z
    phase_img.setSlice(csz + 1);

    Calculate_Magnetic_Moment_3D.logger.addVariable("csx", csx);
    Calculate_Magnetic_Moment_3D.logger.addVariable("csy", csy);
    Calculate_Magnetic_Moment_3D.logger.addVariable("csz", csz);
    // putting phase values into array, following respectively for other 5 loops
    for (int i = csx; i < csx + grid + 1; i++) {
      xPhaseValues_Positive[i] = Math.abs((double) phase_img.getProcessor().getPixelValue(i, csy));
    }
    for (int i = csx; i > csx - grid - 1; i--) {
      xPhaseValues_Negative[i] = Math.abs((double) phase_img.getProcessor().getPixelValue(i, csy));
    }
    for (int j = csy; j < csy + grid + 1; j++) {
      yPhaseValues_Positive[j] = Math.abs((double) phase_img.getProcessor().getPixelValue(csx, j));
    }
    for (int j = csy; j > csy - grid - 1; j--) {
      yPhaseValues_Negative[j] = Math.abs((double) phase_img.getProcessor().getPixelValue(csx, j));
    }
    for (int k = csz; k < csz + grid + 1; k++) {
      phase_img.setSlice(k + 1);
      zPhaseValues_Positive[k] = Math.abs((double) phase_img.getProcessor().getPixelValue(csx, csy));
    }
    for (int k = csz; k > csz - grid - 1; k--) {
      phase_img.setSlice(k + 1);
      zPhaseValues_Negative[k] = Math.abs((double) phase_img.getProcessor().getPixelValue(csx, csy));
    }

    // Getting threshold (should be 50 unless user changes it)
    // M = M_pct;

    // negating all the values that are below the threshold (because phase values
    // can be random within the object)
    negateValues(xPhaseValues_Positive, M, csx, "x", true);
    negateValues(xPhaseValues_Negative, M, csx, "x", false);
    negateValues(yPhaseValues_Positive, M, csy, "y", true);
    negateValues(yPhaseValues_Negative, M, csy, "y", false);
    negateValues(zPhaseValues_Positive, M, csz, "z", true);
    negateValues(zPhaseValues_Negative, M, csz, "z", false);
    Calculate_Magnetic_Moment_3D.logger.addInfo("negated");

    // Putting equitorial plane and MRI field direction into new variables (since we
    // technically don't know which direction is going to be the MRI field
    // direction)
    // neglected will be MRI and PV will be equitorials
    double[] neglectedPVP;
    double[] neglectedPVN;
    int nC;
    // String neglectedAxis;
    switch (findFaultyAxis(grid - 2,
        xPhaseValues_Positive,
        xPhaseValues_Negative,
        yPhaseValues_Positive,
        yPhaseValues_Negative,
        zPhaseValues_Positive,
        zPhaseValues_Negative)) {
      case "x":
        neglectedPVP = xPhaseValues_Positive;
        neglectedPVN = xPhaseValues_Negative;
        neglectedAxis = "x";
        nC = csx;
        // PVP1 = yPhaseValues_Positive;
        // PVN1 = yPhaseValues_Negative;
        // axis1 = "y";
        // P1C = csy;
        // PVP2 = zPhaseValues_Positive;
        // PVN2 = zPhaseValues_Negative;
        // axis2 = "z";
        // P2C = csz;
        break;

      case "y":
        neglectedPVP = yPhaseValues_Positive;
        neglectedPVN = yPhaseValues_Negative;
        neglectedAxis = "y";
        nC = csy;
        // PVP1 = xPhaseValues_Positive;
        // PVN1 = xPhaseValues_Negative;
        // P1C = csx;
        // axis1 = "x";
        // PVP2 = zPhaseValues_Positive;
        // PVN2 = zPhaseValues_Negative;
        // P2C = csz;
        // axis2 = "z";
        break;

      default:
        neglectedPVP = zPhaseValues_Positive;
        neglectedPVN = zPhaseValues_Negative;
        neglectedAxis = "z";
        nC = csz;
        // PVP1 = xPhaseValues_Positive;
        // PVN1 = xPhaseValues_Negative;
        // P1C = csx;
        // axis1 = "x";
        // PVP2 = yPhaseValues_Positive;
        // PVN2 = yPhaseValues_Negative;
        // P2C = csy;
        // axis2 = "y";
        break;
    }

    Calculate_Magnetic_Moment_3D.logger.addVariable("neglectedAxis", neglectedAxis);
    // returning the estimated radius along the MRI field direction
    return findRadiiNeglected(nC, neglectedAxis, neglectedPVP, neglectedPVN);
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

  public void calcCenterL() {
    // ---------- begin to find Center_L, the minimum pixel value in the user-drawn
    // ROI

    // Method to find minimum pixel value in ROI
    float min = 0.0f;
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

    // ---------- end to find Center_L
  }

  public void calcCenterM() {
    // ---------- Begin to find Center_M
    center_m = new Triplet<Double>(0.0, 0.0, 0.0);

    // Making Center_M the small box to start
    center_m.set(0, (double) (roi_mag_belowM_xi + roi_mag_belowM_xi + roi_mag_belowM_Dx));
    center_m.set(1, (double) (roi_mag_belowM_yi + roi_mag_belowM_yi + roi_mag_belowM_Dy));
    center_m.set(2, (double) (roi_mag_belowM_zi + roi_mag_belowM_zi + roi_mag_belowM_Dz));

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

    center_m.set(0, center_m.get(0) / 2.0 + 0.5);
    center_m.set(1, center_m.get(1) / 2.0 + 0.5);
    center_m.set(2, center_m.get(2) / 2.0 + 0.5);

    Calculate_Magnetic_Moment_3D.logger.addVariable("center_m", center_m);

    // ---------- end to find Center_M
  }

  public void calcCenterS() {
    center_s = new Triplet<Double>(0.0, 0.0, 0.0);
    // ------ Begin summing planes in innerbox
    // Initializing arrays for sums of planes in small box
    double[] innerBox_sumOfXPlane = new double[roi_mag_belowM_xi + roi_mag_belowM_Dx + 1];
    double[] innerBox_sumOfYPlane = new double[roi_mag_belowM_yi + roi_mag_belowM_Dy + 1];
    double[] innerBox_sumOfZPlane = new double[roi_mag_belowM_zi + roi_mag_belowM_Dz + 1];

    Calculate_Magnetic_Moment_3D.logger.addVariable("rmb xi", roi_mag_belowM_xi);
    Calculate_Magnetic_Moment_3D.logger.addVariable("rmb yi", roi_mag_belowM_yi);
    Calculate_Magnetic_Moment_3D.logger.addVariable("rmb zi", roi_mag_belowM_zi);
    Calculate_Magnetic_Moment_3D.logger.addVariable("rmb dx", roi_mag_belowM_Dx);
    Calculate_Magnetic_Moment_3D.logger.addVariable("rmb dy", roi_mag_belowM_Dy);
    Calculate_Magnetic_Moment_3D.logger.addVariable("rmb dx", roi_mag_belowM_Dz);

    // Summing z plane and putting it into array
    for (int k = roi_mag_belowM_zi; k <= roi_mag_belowM_zi + roi_mag_belowM_Dz; k++) {
      mag_img.setSlice(k + 1);
      innerBox_sumOfZPlane[k] = 0;
      for (int i = roi_mag_belowM_xi; i <= roi_mag_belowM_xi + roi_mag_belowM_Dx; i++) {
        for (int j = roi_mag_belowM_yi; j <= roi_mag_belowM_yi + roi_mag_belowM_Dy; j++) {
          innerBox_sumOfZPlane[k] += mag_img.getProcessor().getPixelValue(i, j);
        }
      }
    }

    // Summing x plane and putting it into array
    for (int i = roi_mag_belowM_xi; i <= roi_mag_belowM_xi + roi_mag_belowM_Dx; i++) {
      innerBox_sumOfXPlane[i] = 0;
      for (int j = roi_mag_belowM_yi; j <= roi_mag_belowM_yi + roi_mag_belowM_Dy; j++) {
        for (int k = roi_mag_belowM_zi; k <= roi_mag_belowM_zi + roi_mag_belowM_Dz; k++) {
          mag_img.setSlice(k + 1);
          innerBox_sumOfXPlane[i] += mag_img.getProcessor().getPixelValue(i, j);
        }
      }
    }

    // Summing y plane and putting it into array
    for (int j = roi_mag_belowM_yi; j <= roi_mag_belowM_yi + roi_mag_belowM_Dy; j++) {
      innerBox_sumOfYPlane[j] = 0;
      for (int k = roi_mag_belowM_zi; k <= roi_mag_belowM_zi + roi_mag_belowM_Dz; k++) {
        mag_img.setSlice(k + 1);
        for (int i = roi_mag_belowM_xi; i <= roi_mag_belowM_xi + roi_mag_belowM_Dx; i++) {
          innerBox_sumOfYPlane[j] += mag_img.getProcessor().getPixelValue(i, j);
        }
      }
    }

    // for finding the minimum summed plane
    roi_mag_belowM_sumX = innerBox_sumOfXPlane[roi_mag_belowM_xi];
    roi_mag_belowM_sumY = innerBox_sumOfYPlane[roi_mag_belowM_yi];
    roi_mag_belowM_sumZ = innerBox_sumOfZPlane[roi_mag_belowM_zi];

    // ------ End summing planes in innerbox

    // ---------- begin to find Center_S

    center_s.set(0, (double) roi_mag_belowM_xi);
    center_s.set(1, (double) roi_mag_belowM_yi);
    center_s.set(2, (double) roi_mag_belowM_zi);

    // Finding minimum x plane and putting it into center_s.get(0)
    for (int i = roi_mag_belowM_xi; i <= roi_mag_belowM_xi + roi_mag_belowM_Dx; i++) {
      if (innerBox_sumOfXPlane[i] < roi_mag_belowM_sumX) {
        roi_mag_belowM_sumX = innerBox_sumOfXPlane[i];
        center_s.set(0, (double) i);
      }
    }

    // Finding minimum y plane and putting it into center_s.get(1)
    for (int j = roi_mag_belowM_yi; j <= roi_mag_belowM_yi + roi_mag_belowM_Dy; j++) {
      if (innerBox_sumOfYPlane[j] < roi_mag_belowM_sumY) {
        roi_mag_belowM_sumY = innerBox_sumOfYPlane[j];
        center_s.set(1, (double) j);
      }
    }

    // Finding minimum z plane and putting it into center_s.get(2)
    for (int k = roi_mag_belowM_zi; k <= roi_mag_belowM_zi + roi_mag_belowM_Dz; k++) {
      if (innerBox_sumOfZPlane[k] < roi_mag_belowM_sumZ) {
        roi_mag_belowM_sumZ = innerBox_sumOfZPlane[k];
        center_s.set(2, (double) k);
      }
    }
    // ---------- end to find Center_S
  }

  public Triplet<Double> centerL() {
    return center_l;
  }

  public Triplet<Double> centerM() {
    return center_m;
  }

  public Triplet<Double> centerS() {
    return center_s;
  }

  public String neglectedAxis() {
    return neglectedAxis;
  }

  public double RCenter() {
    return RCenter;
  }

  public Double m_R0() {
    return m_R0;
  }

  public void setCenterSX(double x) {
    center_s.set(0, x);
  }

  public void setCenterSY(double y) {
    center_s.set(1, y);
  }

  public void setCenterSZ(double z) {
    center_s.set(2, z);
  }

}
