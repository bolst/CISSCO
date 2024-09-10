import java.util.Arrays;
import ij.ImagePlus;

public class ImageMethods {

  // gets voxel value of a given image
  // specify an echo time if t > 0
  public static float getVoxelValue(ImagePlus image, int x, int y, int z, int t) {
    if (image.isHyperStack()) {
      int currFrame = t > 0 ? t : image.getFrame();
      int T = image.getT();
      int Z = image.getZ();

      image.setT(currFrame);
      image.setZ(z + 1);

      float retval = image.getProcessor().getPixelValue(x, y);

      image.setT(T);
      image.setZ(Z);

      return retval;

    } else {
      int z0 = image.getSlice();
      image.setSlice(z + 1);
      float retval = image.getProcessor().getPixelValue(x, y);
      image.setSlice(z0);
      return retval;
    }
  }

  /*
   * Function to find MRI field direction. Does this by comparing the values in
   * each direction from the estimated center.
   * The two that follow similar trends are going to be on the equitorial plane,
   * while the other will be the MRI field direction.
   *
   * @param accuracy The amount of times the code will iterate along each
   * direction
   *
   * @return A string value ("x", "y", or "z") of the MRI field direction
   */
  public static Axis calculateMRIAxis(double[] XP, double[] XN, double[] YP, double[] YN, double[] ZP, double[] ZN) {
    Calculate_Magnetic_Moment_3D.logger.addVariable("XP", Arrays.toString(XP));
    Calculate_Magnetic_Moment_3D.logger.addVariable("XN", Arrays.toString(XN));
    Calculate_Magnetic_Moment_3D.logger.addVariable("YP", Arrays.toString(YP));
    Calculate_Magnetic_Moment_3D.logger.addVariable("YN", Arrays.toString(YN));
    Calculate_Magnetic_Moment_3D.logger.addVariable("ZP", Arrays.toString(ZP));
    Calculate_Magnetic_Moment_3D.logger.addVariable("ZN", Arrays.toString(ZN));

    // to hold distances from center_s where phase is non-zero
    // indexed xp,xn,yp,yn,zp,zn
    final int nDims = 6;
    int[] r0s = new int[] { 0, 0, 0, 0, 0, 0 };
    double[][] Ps = new double[][] { XP, XN, YP, YN, ZP, ZN };

    // for each direction
    for (int i = 0; i < nDims; i++) {
      // get direction
      double[] P = Ps[i];

      // TODO: what if peak value is found at last index?

      // for each phase value in phase direction from center_s
      for (int r = 0; r < P.length; r++) {
        double phase = P[r];
        // find peak phase value
        if (phase != 0.0) {
          if (phase > P[r0s[i]]) {
            r0s[i] = r;
          }
        }
      }
    }

    Calculate_Magnetic_Moment_3D.logger.addVariable("dynamic box", Arrays.toString(r0s));

    // TODO: fix potential out of bounds
    double[] averages = new double[3];
    for (int i = 0; i < nDims; i += 2) {
      int j = i + 1;

      int r0p = r0s[i];
      int r0n = r0s[j];

      double pos0 = (r0p + 1 < Ps.length) ? Ps[i][r0p + 1] : 0;
      double neg0 = (r0n + 1 < Ps.length) ? Ps[j][r0n + 1] : 0;

      double pos = Math.abs(Ps[i][r0p] - pos0);
      double neg = Math.abs(Ps[j][r0n] - neg0);
      averages[(int) (i / 2)] = (pos + neg) / 2.0;
    }

    Calculate_Magnetic_Moment_3D.logger.addVariable("average x slope", averages[0]);
    Calculate_Magnetic_Moment_3D.logger.addVariable("average y slope", averages[1]);
    Calculate_Magnetic_Moment_3D.logger.addVariable("average z slope", averages[2]);

    double minAvg = Math.min(Math.min(averages[0], averages[1]), averages[2]);
    Axis retval = Axis.X;
    if (minAvg == averages[2])
      retval = Axis.Z;
    else if (minAvg == averages[1])
      retval = Axis.Y;
    else
      retval = Axis.X;

    Calculate_Magnetic_Moment_3D.logger.addVariable("mri axis", retval);
    return retval;
  }

  /*
   * Function used to estimate radius on MRI field axis
   * Uses the theory that on this axis, phase values follow a 2/r^3 trend
   *
   * @return averaged radius in both directions
   */
  public static double[] calculateRC(int mriAxisCenter,
      Axis mriAxis,
      double[] phaseValsPos,
      double[] phaseValsNeg,
      double phaseValue) {

    // radius in each direction
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

    for (int i = 0; i < phaseValsPos.length - 1; i++) {
      // this condition is if the current voxel is nested between two voxels greater
      // and less than the user-defined phase value
      if ((phaseValsPos[i] > phaseValue) && (phaseValsPos[i + 1] < phaseValue)) {
        Calculate_Magnetic_Moment_3D.logger.addVariable("Found pos phase 1", phaseValsPos[i]);
        Calculate_Magnetic_Moment_3D.logger.addVariable("Found pos phase 2", phaseValsPos[i + 1]);
        // this condition is if the interpolation value is greater than 1.6, because the
        // object radius has to be greater than 1.6
        double interpolatedPhase = Utilities.interpolation(phaseValue, phaseValsPos[i], phaseValsPos[i + 1], i, i + 1);
        if (interpolatedPhase > 1.6) {
          r_pos = interpolatedPhase;
        }
        break;
      }
    }

    // Same thing as above just the opposite direction
    for (int i = 0; i < phaseValsNeg.length - 1; i++) {
      if ((phaseValsNeg[i] > phaseValue) && (phaseValsNeg[i + 1] < phaseValue)) {
        Calculate_Magnetic_Moment_3D.logger.addVariable("Found neg phase 1", phaseValsNeg[i]);
        Calculate_Magnetic_Moment_3D.logger.addVariable("Found neg phase 2", phaseValsNeg[i + 1]);
        double interpolatedPhase = Utilities.interpolation(phaseValue, phaseValsNeg[i], phaseValsNeg[i + 1], i, i + 1);
        if (interpolatedPhase > 1.6) {
          r_neg = interpolatedPhase;
        }
        break;
      }
    }

    Calculate_Magnetic_Moment_3D.logger.addVariable("r_pos", r_pos);
    Calculate_Magnetic_Moment_3D.logger.addVariable("r_neg", r_neg);
    Calculate_Magnetic_Moment_3D.logger.addVariable("MRI axis center", mriAxisCenter);

    // RCenter to return
    double RCenter = Math.abs(r_pos + r_neg) / 2.0 / Math.cbrt(2);

    // new center to return
    double newCenter = ((mriAxisCenter + r_pos) + (mriAxisCenter - r_neg)) / 2.0;

    double[] retval = new double[] { RCenter, newCenter };

    return retval;
  }

  /*
   * Function to void all values below the param M.
   *
   * @param values array of values in one direction
   *
   * @param M the % to negate, default is 50
   *
   * @param r0 the center coordinate of the direction inputted
   *
   * @param axis The direction
   *
   * @param direction Flag to specify if function should iterate in the + or -
   * direction
   */
  public static void removeValuesBelow(
      double[] values,
      ImagePlus magImg,
      Vec3<Double> centerS,
      Vec3<Integer> roiCorner,
      Vec3<Integer> roiSize,
      int echoImageIndex,
      int M,
      Axis axis,
      boolean direction) {

    int roi_xi = roiCorner.get(0);
    int roi_yi = roiCorner.get(1);
    int roi_zi = roiCorner.get(2);
    int roi_dx = roiSize.get(0);
    int roi_dy = roiSize.get(1);
    int roi_dz = roiSize.get(2);

    // getting average of ROI
    double avgOfCorners = 0.0;
    for (int x : new int[] { roi_xi, roi_xi + roi_dx }) {
      for (int y : new int[] { roi_yi, roi_yi + roi_dy }) {
        for (int z : new int[] { roi_zi, roi_zi + roi_dz }) {
          avgOfCorners += ImageMethods.getVoxelValue(magImg, x, y, z, echoImageIndex);
        }
      }
    }
    avgOfCorners /= 8.0;

    // all values below this variable will be negated
    double maxMagValue = avgOfCorners * (double) M / 100.0;

    int intCsx = centerS.get(0).intValue();
    int intCsy = centerS.get(1).intValue();
    int intCsz = centerS.get(2).intValue();

    // setting all values in array below maxMagValue to 0
    switch (axis) {
      // x axis
      case X:
        int cx0 = centerS.get(0).intValue();
        // positive x direction
        if (direction) {
          for (int i = 0; i < values.length; i++) {
            if (Math.abs(
                ImageMethods.getVoxelValue(magImg, i + cx0, intCsy, intCsz, echoImageIndex)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        }
        // negative x direction
        else {
          for (int i = 0; i < values.length; i++) {
            if (Math.abs(ImageMethods.getVoxelValue(magImg, cx0 - i, intCsy, intCsz, echoImageIndex)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        }
        // y axis
      case Y:
        int cy0 = centerS.get(1).intValue();
        // positive y direction
        if (direction) {
          for (int i = 0; i < values.length; i++) {
            if (Math.abs(ImageMethods.getVoxelValue(magImg, intCsx, i + cy0, intCsz, echoImageIndex)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        }
        // negative y direction
        else {
          for (int i = 0; i < values.length; i++) {
            if (Math.abs(ImageMethods.getVoxelValue(magImg, intCsx, cy0 - i, intCsz, echoImageIndex)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        }
        // z axis
      case Z:
        int cz0 = centerS.get(2).intValue();
        // positive z direction
        if (direction) {
          for (int i = 0; i < values.length; i++) {
            if (Math.abs(ImageMethods.getVoxelValue(magImg, intCsx, intCsy, i + cz0, echoImageIndex)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        }
        // negative z direction
        else {
          for (int i = 0; i < values.length; i++) {
            if (Math.abs(ImageMethods.getVoxelValue(magImg, intCsx, intCsy, cz0 - i, echoImageIndex)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          return;
        }
      default:
        break;
    }
  }

}
