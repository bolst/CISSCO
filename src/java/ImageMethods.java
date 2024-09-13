import java.util.Arrays;
import ij.ImagePlus;

public class ImageMethods {

  // gets voxel value of a given image
  // specify an echo time if t > 0
  public static float getVoxelValue(ImagePlus image, int x, int y, int z, int t) {
    if (image.isHyperStack()) {
      int currFrame = t > 0 ? t : image.getFrame();
      final int T = image.getT();
      final int Z = image.getZ();

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

  public static float sumBoxCorners(ImagePlus image, Vec3<Integer> corner, Vec3<Integer> size, int echoImageIndex) {
    float retval = 0f;
    for (int x : new int[] { corner.get(0), size.get(0) }) {
      for (int y : new int[] { corner.get(1), size.get(1) }) {
        for (int z : new int[] { corner.get(2), size.get(2) }) {
          retval += ImageMethods.getVoxelValue(image, x, y, z, echoImageIndex);
        }
      }
    }
    return retval;
  }

  /**
   * Function to find MRI field direction. Does this by comparing the values in
   * each direction from the estimated center.
   * The two that follow similar trends are going to be on the equitorial plane,
   * while the other will be the MRI field direction.
   * 
   * @param XP phase values in positive x direction
   * @param XN phase values in negative x direction
   * @param YP phase values in positive y direction
   * @param YN phase values in negative y direction
   * @param ZP phase values in positive z direction
   * @param ZN phase values in negative z direction
   * @return Axis of MRI field direction
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
    double[][] Ps = new double[][] { XP, XN, YP, YN, ZP, ZN };

    double[] weightedAvgs = new double[] { 0, 0, 0 };
    // for each direction X,Y,Z
    for (int iP = 0; iP < Ps.length; iP += 2) {
      // get phase arrays in both directions
      double[] phasesPos = Ps[iP];
      double[] phasesNeg = Ps[iP + 1];

      double weightedSum = 0.0;
      int norm = 0;

      // weighted sum along positive axis
      for (int r = 0; r < phasesPos.length; r++) {
        double phase = phasesPos[r];
        if (phase != 0) {
          weightedSum += phase * Math.pow(r, 3);
          norm++;
        }
      }
      // weighted sum along negative axis
      for (int r = 0; r < phasesNeg.length; r++) {
        double phase = phasesNeg[r];
        if (phase != 0) {
          weightedSum += phase * Math.pow(r, 3);
          norm++;
        }
      }

      Calculate_Magnetic_Moment_3D.logger.addVariable("weighted sum " + String.valueOf((int) (iP / 2)), weightedSum);
      Calculate_Magnetic_Moment_3D.logger.addVariable("norm", norm);

      weightedAvgs[(int) (iP / 2)] = weightedSum / norm;
    }

    double xAvg = weightedAvgs[0];
    double yAvg = weightedAvgs[1];
    double zAvg = weightedAvgs[2];

    // max average is the corresponding direction
    Axis mri = xAvg > yAvg ? (xAvg > zAvg ? Axis.X : Axis.Z) : (yAvg > zAvg ? Axis.Y : Axis.Z);

    return mri;
  }

  /*
   * Function used to estimate radius on MRI field axis
   * Uses the theory that on this axis, phase values follow a 2/r^3 trend
   *
   * @return averaged radius in both directions
   */
  public static double[] calculateRCenter(int mriAxisCenter,
      Axis mriAxis,
      double[] phaseValsPos,
      double[] phaseValsNeg,
      double phaseValue) {

    // radius in each direction (1.6)
    double r_pos = Calculate_Magnetic_Moment_3D.MIN_RCENTER;
    double r_neg = Calculate_Magnetic_Moment_3D.MIN_RCENTER;

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
        double r0 = Utilities.interpolation(phaseValue, phaseValsPos[i], phaseValsPos[i + 1], i, i + 1);
        if (r0 > Calculate_Magnetic_Moment_3D.MIN_RCENTER) {
          r_pos = r0;
        }
        break;
      }
    }

    // Same thing as above just the opposite direction
    for (int i = 0; i < phaseValsNeg.length - 1; i++) {
      if ((phaseValsNeg[i] > phaseValue) && (phaseValsNeg[i + 1] < phaseValue)) {
        Calculate_Magnetic_Moment_3D.logger.addVariable("Found neg phase 1", phaseValsNeg[i]);
        Calculate_Magnetic_Moment_3D.logger.addVariable("Found neg phase 2", phaseValsNeg[i + 1]);
        double r0 = Utilities.interpolation(phaseValue, phaseValsNeg[i], phaseValsNeg[i + 1], i, i + 1);
        if (r0 > Calculate_Magnetic_Moment_3D.MIN_RCENTER) {
          r_neg = r0;
        }
        break;
      }
    }

    Calculate_Magnetic_Moment_3D.logger.addVariable("r_pos", r_pos);
    Calculate_Magnetic_Moment_3D.logger.addVariable("r_neg", r_neg);
    Calculate_Magnetic_Moment_3D.logger.addVariable("MRI axis center", mriAxisCenter);

    // RCenter to return
    double RCenter = Math.abs(r_pos + r_neg) / 2.0 / Math.cbrt(2);
    RCenter = Math.max(RCenter, Calculate_Magnetic_Moment_3D.MIN_RCENTER);

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

    // getting average of ROI
    double avgOfCorners = sumBoxCorners(magImg, roiCorner, roiSize, echoImageIndex) / 8.0;

    Calculate_Magnetic_Moment_3D.logger.addVariable("8 corner average", avgOfCorners);

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
          break;
        }
        // negative x direction
        else {
          for (int i = 0; i < values.length; i++) {
            if (Math.abs(ImageMethods.getVoxelValue(magImg, cx0 - i, intCsy, intCsz, echoImageIndex)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          break;
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
          break;
        }
        // negative y direction
        else {
          for (int i = 0; i < values.length; i++) {
            if (Math.abs(ImageMethods.getVoxelValue(magImg, intCsx, cy0 - i, intCsz, echoImageIndex)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          break;
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
          break;
        }
        // negative z direction
        else {
          for (int i = 0; i < values.length; i++) {
            if (Math.abs(ImageMethods.getVoxelValue(magImg, intCsx, intCsy, cz0 - i, echoImageIndex)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
          break;
        }
      default:
        break;
    }

    // now go back and make sure there are no nonzero values before any zero values
    // if there are, set them to zero
    boolean zeroFound = false;
    for (int i = values.length - 1; i >= 0; i--) {
      if (zeroFound) {
        values[i] = 0.0;
      } else if (values[i] == 0.0) {
        zeroFound = true;
      }
    }

    // also inspect the first 3 nonzero values and make sure they are all descending
    // if there are values that ascend, set them to 0
    // ex. [4,3,7,6,...] -> [0,0,7,6,...]

    // find first nonzero value
    int nonZeroIndex = 0;
    for (int i = 0; i < values.length; i++) {
      if (values[i] != 0.0) {
        nonZeroIndex = i;
        break;
      }
    }

    // go through first three nonzero values (or end of array)
    for (int i = nonZeroIndex; i < Math.min(nonZeroIndex + 3, values.length - 1); i++) {
      // if value is smaller than next, value is ascending, zero all values before and
      // including this one
      if (values[i] < values[i + 1]) {
        for (int j = nonZeroIndex; j <= i; j++) {
          values[j] = 0.0;
        }
      }
    }

  }

}
