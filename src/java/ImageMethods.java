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

    // indexed xp,xn,yp,yn,zp,zn
    double[][] Ps = new double[][] { XP, XN, YP, YN, ZP, ZN };

    // to hold moment in each direction x,y,z
    double[] moments = new double[3];

    // for each direction X,Y,Z
    for (int iP = 0; iP < Ps.length; iP += 2) {
      // get phase arrays in both + and - directions
      double[] phasesPos = Ps[iP];
      double[] phasesNeg = Ps[iP + 1];

      double moment = 0.0;

      // at the first non-zero phase value, calculate the moment
      for (int r = 0; r < phasesPos.length; r++) {
        double phase = phasesPos[r];
        if (phase != 0) {
          moment += phase * Math.pow(r, 3);
          break;
        }
      }
      // repeat in negative direction
      for (int r = 0; r < phasesNeg.length; r++) {
        double phase = phasesNeg[r];
        if (phase != 0) {
          moment += phase * Math.pow(r, 3);
          break;
        }
      }

      Calculate_Magnetic_Moment_3D.logger.addVariable("moment " + String.valueOf((int) (iP / 2)), moment);

      // divide by 2 for average
      moments[(int) (iP / 2)] = moment / 2.0;
    }

    double xAvg = moments[0];
    double yAvg = moments[1];
    double zAvg = moments[2];

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
    double rPos = Calculate_Magnetic_Moment_3D.MIN_RCENTER;
    double maxPhasePos = 0.0;
    for (int r = 0; r < phaseValsPos.length; r++) {
      double phase = phaseValsPos[r];
      if (phase > maxPhasePos) {
        maxPhasePos = phase;
        rPos = r;
      }
    }

    double rNeg = Calculate_Magnetic_Moment_3D.MIN_RCENTER;
    double maxPhaseNeg = 0.0;
    for (int r = 0; r < phaseValsNeg.length; r++) {
      double phase = phaseValsNeg[r];
      if (phase > maxPhaseNeg) {
        maxPhaseNeg = phase;
        rNeg = r;
      }
    }

    // half distance between max phase values
    double halfMaxDistance = (rPos + rNeg) / 2.0;
    double avgMaxPhase = (maxPhasePos + maxPhaseNeg) / 2.0;
    double RCenter = halfMaxDistance * Math.cbrt(avgMaxPhase / 2.0);
    RCenter = Math.max(RCenter, Calculate_Magnetic_Moment_3D.MIN_RCENTER);
    // new center to return
    double newCenter = ((mriAxisCenter + rPos) + (mriAxisCenter - rNeg)) / 2.0;

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
      case X:
        // positive x direction
        if (direction) {
          for (int i = 0; i < values.length; i++) {
            if (Math
                .abs(ImageMethods.getVoxelValue(magImg, i + intCsx, intCsy, intCsz, echoImageIndex)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
        }
        // negative x direction
        else {
          for (int i = 0; i < values.length; i++) {
            if (Math
                .abs(ImageMethods.getVoxelValue(magImg, intCsx - i, intCsy, intCsz, echoImageIndex)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
        }
        break;
      case Y:
        // positive y direction
        if (direction) {
          for (int i = 0; i < values.length; i++) {
            if (Math
                .abs(ImageMethods.getVoxelValue(magImg, intCsx, i + intCsy, intCsz, echoImageIndex)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
        }
        // negative y direction
        else {
          for (int i = 0; i < values.length; i++) {
            if (Math
                .abs(ImageMethods.getVoxelValue(magImg, intCsx, intCsy - i, intCsz, echoImageIndex)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
        }
        break;
      case Z:
        // positive z direction
        if (direction) {
          for (int i = 0; i < values.length; i++) {
            if (Math
                .abs(ImageMethods.getVoxelValue(magImg, intCsx, intCsy, i + intCsz, echoImageIndex)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
        }
        // negative z direction
        else {
          for (int i = 0; i < values.length; i++) {
            if (Math
                .abs(ImageMethods.getVoxelValue(magImg, intCsx, intCsy, intCsz - i, echoImageIndex)) < maxMagValue) {
              values[i] = 0.0;
            }
          }
        }
        break;
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
