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

  public static float sumBoxCorners(ImagePlus image, Vec3<Integer> corner1, Vec3<Integer> corner2, int echoImageIndex) {
    float retval = 0f;
    int[] xs = new int[] { corner1.get(0), corner2.get(0) };
    int[] ys = new int[] { corner1.get(1), corner2.get(1) };
    int[] zs = new int[] { corner1.get(2), corner2.get(2) };
    for (int x : xs) {
      for (int y : ys) {
        for (int z : zs) {
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
    int rPos = 0;
    double maxPhasePos = 0.0;
    for (int r = 0; r < phaseValsPos.length; r++) {
      double phase = phaseValsPos[r];
      if (phase > maxPhasePos) {
        maxPhasePos = phase;
        rPos = r;
      }
    }

    int rNeg = 0;
    double maxPhaseNeg = 0.0;
    for (int r = 0; r < phaseValsNeg.length; r++) {
      double phase = phaseValsNeg[r];
      if (phase > maxPhaseNeg) {
        maxPhaseNeg = phase;
        rNeg = r;
      }
    }

    double smallestDistanceToOne = Double.POSITIVE_INFINITY;

    double eqPhasePos = maxPhasePos;
    for (int r = rPos; r < phaseValsPos.length; r++) {
      double phase = phaseValsPos[r];
      double dist = Math.abs(phase - 1);
      if (dist < smallestDistanceToOne) {
        smallestDistanceToOne = dist;
        eqPhasePos = phase;
        rPos = r;
      }
    }

    smallestDistanceToOne = Double.POSITIVE_INFINITY;
    double eqPhaseNeg = maxPhaseNeg;
    for (int r = rNeg; r < phaseValsNeg.length; r++) {
      double phase = phaseValsNeg[r];
      double dist = Math.abs(phase - 1);
      if (dist < smallestDistanceToOne) {
        smallestDistanceToOne = dist;
        eqPhaseNeg = phase;
        rNeg = r;
      }
    }

    double avgDistance = (rPos + rNeg) / 2.0;
    double avgEqPhase = (eqPhasePos + eqPhaseNeg) / 2.0;
    double RCenter = avgDistance * Math.cbrt(avgEqPhase / 2.0);
    RCenter = Math.max(RCenter, Calculate_Magnetic_Moment_3D.MIN_RCENTER);
    // new center to return
    double newCenter = mriAxisCenter;

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

  }

}
