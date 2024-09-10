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

    Calculate_Magnetic_Moment_3D.logger.addInfo("HI");

    // for each direction
    for (int i = 0; i < nDims; i++) {
      // get direction
      double[] P = Ps[i];

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

    Calculate_Magnetic_Moment_3D.logger.addInfo("HI");

    // TODO: fix potential out of bounds
    double[] averages = new double[3];
    for (int i = 0; i < nDims; i += 2) {
      int j = i + 1;

      int r0p = r0s[i];
      int r0n = r0s[j];

      double pos0 = (r0p < Ps.length - 1) ? Ps[i][r0p + 1] : 0;
      double neg0 = (r0n < Ps.length - 1) ? Ps[j][r0n + 1] : 0;

      double pos = Math.abs(Ps[i][r0p] - pos0);
      double neg = Math.abs(Ps[j][r0n] - neg0);
      averages[(int) (i / 2)] = (pos + neg) / 2.0;
    }

    Calculate_Magnetic_Moment_3D.logger.addVariable("average x", averages[0]);
    Calculate_Magnetic_Moment_3D.logger.addVariable("average y", averages[1]);
    Calculate_Magnetic_Moment_3D.logger.addVariable("average z", averages[2]);

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
        if (Utilities.interpolation(phaseValue, phaseValsPos[i], phaseValsPos[i + 1], i, i + 1) > 1.6) {
          r_pos = Utilities.interpolation(phaseValue, phaseValsPos[i], phaseValsPos[i + 1], i, i + 1);
        }
        break;
      }
    }

    // Same thing as above just the opposite direction
    for (int i = 0; i < phaseValsNeg.length - 1; i++) {
      if ((phaseValsNeg[i] > phaseValue) && (phaseValsNeg[i + 1] < phaseValue)) {
        Calculate_Magnetic_Moment_3D.logger.addVariable("Found neg phase 1", phaseValsNeg[i]);
        Calculate_Magnetic_Moment_3D.logger.addVariable("Found neg phase 2", phaseValsNeg[i + 1]);
        if (Utilities.interpolation(phaseValue, phaseValsNeg[i], phaseValsNeg[i + 1], i, i + 1) > 1.6) {
          r_neg = Utilities.interpolation(phaseValue, phaseValsNeg[i], phaseValsNeg[i + 1], i, i + 1);
        }
        break;
      }
    }

    Calculate_Magnetic_Moment_3D.logger.addVariable("r_pos", r_pos);
    Calculate_Magnetic_Moment_3D.logger.addVariable("r_neg", r_neg);

    // RCenter to return
    double RCenter = Math.abs(r_pos + r_neg) / 2.0 / Math.cbrt(2);

    // new center to return
    double newCenter = ((mriAxisCenter + r_pos) + (mriAxisCenter - r_neg)) / 2.0;

    double[] retval = new double[] { RCenter, newCenter };

    return retval;
  }

}
