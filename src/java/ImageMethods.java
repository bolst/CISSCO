import ij.ImagePlus;

public class ImageMethods {

    public static float getVoxelValue(ImagePlus image, int x, int y, int z) {
        if (image.isHyperStack()) {
            int currFrame = image.getFrame();
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
}
