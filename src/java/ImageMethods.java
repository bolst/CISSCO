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
}
