public class Utilities {
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
    public static double interpolation(double x, double x1, double x2, double y1, double y2) {
        return y1 + ((x - x1) / (x2 - x1)) * (y2 - y1);
    }
}
