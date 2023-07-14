
/**
 * Class for handling multiple ROIs in ImageJ - specifically for Calculate_Magnetic_Moment_3D
 *
 * @author Nicholas Bolton (bolton21@uwindsor.ca)
 * @version 2.0
 */

import java.util.ArrayList;
import ij.WindowManager;
import ij.gui.OvalRoi;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.IJ;
import ij.ImagePlus;
import java.awt.Rectangle;

public class ROIS {

    private ArrayList<Roi> rois;
    private String imageTitle;

    public ROIS(String image_title) {
        rois = new ArrayList<Roi>();
        imageTitle = image_title;

        ImagePlus img = WindowManager.getImage(imageTitle);
        if (img != null)
            IJ.run(img, "Select None", "");
    }

    /**
     * Adds an ROI to the list
     * 
     * @param roi the ROI to add
     */
    public void addROI(Roi roi) {
        rois.add(roi);
    }

    /**
     * Adds a circle ROI to the list
     * 
     * @param x      the x coordinate of the center of the circle
     * @param y      the y coordinate of the center of the circle
     * @param radius
     */
    public void addCircleROI(int x, int y, double radius) {
        OvalRoi circle = new OvalRoi((double) x - radius, (double) y - radius, radius * 2.0, radius * 2.0);
        rois.add(circle);
    }

    /**
     * Adds a point ROI to the list
     * 
     * @param x the x coordinate of the point
     * @param y the y coordinate of the point
     */
    public void addPointROI(int x, int y) {
        PointRoi point = new PointRoi(x, y);
        // ShapeRoi p = new ShapeRoi(point);
        rois.add(point);
    }

    /**
     * Adds a rectangle ROI to the list
     * 
     * @param x      the x coordinate of the top-left point of the rectangle
     * @param y      the y coordinate of the top-left point of the rectangle
     * @param width  the width of the rectangle
     * @param height the height of the rectangle
     */
    public void addRectangle(int x, int y, int width, int height) {
        Roi roi = new Roi(new Rectangle(x, y, width, height));
        rois.add(roi);
    }

    /**
     * Displays the ROI on the image
     */
    public void displayROIS() {
        ImagePlus img = WindowManager.getImage(imageTitle);
        for (Roi r : rois) {
            img.setRoi(r);
            IJ.run(img, "Add Selection...", "");
        }
        IJ.run(img, "Select None", "");
    }

    /**
     * Clears all the ROIs on the image and from the list
     */
    public void clear() {
        ImagePlus img = WindowManager.getImage(imageTitle);
        IJ.run(img, "Remove Overlay", "");
        rois.clear();
    }
}
