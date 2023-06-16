
/**
 * Class for handling multiple ROIs in ImageJ - specifically for Calculate_Magnetic_Moment_3D
 *
 * @author Nicholas Bolton (bolton21@uwindsor.ca)
 * @version 1.0
 */

import java.util.HashMap;
import ij.gui.OvalRoi;
import ij.gui.PointRoi;
import ij.gui.Roi;
//import ij.gui.ShapeRoi;
import ij.plugin.frame.RoiManager;
import java.awt.Rectangle;

public class ROIS {

    private HashMap<String, Roi> rois;
    private RoiManager roiManager;
    private String axis;

    /**
     * Default constructor
     * 
     * @param dimension flag for what image this object will be used for: "MXY" for
     *                  magnitude XY, "MXZ" for magnitude XZ, "PXY" for phase XY,
     *                  "PXZ" for phase XZ, "MAG2XY" for V1SE XY, "MAG2XZ" for V1SE
     *                  XZ
     */
    public ROIS(String dimension) {
        rois = new HashMap<String, Roi>();
        roiManager = new RoiManager(true);
        axis = dimension;

        if (axis.compareTo("MXY") == 0) {
            Calculate_Magnetic_Moment_3D.subpixelMagImage.killRoi();
        }
        if (axis.compareTo("MXZ") == 0) {
            Calculate_Magnetic_Moment_3D.subpixelMagImageXZ.killRoi();
        }
        if (axis.compareTo("PXY") == 0) {
            Calculate_Magnetic_Moment_3D.subpixelPhaseImage.killRoi();
        }
        if (axis.compareTo("PXZ") == 0) {
            Calculate_Magnetic_Moment_3D.subpixelPhaseImageXZ.killRoi();
        }
        if (axis.compareTo("MAG2XY") == 0) {
            Calculate_Magnetic_Moment_3D.V1SE_XYImage.killRoi();
        }
        if (axis.compareTo("MAG2XZ") == 0) {
            Calculate_Magnetic_Moment_3D.V1SE_XZImage.killRoi();
        }
    }

    /**
     * Adds an ROI to the list
     * 
     * @param name the name of the ROI
     * @param roi  the ROI to add
     */
    public void addROI(String name, Roi roi) {
        rois.put(name, roi);
    }

    /**
     * Adds a circle ROI to the list
     * 
     * @param name   the name of the ROI
     * @param x      the x coordinate of the center of the circle
     * @param y      the y coordinate of the center of the circle
     * @param radius
     */
    public void addCircleROI(String name, int x, int y, double radius) {
        OvalRoi circle = new OvalRoi((double) x - radius, (double) y - radius, radius * 2.0, radius * 2.0);
        rois.put(name, circle);
    }

    /**
     * Adds a point ROI to the list
     * 
     * @param name the name of the ROI
     * @param x    the x coordinate of the point
     * @param y    the y coordinate of the point
     */
    public void addPointROI(String name, int x, int y) {
        PointRoi point = new PointRoi(x, y);
        // ShapeRoi p = new ShapeRoi(point);
        rois.put(name, point);
    }

    /**
     * Adds a rectangle ROI to the list
     * 
     * @param name   the name of the ROI
     * @param x      the x coordinate of the top-left point of the rectangle
     * @param y      the y coordinate of the top-left point of the rectangle
     * @param width  the width of the rectangle
     * @param height the height of the rectangle
     */
    public void addRectangle(String name, int x, int y, int width, int height) {
        Roi roi = new Roi(new Rectangle(x, y, width, height));
        rois.put(name, roi);
    }

    /**
     * Removes an ROI from the list
     * 
     * @param name the name of the ROI to remove
     */
    public void removeROI(String name) {
        rois.remove(name);
    }

    /**
     * Gets an ROI from the list
     * 
     * @param name the name of the ROI to get
     * @return the ROI
     */
    public Roi getROI(String name) {
        return rois.get(name);
    }

    /**
     * Displays the ROI on the image
     */
    public void displayROIS() {

        switch (axis) {
            case "MXY":
                for (Roi r : rois.values()) {
                    Calculate_Magnetic_Moment_3D.subpixelMagImage.setRoi(r);
                    r.setPosition(1);
                    roiManager.addRoi(r);
                    Calculate_Magnetic_Moment_3D.logger
                            .addVariable("MXY" + Calculate_Magnetic_Moment_3D.subpixelMagImage.getTitle(),
                                    r.toString());
                }
                roiManager.runCommand(Calculate_Magnetic_Moment_3D.subpixelMagImage, "Show All");
                break;

            case "MXZ":
                for (Roi r : rois.values()) {
                    Calculate_Magnetic_Moment_3D.subpixelMagImageXZ.setRoi(r);
                    r.setPosition(1);
                    roiManager.addRoi(r);
                    Calculate_Magnetic_Moment_3D.logger
                            .addVariable("MXZ" + Calculate_Magnetic_Moment_3D.subpixelMagImageXZ.getTitle(),
                                    r.toString());
                }
                roiManager.runCommand(Calculate_Magnetic_Moment_3D.subpixelMagImageXZ, "Show All");
                break;

            case "PXY":
                for (Roi r : rois.values()) {
                    Calculate_Magnetic_Moment_3D.subpixelPhaseImage.setRoi(r);
                    r.setPosition(1);
                    roiManager.addRoi(r);
                    Calculate_Magnetic_Moment_3D.logger
                            .addVariable("PXY" + Calculate_Magnetic_Moment_3D.subpixelPhaseImage.getTitle(),
                                    r.toString());
                }
                roiManager.runCommand(Calculate_Magnetic_Moment_3D.subpixelPhaseImage, "Show All");
                break;

            case "PXZ":
                for (Roi r : rois.values()) {
                    Calculate_Magnetic_Moment_3D.subpixelPhaseImageXZ.setRoi(r);
                    r.setPosition(1);
                    roiManager.addRoi(r);
                    Calculate_Magnetic_Moment_3D.logger
                            .addVariable("PXZ" + Calculate_Magnetic_Moment_3D.subpixelPhaseImageXZ.getTitle(),
                                    r.toString());
                }
                roiManager.runCommand(Calculate_Magnetic_Moment_3D.subpixelPhaseImageXZ, "Show All");
                break;

            case "MAG2XY":
                for (Roi r : rois.values()) {
                    Calculate_Magnetic_Moment_3D.V1SE_XYImage.setRoi(r);
                    r.setPosition(1);
                    roiManager.addRoi(r);
                }
                roiManager.runCommand(Calculate_Magnetic_Moment_3D.V1SE_XYImage, "Show All");
                break;

            case "MAG2XZ":
                for (Roi r : rois.values()) {
                    Calculate_Magnetic_Moment_3D.V1SE_XZImage.setRoi(r);
                    r.setPosition(1);
                    roiManager.addRoi(r);
                }
                roiManager.runCommand(Calculate_Magnetic_Moment_3D.V1SE_XZImage, "Show All");
                break;

            default:
                break;
        }

    }

    /**
     * Clears all the ROIs on the image and from the list
     */
    public void clear() {
        if (axis.compareTo("MXY") == 0) {
            Calculate_Magnetic_Moment_3D.subpixelMagImage.killRoi();
        }
        if (axis.compareTo("MXZ") == 0) {
            Calculate_Magnetic_Moment_3D.subpixelMagImageXZ.killRoi();
        }
        if (axis.compareTo("PXY") == 0) {
            Calculate_Magnetic_Moment_3D.subpixelPhaseImage.killRoi();
        }
        if (axis.compareTo("PXZ") == 0) {
            Calculate_Magnetic_Moment_3D.subpixelPhaseImageXZ.killRoi();
        }
        if (axis.compareTo("MAG2XY") == 0) {
            Calculate_Magnetic_Moment_3D.V1SE_XYImage.killRoi();
        }
        if (axis.compareTo("MAG2XZ") == 0) {
            Calculate_Magnetic_Moment_3D.V1SE_XZImage.killRoi();
        }
        rois.clear();
        roiManager.reset();
    }
}
