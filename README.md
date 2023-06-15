# __Calculate Magnetic Moment 3D__

Calculate Magnetic Moment 3D (CMM3D) is an ImageJ plug-in developed as an implementation of __Magnetic moment quantifications of small spherical objects in MRI__ by Dr. Yu-Chung Cheng.

## __Developing with this repo__

If one wishes to develop/modify the code in `src/` and test within ImageJ, I have written build+run scripts in cmd for convenience.

Java 8 must be used. ImageJ is written in Java 8 and thus CISSCO is also. My C++ compiler is g++ 12.2.0 (MinGW-W64 x86_64-ucrt-posix-seh) and is also required.

To test CISSCO within ImageJ, perform the following steps:

1. Clone this repo to a folder
2. Create a folder called `ext`. It should be on the same level as `src`,`lib`,etc
3. Place a version of ImageJ within `ext`. The path should look like `CISSCO/ext/ImageJ/...`. Create a folder `CISSCO` within `ImageJ/plugins/`
4. Execute `build.*` from its location (I suggest using a terminal, `cd` to its location and execute it). ImageJ should open with the CISSCO plug-in included

## __Dependencies__

This plug-in uses MigLayout. It is a layout manager for Java Swing and is used to build the GUI. MigLayout produces flowing, grid based, absolute, grouped and docking layouts. It can also be configured with Eclipse so that you can edit the GUI interactively instead of by source code.

You can download MigLayout from https://www.miglayout.com/ but it is not required as it is already in this repository within ``lib/``

## __Usage__

- Step 1:

Here you can load the magnitude and phase images to study. The program will prompt you to first upload the magnitude, then the phase images. You will know if you successfully did this if you can see them displayed.

- Step 2:

You can draw an ROI using ImageJ's ROI tools around the object you wish to examine. The entire object should be inside this ROI. Once you do this you can click the "Estimate Center/Radii" button and the program will estimate the center coordinate of the object. You will notice the "-1" next to the z-coordinate. This is because the z coordinate will display the slice number, but the true coordinate is the displayed number minus one. This follows for all XYZ coordinates in this program.

- Step 3:

You can click the "Generate Subpixel Grid/Data" to generate the subpixel images associated with the object. You will see four images, XY and XZ magnitude images as well as XY and XZ phase images.
You will also notice the "Remove Bkg" button in this step. This will remove the background phase off of the phase matrices for the calculations. You can click this button at any time in the program as long as there is a background phase to remove.

- Step 4:

You can click the "Estimate Subpixel Center" after generating the subpixel images to estimate the center of the object with more precision. It should be something similar to the one calculated in Step 2.
Now you can change this center and click "Redraw Center" if you wish. You can also click "Verify Radii" to check the calculated R1, R2, and R3 from step 1. These radii should be at least 1 pixel apart.
You will also notice the "Plot Phase Profiles" button for each axis. This just displays a graph of the phase values along the respective axis through the center.
Once you find suitable R1, R2, and R3 values, you can click "Estimate Bkg & rho_0". This will calculate the spin density and the background phase again. The background phase should be updated with a more precise result.

- Step 5 (incomplete):

You can click "Calculate Magnetic Moment" after completing the above to calculate the magnetic moment. If you do not have values for e12 and e23 you will not see the uncertainty. You can calculate these values by uploading simulated images. These images are generated with the Python code. Loading these images will have the program automatically calculate e12 and e23 as well as the uncertainty. Make sure to input the appropriate SNR, which is defaulted to 1.
You can also calculate the sum of the real and imaginary values inside an inputted radius Ri and have it display.

- Step 6 (incomplete):

If you have a second set of images that have different echo times, you can use this to approximate deltaChi and a. Make sure to input the correct TE times as well as your B0, and R_deltaChi. The button to perform the calculations is currently unnamed.

- Step 7 (incomplete):

If you have a spin echo image, you can upload it to determine V0, rho_0SE, deltaChi and a. Once you upload the image you can input the center of this image as well as the two boxes for V1SE and V2SE. Clicking the "Estimate Object Radius From Spin Echo" button will perform these calculations.
