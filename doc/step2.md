
# Step 2: 
Use magnitude images to estimate the center of the object at the pixel level. First, a user draws a rectangle which encloses the object. Second, our program adds intensities along each direction within the rectangle. When the sum of the intensities along each direction is minimal, that coordinate provides us the center of the object. This idea assumes that the intensities inside the object are zero, while the intensities outside the object are a constant. (This procedure may become insensitive if a box is too large. We are fine-tuning this procedure below.)

# Step 2:
a. The user draws a rectangle enclosing the spherical object displayed on one slice. The Java code will take the larger of the width and length of the rectangle and use that size to place a 3D box enclosing the object.

b. Record corner coordinates of the 3D box around the object. Find the averaged magnitude intensity based on the coordinates of 8 corners of the box.

c. Find the voxel with the lowest magnitude intensity inside the box. Label this voxel as Center_L.

d. Consider a fraction (around 50%) of the averaged magnitude intensity from Step 2b. Allow the user to change this fraction through the GUI. Identify all voxels with magnitude intensities lower than that fraction of the averaged intensity inside the box.

e. Determine the smallest box that contains all voxels identified from the previous Step 2d. Record the coordinates of 8 corners of this box. Label the center of this box as Center_M. If the difference between the two coordinates is an even number, the middle coordinate (or pixel) between the two corner coordinates is the coordinate for Center_M. However, if the difference between the two coordinates is an odd number, 0.5 pixel should be added to the middle coordinate (or pixel).

At this stage, we can see whether Center_M agrees well with the actual center of the object. Typically, due to the presence of Gibbs ringing, Center_L is not Center_M. In addition, we can further try the following:

f. At each coordinate along each axis, add all magnitude intensities inside the smallest box determined from the previous Step 2e. These additions are performed over a plane which is perpendicular to the axis. The minimal value among these sums along each axis leads to the coordinate of a center. Label this center as Center_S, which is not identical to Center_M.

Because the first element of an array in Java or C++ is labeled to be the zeroth element, the GUI display along the slice direction needs to subtract one. Thus we add “-1” next to the z-coordinate in Step 2.

The y-axis of images is inverted between ImageJ and SPIN. Here is an example. Assume that the dimension of the y-axis is 32, with indices running from 0 to 31. Also assume that we have y = 17.45 in ImageJ. Because the images shown in SPIN are inverted from images shown in ImageJ, we have y = 31 - 17.45 = 13.55 for SPIN.

g. Average phase values from the 8 corner voxels of the large box determined from Step 2a. Subtract this average phase value from phase images of all voxels inside the large box.

h. From Step 2d, assign phase values to be zero for those voxels with low magnitude intensities. Take absolute value of phase inside the large box determined from Step 2a. The procedures here are limited to this Step 2.

i. Because phase values along the MRI field direction are twice of the phase values on the equatorial plane at the same distance from the center of the spherical object, based on the three coordinates of Center_S, identify the MRI field direction, which is considered to be the MRI z-direction.

(Note: Find out the highest absolute phase value and its corresponding voxel index along all three axes through Center_S.)

j. Along each of the positive and negative MRI field directions, identify two adjacent voxels that have absolute phase values covering $\phi_{R\text{center}}$, which is chosen by the user and is 1 radian by default. Interpolate phase values between those two adjacent voxels such that at the distance away from Center_S the interpolated subvoxel has a phase value of roughly $\phi_{R\text{center}}$. Label the two distances away from Center_S along the positive and negative MRI field directions to be $r_{z+}$ and $r_{z-}
$, and determine both values at the subvoxel level.

(Notes: If we choose $R_\text{center}$ such that its corresponding phase value on the equatorial plane is $\phi_{R\text{center}}$, then

$$ p/(R_\text{center})^3 = \phi_{R\text{center}} $$
This leads to $p = (R_\text{center})^3  \phi_{R\text{center}}$
Along the z-direction, we have
$$ 2p/r^3 = \phi $$
Now, if we choose the same $\phi_{R\text{center}}$ along the field direction, then we have:
$$ \phi_{R\text{center}} = 2p/r^3 = 2\phi_{R\text{center}} (R_\text{center}^3) /r^3 $$

Then we still solve $ 2(R_\text{center}/r)^3 = 1 $ for Step.2j, except that the corresponding phase value is $\phi_{R\text{center}}$. The solution is $$ R_\text{center} = r/2^{1/3} $$

k. Update the coordinate of Center_S along the MRI field direction to be $(r_{z+} + r_{z-}
)/2$. Assign 
$$ 
R_\text{center} = \frac{|r_{z+} - r_{z-}|}{2\cdot 2^{1/3}}
$$

Step 2a: For the first time running the software, the initial background phase should be assigned to zero until it is calculated later.

Step 2b: How to estimate Rcenter with the phase value ($\phi_{R\text{center}}$) entered by a user in the GUI?

We first need to null phase values inside the sphere. We need to use magnitude images for this purpose. From averaging the 8 corners of the box that a user draws in the GUI, we know the averaged intensity value outside the sphere. Then consider that the magnitude intensity inside the sphere is to be less than 50% or 60%. Thus you can write a function to exclude pixels with magnitude intensities less than a specified value. The overall pixels you will exclude reflect the volume of the sphere. Use that volume and estimate the radius of the spherical object. This radius will be needed at the end of the GUI.