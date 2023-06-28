CISSCO uses JNI to work with both Java and C++. This is a documentation for which buttons require what variables/values to be passed from Java to C++. Since C++ does not have direct access to Java variables and vice versa, the only way of passing variables is through functions (which gets annoying). So we need to call these functions before calling the functions that actually "do the work".

## Load Magnitude and Phase Images
- none

## Estimate Center/Radii
- none

## Generate Subpixel Grid/Data
| C++         | Java        |
| ----------- | ----------- |
| **OnBnClickedGenerateSubpixel** | **jni.generateSubpixelArray()** |
| m_R0 | item.m_R0|
| m_SubPixels | grid|
| m_RCenter | gui->RCenter|
| RealNumbers[][][] | croppedRealNumbers3D[][][]|
| ImagNumbers[][][] | croppedImaginaryNumbers3D[][][]| 
| BackPhase | item.bkgPhase | 


## Remove Bkg
| C++         | Java        |
| ----------- | ----------- |
| **removeBGPhaseAndInterpolateVoxels(double)** | **jni.removeBackgroundPhase(double)** |
| m_SubPixels | grid |
| RealNumbers[][][] | croppedRealNumbers3D[][][] |
| ImagNumbers[][][] | croppedImaginaryNumbers3D[][][] |

However for this button, "Generate Subpixel Grid/Data" must be clicked first, which already passes the above values.


## Estimate Subpixel Center
| C++         | Java        |
| ----------- | ----------- |
| **OnBnClickedEstimatecenter()** | **jni.estimateSubpixelCenter()** |
| m_SubPixels | grid |
| smallBox_X,Y,Z | item.roi_mag_belowM_x,y,zi|
| smallBox_X,Y,ZSize | item.roi_mag_belowM_dx,y,z|
| m_RCenter | gui -> rc |
| m_R0 | item.m_R0 |
| centerL_x,y,z | item.CenterL() |
| centerM_x,y,z | item.CenterM() |
| centerS_x,y,z | gui -> ltf_rcx,y,z |
| m_CenterX,Y,Z2 | (int) Double.parseDouble(gui.ltf_rcx,y,z.getValue()) |


## Estimate Bkg & Spin Density
| C++         | Java        |
| ----------- | ----------- |
| __CalculateSpinDensity(double, ... ,double, double*, double*)__ | **jni.estBkgAndSpinDensity()** |
| m_CenterX,Y,Z | gui -> ltf_spx,y,z |
| m_R1,2,3 | gui -> ltf_r1,2,3 |
| BkgPhase | item.bkgPhase |
| m_MagMoment | gui -> ltf_magMom |

## Calculate Magnetic Moment
| C++         | Java        |
| ----------- | ----------- |
| **OnBnClickedCalcmagmoment()** | **jni.calculateMagneticMoment()** |
| m_R1,2,3 | gui -> ltf_r1,2,3 |
| m_CenterX,Y,Z | gui -> ltf_spx,y,z |
| m_SubPixels | grid |
| m_R0 | item.m_R0 | 
| BkgPhase | item.bkgPhase |
| m_RChi | gui -> ltf_RChi|
| m_B0 | gui -> ltf_B0|
| m_TE_first | gui -> ltf_TEFirst |
| m_SNR | gui -> ltf_snr |



## Misc
| C++         | Java        |
| ----------- | ----------- |