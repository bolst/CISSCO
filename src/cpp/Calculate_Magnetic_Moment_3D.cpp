// MagneticMomentDlg.cpp : implementation file
//

/*
The following file is used in parallel with the corresponding Java file. It uses JNI to exchange data between both files.
Steps for compiling and running the files are in the README.

You will see many functions below prefixed "Java_JNIMethods_". These are the JNI functions that will be called from Java.
You can see the return type in the function declaration as well, just like any other function.
Only primitive data types can be passed between programs (int, double, char, etc). Arrays can be passed too but are a bit tricky since
C++ technically does not have arrays but pointers, unlike Java. So the firstlevel, secondlevel, and thirdlevel functions assist with this.

*/

using namespace std;

#include <cmath>
#include <limits>
#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <jni.h>
#include "JNIMethods.h"

#define PI 3.1415926535
#define MAX_SUBPIXEL_DIM 640

#include <ctime>

int m_SubPixels;
double m_R0, m_R1, m_R2, m_R3, m_radians, m_RCenter, m_CenterX, m_CenterY, m_CenterZ, m_CenterX2, m_CenterY2, m_CenterZ2, m_CenterX3, m_CenterY3, m_CenterZ3;
float ***RealNumbers, ***ImagNumbers;
float ***RealNumbers_S5, ***ImagNumbers_S5;
float **subPhaseMatrix, **subPhaseMatrixXZ, **subMagMatrix, **subMagMatrixXZ;
float ***SimRealNumbers;
int m_RCenterPhase;
int OBcount;
// int Xfirst, Yfirst, Zfirst;
int halfreal, Zhalfreal, halfdisplay, displaydim, realdim, Zrealdim, subpixeldisplay, subpixelreal, Zsubpixelreal;
double BackPhase;
double ZoomedX, ZoomedY, ZoomedZ, lastValueSlice;
// double REALXfirst, REALYfirst, REALSlice;
// double Xinitial, Yinitial, Zinitial;
string magFileName, phaseFileName, errorMessage, errorMessage_magMom, errorMessage_sums;
bool m_Rcentercheck;
vector<vector<vector<float>>> SubpixelRealMatrix3D;
vector<vector<vector<float>>> SubpixelImagMatrix3D;
vector<vector<vector<float>>> SubpixelRealMatrix3D_S5;
vector<vector<vector<float>>> SubpixelImagMatrix3D_S5;
vector<vector<vector<float>>> SubpixelSimulatedRealMatrix3D;
vector<vector<float>> SubpixelPhaseMatrix;
vector<vector<float>> SubpixelMagMatrix;
vector<vector<float>> SubpixelPhaseMatrixXZ;
vector<vector<float>> SubpixelMagMatrixXZ;
vector<vector<vector<float>>> SmallReal(realdim, vector<vector<float>>(realdim, vector<float>(Zrealdim, 0)));
vector<vector<vector<float>>> SmallImag(realdim, vector<vector<float>>(realdim, vector<float>(Zrealdim, 0)));
float ***tempReal_BG, ***tempImag_BG;
int smallBox_X, smallBox_Y, smallBox_Z, smallBox_XSize, smallBox_YSize, smallBox_ZSize;
double centerL_x, centerL_y, centerL_z, centerM_x, centerM_y, centerM_z, centerS_x, centerS_y, centerS_z;
double subCenter;
double m_MagMoment, m_p_first, m_p_last, m_RInnerPhase, m_RMiddlePhase, m_ROuterPhase, m_rho, m_p0, m_BackPhase, m_RChi, m_B0, m_TE_first, m_TE_last, m_g, m_a, m_Chi, m_e12, m_e23, m_Uncertainty, m_Si, m_Si2;

double m_resx = 1.0, m_resy = 1.0, m_resz = 1.0;

double m_RInnerFrom = 2.5;
double m_RMiddleFrom = 0.9;
double m_ROuterFrom = 0.2;
double m_SNR, m_SNRoc;
double rho, BkgPhase;
double m_Ri;

int counter = 0;

int Xup, Yup;
int halfCONSTmat = 26;
// int halfsubCONSTmat = halfCONSTmat * 10;
int CONSTmat = halfCONSTmat * 2;
int subCONSTmat = CONSTmat * 10;
// subCONSTmat left in for handling error messages

// int halfCONSTmat3D = 15;
// int halfsubCONSTmat3D = halfCONSTmat3D * 10;
// int CONSTmat3D = halfCONSTmat3D * 2;
// int subCONSTmat3D = CONSTmat3D * 10;

// firstlevel, secondlevel, and thirdlevel are helper functions used for converting a jobjectArray to an nD array of doubles so that it can be easily used in the functions above

// Use this for converting 1D arrays
float *thirdLevel(JNIEnv *env, jfloatArray arr)
{
    jsize len = env->GetArrayLength(arr);
    float *ret = new float[len];
    env->GetFloatArrayRegion(arr, 0, len, ret);
    return ret;
}

// Use this for converting 2D arrays
float **secondLevel(JNIEnv *env, jobjectArray arr)
{
    jsize len = env->GetArrayLength(arr);
    float **ret = new float *[len];
    for (int i = 0; i < len; i++)
    {
        jobject item = env->GetObjectArrayElement(arr, i);
        ret[i] = thirdLevel(env, (jfloatArray)item);
        env->DeleteLocalRef(item);
    }
    return ret;
}

// Use this for converting 3D arrays
float ***firstLevel(JNIEnv *env, jobjectArray arr)
{
    jsize len = env->GetArrayLength(arr);
    float ***ret = new float **[len];
    for (int i = 0; i < len; i++)
    {
        jobject item = env->GetObjectArrayElement(arr, i);
        ret[i] = secondLevel(env, (jobjectArray)item);
        env->DeleteLocalRef(item);
    }
    return ret;
}

// These JNI functions are used to set variables from Java to C++, so Java is passing the data through these functions
/******************************************************************************************************************/
/******************************************************************************************************************/
/**********************************************    BEGIN    *******************************************************/
/**********************************************JNI FUNCTIONS*******************************************************/
/**********************************************             *******************************************************/
/******************************************************************************************************************/
/******************************************************************************************************************/

JNIEXPORT void JNICALL Java_JNIMethods_passGenSubpixelValues(JNIEnv *env, jobject thisObj,
                                                             jdouble jm_R0,
                                                             jint jm_SubPixels,
                                                             jdouble jm_RCenter,
                                                             jdouble jm_BackPhase)
{
    m_R0 = jm_R0;
    m_SubPixels = jm_SubPixels;
    m_RCenter = jm_RCenter;
    m_BackPhase = jm_BackPhase;
}

JNIEXPORT void JNICALL Java_JNIMethods_passCalcSubCenterValues(JNIEnv *env, jobject thisObj,
                                                               jint jsmallx, jint jsmalldx, jint jsmally, jint jsmalldy, jint jsmallz, jint jsmalldz,
                                                               jdouble jm_RCenter, jdouble jm_R0,
                                                               jdouble jclx, jdouble jcly, jdouble jclz,
                                                               jdouble jcmx, jdouble jcmy, jdouble jcmz,
                                                               jdouble jcsx, jdouble jcsy, jdouble jcsz,
                                                               jdouble jcx2, jdouble jcy2, jdouble jcz2)
{
    smallBox_X = jsmallx;
    smallBox_XSize = jsmalldx;
    smallBox_Y = jsmally;
    smallBox_YSize = jsmalldy;
    smallBox_Z = jsmallz;
    smallBox_ZSize = jsmalldz;

    m_RCenter = jm_RCenter;
    m_R0 = jm_R0;

    centerL_x = jclx;
    centerL_y = jcly;
    centerL_z = jclz;

    centerM_x = jcmx;
    centerM_y = jcmy;
    centerM_z = jcmz;

    centerS_x = jcsx;
    centerS_y = jcsy;
    centerS_z = jcsz;

    m_CenterX2 = jcx2;
    m_CenterY2 = jcy2;
    m_CenterZ2 = jcz2;

    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_passMagMomValues(JNIEnv *env, jobject thisObj,
                                                        jdouble jm_r1, jdouble jm_r2, jdouble jm_r3,
                                                        jdouble jm_csx, jdouble jm_csy, jdouble jm_csz,
                                                        jdouble jm_R0, jdouble jm_bkg,
                                                        jdouble jm_RChi, jdouble jm_B0, jdouble jm_TEFirst, jdouble jm_snr,
                                                        jdouble jm_in, jdouble jm_mid, jdouble jm_out,
                                                        jdouble jm_eps12, jdouble jm_eps23)
{
    m_R1 = jm_r1;
    m_R2 = jm_r2;
    m_R3 = jm_r3;

    m_CenterX = jm_csx;
    m_CenterY = jm_csy;
    m_CenterZ = jm_csz;

    m_R0 = jm_R0;
    BkgPhase = jm_bkg;

    m_RChi = jm_RChi;
    m_B0 = jm_B0;
    m_TE_first = jm_TEFirst;
    m_SNR = jm_snr;

    m_RInnerPhase = jm_in;
    m_RMiddlePhase = jm_mid;
    m_ROuterPhase = jm_out;

    m_e12 = jm_eps12;
    m_e23 = jm_eps23;

    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_passSpinDensValues(JNIEnv *env, jobject thisObj, jdouble jm_cx, jdouble jm_cy, jdouble jm_cz,
                                                          jdouble jm_r1, jdouble jm_r2, jdouble jm_r3,
                                                          jdouble jm_bkgphase, jdouble jm_magmom)
{
    m_CenterX = jm_cx;
    m_CenterY = jm_cy;
    m_CenterZ = jm_cz;

    m_R1 = jm_r1;
    m_R2 = jm_r2;
    m_R3 = jm_r3;

    BkgPhase = jm_bkgphase;
    m_MagMoment = jm_magmom;

    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_passSumValues(JNIEnv *env, jobject thisObj, jdouble jm_ri, jdouble jm_cx, jdouble jm_cy, jdouble jm_cz)
{
    m_Ri = jm_ri;

    m_CenterX = jm_cx;
    m_CenterY = jm_cy;
    m_CenterZ = jm_cz;

    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_setmVariables(JNIEnv *env, jobject thisObj, jint nm_SubPixels, jdouble nm_R0, jdouble nm_RCenter, jdouble nm_CenterX2, jdouble nm_CenterY2, jdouble nm_CenterZ2, jdouble nmphasevalue)
{
    m_SubPixels = nm_SubPixels;
    m_R0 = nm_R0;
    m_RCenter = nm_RCenter;
    m_CenterX2 = nm_CenterX2;
    m_CenterY2 = nm_CenterY2;
    m_CenterZ2 = nm_CenterZ2;
    m_RCenterPhase = nmphasevalue;
    // m_MagMoment = m_RCenterPhase * pow(m_RCenter, 3);
    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_setMagMoment(JNIEnv *env, jobject thisObj, jdouble nmoment)
{
    m_MagMoment = nmoment;
    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_setBackPhase(JNIEnv *env, jobject thisObj, jdouble nBackPhase)
{
    BackPhase = nBackPhase;
    BkgPhase = nBackPhase;
    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_setRealImagNumbers(JNIEnv *env, jobject thisObj, jobjectArray _real, jobjectArray _imag)
{
    /*
    RealNumbers = new float **[(int)(m_CenterX2 + m_R0) + 1];
    ImagNumbers = new float **[(int)(m_CenterX2 + m_R0) + 1];
    for (int i = 0; i <= (int)(m_CenterX2 + m_R0); i++)
    {
        RealNumbers[i] = new float *[(int)(m_CenterY2 + m_R0) + 1];
        ImagNumbers[i] = new float *[(int)(m_CenterY2 + m_R0) + 1];
        for (int j = 0; j <= (int)(m_CenterY2 + m_R0); j++)
        {
            RealNumbers[i][j] = new float[(int)(m_CenterZ2 + m_R0) + 1];
            ImagNumbers[i][j] = new float[(int)(m_CenterZ2 + m_R0) + 1];
        }
    }
    */

    RealNumbers = new float **[subpixeldisplay + 1];
    ImagNumbers = new float **[subpixeldisplay + 1];

    for (int i = 0; i <= subpixeldisplay; i++)
    {
        RealNumbers[i] = new float *[subpixeldisplay + 1];
        ImagNumbers[i] = new float *[subpixeldisplay + 1];
        for (int j = 0; j <= subpixeldisplay; j++)
        {
            RealNumbers[i][j] = new float[subpixeldisplay + 1];
            ImagNumbers[i][j] = new float[subpixeldisplay + 1];
        }
    }

    RealNumbers = firstLevel(env, _real);
    ImagNumbers = firstLevel(env, _imag);

    // --------------- Setting up matrices for when the background phase is estimated for the second time

    RealNumbers_S5 = new float **[subpixeldisplay + 1];
    ImagNumbers_S5 = new float **[subpixeldisplay + 1];

    for (int i = 0; i <= subpixeldisplay; i++)
    {
        RealNumbers_S5[i] = new float *[subpixeldisplay + 1];
        ImagNumbers_S5[i] = new float *[subpixeldisplay + 1];
        for (int j = 0; j <= subpixeldisplay; j++)
        {
            RealNumbers_S5[i][j] = new float[subpixeldisplay + 1];
            ImagNumbers_S5[i][j] = new float[subpixeldisplay + 1];
        }
    }

    RealNumbers_S5 = firstLevel(env, _real);
    ImagNumbers_S5 = firstLevel(env, _imag);

    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_setXYZ(JNIEnv *env, jobject thisObj, jdouble nx, jdouble ny, jdouble nz)
{
    m_CenterX = nx;
    m_CenterY = ny;
    m_CenterZ = nz;

    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_setPhaseXYMatrix(JNIEnv *env, jobject thisObj, jobjectArray jarr)
{
    subPhaseMatrix = new float *[subpixeldisplay];
    for (int i = 0; i < subpixeldisplay; i++)
    {
        subPhaseMatrix[i] = new float[subpixeldisplay];
    }

    subPhaseMatrix = secondLevel(env, jarr);

    SubpixelPhaseMatrix.clear();
    SubpixelPhaseMatrix.resize(subpixeldisplay, vector<float>(subpixeldisplay, 0));

    for (int i = 0; i < subpixeldisplay; i++)
    {
        for (int j = 0; j < subpixeldisplay; j++)
        {
            SubpixelPhaseMatrix[j][i] = subPhaseMatrix[j][i];
        }
    }

    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_setPhaseXZMatrix(JNIEnv *env, jobject thisObj, jobjectArray jarr)
{
    subPhaseMatrixXZ = new float *[subpixeldisplay];
    for (int i = 0; i < subpixeldisplay; i++)
    {
        subPhaseMatrixXZ[i] = new float[subpixeldisplay];
    }

    subPhaseMatrixXZ = secondLevel(env, jarr);

    SubpixelPhaseMatrixXZ.clear();
    SubpixelPhaseMatrixXZ.resize(subpixeldisplay, vector<float>(subpixeldisplay, 0));

    for (int i = 0; i < subpixeldisplay; i++)
    {
        for (int j = 0; j < subpixeldisplay; j++)
        {
            SubpixelPhaseMatrixXZ[j][i] = subPhaseMatrixXZ[j][i];
        }
    }

    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_setMagXYMatrix(JNIEnv *env, jobject thisObj, jobjectArray jarr)
{
    subMagMatrix = new float *[subpixeldisplay];
    for (int i = 0; i < subpixeldisplay; i++)
    {
        subMagMatrix[i] = new float[subpixeldisplay];
    }

    subMagMatrix = secondLevel(env, jarr);

    SubpixelMagMatrix.clear();
    SubpixelMagMatrix.resize(subpixeldisplay, vector<float>(subpixeldisplay, 0));

    for (int i = 0; i < subpixeldisplay; i++)
    {
        for (int j = 0; j < subpixeldisplay; j++)
        {
            SubpixelMagMatrix[j][i] = subMagMatrix[j][i];
        }
    }

    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_setMagXZMatrix(JNIEnv *env, jobject thisObj, jobjectArray jarr)
{
    subMagMatrixXZ = new float *[subpixeldisplay];
    for (int i = 0; i < subpixeldisplay; i++)
    {
        subMagMatrixXZ[i] = new float[subpixeldisplay];
    }

    subMagMatrixXZ = secondLevel(env, jarr);

    SubpixelMagMatrixXZ.clear();
    SubpixelMagMatrixXZ.resize(subpixeldisplay, vector<float>(subpixeldisplay, 0));

    for (int i = 0; i < subpixeldisplay; i++)
    {
        for (int j = 0; j < subpixeldisplay; j++)
        {
            SubpixelMagMatrixXZ[j][i] = subMagMatrixXZ[j][i];
        }
    }

    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_setmR123(JNIEnv *env, jobject thisObj, jdouble nm_R1, jdouble nm_R2, jdouble nm_R3)
{
    m_R1 = nm_R1;
    m_R2 = nm_R2;
    m_R3 = nm_R3;
    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_setR123PhaseCalc(JNIEnv *env, jobject thisObj, jdouble pr1, jdouble pr2, jdouble pr3)
{
    m_ROuterPhase = pr1;
    m_RMiddlePhase = pr2;
    m_RInnerPhase = pr3;
    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_setR123PhaseActual(JNIEnv *env, jobject thisObj, jdouble pr1, jdouble pr2, jdouble pr3)
{
    m_R1 = pr1;
    m_R2 = pr2;
    m_R3 = pr3;
}

JNIEXPORT void JNICALL Java_JNIMethods_setSmallBox(JNIEnv *env, jobject thisObj, jint xi, jint yi, jint zi, jint xsize, jint ysize, jint zsize)
{

    smallBox_X = xi;
    smallBox_Y = yi;
    smallBox_Z = zi;
    smallBox_XSize = xsize;
    smallBox_YSize = ysize;
    smallBox_ZSize = zsize;

    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_setCenterL(JNIEnv *env, jobject thisObj, jdouble x, jdouble y, jdouble z)
{
    centerL_x = x;
    centerL_y = y;
    centerL_z = z;

    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_setCenterM(JNIEnv *env, jobject thisObj, jdouble x, jdouble y, jdouble z)
{
    centerM_x = x;
    centerM_y = y;
    centerM_z = z;

    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_setCenterS(JNIEnv *env, jobject thisObj, jdouble x, jdouble y, jdouble z)
{
    centerS_x = x;
    centerS_y = y;
    centerS_z = z;

    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_setMagMomentVariables(JNIEnv *env, jobject thisObj, jdouble snr, jdouble e12, jdouble e23, jdouble b0, jdouble te, jdouble rchi)
{
    m_SNR = snr;
    m_e12 = e12;
    m_e23 = e23;
    m_B0 = b0;
    m_TE_first = te;
    m_RChi = rchi;
    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_setRi(JNIEnv *env, jobject thisObj, jdouble mri)
{
    m_Ri = mri;
    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_setStep6Variables(JNIEnv *env, jobject thisObj, jdouble jtfirst, jdouble jtlast, jdouble jb0, jdouble jrx)
{
    m_TE_first = jtfirst;
    m_TE_last = jtlast;
    m_B0 = jb0;
    m_RChi = jrx;

    m_p_first = m_p_last * m_TE_first / m_TE_last;
}

JNIEXPORT void JNICALL Java_JNIMethods_setSimulatedMatrices(JNIEnv *env, jobject thisObj, jobjectArray jreal, jint jsize)
{
    SimRealNumbers = new float **[jsize];

    for (int i = 0; i < jsize; i++)
    {
        SimRealNumbers[i] = new float *[jsize];

        for (int j = 0; j < jsize; j++)
        {
            SimRealNumbers[i][j] = new float[jsize];
        }
    }

    SimRealNumbers = firstLevel(env, jreal);

    return;
}

/******************************************************************************************************************/
/******************************************************************************************************************/
/**********************************************     END     *******************************************************/
/**********************************************JNI FUNCTIONS*******************************************************/
/**********************************************             *******************************************************/
/******************************************************************************************************************/
/******************************************************************************************************************/

// Helper function to convert original coordinates to subpixel coordinates
double pixelToSubpixel(double coordinate, int axisFlag)
{
    // Calculated by taking total size of cropped image and dividing it by 2 - this will give the center which is also the estimated center of the image in step 2

    subCenter = (2 * m_R0 + 1) * (10.0 / 2.0);
    double subpixelCoordinate = 0.0;

    // if x axis
    if (axisFlag == 0)
    {
        subpixelCoordinate = subCenter + (coordinate - m_CenterX2) * 10.0;
    }
    // if y axis
    if (axisFlag == 1)
    {
        subpixelCoordinate = subCenter + (coordinate - m_CenterY2) * 10.0;
    }
    // if z axis
    if (axisFlag == 2)
    {
        subpixelCoordinate = subCenter + (coordinate - m_CenterZ2) * 10.0;
    }

    return subpixelCoordinate;
}

// Helper function to convert subpixel coordinates to original coordinates
double subpixelToPixel(double coordinate, int axisFlag)
{
    subCenter = (2 * m_R0 + 1) * (10.0 / 2.0);
    double voxelCoordinate = 0.0;

    // if x axis
    if (axisFlag == 0)
    {
        voxelCoordinate = ((coordinate - subCenter) / 10.0) + m_CenterX2;
    }
    // if y axis
    if (axisFlag == 1)
    {
        voxelCoordinate = ((coordinate - subCenter) / 10.0) + m_CenterY2;
    }
    // if z axis
    if (axisFlag == 2)
    {
        voxelCoordinate = ((coordinate - subCenter) / 10.0) + m_CenterZ2;
    }

    return voxelCoordinate;
}
// Function to remove BG phase and interpolate real and imag matrices
void removeBGPhaseAndInterpolateVoxels(double bPhase)
{

    // ---------- begin to remove background phase  ---------- Do NOT change without discussion

    // double tmpreal, tmpimag, tmpPhase = BkgPhase;
    float tmpreal, tmpimag;
    double tmpPhase = bPhase;
    int Nfinal = m_SubPixels;
    int ratio = Nfinal * Nfinal * Nfinal;
    int m, k, t, p, i, j, pfrom, pto, ifrom, ito, jfrom, jto;

    for (m = 0; m < Zrealdim; m++)
    {
        pfrom = m * Nfinal;
        pto = pfrom + Nfinal;

        for (k = 0; k < realdim; k++)
        {
            ifrom = k * Nfinal;
            ito = ifrom + Nfinal;

            for (t = 0; t < realdim; t++)
            {
                jfrom = t * Nfinal;
                jto = jfrom + Nfinal;
                /*
                tmpreal = SmallReal[t][k][m];
                tmpimag = SmallImag[t][k][m];
                SmallReal[t][k][m] = tmpreal * cos(tmpPhase) + tmpimag * sin(tmpPhase);
                SmallImag[t][k][m] = tmpimag * cos(tmpPhase) - tmpreal * sin(tmpPhase);
                tmpreal = (SmallReal[t][k][m]) / ratio;
                tmpimag = (SmallImag[t][k][m]) / ratio;
                */
                tmpreal = RealNumbers[t][k][m];
                tmpimag = ImagNumbers[t][k][m];
                RealNumbers[t][k][m] = tmpreal * cos(tmpPhase) + tmpimag * sin(tmpPhase);
                ImagNumbers[t][k][m] = tmpimag * cos(tmpPhase) - tmpreal * sin(tmpPhase);
                tmpreal = (RealNumbers[t][k][m]) / ratio;
                tmpimag = (ImagNumbers[t][k][m]) / ratio;

                for (p = pfrom; p < pto; p++)
                {
                    for (i = ifrom; i < ito; i++)
                    {
                        for (j = jfrom; j < jto; j++)
                        {
                            SubpixelRealMatrix3D[j][i][p] = tmpreal;
                            SubpixelImagMatrix3D[j][i][p] = tmpimag;
                        }
                    }
                }
            }
        }
    }
    // ---------- end to remove background phase

    return;
}

// Same as above function, but for the simulated images. Simulated images have no background phase so no need to include
void interpolateVoxels_SIM(int msize)
{

    // ---------- begin to interpolate voxels for simmed images  ---------- Do NOT change without discussion

    // double tmpreal, tmpimag, tmpPhase = BkgPhase;
    float tmpreal, tmpimag;
    int Nfinal = m_SubPixels;
    int ratio = Nfinal * Nfinal * Nfinal;
    int m, k, t, p, i, j, pfrom, pto, ifrom, ito, jfrom, jto;

    for (m = 0; m < msize; m++)
    {
        pfrom = m * Nfinal;
        pto = pfrom + Nfinal;

        for (k = 0; k < msize; k++)
        {
            ifrom = k * Nfinal;
            ito = ifrom + Nfinal;

            for (t = 0; t < msize; t++)
            {
                jfrom = t * Nfinal;
                jto = jfrom + Nfinal;

                tmpreal = SimRealNumbers[t][k][m];
                SimRealNumbers[t][k][m] = tmpreal;
                tmpreal = (SimRealNumbers[t][k][m]) / ratio;

                for (p = pfrom; p < pto; p++)
                {
                    for (i = ifrom; i < ito; i++)
                    {
                        for (j = jfrom; j < jto; j++)
                        {
                            SubpixelSimulatedRealMatrix3D[j][i][p] = tmpreal;
                        }
                    }
                }
            }
        }
    }
    // ---------- end to interpolate voxels for simmed images

    return;
}

/* This deals with the same values as removeBGPhaseAndInterpolateVoxels(double bPhase) but the background phase is not removed.
Only interpolation happens here, this is so that when the background phase is calculated again in step 5 it does not perform calculations
on data with an already removed background phase.
*/
void interpolateVoxels_S5()
{

    // ---------- begin to interpolate voxels for estimation of background phase  ---------- Do NOT change without discussion

    // double tmpreal, tmpimag, tmpPhase = BkgPhase;
    float tmpreal, tmpimag;
    int Nfinal = m_SubPixels;
    int ratio = Nfinal * Nfinal * Nfinal;
    int m, k, t, p, i, j, pfrom, pto, ifrom, ito, jfrom, jto;

    for (m = 0; m < Zrealdim; m++)
    {
        pfrom = m * Nfinal;
        pto = pfrom + Nfinal;

        for (k = 0; k < realdim; k++)
        {
            ifrom = k * Nfinal;
            ito = ifrom + Nfinal;

            for (t = 0; t < realdim; t++)
            {
                jfrom = t * Nfinal;
                jto = jfrom + Nfinal;

                tmpreal = RealNumbers_S5[t][k][m];
                tmpimag = ImagNumbers_S5[t][k][m];
                RealNumbers_S5[t][k][m] = tmpreal;
                ImagNumbers_S5[t][k][m] = tmpimag;
                tmpreal = (RealNumbers_S5[t][k][m]) / ratio;
                tmpimag = (ImagNumbers_S5[t][k][m]) / ratio;

                for (p = pfrom; p < pto; p++)
                {
                    for (i = ifrom; i < ito; i++)
                    {
                        for (j = jfrom; j < jto; j++)
                        {
                            SubpixelRealMatrix3D_S5[j][i][p] = tmpreal;
                            SubpixelImagMatrix3D_S5[j][i][p] = tmpimag;
                        }
                    }
                }
            }
        }
    }
    // ---------- end to interpolate voxels for estimation of background phase

    return;
}

/*
This is the heavily modified code to generate the subpixel images. We generate the images on the Java side but handle
the interpolation and background removal of the real and imaginary numbers here. These real and imaginary matrices
are used in further calculations.
*/
void OnBnClickedGenerateSubpixel()
{
    // Define pDoc. Get title to check if its subpixel or not.
    // CMDIFrameWnd *pFrame = (CMDIFrameWnd *)AfxGetApp()->m_pMainWnd;
    // CMDIChildWnd *pChild = (CMDIChildWnd *)pFrame->GetActiveFrame();
    // CHipoView *pView = (CHipoView *)pChild->GetActiveView();
    // CHipoDoc *pDoc = (CHipoDoc *)pView->GetDocument();
    // CString TITLE = pDoc->GetTitle();
    // m_R0=12;

    errorMessage = "";
    halfdisplay = m_R0;
    halfreal = m_R0;
    Zhalfreal = m_R0;
    displaydim = 2 * m_R0;
    realdim = 2 * m_R0;
    Zrealdim = 2 * Zhalfreal;
    subpixeldisplay = 2 * m_R0 * m_SubPixels;
    subpixelreal = 2 * m_R0 * m_SubPixels;
    Zsubpixelreal = 2 * Zhalfreal * m_SubPixels;

    double RCenter = m_RCenter;
    int Nfinal = m_SubPixels;

    // Xfirst = m_CenterX2;
    // Yfirst = m_CenterY2;
    // Zfirst = m_CenterZ2;

    // vector<vector<vector<float>>> SmallImag(realdim, vector<vector<float>>(realdim, vector<float>(Zrealdim, 0)));
    // vector<vector<vector<float>>> SmallReal(realdim, vector<vector<float>>(realdim, vector<float>(Zrealdim, 0)));
    SmallReal.clear();
    SmallReal.resize(realdim, vector<vector<float>>(realdim, vector<float>(Zrealdim, 0)));
    SmallImag.clear();
    SmallImag.resize(realdim, vector<vector<float>>(realdim, vector<float>(Zrealdim, 0)));
    SubpixelRealMatrix3D.clear();
    SubpixelRealMatrix3D.resize(subpixelreal, vector<vector<float>>(subpixelreal, vector<float>(Zsubpixelreal, 0)));
    SubpixelImagMatrix3D.clear();
    SubpixelImagMatrix3D.resize(subpixelreal, vector<vector<float>>(subpixelreal, vector<float>(Zsubpixelreal, 0)));
    SubpixelRealMatrix3D_S5.clear();
    SubpixelRealMatrix3D_S5.resize(subpixelreal, vector<vector<float>>(subpixelreal, vector<float>(Zsubpixelreal, 0)));
    SubpixelImagMatrix3D_S5.clear();
    SubpixelImagMatrix3D_S5.resize(subpixelreal, vector<vector<float>>(subpixelreal, vector<float>(Zsubpixelreal, 0)));

    for (int k = 0; k < Zrealdim; k++)
    {
        for (int i = 0; i < realdim; i++)
        {
            for (int j = 0; j < realdim; j++)
            {
                // SmallReal[j][i][k] = (double)RealNumbers[(Zfirst - Zhalfreal + k) * dimx * dimy + (Yfirst - halfreal + i) * dimx + (Xfirst - halfreal + j)];
                // SmallImag[j][i][k] = (double)ImagNumbers[(Zfirst - Zhalfreal + k) * dimx * dimy + (Yfirst - halfreal + i) * dimx + (Xfirst - halfreal + j)];
                // SmallReal[j][i][k] = RealNumbers[(int)m_CenterX2 - halfreal + j][(int)m_CenterY2 - halfreal + i][(int)m_CenterZ2 - Zhalfreal + k];
                // SmallImag[j][i][k] = ImagNumbers[(int)m_CenterX2 - halfreal + j][(int)m_CenterY2 - halfreal + i][(int)m_CenterZ2 - Zhalfreal + k];
                SmallReal[j][i][k] = RealNumbers[j][i][k];
                SmallImag[j][i][k] = ImagNumbers[j][i][k];
            }
        }
    }

    removeBGPhaseAndInterpolateVoxels(BackPhase);

    // interpolating same matrices but not removing background phase - for calculating bkg phase in step 5
    interpolateVoxels_S5();

    return;

    // The remaining lines under this function have been commented out

    /*
        // make sure we have a center radius
        if (RCenter == 0)
        {
            //::MessageBox(NULL, "Must estimate center first", "Message", MB_OK);
            errorMessage = "Must estimate center first";
            return;
        }
    if ((phaseFileName != "PhaseSubpixel.dcm") && (magFileName != "MagSubpixel.dcm")
        fileNameFlag == 1)
    {
        // Get Real Numbers of entire dataset and Phase values. Store in 1D array. dim globally defined.
        //  GetRealNumbers3D is only called when you generate subpixel data from scratch.
        // GetRealNumbers3D(RealNumbers, PhaseValues, ImagNumbers, MagValues, dimx, dimy, dimz, 0);
                if ((RealNumbers == NULL))
                {
                    //::MessageBox(NULL, "Size of data matrix is 0. Check to make sure image is part of an entire 3D dataset.", "Message", MB_OK);
                    errorMessage = "Size of data matrix is 0. Check to make sure image is part of an entire 3D dataset.";
                    return;
                }

                if ((Xfirst + halfdisplay > dimx) || (Yfirst + halfdisplay > dimy))
                {
                    //::MessageBox(NULL, "Object is too close to edge of image", "Message", MB_OK);
                    errorMessage = "Object is too close to edge of image (+)";
                    return;
                }
                if ((Xfirst - halfdisplay / 2 < 0) || (Yfirst - halfdisplay / 2 < 0))
                {
                    //::MessageBox(NULL, "Object is too close to edge of image", "Message", MB_OK);
                    errorMessage = "Object is too close to edge of image";
                    return;
                }

                // Create Subpixel for first zoom********************************************************************

                // Phase arrays for display. Real Number arrays for center and magnetic moment calculations
                // The reason Paul clears and resizes was in case we had to regenerate the subpixel data by clicking the button again.
                vector<vector<float>> SmallPhase(displaydim, vector<float>(displaydim, 0));
                vector<vector<float>> SmallMag(displaydim, vector<float>(displaydim, 0));
                SmallPhase.clear();
                SmallPhase.resize(displaydim, vector<float>(displaydim, 0));
                SmallMag.clear();
                SmallMag.resize(displaydim, vector<float>(displaydim, 0));
                SubpixelPhaseMatrix.clear();
                SubpixelPhaseMatrix.resize(subpixeldisplay, vector<float>(subpixeldisplay, 0));
                SubpixelMagMatrix.clear();
                SubpixelMagMatrix.resize(subpixeldisplay, vector<float>(subpixeldisplay, 0));

                int i, j, k, m, t, p;
                double ratio;

                vector<vector<float>> SmallPhaseXZ(displaydim, vector<float>(displaydim, 0));
                vector<vector<float>> SmallMagXZ(displaydim, vector<float>(displaydim, 0));
                SmallPhaseXZ.clear();
                SmallPhaseXZ.resize(displaydim, vector<float>(displaydim, 0));
                SmallMagXZ.clear();
                SmallMagXZ.resize(displaydim, vector<float>(displaydim, 0));
                SubpixelPhaseMatrixXZ.clear();
                SubpixelPhaseMatrixXZ.resize(subpixeldisplay, vector<float>(subpixeldisplay, 0));
                SubpixelMagMatrixXZ.clear();
                SubpixelMagMatrixXZ.resize(subpixeldisplay, vector<float>(subpixeldisplay, 0));

                // vector<vector<vector<float>>> SmallReal(realdim,vector<vector<float>>(realdim,vector<float>(Zrealdim,0)));
                // vector<vector<vector<float>>> SmallImag(realdim,vector<vector<float>>(realdim,vector<float>(Zrealdim,0)));
                SmallReal.clear();
                SmallReal.resize(realdim, vector<vector<float>>(realdim, vector<float>(Zrealdim, 0)));
                SmallImag.clear();
                SmallImag.resize(realdim, vector<vector<float>>(realdim, vector<float>(Zrealdim, 0)));
                SubpixelRealMatrix3D.clear();
                SubpixelRealMatrix3D.resize(subpixelreal, vector<vector<float>>(subpixelreal, vector<float>(Zsubpixelreal, 0)));
                SubpixelImagMatrix3D.clear();
                SubpixelImagMatrix3D.resize(subpixelreal, vector<vector<float>>(subpixelreal, vector<float>(Zsubpixelreal, 0)));
                ratio = Nfinal * Nfinal * Nfinal;

                // initializing 2d arrays
                spMagMatrix = new float *[MAX_SUBPIXEL_DIM];
                for (int i = 0; i < MAX_SUBPIXEL_DIM; i++)
                {
                    spMagMatrix[i] = new float[MAX_SUBPIXEL_DIM];
                }

                spMagMatrixXZ = new float *[MAX_SUBPIXEL_DIM];
                for (int i = 0; i < MAX_SUBPIXEL_DIM; i++)
                {
                    spMagMatrixXZ[i] = new float[MAX_SUBPIXEL_DIM];
                }

                spPhaseMatrix = new float *[MAX_SUBPIXEL_DIM];
                for (int i = 0; i < MAX_SUBPIXEL_DIM; i++)
                {
                    spPhaseMatrix[i] = new float[MAX_SUBPIXEL_DIM];
                }

                spPhaseMatrixXZ = new float *[MAX_SUBPIXEL_DIM];
                for (int i = 0; i < MAX_SUBPIXEL_DIM; i++)
                {
                    spPhaseMatrixXZ[i] = new float[MAX_SUBPIXEL_DIM];
                }

                // initializing 3d arrays

                spRealMatrix3D = new float **[MAX_SUBPIXEL_DIM];
                for (int i = 0; i < MAX_SUBPIXEL_DIM; i++)
                {
                    spRealMatrix3D[i] = new float *[MAX_SUBPIXEL_DIM];
                    for (int j = 0; j < MAX_SUBPIXEL_DIM; j++)
                    {
                        spRealMatrix3D[i][j] = new float[MAX_SUBPIXEL_DIM];
                    }
                }

                spImagMatrix3D = new float **[MAX_SUBPIXEL_DIM];
                for (int i = 0; i < MAX_SUBPIXEL_DIM; i++)
                {
                    spImagMatrix3D[i] = new float *[MAX_SUBPIXEL_DIM];
                    for (int j = 0; j < MAX_SUBPIXEL_DIM; j++)
                    {
                        spImagMatrix3D[i][j] = new float[MAX_SUBPIXEL_DIM];
                    }
                }

                for (i = 0; i < displaydim; i++)
                {
                    for (j = 0; j < displaydim; j++)
                    {
                        // SmallPhase[j][i] = PhaseValues[Zfirst * dimx * dimy + (Yfirst - halfdisplay + i) * dimx + (Xfirst - halfdisplay + j)];
                        SmallPhase[j][i] = PhaseValues[Xfirst - halfdisplay + j][Yfirst - halfdisplay + i][Zfirst];
                    }
                }

                for (i = 0; i < displaydim; i++)
                {
                    for (j = 0; j < displaydim; j++)
                    {
                        // SmallMag[j][i] = MagValues[Zfirst * dimx * dimy + (Yfirst - halfdisplay + i) * dimx + (Xfirst - halfdisplay + j)];
                        SmallMag[j][i] = MagValues[Xfirst - halfdisplay + j][Yfirst - halfdisplay + i][Zfirst];
                    }
                }

                for (i = 0; i < displaydim; i++)
                {
                    for (j = 0; j < displaydim; j++)
                    {
                        // SmallPhaseXZ[j][i] = PhaseValues[(Zfirst - halfdisplay + i) * dimx * dimy + Yfirst * dimx + (Xfirst - halfdisplay + j)];
                        SmallPhaseXZ[j][i] = PhaseValues[Xfirst - halfdisplay + j][Yfirst][Zfirst - halfdisplay + i];
                    }
                }

                for (i = 0; i < displaydim; i++)
                {
                    for (j = 0; j < displaydim; j++)
                    {
                        // SmallMagXZ[j][i] = MagValues[(Zfirst - halfdisplay + i) * dimx * dimy + Yfirst * dimx + (Xfirst - halfdisplay + j)];
                        SmallMagXZ[j][i] = MagValues[Xfirst - halfdisplay + j][Yfirst][Zfirst - halfdisplay + i];
                    }
                }

        for (k = 0; k < Zrealdim; k++)
        {
            for (i = 0; i < realdim; i++)
            {
                for (j = 0; j < realdim; j++)
                {
                    SmallReal[j][i][k] = (double)RealNumbers[(Zfirst - Zhalfreal + k) * dimx * dimy + (Yfirst - halfreal + i) * dimx + (Xfirst - halfreal + j)];
                    SmallImag[j][i][k] = (double)ImagNumbers[(Zfirst - Zhalfreal + k) * dimx * dimy + (Yfirst - halfreal + i) * dimx + (Xfirst - halfreal + j)];
                }
            }
        }
                for (m = 0; m < displaydim; m++)
                {
                    for (k = 0; k < displaydim; k++)
                    {
                        for (i = m * Nfinal; i < (m + 1) * Nfinal; i++)
                        {
                            for (j = k * Nfinal; j < (k + 1) * Nfinal; j++)
                            {
                                SubpixelPhaseMatrix[j][i] = (SmallPhase[k][m]);
                                spPhaseMatrix[j][i] = SubpixelPhaseMatrix[j][i];
                            }
                        }
                    }
                }

                for (m = 0; m < displaydim; m++)
                {
                    for (k = 0; k < displaydim; k++)
                    {
                        for (i = m * Nfinal; i < (m + 1) * Nfinal; i++)
                        {
                            for (j = k * Nfinal; j < (k + 1) * Nfinal; j++)
                            {
                                SubpixelMagMatrix[j][i] = (SmallMag[k][m]);
                                spMagMatrix[j][i] = SubpixelMagMatrix[j][i];
                            }
                        }
                    }
                }

                for (m = 0; m < displaydim; m++)
                {
                    for (k = 0; k < displaydim; k++)
                    {
                        for (i = m * Nfinal; i < (m + 1) * Nfinal; i++)
                        {
                            for (j = k * Nfinal; j < (k + 1) * Nfinal; j++)
                            {
                                SubpixelPhaseMatrixXZ[j][i] = (SmallPhaseXZ[k][m]);
                                spPhaseMatrixXZ[j][i] = SubpixelPhaseMatrixXZ[j][i];
                            }
                        }
                    }
                }

                for (m = 0; m < displaydim; m++)
                {
                    for (k = 0; k < displaydim; k++)
                    {
                        for (i = m * Nfinal; i < (m + 1) * Nfinal; i++)
                        {
                            for (j = k * Nfinal; j < (k + 1) * Nfinal; j++)
                            {
                                SubpixelMagMatrixXZ[j][i] = (SmallMagXZ[k][m]);
                                spMagMatrixXZ[j][i] = SubpixelMagMatrixXZ[j][i];
                            }
                        }
                    }
                }
        */
    /*
            double val, ival;
            int pfrom, pto, ifrom, ito, jfrom, jto;

            for (m = 0; m < Zrealdim; m++)
            {
                for (k = 0; k < realdim; k++)
                {
                    for (t = 0; t < realdim; t++)
                    {
                        val = (SmallReal[t][k][m]) / ratio;
                        ival = (SmallImag[t][k][m]) / ratio;
                        pfrom = m * Nfinal;
                        pto = (m + 1) * Nfinal;
                        ifrom = k * Nfinal;
                        ito = (k + 1) * Nfinal;
                        jfrom = t * Nfinal;
                        jto = (t + 1) * Nfinal;

                        for (p = pfrom; p < pto; p++)
                        {
                            for (i = ifrom; i < ito; i++)
                            {
                                for (j = jfrom; j < jto; j++)
                                {
                                    SubpixelRealMatrix3D[j][i][p] = val;
                                    SubpixelImagMatrix3D[j][i][p] = ival;
                                    spRealMatrix3D[j][i][p] = SubpixelRealMatrix3D[j][i][p];
                                    spImagMatrix3D[j][i][p] = SubpixelImagMatrix3D[j][i][p];
                                }
                            }
                        }
                    }
                }
            }
    */
    // SubPixelImage(halfdisplay);
    // SubPixelMagImage(halfdisplay);

    // issubpix = TRUE;

    // CMDIFrameWnd *pFrame = (CMDIFrameWnd *)AfxGetApp()->m_pMainWnd;
    // CMDIChildWnd *pChild = (CMDIChildWnd *)pFrame->GetActiveFrame();
    // CHipoView *pView = (CHipoView *)pChild->GetActiveView();

    /*
            if (m_Rcentercheck)
            {

                CCtrlPoint RCStartPoint, RCEndPoint;

                RCStartPoint.x = (subpixeldisplay / 2 + 5) - (m_RCenter * m_SubPixels) + 5;
                RCStartPoint.y = (subpixeldisplay / 2 + 5) - (m_RCenter * m_SubPixels) + 5;
                RCEndPoint.x = (subpixeldisplay / 2 + 5) + (m_RCenter * m_SubPixels) + 5;
                RCEndPoint.y = (subpixeldisplay / 2 + 5) + (m_RCenter * m_SubPixels) + 5;

                pView->mROIContainer.RemoveAll();

                pView->psr->Remove();
                pView->psr->m_stPropertyOfROI.m_iROIType = ROI_ELLIPSE;
                pView->psr->AddTail(RCStartPoint);
                pView->psr->AddTail(RCEndPoint);
                pView->mROIContainer.AddTail(pView->psr);

                pView->Invalidate(false);

                pFrame->MDINext();
                pChild = (CMDIChildWnd *)pFrame->GetActiveFrame();
                pView = (CHipoView *)pChild->GetActiveView();

                pView->mROIContainer.RemoveAll();

                pView->psr->Remove();
                pView->psr->m_stPropertyOfROI.m_iROIType = ROI_ELLIPSE;
                pView->psr->AddTail(RCStartPoint);
                pView->psr->AddTail(RCEndPoint);
                pView->mROIContainer.AddTail(pView->psr);

                pView->Invalidate(false);
            }

}
else
{
    //::MessageBox(NULL, "Subpixel grid already generated.", "Message", MB_OK);
    errorMessage = "Subpixel grid already generated.";
    return;
}
*/
}

double SumSphericalMask(int radius, int Scan1, int Scan2, int Scan3, vector<vector<vector<int>>> CubeMask)
{
    int diameter = 2 * radius;
    int ktmp, itmp, jtmp;
    double sum = 0, sumi = 0;

    Scan1 = Scan1 - radius;
    Scan2 = Scan2 - radius;
    Scan3 = Scan3 - radius;

    for (int k = 0; k <= diameter; k++)
    {
        ktmp = Scan3 + k;
        for (int i = 0; i <= diameter; i++)
        {
            itmp = Scan2 + i;
            for (int j = 0; j <= diameter; j++)
            {
                jtmp = Scan1 + j;
                sumi = SubpixelRealMatrix3D[jtmp][itmp][ktmp] * CubeMask[j][i][k];
                sum += sumi;
            }
        }
    }
    return sum;
}

// ---------- begin GetSum to be used in Amoeba ---------- Do NOT change without discussion

void GetSum(vector<vector<double>> &p, vector<double> &psum)
{
    int i, j;
    double sum;

    int mpts = p.size();
    int ndim = p[0].size();
    for (j = 0; j < ndim; j++)
    {
        for (sum = 0.0, i = 0; i < mpts; i++)
            sum += p[i][j];
        psum[j] = sum;
    }
}
// ---------- end GetSum

// ---------- begin to use Amotry ---------- Do NOT change without discussion
// This function is from Numerical Recipes with changes in the argument, checking in-bound and ytry below.

double amotry(vector<vector<double>> &p, vector<double> &y, vector<double> &psum,
              const int ihi, const double fac, int &radiusCircle, double &GRID, int constgrid, vector<vector<vector<int>>> CubeMask)
{
    int j;
    double fac1, fac2, ytry;
    // MagMomentFunctions *Functions = NULL;
    // CString messagea = TEXT("");
    /*
        double Xmin = subCenter + (smallBox_X - m_CenterX2) * 10.0;
        double Ymin = subCenter + (smallBox_Y - m_CenterY2) * 10.0;
        double Zmin = subCenter + (smallBox_Z - m_CenterZ2) * 10.0;
        double Xmax = subCenter + (smallBox_X + smallBox_XSize - m_CenterX2) * 10.0;
        double Ymax = subCenter + (smallBox_Y + smallBox_YSize - m_CenterY2) * 10.0;
        double Zmax = subCenter + (smallBox_Z + smallBox_ZSize - m_CenterZ2) * 10.0;
    */

    double Xmin = pixelToSubpixel((double)(smallBox_X), 0);
    double Xmax = pixelToSubpixel((double)(smallBox_X + smallBox_XSize), 0);
    double Ymin = pixelToSubpixel((double)(smallBox_Y), 1);
    double Ymax = pixelToSubpixel((double)(smallBox_Y + smallBox_YSize), 1);
    double Zmin = pixelToSubpixel((double)(smallBox_Z), 2);
    double Zmax = pixelToSubpixel((double)(smallBox_Z + smallBox_ZSize), 2);

    int ndim = p[0].size();
    vector<double> ptry(ndim);
    fac1 = (1.0 - fac) / ndim;
    fac2 = fac1 - fac;
    for (j = 0; j < ndim; j++)
        ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2;

    // check to keep in bounds and restrict to integer

    /*
    if (ptry[0] < (Xinitial - (GRID / 2)))
    {
        ptry[0] = Xinitial - (GRID / 2) + 1;
        if (constgrid == 0)
            GRID = GRID + 1;
        OBcount = OBcount + 1;
        ptry[0] = (int)ptry[0];
    }
    if (ptry[0] > (Xinitial + (GRID / 2)))
    {
        ptry[0] = Xinitial + (GRID / 2) - 1;
        if (constgrid == 0)
            GRID = GRID + 1;
        OBcount = OBcount + 1;
        ptry[0] = (int)ptry[0];
    }

    if (ptry[1] < (Yinitial - (GRID / 2)))
    {
        ptry[1] = Yinitial - (GRID / 2) + 1;
        if (constgrid == 0)
            GRID = GRID + 1;
        OBcount = OBcount + 1;
        ptry[1] = (int)ptry[1];
    }
    if (ptry[1] > (Yinitial + (GRID / 2)))
    {
        ptry[1] = Yinitial + (GRID / 2) - 1;
        if (constgrid == 0)
            GRID = GRID + 1;
        OBcount = OBcount + 1;
        ptry[1] = (int)ptry[1];
    }

    if (ptry[2] < (Zinitial - (GRID / 2)))
    {
        ptry[2] = Zinitial - (GRID / 2) + 1;
        if (constgrid == 0)
            GRID = GRID + 1;
        OBcount = OBcount + 1;
        ptry[2] = (int)ptry[2];
    }
    if (ptry[2] > (Zinitial + (GRID / 2)))
    {
        ptry[2] = Zinitial + (GRID / 2) - 1;
        if (constgrid == 0)
            GRID = GRID + 1;
        OBcount = OBcount + 1;
        ptry[2] = (int)ptry[2];
    }
    */

    if (ptry[0] < Xmin)
    {
        ptry[0] = (int)Xmin;
    }
    if (ptry[0] > Xmax)
    {
        ptry[0] = (int)Xmax;
    }
    if (ptry[1] < Ymin)
    {
        ptry[1] = (int)Ymin;
    }
    if (ptry[1] > Ymax)
    {
        ptry[1] = (int)Ymax;
    }
    if (ptry[2] < Zmin)
    {
        ptry[2] = (int)Zmin;
    }
    if (ptry[2] > Zmax)
    {
        ptry[2] = (int)Zmax;
    }

    ptry[2] = (int)ptry[2];
    ptry[0] = (int)ptry[0];
    ptry[1] = (int)ptry[1];

    ytry = SumSphericalMask(radiusCircle, (int)(ptry[0]), (int)(ptry[1]), (int)(ptry[2]), CubeMask);

    if (ytry < y[ihi])
    {
        y[ihi] = ytry;
        for (j = 0; j < ndim; j++)
        {
            psum[j] += ptry[j] - p[ihi][j];
            p[ihi][j] = ptry[j];
        }
    }

    return ytry;
}
// ---------- end to use Amotry

// ---------- begin Swap to be used in Amoeba ---------- Do NOT change without discussion

void SWAP(double &A, double &B)
{
    double swap = A;
    A = B;
    B = swap;
}
// ---------- end Swap

// ---------- begin to use Amoeba ---------- Do NOT change without discussion
// This function is from Numerical Recipes with changes in the argument, amotry, and y[i] below.

void Amoeba(vector<vector<double>> &p, vector<double> &y, const double ftol, int &nfunk, int &radiusCircle, double GRID, int constgrid, vector<vector<vector<int>>> CubeMask)
{
    const int NMAX = 5000;
    const double TINY = 1.0e-10;
    int i, ihi, ilo, inhi, j;
    double rtol, ysave, ytry;
    int mpts = y.size();
    int ndim = p[0].size();
    vector<double> psum(ndim);
    nfunk = 0;
    GetSum(p, psum);
    OBcount = 1;

    for (;;)
    {

        ilo = 0;
        ihi = y[0] > y[1] ? (inhi = 1, 0) : (inhi = 0, 1);
        for (i = 0; i < mpts; i++)
        {
            if (y[i] <= y[ilo])
                ilo = i;
            if (y[i] > y[ihi])
            {
                inhi = ihi;
                ihi = i;
            }
            else if (y[i] > y[inhi] && i != ihi)
                inhi = i;
        }
        rtol = 2.0 * abs(y[ihi] - y[ilo]) / (abs(y[ihi]) + abs(y[ilo]) + TINY);
        if (rtol < ftol)
        {
            SWAP(y[0], y[ilo]);
            for (i = 0; i < ndim; i++)
                SWAP(p[0][i], p[ilo][i]);
            break;
        }
        if (nfunk >= NMAX)
        {
            //::MessageBox(NULL, "NMAX exceeded", "Message", MB_OK);
            errorMessage = "NMAX exceeded.";
            return;
        }
        nfunk += 2;
        ytry = amotry(p, y, psum, ihi, -1.0, radiusCircle, GRID, constgrid, CubeMask);
        if (ytry <= y[ilo])
            ytry = amotry(p, y, psum, ihi, 2.0, radiusCircle, GRID, constgrid, CubeMask);
        else if (ytry >= y[inhi])
        {
            ysave = y[ihi];
            ytry = amotry(p, y, psum, ihi, 0.5, radiusCircle, GRID, constgrid, CubeMask);
            if (ytry >= ysave)
            {
                for (i = 0; i < mpts; i++)
                {
                    if (i != ilo)
                    {
                        for (j = 0; j < ndim; j++)
                            p[i][j] = psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
                        y[i] = SumSphericalMask(radiusCircle, (int)(psum[0]), (int)(psum[1]), (int)(psum[2]), CubeMask);
                    }
                }
                nfunk += ndim;
                GetSum(p, psum);
            }
        }
        else
            --nfunk;
    }
}
// ---------- end to use Amoeba

// This function is used to estimate the subpixel center. It uses various tools from numerical recipes that can be found online
void OnBnClickedEstimatecenter()
{
    // UpdateData(true);
    //  Define pDoc. Get title to check if its subpixel or not.
    // CMDIFrameWnd *pFrame = (CMDIFrameWnd *)AfxGetApp()->m_pMainWnd;
    // CMDIChildWnd *pChild = (CMDIChildWnd *)pFrame->GetActiveFrame();
    // CHipoView *pView = (CHipoView *)pChild->GetActiveView();
    // CHipoDoc *pDoc = (CHipoDoc *)pView->GetDocument();
    // CString TITLE = pDoc->GetTitle();

    // pDoc->OnButtonRestore();

    errorMessage = "";

    if (SubpixelPhaseMatrix.empty())
    {
        //::MessageBox(NULL, "Must generate Subpixel grid first.", "Message", MB_OK);
        errorMessage = "Must generate Subpixel grid first.";
        return;
    }

    double Xcenter, Ycenter, Zcenter;

    int Grid, GridZ;
    int Nfinal = m_SubPixels;
    int RCenterSubpixels;

    double Xmin = pixelToSubpixel((double)(smallBox_X), 0);
    double Xmax = pixelToSubpixel((double)(smallBox_X + smallBox_XSize), 0);
    double Ymin = pixelToSubpixel((double)(smallBox_Y), 1);
    double Ymax = pixelToSubpixel((double)(smallBox_Y + smallBox_YSize), 1);
    double Zmin = pixelToSubpixel((double)(smallBox_Z), 2);
    double Zmax = pixelToSubpixel((double)(smallBox_Z + smallBox_ZSize), 2);

    // Dont think this is necessary
    double RCenter = m_RCenter;
    // double R1 = m_R1;
    double R0 = m_R0;
    // double factor = m_RCenter * m_RCenter * m_RCenterPhase;

    halfdisplay = m_R0;
    halfreal = m_R0;
    Zhalfreal = m_R0;
    displaydim = 2 * m_R0;
    realdim = 2 * m_R0;
    Zrealdim = 2 * Zhalfreal;
    subpixeldisplay = 2 * m_R0 * m_SubPixels;
    subpixelreal = 2 * m_R0 * m_SubPixels;
    Zsubpixelreal = 2 * Zhalfreal * m_SubPixels;

    // First time drawing rectangle on original image
    /*/Draw Rectangle on subpixel image NOT ANYMORE
        ROIcheck=0;
        while (ROIcheck==0){
            GetRectangle(Xinitial, Yinitial, Grid); //center point now referenced to subpixel displayed image
        }
    */
    // center of subpixel grid comes from Find_RCenter
    // Shift 5 subpixels = 0.5 pixel below so the iso-center is at the center of a pixel

    // Commented out because why is this even here?
    /*
        if ((REALXfirst - Xfirst) == 0)
        {
            Xinitial = subpixelreal / 2 + 5;
        }
        else
        {
            Xinitial = subpixelreal / 2; // Unclear about this condition;
        }

        if ((REALYfirst - Yfirst) == 0)
        {
            Yinitial = subpixelreal / 2 + 5;
        }
        else
        {
            Yinitial = subpixelreal / 2;
        }
    */

    // Xinitial = (int)(((m_CenterX - Xfirst + halfdisplay + 0.46) / 0.1));
    // Yinitial = (int)(((m_CenterY - Yfirst + halfdisplay + 0.46) / 0.1));

    // Grid = 30; // in subpixels
    //         UpdateData(true); // in case Rcenter is changes after first rectangle
    // double RcenterMIN;
    //        CString messagea = TEXT("");

    // Removed because why is this here?
    /*
        if ((REALSlice - Zfirst) == 0)
        {
            Zinitial = Zsubpixelreal / 2 + 5;
        }
        else
        {
            Zinitial = Zsubpixelreal / 2;
        }

    Zinitial = (int)(((m_CenterZ - m_CenterZ2 + halfdisplay + 0.46) / 0.1));
*/
    /*
        GridZ = 30; // maybe just try without this first.
        RcenterMIN = ceil(Grid * sqrt(3.0)) / 10;

        if (Grid > (Nfinal * RCenter / sqrt(3.0)))
        {
            RCenter = RcenterMIN; // make sure RCenter contains predefined GRID if not already.
            // messagea.Format(TEXT("Search grid is too large. Changing RCenter to %3.2f to contain this grid."), RcenterMIN);
            // AfxMessageBox(messagea, MB_OK | MB_ICONINFORMATION);
            //  return ;
        }
    */
    // if ((Xinitial + Grid / 2 + ceil(Nfinal * RCenter) > subCONSTmat) || (Yinitial + Grid / 2 + ceil(Nfinal * RCenter) > subCONSTmat))
    if ((Xmax + ceil(Nfinal * RCenter) > subCONSTmat) || (Ymax + ceil(Nfinal * RCenter) > subCONSTmat))
    {
        //::MessageBox(NULL, "Initial center is too close to edge of image", "Message", MB_OK);
        errorMessage = "Initial center is too close to edge of image.";
        return;
    }
    // if ((Xinitial - Grid / 2 - ceil(Nfinal * RCenter) < 0) || (Yinitial - Grid / 2 - ceil(Nfinal * RCenter) < 0))
    if ((Xmin - ceil(Nfinal * RCenter) > subCONSTmat) || (Ymin - ceil(Nfinal * RCenter) > subCONSTmat))
    {
        //::MessageBox(NULL, "Initial center is too close to edge of image", "Message", MB_OK);
        errorMessage = "Initial center is too close to edge of image.";
        return;
    }

    // Condition if the box has same coordinate in respective directions
    // If box in one direction lies in the same subvoxel it will be expanded to the edges of that subvoxel
    if (floor(Xmin / 10.0) == floor(Xmax / 10.0))
    {
        Xmin = floor(Xmin / 10.0) * 10.0;
        Xmax = ceil(Xmax / 10.0) * 10.0 - 1;
    }
    if (floor(Ymin / 10.0) == floor(Ymax / 10.0))
    {
        Ymin = floor(Ymin / 10.0) * 10.0;
        Ymax = ceil(Ymax / 10.0) * 10.0 - 1;
    }
    if (floor(Zmin / 10.0) == floor(Zmax / 10.0))
    {
        Zmin = floor(Zmin / 10.0) * 10.0;
        Zmax = ceil(Zmax / 10.0) * 10.0 - 1;
    }

    // Find center
    //  ---------- begin to set up a unit sphere inside a tight cube---------- Do NOT change without discussion
    //  The center of the sphere is located at (radius, radius, radius). The cube has a size from 0 to 2*radius.
    //  This function should be executed when RCenter is given or changed.
    //  This 3D matrix, CubeMask, it needs to be stored in the memory and be used by the amoeba function.

    RCenterSubpixels = round(RCenter * Nfinal);
    int Dsubpixels = 2 * RCenterSubpixels;
    int diam1 = Dsubpixels + 1;
    double distance;
    int Xdiff, Ydiff, Zdiff;

    vector<vector<vector<int>>> CubeMask(diam1, vector<vector<int>>(diam1, vector<int>(diam1, 0))); // initiate a 3D zero matrix
    CubeMask.clear();                                                                               // in case we need to regenerate this matrix again by re-clicking
    CubeMask.resize(diam1, vector<vector<int>>(diam1, vector<int>(diam1, 0)));

    for (int k = 0; k <= Dsubpixels; k++)
    {
        Zdiff = k - RCenterSubpixels;
        for (int i = 0; i <= Dsubpixels; i++)
        {
            Ydiff = i - RCenterSubpixels;
            for (int j = 0; j <= Dsubpixels; j++)
            {
                CubeMask[j][i][k] = 0;
                Xdiff = j - RCenterSubpixels;
                distance = sqrt((double)Xdiff * Xdiff + Ydiff * Ydiff + Zdiff * Zdiff);
                if (distance <= RCenterSubpixels)
                {
                    CubeMask[j][i][k] = 1;
                }
            }
        }
    }
    //  ---------- end to set up a unit sphere inside a tight cube

    // ---------- begin to set up corner points for the first amoeba search

    const double ftol = 0.00001;
    int nfunk;

    // From Nic: Had to change this - not sure how this vector was initialized like this in the first place?
    // vector<vector<double>> Simplex(4, 3);

    // vector<vector<double>> Simplex(4, vector<double>(3, 0));
    // vector<double> SimplexY(4);
    vector<vector<double>> Simplex(11, vector<double>(3, 0));
    vector<double> SimplexY(11);

    RCenterSubpixels = round(RCenter * Nfinal);

    // Order does not matter - check NUMERICAL RECIPES for clarification
    /*
        Simplex[0][0] = Xinitial;
        Simplex[0][1] = Yinitial;
        Simplex[0][2] = Zinitial;
        Simplex[1][0] = Xinitial + 4;
        Simplex[1][1] = Yinitial;
        Simplex[1][2] = Zinitial;
        Simplex[2][0] = Xinitial;
        Simplex[2][1] = Yinitial + 4;
        Simplex[2][2] = Zinitial;
        Simplex[3][0] = Xinitial;
        Simplex[3][1] = Yinitial;
        Simplex[3][2] = Zinitial + 4;

        SimplexY[0] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]), (int)(Simplex[0][1]), (int)(Simplex[0][2]), CubeMask);
        SimplexY[1] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[1][0]), (int)(Simplex[1][1]), (int)(Simplex[1][2]), CubeMask);
        SimplexY[2] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[2][0]), (int)(Simplex[2][1]), (int)(Simplex[2][2]), CubeMask);
        SimplexY[3] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[3][0]), (int)(Simplex[3][1]), (int)(Simplex[3][2]), CubeMask);
    */

    Simplex[0][0] = Xmin;
    Simplex[0][1] = Ymin;
    Simplex[0][2] = Zmin;

    Simplex[1][0] = Xmax;
    Simplex[1][1] = Ymin;
    Simplex[1][2] = Zmin;

    Simplex[2][0] = Xmin;
    Simplex[2][1] = Ymax;
    Simplex[2][2] = Zmin;

    Simplex[3][0] = Xmin;
    Simplex[3][1] = Ymin;
    Simplex[3][2] = Zmax;

    Simplex[4][0] = Xmax;
    Simplex[4][1] = Ymax;
    Simplex[4][2] = Zmin;

    Simplex[5][0] = Xmax;
    Simplex[5][1] = Ymin;
    Simplex[5][2] = Zmax;

    Simplex[6][0] = Xmin;
    Simplex[6][1] = Ymax;
    Simplex[6][2] = Zmax;

    Simplex[7][0] = Xmax;
    Simplex[7][1] = Ymax;
    Simplex[7][2] = Zmax;

    Simplex[8][0] = pixelToSubpixel((double)(centerL_x), 0);
    Simplex[8][1] = pixelToSubpixel((double)(centerL_y), 1);
    Simplex[8][2] = pixelToSubpixel((double)(centerL_z), 2);

    Simplex[9][0] = pixelToSubpixel((double)(centerM_x), 0);
    Simplex[9][1] = pixelToSubpixel((double)(centerM_y), 1);
    Simplex[9][2] = pixelToSubpixel((double)(centerM_z), 2);

    Simplex[10][0] = pixelToSubpixel((double)(centerS_x), 0);
    Simplex[10][1] = pixelToSubpixel((double)(centerS_y), 1);
    Simplex[10][2] = pixelToSubpixel((double)(centerS_z), 2);

    SimplexY[0] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]), (int)(Simplex[0][1]), (int)(Simplex[0][2]), CubeMask);
    SimplexY[1] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[1][0]), (int)(Simplex[1][1]), (int)(Simplex[1][2]), CubeMask);
    SimplexY[2] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[2][0]), (int)(Simplex[2][1]), (int)(Simplex[2][2]), CubeMask);
    SimplexY[3] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[3][0]), (int)(Simplex[3][1]), (int)(Simplex[3][2]), CubeMask);
    SimplexY[4] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[4][0]), (int)(Simplex[4][1]), (int)(Simplex[4][2]), CubeMask);
    SimplexY[5] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[5][0]), (int)(Simplex[5][1]), (int)(Simplex[5][2]), CubeMask);
    SimplexY[6] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[6][0]), (int)(Simplex[6][1]), (int)(Simplex[6][2]), CubeMask);
    SimplexY[7] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[7][0]), (int)(Simplex[7][1]), (int)(Simplex[7][2]), CubeMask);
    SimplexY[8] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[8][0]), (int)(Simplex[8][1]), (int)(Simplex[8][2]), CubeMask);
    SimplexY[9] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[9][0]), (int)(Simplex[9][1]), (int)(Simplex[9][2]), CubeMask);
    SimplexY[10] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[10][0]), (int)(Simplex[10][1]), (int)(Simplex[10][2]), CubeMask);

    //  ---------- end to set up corner points for the first amoeba search

    Grid = 0;
    Amoeba(Simplex, SimplexY, ftol, nfunk, RCenterSubpixels, Grid, 1, CubeMask);

    ZoomedY = (int)Simplex[0][1];
    ZoomedX = (int)Simplex[0][0];
    lastValueSlice = (int)Simplex[0][2];
    m_CenterX = subpixelToPixel(ZoomedX, 0);
    m_CenterY = subpixelToPixel(ZoomedY, 1);
    m_CenterZ = subpixelToPixel(lastValueSlice, 2);

    /*
    // int dx = abs(ZoomedX - Xinitial);
    // int dy = abs(ZoomedY - Yinitial);
    // int dz = abs(lastValueSlice - Zinitial);
    subCenter = (2 * m_R0 + 1) * (10.0 / 2.0);
    int dx = abs(ZoomedX - subCenter);
    int dy = abs(ZoomedY - subCenter);
    int dz = abs(lastValueSlice - subCenter);

    if ((dx >= dy) && (dx >= dz) && (dx > 10))
        Grid = dx;
    else if ((dy >= dx) && (dy >= dz) && (dy > 10))
        Grid = dy;
    else if ((dz >= dy) && (dz >= dx) && (dz > 10))
        Grid = dz;
    else
        Grid = 10;

    RCenter = m_RCenter;
    if (Grid > (Nfinal * RCenter / sqrt(3.0)))
    {
        m_RCenter = ceil(Grid * sqrt(3.0)) / 10;
    }
        RCenterSubpixels = round(RCenter * m_SubPixels);

        // ---------- begin to set up corner points for the second amoeba search

        Simplex[1][0] = Simplex[0][0] + 4;
        Simplex[1][1] = Simplex[0][1];
        Simplex[1][2] = Simplex[0][2];
        Simplex[2][0] = Simplex[0][0];
        Simplex[2][1] = Simplex[0][1] + 4;
        Simplex[2][2] = Simplex[0][2];
        Simplex[3][0] = Simplex[0][0];
        Simplex[3][1] = Simplex[0][1];
        Simplex[3][2] = Simplex[0][2] + 4;

        SimplexY[0] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]), (int)(Simplex[0][1]), (int)(Simplex[0][2]), CubeMask);
        SimplexY[1] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[1][0]), (int)(Simplex[1][1]), (int)(Simplex[1][2]), CubeMask);
        SimplexY[2] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[2][0]), (int)(Simplex[2][1]), (int)(Simplex[2][2]), CubeMask);
        SimplexY[3] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[3][0]), (int)(Simplex[3][1]), (int)(Simplex[3][2]), CubeMask);

        // ---------- end to set up corner points for the second amoeba search

        Grid = 0;
        Amoeba(Simplex, SimplexY, ftol, nfunk, RCenterSubpixels, Grid, 1, CubeMask);

        ZoomedY = (int)Simplex[0][1];
        ZoomedX = (int)Simplex[0][0];
        lastValueSlice = (int)Simplex[0][2];

        Xcenter = m_CenterX2 - halfdisplay - 0.45 + 0.1 * ZoomedX; // coordinate on original image
        Ycenter = m_CenterY2 - halfdisplay - 0.45 + 0.1 * ZoomedY;
        Zcenter = m_CenterZ2 - Zhalfreal - 0.45 + 0.1 * lastValueSlice;
    */
    // m_CenterX = Xcenter;
    // m_CenterY = Ycenter;
    // m_CenterZ = Zcenter;

    // UpdateData(false);

    // ---------- begin to check neighboring pixels with original Rcenter

    RCenterSubpixels = round(RCenter * m_SubPixels);

    SimplexY.resize(27);
    int i = 1;
    int ilo, max_it;
    max_it = 1;
    /*
        for (;;)
        {

            if (max_it == 10)
            {
                //::MessageBox(NULL, "Too many minimums found. Try reducing RCenter.", "Message", MB_OK);
                errorMessage = "Too many minimums found. Try reducing RCenter.";
                break;
            }
    */
    /*
     if (Simplex[0][0] > Xmax)
     {
         Simplex[0][0] = Xmax;
     }
     if (Simplex[0][0] < Xmin)
     {
         Simplex[0][0] = Xmin;
     }
     if (Simplex[0][1] > Ymax)
     {
         Simplex[0][1] = Ymax;
     }
     if (Simplex[0][1] < Ymin)
     {
         Simplex[0][1] = Ymin;
     }
     if (Simplex[0][2] > Zmax)
     {
         Simplex[0][2] = Zmax;
     }
     if (Simplex[0][2] < Zmin)
     {
         Simplex[0][2] = Zmin;
     }
     */

    // making sure simplex coordinates are in small box
    int neighborXMAX = 3;
    int neighborXMIN = 3;
    int neighborYMAX = 3;
    int neighborYMIN = 3;
    int neighborZMAX = 3;
    int neighborZMIN = 3;

    if (Xmax - Simplex[0][0] < neighborXMAX)
    {
        neighborXMAX = Xmax - Simplex[0][0];
    }
    if (Simplex[0][0] - Xmin < neighborXMIN)
    {
        neighborXMIN = Simplex[0][0] - Xmin;
    }
    if (Ymax - Simplex[0][1] < neighborYMAX)
    {
        neighborYMAX = Ymax - Simplex[0][1];
    }
    if (Simplex[0][1] - Ymin < neighborYMIN)
    {
        neighborYMIN = Simplex[0][1] - Ymin;
    }
    if (Zmax - Simplex[0][2] < neighborZMAX)
    {
        neighborZMAX = Zmax - Simplex[0][2];
    }
    if (Simplex[0][2] - Zmin < neighborZMIN)
    {
        neighborZMIN = Simplex[0][2] - Zmin;
    }
    /*
     SimplexY[0] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]), (int)(Simplex[0][1]), (int)(Simplex[0][2]), CubeMask);
            SimplexY[1] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]), (int)(Simplex[0][1]) + 3, (int)(Simplex[0][2]) + 3, CubeMask);
            SimplexY[2] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]), (int)(Simplex[0][1]) + 3, (int)(Simplex[0][2]), CubeMask);
            SimplexY[3] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]), (int)(Simplex[0][1]) + 3, (int)(Simplex[0][2]) - 3, CubeMask);
            SimplexY[4] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]), (int)(Simplex[0][1]) - 3, (int)(Simplex[0][2]) + 3, CubeMask);
            SimplexY[5] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]), (int)(Simplex[0][1]) - 3, (int)(Simplex[0][2]), CubeMask);
            SimplexY[6] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]), (int)(Simplex[0][1]) - 3, (int)(Simplex[0][2]) - 3, CubeMask);
            SimplexY[7] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]), (int)(Simplex[0][1]), (int)(Simplex[0][2]) + 3, CubeMask);
            SimplexY[8] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]), (int)(Simplex[0][1]), (int)(Simplex[0][2]) - 3, CubeMask);
            SimplexY[9] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]) + 3, (int)(Simplex[0][1]) + 3, (int)(Simplex[0][2]) + 3, CubeMask);
            SimplexY[10] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]) + 3, (int)(Simplex[0][1]) + 3, (int)(Simplex[0][2]), CubeMask);
            SimplexY[11] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]) + 3, (int)(Simplex[0][1]) + 3, (int)(Simplex[0][2]) - 3, CubeMask);
            SimplexY[12] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]) + 3, (int)(Simplex[0][1]), (int)(Simplex[0][2]) + 3, CubeMask);
            SimplexY[13] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]) + 3, (int)(Simplex[0][1]), (int)(Simplex[0][2]), CubeMask);
            SimplexY[14] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]) + 3, (int)(Simplex[0][1]), (int)(Simplex[0][2]) - 3, CubeMask);
            SimplexY[15] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]) + 3, (int)(Simplex[0][1]) - 3, (int)(Simplex[0][2]) + 3, CubeMask);
            SimplexY[16] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]) + 3, (int)(Simplex[0][1]) - 3, (int)(Simplex[0][2]), CubeMask);
            SimplexY[17] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]) + 3, (int)(Simplex[0][1]) - 3, (int)(Simplex[0][2]) - 3, CubeMask);
            SimplexY[18] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]) - 3, (int)(Simplex[0][1]) + 3, (int)(Simplex[0][2]) + 3, CubeMask);
            SimplexY[19] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]) - 3, (int)(Simplex[0][1]) + 3, (int)(Simplex[0][2]), CubeMask);
            SimplexY[20] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]) - 3, (int)(Simplex[0][1]) + 3, (int)(Simplex[0][2]) - 3, CubeMask);
            SimplexY[21] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]) - 3, (int)(Simplex[0][1]), (int)(Simplex[0][2]) + 3, CubeMask);
            SimplexY[22] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]) - 3, (int)(Simplex[0][1]), (int)(Simplex[0][2]), CubeMask);
            SimplexY[23] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]) - 3, (int)(Simplex[0][1]), (int)(Simplex[0][2]) - 3, CubeMask);
            SimplexY[24] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]) - 3, (int)(Simplex[0][1]) - 3, (int)(Simplex[0][2]) + 3, CubeMask);
            SimplexY[25] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]) - 3, (int)(Simplex[0][1]) - 3, (int)(Simplex[0][2]), CubeMask);
            SimplexY[26] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]) - 3, (int)(Simplex[0][1]) - 3, (int)(Simplex[0][2]) - 3, CubeMask);
    */
    Simplex.resize(27, vector<double>(3, 0));

    // fst << "resized\n";
    Simplex[1][0] = Simplex[0][0];
    Simplex[1][1] = Simplex[0][1] + neighborYMAX;
    Simplex[1][2] = Simplex[0][2] + neighborZMAX;

    Simplex[2][0] = Simplex[0][0];
    Simplex[2][1] = Simplex[0][1] + neighborYMAX;
    Simplex[2][2] = Simplex[0][2];

    Simplex[3][0] = Simplex[0][0];
    Simplex[3][1] = Simplex[0][1] + neighborYMAX;
    Simplex[3][2] = Simplex[0][2] - neighborZMIN;

    Simplex[4][0] = Simplex[0][0];
    Simplex[4][1] = Simplex[0][1] - neighborYMIN;
    Simplex[4][2] = Simplex[0][2] + neighborZMAX;

    Simplex[5][0] = Simplex[0][0];
    Simplex[5][1] = Simplex[0][1] - neighborYMIN;
    Simplex[5][2] = Simplex[0][2];

    Simplex[6][0] = Simplex[0][0];
    Simplex[6][1] = Simplex[0][1] - neighborYMIN;
    Simplex[6][2] = Simplex[0][2] - neighborZMIN;

    Simplex[7][0] = Simplex[0][0];
    Simplex[7][1] = Simplex[0][1];
    Simplex[7][2] = Simplex[0][2] + neighborZMAX;

    Simplex[8][0] = Simplex[0][0];
    Simplex[8][1] = Simplex[0][1];
    Simplex[8][2] = Simplex[0][2] - neighborZMIN;

    Simplex[9][0] = Simplex[0][0] + neighborXMAX;
    Simplex[9][1] = Simplex[0][1] + neighborYMAX;
    Simplex[9][2] = Simplex[0][2] + neighborZMAX;

    Simplex[10][0] = Simplex[0][0] + neighborXMAX;
    Simplex[10][1] = Simplex[0][1] + neighborYMAX;
    Simplex[10][2] = Simplex[0][2];

    Simplex[11][0] = Simplex[0][0] + neighborXMAX;
    Simplex[11][1] = Simplex[0][1] + neighborYMAX;
    Simplex[11][2] = Simplex[0][2] - neighborZMIN;

    Simplex[12][0] = Simplex[0][0] + neighborXMAX;
    Simplex[12][1] = Simplex[0][1];
    Simplex[12][2] = Simplex[0][2] + neighborZMAX;

    Simplex[13][0] = Simplex[0][0] + neighborXMAX;
    Simplex[13][1] = Simplex[0][1];
    Simplex[13][2] = Simplex[0][2];

    Simplex[14][0] = Simplex[0][0] + neighborXMAX;
    Simplex[14][1] = Simplex[0][1];
    Simplex[14][2] = Simplex[0][2] - neighborZMIN;

    Simplex[15][0] = Simplex[0][0] + neighborXMAX;
    Simplex[15][1] = Simplex[0][1] - neighborYMIN;
    Simplex[15][2] = Simplex[0][2] + neighborZMAX;

    Simplex[16][0] = Simplex[0][0] + neighborXMAX;
    Simplex[16][1] = Simplex[0][1] - neighborYMIN;
    Simplex[16][2] = Simplex[0][2];

    Simplex[17][0] = Simplex[0][0] + neighborXMAX;
    Simplex[17][1] = Simplex[0][1] - neighborYMIN;
    Simplex[17][2] = Simplex[0][2] - neighborZMIN;

    Simplex[18][0] = Simplex[0][0] - neighborXMIN;
    Simplex[18][1] = Simplex[0][1] + neighborYMAX;
    Simplex[18][2] = Simplex[0][2] + neighborZMAX;

    Simplex[19][0] = Simplex[0][0] - neighborXMIN;
    Simplex[19][1] = Simplex[0][1] + neighborYMAX;
    Simplex[19][2] = Simplex[0][2];

    Simplex[20][0] = Simplex[0][0] - neighborXMIN;
    Simplex[20][1] = Simplex[0][1] + neighborYMAX;
    Simplex[20][2] = Simplex[0][2] - neighborZMIN;

    Simplex[21][0] = Simplex[0][0] - neighborXMIN;
    Simplex[21][1] = Simplex[0][1];
    Simplex[21][2] = Simplex[0][2] + neighborZMAX;

    Simplex[22][0] = Simplex[0][0] - neighborXMIN;
    Simplex[22][1] = Simplex[0][1];
    Simplex[22][2] = Simplex[0][2];

    Simplex[23][0] = Simplex[0][0] - neighborXMIN;
    Simplex[23][1] = Simplex[0][1];
    Simplex[23][2] = Simplex[0][2] - neighborZMIN;

    Simplex[24][0] = Simplex[0][0] - neighborXMIN;
    Simplex[24][1] = Simplex[0][1] - neighborYMIN;
    Simplex[24][2] = Simplex[0][2] + neighborZMAX;

    Simplex[25][0] = Simplex[0][0] - neighborXMIN;
    Simplex[25][1] = Simplex[0][1] - neighborYMIN;
    Simplex[25][2] = Simplex[0][2];

    Simplex[26][0] = Simplex[0][0] - neighborXMIN;
    Simplex[26][1] = Simplex[0][1] - neighborYMIN;
    Simplex[26][2] = Simplex[0][2] - neighborZMIN;

    SimplexY[0] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[0][0]), (int)(Simplex[0][1]), (int)(Simplex[0][2]), CubeMask);
    SimplexY[1] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[1][0]), (int)(Simplex[1][1]), (int)(Simplex[1][2]), CubeMask);
    SimplexY[2] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[2][0]), (int)(Simplex[2][1]), (int)(Simplex[2][2]), CubeMask);
    SimplexY[3] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[3][0]), (int)(Simplex[3][1]), (int)(Simplex[3][2]), CubeMask);
    SimplexY[4] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[4][0]), (int)(Simplex[4][1]), (int)(Simplex[4][2]), CubeMask);
    SimplexY[5] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[5][0]), (int)(Simplex[5][1]), (int)(Simplex[5][2]), CubeMask);
    SimplexY[6] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[6][0]), (int)(Simplex[6][1]), (int)(Simplex[6][2]), CubeMask);
    SimplexY[7] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[7][0]), (int)(Simplex[7][1]), (int)(Simplex[7][2]), CubeMask);
    SimplexY[8] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[8][0]), (int)(Simplex[8][1]), (int)(Simplex[8][2]), CubeMask);
    SimplexY[9] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[9][0]), (int)(Simplex[9][1]), (int)(Simplex[9][2]), CubeMask);
    SimplexY[10] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[10][0]), (int)(Simplex[10][1]), (int)(Simplex[10][2]), CubeMask);
    SimplexY[11] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[11][0]), (int)(Simplex[11][1]), (int)(Simplex[11][2]), CubeMask);
    SimplexY[12] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[12][0]), (int)(Simplex[12][1]), (int)(Simplex[12][2]), CubeMask);
    SimplexY[13] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[13][0]), (int)(Simplex[13][1]), (int)(Simplex[13][2]), CubeMask);
    SimplexY[14] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[14][0]), (int)(Simplex[14][1]), (int)(Simplex[14][2]), CubeMask);
    SimplexY[15] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[15][0]), (int)(Simplex[15][1]), (int)(Simplex[15][2]), CubeMask);
    SimplexY[16] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[16][0]), (int)(Simplex[16][1]), (int)(Simplex[16][2]), CubeMask);
    SimplexY[17] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[17][0]), (int)(Simplex[17][1]), (int)(Simplex[17][2]), CubeMask);
    SimplexY[18] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[18][0]), (int)(Simplex[18][1]), (int)(Simplex[18][2]), CubeMask);
    SimplexY[19] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[19][0]), (int)(Simplex[19][1]), (int)(Simplex[19][2]), CubeMask);
    SimplexY[20] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[20][0]), (int)(Simplex[20][1]), (int)(Simplex[20][2]), CubeMask);
    SimplexY[21] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[21][0]), (int)(Simplex[21][1]), (int)(Simplex[21][2]), CubeMask);
    SimplexY[22] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[22][0]), (int)(Simplex[22][1]), (int)(Simplex[22][2]), CubeMask);
    SimplexY[23] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[23][0]), (int)(Simplex[23][1]), (int)(Simplex[23][2]), CubeMask);
    SimplexY[24] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[24][0]), (int)(Simplex[24][1]), (int)(Simplex[24][2]), CubeMask);
    SimplexY[25] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[25][0]), (int)(Simplex[25][1]), (int)(Simplex[25][2]), CubeMask);
    SimplexY[26] = SumSphericalMask(RCenterSubpixels, (int)(Simplex[26][0]), (int)(Simplex[26][1]), (int)(Simplex[26][2]), CubeMask);

    /*
        ilo = 0;
        i = 1;

        while (i < SimplexY.size())
        {
            if (SimplexY[i] < SimplexY[0])
                ilo = i;
            i++;
        }
    */
    Amoeba(Simplex, SimplexY, ftol, nfunk, RCenterSubpixels, Grid, 1, CubeMask);
    // fst << "amoeba\n";

    // if (ilo == 0)
    //     break;
    /*
        switch (ilo)
        {
            case 1:
            Simplex[0][0] = Simplex[0][0];
            Simplex[0][1] = Simplex[0][1] + 3;
            Simplex[0][2] = Simplex[0][2] + 3;
            break;
        case 2:
            Simplex[0][0] = Simplex[0][0];
            Simplex[0][1] = Simplex[0][1] + 3;
            Simplex[0][2] = Simplex[0][2];
            break;
        case 3:
            Simplex[0][0] = Simplex[0][0];
            Simplex[0][1] = Simplex[0][1] + 3;
            Simplex[0][2] = Simplex[0][2] - 3;
            break;
        case 4:
            Simplex[0][0] = Simplex[0][0];
            Simplex[0][1] = Simplex[0][1] - 3;
            Simplex[0][2] = Simplex[0][2] + 3;
            break;
        case 5:
            Simplex[0][0] = Simplex[0][0];
            Simplex[0][1] = Simplex[0][1] - 3;
            Simplex[0][2] = Simplex[0][2];
            break;
        case 6:
            Simplex[0][0] = Simplex[0][0];
            Simplex[0][1] = Simplex[0][1] - 3;
            Simplex[0][2] = Simplex[0][2] - 3;
            break;
        case 7:
            Simplex[0][0] = Simplex[0][0];
            Simplex[0][1] = Simplex[0][1];
            Simplex[0][2] = Simplex[0][2] + 3;
            break;
        case 8:
            Simplex[0][0] = Simplex[0][0];
            Simplex[0][1] = Simplex[0][1];
            Simplex[0][2] = Simplex[0][2] - 3;
            break;
        case 9:
            Simplex[0][0] = Simplex[0][0] + 3;
            Simplex[0][1] = Simplex[0][1] + 3;
            Simplex[0][2] = Simplex[0][2] + 3;
            break;
        case 10:
            Simplex[0][0] = Simplex[0][0] + 3;
            Simplex[0][1] = Simplex[0][1] + 3;
            Simplex[0][2] = Simplex[0][2];
            break;
        case 11:
            Simplex[0][0] = Simplex[0][0] + 3;
            Simplex[0][1] = Simplex[0][1] + 3;
            Simplex[0][2] = Simplex[0][2] - 3;
            break;
        case 12:
            Simplex[0][0] = Simplex[0][0] + 3;
            Simplex[0][1] = Simplex[0][1];
            Simplex[0][2] = Simplex[0][2] + 3;
            break;
        case 13:
            Simplex[0][0] = Simplex[0][0] + 3;
            Simplex[0][1] = Simplex[0][1];
            Simplex[0][2] = Simplex[0][2];
            break;
        case 14:
            Simplex[0][0] = Simplex[0][0] + 3;
            Simplex[0][1] = Simplex[0][1];
            Simplex[0][2] = Simplex[0][2] - 3;
            break;
        case 15:
            Simplex[0][0] = Simplex[0][0] + 3;
            Simplex[0][1] = Simplex[0][1] - 3;
            Simplex[0][2] = Simplex[0][2] + 3;
            break;
        case 16:
            Simplex[0][0] = Simplex[0][0] + 3;
            Simplex[0][1] = Simplex[0][1] - 3;
            Simplex[0][2] = Simplex[0][2];
            break;
        case 17:
            Simplex[0][0] = Simplex[0][0] + 3;
            Simplex[0][1] = Simplex[0][1] - 3;
            Simplex[0][2] = Simplex[0][2] - 3;
            break;
        case 18:
            Simplex[0][0] = Simplex[0][0] - 3;
            Simplex[0][1] = Simplex[0][1] + 3;
            Simplex[0][2] = Simplex[0][2] + 3;
            break;
        case 19:
            Simplex[0][0] = Simplex[0][0] - 3;
            Simplex[0][1] = Simplex[0][1] + 3;
            Simplex[0][2] = Simplex[0][2];
            break;
        case 20:
            Simplex[0][0] = Simplex[0][0] - 3;
            Simplex[0][1] = Simplex[0][1] + 3;
            Simplex[0][2] = Simplex[0][2] - 3;
            break;
        case 21:
            Simplex[0][0] = Simplex[0][0] - 3;
            Simplex[0][1] = Simplex[0][1];
            Simplex[0][2] = Simplex[0][2] + 3;
            break;
        case 22:
            Simplex[0][0] = Simplex[0][0] - 3;
            Simplex[0][1] = Simplex[0][1];
            Simplex[0][2] = Simplex[0][2];
            break;
        case 23:
            Simplex[0][0] = Simplex[0][0] - 3;
            Simplex[0][1] = Simplex[0][1];
            Simplex[0][2] = Simplex[0][2] - 3;
            break;
        case 24:
            Simplex[0][0] = Simplex[0][0] - 3;
            Simplex[0][1] = Simplex[0][1] - 3;
            Simplex[0][2] = Simplex[0][2] + 3;
            break;
        case 25:
            Simplex[0][0] = Simplex[0][0] - 3;
            Simplex[0][1] = Simplex[0][1] - 3;
            Simplex[0][2] = Simplex[0][2];
            break;
        case 26:
            Simplex[0][0] = Simplex[0][0] - 3;
            Simplex[0][1] = Simplex[0][1] - 3;
            Simplex[0][2] = Simplex[0][2] - 3;
            break;


        case 1:
            Simplex[0][0] = Simplex[1][0];
            Simplex[0][1] = Simplex[1][1];
            Simplex[0][2] = Simplex[1][2];
            break;
        case 2:
            Simplex[0][0] = Simplex[2][0];
            Simplex[0][1] = Simplex[2][1];
            Simplex[0][2] = Simplex[2][2];
            break;
        case 3:
            Simplex[0][0] = Simplex[3][0];
            Simplex[0][1] = Simplex[3][1];
            Simplex[0][2] = Simplex[3][2];
            break;
        case 4:
            Simplex[0][0] = Simplex[4][0];
            Simplex[0][1] = Simplex[4][1];
            Simplex[0][2] = Simplex[4][2];
            break;
        case 5:
            Simplex[0][0] = Simplex[5][0];
            Simplex[0][1] = Simplex[5][1];
            Simplex[0][2] = Simplex[5][2];
            break;
        case 6:
            Simplex[0][0] = Simplex[6][0];
            Simplex[0][1] = Simplex[6][1];
            Simplex[0][2] = Simplex[6][2];
            break;
        case 7:
            Simplex[0][0] = Simplex[7][0];
            Simplex[0][1] = Simplex[7][1];
            Simplex[0][2] = Simplex[7][2];
            break;
        case 8:
            Simplex[0][0] = Simplex[8][0];
            Simplex[0][1] = Simplex[8][1];
            Simplex[0][2] = Simplex[8][2];
            break;
        case 9:
            Simplex[0][0] = Simplex[9][0];
            Simplex[0][1] = Simplex[9][1];
            Simplex[0][2] = Simplex[9][2];
            break;
        case 10:
            Simplex[0][0] = Simplex[10][0];
            Simplex[0][1] = Simplex[10][1];
            Simplex[0][2] = Simplex[10][2];
            break;
        case 11:
            Simplex[0][0] = Simplex[11][0];
            Simplex[0][1] = Simplex[11][1];
            Simplex[0][2] = Simplex[11][2];
            break;
        case 12:
            Simplex[0][0] = Simplex[12][0];
            Simplex[0][1] = Simplex[12][1];
            Simplex[0][2] = Simplex[12][2];
            break;
        case 13:
            Simplex[0][0] = Simplex[13][0];
            Simplex[0][1] = Simplex[13][1];
            Simplex[0][2] = Simplex[13][2];
            break;
        case 14:
            Simplex[0][0] = Simplex[14][0];
            Simplex[0][1] = Simplex[14][1];
            Simplex[0][2] = Simplex[14][2];
            break;
        case 15:
            Simplex[0][0] = Simplex[15][0];
            Simplex[0][1] = Simplex[15][1];
            Simplex[0][2] = Simplex[15][2];
            break;
        case 16:
            Simplex[0][0] = Simplex[16][0];
            Simplex[0][1] = Simplex[16][1];
            Simplex[0][2] = Simplex[16][2];
            break;
        case 17:
            Simplex[0][0] = Simplex[17][0];
            Simplex[0][1] = Simplex[17][1];
            Simplex[0][2] = Simplex[17][2];
            break;
        case 18:
            Simplex[0][0] = Simplex[18][0];
            Simplex[0][1] = Simplex[18][1];
            Simplex[0][2] = Simplex[18][2];
            break;
        case 19:
            Simplex[0][0] = Simplex[19][0];
            Simplex[0][1] = Simplex[19][1];
            Simplex[0][2] = Simplex[19][2];
            break;
        case 20:
            Simplex[0][0] = Simplex[20][0];
            Simplex[0][1] = Simplex[20][1];
            Simplex[0][2] = Simplex[20][2];
            break;
        case 21:
            Simplex[0][0] = Simplex[21][0];
            Simplex[0][1] = Simplex[21][1];
            Simplex[0][2] = Simplex[21][2];
            break;
        case 22:
            Simplex[0][0] = Simplex[22][0];
            Simplex[0][1] = Simplex[22][1];
            Simplex[0][2] = Simplex[22][2];
            break;
        case 23:
            Simplex[0][0] = Simplex[23][0];
            Simplex[0][1] = Simplex[23][1];
            Simplex[0][2] = Simplex[23][2];
            break;
        case 24:
            Simplex[0][0] = Simplex[24][0];
            Simplex[0][1] = Simplex[24][1];
            Simplex[0][2] = Simplex[24][2];
            break;
        case 25:
            Simplex[0][0] = Simplex[25][0];
            Simplex[0][1] = Simplex[25][1];
            Simplex[0][2] = Simplex[25][2];
            break;
        case 26:
            Simplex[0][0] = Simplex[26][0];
            Simplex[0][1] = Simplex[26][1];
            Simplex[0][2] = Simplex[26][2];
            break;
        // From Nic: this is proper for a switch case, you should always have a default
        default:
            break;
        }*/

    // From Nic: Added this since it has the same logic that the switch() condition had. This is obviously more efficient
    // Simplex[0][0] = Simplex[ilo][0];
    // Simplex[0][1] = Simplex[ilo][1];
    // Simplex[0][2] = Simplex[ilo][2];

    /*
                    if (ilo != 0)
                    {

                        double CenterX3 = m_CenterX2 - halfdisplay - 0.45 + 0.1 * (int)Simplex[0][0];
                        double CenterY3 = m_CenterY2 - halfdisplay - 0.45 + 0.1 * (int)Simplex[0][1];
                        double CenterZ3 = m_CenterZ2 - halfdisplay - 0.45 + 0.1 * (int)Simplex[0][2];

                        // From Nic: I think this just set the other subpixel coordinate values in the GUI, so instead I'm passing these to Java

                        m_CenterX3.Format(TEXT("%3.2f"), CenterX3);
                        m_CenterY3.Format(TEXT("%3.2f"), CenterY3);
                        m_CenterZ3.Format(TEXT("%3.2f"), CenterZ3);

            m_CenterX3 = CenterX3;
            m_CenterY3 = CenterY3;
            m_CenterZ3 = CenterZ3;
        }
        */

    // max_it = max_it + 1;
    // ---------- end to check neighboring pixels with original Rcenter

    m_CenterX = subpixelToPixel((int)Simplex[0][0] + 0.5, 0);
    m_CenterY = subpixelToPixel((int)Simplex[0][1] + 0.5, 1);
    m_CenterZ = subpixelToPixel((int)Simplex[0][2] + 0.5, 2);
    // fst << m_CenterX << " " << m_CenterY << " " << m_CenterZ << '\n';

    return;
}

// ---------- begin to add up real values inside a sphere for simulated images ---------- Do NOT change without discussion
double SumCircleElementsReal3D_SIM(int radius, int Scan1, int Scan2, int Scan3)
{

    int diameter = 2 * radius;
    int newx, newy, newz, Xdiff, Ydiff, Zdiff;
    double sum = 0, sumi = 0;
    double distance;

    for (int k = 0; k <= diameter; k++)
    {
        newz = Scan3 - radius + k;
        Zdiff = newz - Scan3;
        for (int i = 0; i <= diameter; i++)
        {
            newy = Scan2 - radius + i;
            Ydiff = newy - Scan2;
            for (int j = 0; j <= diameter; j++)
            {
                newx = Scan1 - radius + j;
                Xdiff = newx - Scan1;

                distance = sqrt((double)Xdiff * Xdiff + Ydiff * Ydiff + Zdiff * Zdiff);
                if (distance <= radius)
                {
                    // sumi = SubpixelRealMatrix3D[(int)newx][(int)newy][(int)newz];
                    sumi = SubpixelSimulatedRealMatrix3D[(int)newx][(int)newy][(int)newz];
                    sum += sumi;
                }
            }
        }
    }
    return sum;
}
// ---------- end to add up real values inside a sphere for simulated images

// ---------- begin to add up real values inside a sphere
double SumCircleElementsReal3D(int radius, int Scan1, int Scan2, int Scan3)
{

    int diameter = 2 * radius;
    int newx, newy, newz, Xdiff, Ydiff, Zdiff;
    double sum = 0, sumi = 0;
    double distance;

    for (int k = 0; k <= diameter; k++)
    {
        newz = Scan3 - radius + k;
        Zdiff = newz - Scan3;
        for (int i = 0; i <= diameter; i++)
        {
            newy = Scan2 - radius + i;
            Ydiff = newy - Scan2;
            for (int j = 0; j <= diameter; j++)
            {
                newx = Scan1 - radius + j;
                Xdiff = newx - Scan1;

                distance = sqrt((double)Xdiff * Xdiff + Ydiff * Ydiff + Zdiff * Zdiff);
                if (distance <= radius)
                {
                    sumi = SubpixelRealMatrix3D[(int)newx][(int)newy][(int)newz];
                    sum += sumi;
                }
            }
        }
    }
    return sum;
}
// ---------- end to add up real values inside a sphere

// ---------- begin to add up imaginary values inside a sphere ---------- Do NOT change without discussion

double SumCircleElementsImag3D(int radius, int Scan1, int Scan2, int Scan3)
{

    int diameter = 2 * radius;
    int newx, newy, newz, Xdiff, Ydiff, Zdiff;
    double sum = 0, sumi = 0;
    double distance;

    for (int k = 0; k <= diameter; k++)
    {
        newz = Scan3 - radius + k;
        Zdiff = newz - Scan3;
        for (int i = 0; i <= diameter; i++)
        {
            newy = Scan2 - radius + i;
            Ydiff = newy - Scan2;
            for (int j = 0; j <= diameter; j++)
            {
                newx = Scan1 - radius + j;
                Xdiff = newx - Scan1;

                distance = sqrt((double)Xdiff * Xdiff + Ydiff * Ydiff + Zdiff * Zdiff);
                if (distance <= radius)
                {
                    sumi = SubpixelImagMatrix3D[(int)newx][(int)newy][(int)newz];
                    sum += sumi;
                }
            }
        }
    }
    return sum;
}
// ---------- end to add up imaginary values inside a sphere

// ---------- begin to add up real values inside a sphere for estimating bkg and spin density ---------- Do NOT change without discussion
double SumCircleElementsReal3D_S5(int radius, int Scan1, int Scan2, int Scan3)
{

    int diameter = 2 * radius;
    int newx, newy, newz, Xdiff, Ydiff, Zdiff;
    double sum = 0, sumi = 0;
    double distance;

    for (int k = 0; k <= diameter; k++)
    {
        newz = Scan3 - radius + k;
        Zdiff = newz - Scan3;
        for (int i = 0; i <= diameter; i++)
        {
            newy = Scan2 - radius + i;
            Ydiff = newy - Scan2;
            for (int j = 0; j <= diameter; j++)
            {
                newx = Scan1 - radius + j;
                Xdiff = newx - Scan1;

                distance = sqrt((double)Xdiff * Xdiff + Ydiff * Ydiff + Zdiff * Zdiff);
                if (distance <= radius)
                {
                    sumi = SubpixelRealMatrix3D_S5[(int)newx][(int)newy][(int)newz];
                    sum += sumi;
                }
            }
        }
    }
    return sum;
}
// ---------- end to add up real values inside a sphere

// ---------- begin to add up imaginary values inside a sphere for estimating bkg and spin density ---------- Do NOT change without discussion

double SumCircleElementsImag3D_S5(int radius, int Scan1, int Scan2, int Scan3)
{

    int diameter = 2 * radius;
    int newx, newy, newz, Xdiff, Ydiff, Zdiff;
    double sum = 0, sumi = 0;
    double distance;

    for (int k = 0; k <= diameter; k++)
    {
        newz = Scan3 - radius + k;
        Zdiff = newz - Scan3;
        for (int i = 0; i <= diameter; i++)
        {
            newy = Scan2 - radius + i;
            Ydiff = newy - Scan2;
            for (int j = 0; j <= diameter; j++)
            {
                newx = Scan1 - radius + j;
                Xdiff = newx - Scan1;

                distance = sqrt((double)Xdiff * Xdiff + Ydiff * Ydiff + Zdiff * Zdiff);
                if (distance <= radius)
                {
                    sumi = SubpixelImagMatrix3D_S5[(int)newx][(int)newy][(int)newz];
                    sum += sumi;
                }
            }
        }
    }
    return sum;
}
// ---------- end to add up imaginary values inside a sphere

void polint(double xa[5], double ya[5], const double x, double &y, double &dy)
{
    int i, m, ns = 0;
    double den, dif, dift, ho, hp, w;

    // int n=xa.size();
    // double c[n],d[n]; changed to 5, which is constant set in qromb function
    int n = 5;
    double c[5], d[5];
    dif = fabs(x - xa[0]);
    for (i = 0; i < n; i++)
    {
        if ((dift = fabs(x - xa[i])) < dif)
        {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }
    y = ya[ns--];
    for (m = 1; m < n; m++)
    {
        for (i = 0; i < n - m; i++)
        {
            ho = xa[i] - x;
            hp = xa[i + m] - x;
            w = c[i + 1] - d[i];
            if ((den = ho - hp) == 0.0)
                //::MessageBox(NULL, "Error in routine polint", "Message", MB_OK);
                errorMessage_magMom = "Error in routine polint";
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        y += (dy = (2 * (ns + 1) < (n - m) ? c[ns + 1] : d[ns--]));
    }
}

// ---------- begin to calculate the integrands of f_ij ---------- Do NOT change without discussion

double innerfunction3D_REa(double x, double ph1)
{
    return ((1 / (x * x)) * (2 * cos(x * ph1) + cos(2 * x * ph1)));
}
double innerfunction3D_REb(double x, double p, double ph1, double ph2)
{
    if (x == 0)
        return (3 / 4) * ((p / ph1) - (p / ph2));

    return ((1 / (x * x)) * (2 - (2 - x) * sqrt(1 + x)) * ((p / ph1) * cos(x * ph1) - (p / ph2) * cos(x * ph2)));
}

double innerfunction3D_IMa(double x, double ph1)
{
    return ((1 / (x * x)) * (2 * sin(x * ph1) - sin(2 * x * ph1)));
}
double innerfunction3D_IMb(double x, double p, double ph1, double ph2)
{
    if (x == 0)
        return 0;

    return ((1 / (x * x)) * (2 - (2 - x) * sqrt(1 + x)) * ((p / ph1) * (-1) * sin(x * ph1) + (p / ph2) * sin(x * ph2)));
}
// ---------- end to calculate the integrands of f_ij

double innerfunction3D_dREa(double x, double ph1)
{
    return ((1 / (x)) * (sin(x * ph1) + sin(2 * x * ph1)));
}

// ---------- begin to set up integrands in Eq. 1 in CISSCO paper ---------- Do NOT change without discussion
double innerfunction3D_REa_NoDiff(double x, double g)
{
    return ((1 / (x * x)) * (2 - x) * sqrt(1 + x)) * cos(-1 * x * g);
}

double innerfunction3D_REb_NoDiff(double x, double lambda, double g)
{
    if (x == 0)
    {
        return (3 / 4) * (1 / (lambda * lambda) - 1);
    }

    return ((1 / (x * x)) * (((2 - x) * sqrt(1 + x)) - ((2 - x / lambda) * sqrt(1 + x / lambda)))) * cos(-1 * x * g);
}
// ---------- end to set up integrands in Eq. 1 in CISSCO paper

double innerfunction3D_dREb(double x, double ph1, double ph2)
{
    if (x == 0)
    {
        return 0;
    }

    return ((1 / (x)) * (2 - (2 - x) * sqrt(1 + x)) * (sin(x * ph2) - sin(x * ph1)));
}

double innerfunction3D_dIMa(double x, double ph1)
{
    return ((1 / (x)) * (cos(x * ph1) - cos(2 * x * ph1)));
}
double innerfunction3D_dIMb(double x, double ph1, double ph2)
{

    return ((1 / (x)) * (2 - (2 - x) * sqrt(1 + x)) * (cos(x * ph2) - cos(x * ph1)));
}

// he functions below are used for calculating the integrals of the above equations

// ---------- Begin functions for calculating integrals

double trapzd3D_REa(const double a, const double b, const int n, double ph1)
{
    double x, tnm, sum, del;
    static double s;
    int it, j;

    if (n == 1)
    {
        return (s = 0.5 * (b - a) * (innerfunction3D_REa(a, ph1) + innerfunction3D_REa(b, ph1)));
    }
    else
    {
        for (it = 1, j = 1; j < n - 1; j++)
            it <<= 1;
        tnm = it;
        del = (b - a) / tnm;
        x = a + 0.5 * del;
        for (sum = 0.0, j = 0; j < it; j++, x += del)
            sum += innerfunction3D_REa(x, ph1);
        s = 0.5 * (s + (b - a) * sum / tnm);
        return s;
    }
}
double trapzd3D_REb(const double a, const double b, const int n, double p, double ph1, double ph2)
{
    double x, tnm, sum, del;
    static double s;
    int it, j;

    if (n == 1)
    {
        return (s = 0.5 * (b - a) * (innerfunction3D_REb(a, p, ph1, ph2) + innerfunction3D_REb(b, p, ph1, ph2)));
    }
    else
    {
        for (it = 1, j = 1; j < n - 1; j++)
            it <<= 1;
        tnm = it;
        del = (b - a) / tnm;
        x = a + 0.5 * del;
        for (sum = 0.0, j = 0; j < it; j++, x += del)
            sum += innerfunction3D_REb(x, p, ph1, ph2);
        s = 0.5 * (s + (b - a) * sum / tnm);
        return s;
    }
}

double trapzd3D_REa_NoDiff(const double a, const double b, const int n, double g)
{
    double x, tnm, sum, del;
    static double s;
    int it, j;

    if (n == 1)
    {
        return (s = 0.5 * (b - a) * (innerfunction3D_REa_NoDiff(a, g) + innerfunction3D_REa_NoDiff(b, g)));
    }
    else
    {
        for (it = 1, j = 1; j < n - 1; j++)
            it <<= 1;
        tnm = it;
        del = (b - a) / tnm;
        x = a + 0.5 * del;
        for (sum = 0.0, j = 0; j < it; j++, x += del)
            sum += innerfunction3D_REa_NoDiff(x, g);
        s = 0.5 * (s + (b - a) * sum / tnm);
        return s;
    }
}

double trapzd3D_REb_NoDiff(const double a, const double b, const int n, double lambda, double g)
{
    double x, tnm, sum, del;
    static double s;
    int it, j;

    if (n == 1)
    {
        return (s = 0.5 * (b - a) * (innerfunction3D_REb_NoDiff(a, lambda, g) + innerfunction3D_REb_NoDiff(b, lambda, g)));
    }
    else
    {
        for (it = 1, j = 1; j < n - 1; j++)
            it <<= 1;
        tnm = it;
        del = (b - a) / tnm;
        x = a + 0.5 * del;
        for (sum = 0.0, j = 0; j < it; j++, x += del)
            sum += innerfunction3D_REb_NoDiff(x, lambda, g);
        s = 0.5 * (s + (b - a) * sum / tnm);
        return s;
    }
}

double trapzd3D_REc_NoDiff(const double a, const double b, const int n, double g)
{
    double x, tnm, sum, del;
    static double s;
    int it, j;

    if (n == 1)
    {
        return (s = 0.5 * (b - a) * (innerfunction3D_REa_NoDiff(a, g) + innerfunction3D_REa_NoDiff(b, g)));
    }
    else
    {
        for (it = 1, j = 1; j < n - 1; j++)
            it <<= 1;
        tnm = it;
        del = (b - a) / tnm;
        x = a + 0.5 * del;
        for (sum = 0.0, j = 0; j < it; j++, x += del)
            sum += innerfunction3D_REa_NoDiff(x, g);
        s = 0.5 * (s + (b - a) * sum / tnm);
        return s;
    }
}

double trapzd3D_IMa(const double a, const double b, const int n, double ph1)
{
    double x, tnm, sum, del;
    static double s;
    int it, j;

    if (n == 1)
    {
        return (s = 0.5 * (b - a) * (innerfunction3D_IMa(a, ph1) + innerfunction3D_IMa(b, ph1)));
    }
    else
    {
        for (it = 1, j = 1; j < n - 1; j++)
            it <<= 1;
        tnm = it;
        del = (b - a) / tnm;
        x = a + 0.5 * del;
        for (sum = 0.0, j = 0; j < it; j++, x += del)
            sum += innerfunction3D_IMa(x, ph1);
        s = 0.5 * (s + (b - a) * sum / tnm);
        return s;
    }
}

double trapzd3D_IMb(const double a, const double b, const int n, double p, double ph1, double ph2)
{
    double x, tnm, sum, del;
    static double s;
    int it, j;

    if (n == 1)
    {
        return (s = 0.5 * (b - a) * (innerfunction3D_IMb(a, p, ph1, ph2) + innerfunction3D_IMb(b, p, ph1, ph2)));
    }
    else
    {
        for (it = 1, j = 1; j < n - 1; j++)
            it <<= 1;
        tnm = it;
        del = (b - a) / tnm;
        x = a + 0.5 * del;
        for (sum = 0.0, j = 0; j < it; j++, x += del)
            sum += innerfunction3D_IMb(x, p, ph1, ph2);
        s = 0.5 * (s + (b - a) * sum / tnm);
        return s;
    }
}

double trapzd3D_dREa(const double a, const double b, const int n, double ph1)
{
    double x, tnm, sum, del;
    static double s;
    int it, j;

    if (n == 1)
    {
        return (s = 0.5 * (b - a) * (innerfunction3D_dREa(a, ph1) + innerfunction3D_dREa(b, ph1)));
    }
    else
    {
        for (it = 1, j = 1; j < n - 1; j++)
            it <<= 1;
        tnm = it;
        del = (b - a) / tnm;
        x = a + 0.5 * del;
        for (sum = 0.0, j = 0; j < it; j++, x += del)
            sum += innerfunction3D_dREa(x, ph1);
        s = 0.5 * (s + (b - a) * sum / tnm);
        return s;
    }
}
double trapzd3D_dREb(const double a, const double b, const int n, double ph1, double ph2)
{
    double x, tnm, sum, del;
    static double s;
    int it, j;

    if (n == 1)
    {
        return (s = 0.5 * (b - a) * (innerfunction3D_dREb(a, ph1, ph2) + innerfunction3D_dREb(b, ph1, ph2)));
    }
    else
    {
        for (it = 1, j = 1; j < n - 1; j++)
            it <<= 1;
        tnm = it;
        del = (b - a) / tnm;
        x = a + 0.5 * del;
        for (sum = 0.0, j = 0; j < it; j++, x += del)
            sum += innerfunction3D_dREb(x, ph1, ph2);
        s = 0.5 * (s + (b - a) * sum / tnm);
        return s;
    }
}

double trapzd3D_dIMa(const double a, const double b, const int n, double ph1)
{
    double x, tnm, sum, del;
    static double s;
    int it, j;

    if (n == 1)
    {
        return (s = 0.5 * (b - a) * (innerfunction3D_dIMa(a, ph1) + innerfunction3D_dIMa(b, ph1)));
    }
    else
    {
        for (it = 1, j = 1; j < n - 1; j++)
            it <<= 1;
        tnm = it;
        del = (b - a) / tnm;
        x = a + 0.5 * del;
        for (sum = 0.0, j = 0; j < it; j++, x += del)
            sum += innerfunction3D_dIMa(x, ph1);
        s = 0.5 * (s + (b - a) * sum / tnm);
        return s;
    }
}

double trapzd3D_dIMb(const double a, const double b, const int n, double ph1, double ph2)
{
    double x, tnm, sum, del;
    static double s;
    int it, j;

    if (n == 1)
    {
        return (s = 0.5 * (b - a) * (innerfunction3D_dIMb(a, ph1, ph2) + innerfunction3D_dIMb(b, ph1, ph2)));
    }
    else
    {
        for (it = 1, j = 1; j < n - 1; j++)
            it <<= 1;
        tnm = it;
        del = (b - a) / tnm;
        x = a + 0.5 * del;
        for (sum = 0.0, j = 0; j < it; j++, x += del)
            sum += innerfunction3D_dIMb(x, ph1, ph2);
        s = 0.5 * (s + (b - a) * sum / tnm);
        return s;
    }
}

double qromb3D_REa(double a, double b, double ph1)
{
    const int JMAX = 20, JMAXP = JMAX + 1, K = 5;
    const double EPS = 0.00001;
    double ss, dss;
    double s[JMAX], h[JMAXP], s_t[K], h_t[K];
    int i, j;

    h[0] = 1.0;
    for (j = 1; j <= JMAX; j++)
    {
        s[j - 1] = trapzd3D_REa(a, b, j, ph1);
        if (j >= K)
        {
            for (i = 0; i < K; i++)
            {
                h_t[i] = h[j - K + i];
                s_t[i] = s[j - K + i];
            }
            polint(h_t, s_t, 0.0, ss, dss);
            if (fabs(dss) <= EPS * fabs(ss))
                return ss;
        }
        h[j] = 0.25 * h[j - 1];
    }
    //::MessageBox(NULL, "Too many steps in routine qromb", "Message", MB_OK);
    errorMessage_magMom = "Too many steps in routine qromb";
    return 0.0;
}
double qromb3D_REb(double a, double b, double p, double ph1, double ph2)
{
    const int JMAX = 20, JMAXP = JMAX + 1, K = 5;
    const double EPS = 0.00001;
    double ss, dss;
    double s[JMAX], h[JMAXP], s_t[K], h_t[K];
    int i, j;

    h[0] = 1.0;
    for (j = 1; j <= JMAX; j++)
    {
        s[j - 1] = trapzd3D_REb(a, b, j, p, ph1, ph2);
        if (j >= K)
        {
            for (i = 0; i < K; i++)
            {
                h_t[i] = h[j - K + i];
                s_t[i] = s[j - K + i];
            }
            polint(h_t, s_t, 0.0, ss, dss);
            if (fabs(dss) <= EPS * fabs(ss))
                return ss;
        }
        h[j] = 0.25 * h[j - 1];
    }
    //::MessageBox(NULL, "Too many steps in routine qromb", "Message", MB_OK);
    errorMessage_magMom = "Too many steps in routine qromb";
    return 0.0;
}

double qromb3D_REa_NoDiff(double a, double b, double g)
{
    const int JMAX = 20, JMAXP = JMAX + 1, K = 5;
    const double EPS = 0.00001;
    double ss, dss;
    double s[JMAX], h[JMAXP], s_t[K], h_t[K];
    int i, j;

    h[0] = 1.0;
    for (j = 1; j <= JMAX; j++)
    {
        s[j - 1] = trapzd3D_REa_NoDiff(a, b, j, g);
        if (j >= K)
        {
            for (i = 0; i < K; i++)
            {
                h_t[i] = h[j - K + i];
                s_t[i] = s[j - K + i];
            }
            polint(h_t, s_t, 0.0, ss, dss);
            if (fabs(dss) <= EPS * fabs(ss))
                return ss;
        }
        h[j] = 0.25 * h[j - 1];
    }
    //::MessageBox(NULL, "Too many steps in routine qromb", "Message", MB_OK);
    errorMessage_magMom = "Too many steps in routine qromb";
    return 0.0;
}

double qromb3D_REb_NoDiff(double a, double b, double lambda, double g)
{
    const int JMAX = 20, JMAXP = JMAX + 1, K = 5;
    const double EPS = 0.00001;
    double ss, dss;
    double s[JMAX], h[JMAXP], s_t[K], h_t[K];
    int i, j;

    h[0] = 1.0;
    for (j = 1; j <= JMAX; j++)
    {
        s[j - 1] = trapzd3D_REb_NoDiff(a, b, j, lambda, g);
        if (j >= K)
        {
            for (i = 0; i < K; i++)
            {
                h_t[i] = h[j - K + i];
                s_t[i] = s[j - K + i];
            }
            polint(h_t, s_t, 0.0, ss, dss);
            if (fabs(dss) <= EPS * fabs(ss))
                return ss;
        }
        h[j] = 0.25 * h[j - 1];
    }
    //::MessageBox(NULL, "Too many steps in routine qromb", "Message", MB_OK);
    errorMessage_magMom = "Too many steps in routine qromb";
    return 0.0;
}

double qromb3D_REc_NoDiff(double a, double b, double g)
{
    const int JMAX = 20, JMAXP = JMAX + 1, K = 5;
    const double EPS = 0.00001;
    double ss, dss;
    double s[JMAX], h[JMAXP], s_t[K], h_t[K];
    int i, j;

    h[0] = 1.0;
    for (j = 1; j <= JMAX; j++)
    {
        s[j - 1] = trapzd3D_REc_NoDiff(a, b, j, g);
        if (j >= K)
        {
            for (i = 0; i < K; i++)
            {
                h_t[i] = h[j - K + i];
                s_t[i] = s[j - K + i];
            }
            polint(h_t, s_t, 0.0, ss, dss);
            if (fabs(dss) <= EPS * fabs(ss))
                return ss;
        }
        h[j] = 0.25 * h[j - 1];
    }
    //::MessageBox(NULL, "Too many steps in routine qromb", "Message", MB_OK);
    errorMessage_magMom = "Too many steps in routine qromb";
    return 0.0;
}

double qromb3D_IMa(double a, double b, double ph1)
{
    const int JMAX = 20, JMAXP = JMAX + 1, K = 5;
    const double EPS = 0.00001;
    double ss, dss;
    double s[JMAX], h[JMAXP], s_t[K], h_t[K];
    int i, j;

    h[0] = 1.0;
    for (j = 1; j <= JMAX; j++)
    {
        s[j - 1] = trapzd3D_IMa(a, b, j, ph1);
        if (j >= K)
        {
            for (i = 0; i < K; i++)
            {
                h_t[i] = h[j - K + i];
                s_t[i] = s[j - K + i];
            }
            polint(h_t, s_t, 0.0, ss, dss);
            if (fabs(dss) <= EPS * fabs(ss))
                return ss;
        }
        h[j] = 0.25 * h[j - 1];
    }
    //::MessageBox(NULL, "Too many steps in routine qromb", "Message", MB_OK);
    errorMessage_magMom = "Too many steps in routine qromb";
    return 0.0;
}
double qromb3D_IMb(double a, double b, double p, double ph1, double ph2)
{
    const int JMAX = 20, JMAXP = JMAX + 1, K = 5;
    const double EPS = 0.00001;
    double ss, dss;
    double s[JMAX], h[JMAXP], s_t[K], h_t[K];
    int i, j;

    h[0] = 1.0;
    for (j = 1; j <= JMAX; j++)
    {
        s[j - 1] = trapzd3D_IMb(a, b, j, p, ph1, ph2);
        if (j >= K)
        {
            for (i = 0; i < K; i++)
            {
                h_t[i] = h[j - K + i];
                s_t[i] = s[j - K + i];
            }
            polint(h_t, s_t, 0.0, ss, dss);
            if (fabs(dss) <= EPS * fabs(ss))
                return ss;
        }
        h[j] = 0.25 * h[j - 1];
    }
    //::MessageBox(NULL, "Too many steps in routine qromb", "Message", MB_OK);
    errorMessage_magMom = "Too many steps in routine qromb";
    return 0.0;
}

double qromb3D_dREa(double a, double b, double ph1)
{
    const int JMAX = 20, JMAXP = JMAX + 1, K = 5;
    const double EPS = 0.00001;
    double ss, dss;
    double s[JMAX], h[JMAXP], s_t[K], h_t[K];
    int i, j;

    h[0] = 1.0;
    for (j = 1; j <= JMAX; j++)
    {
        s[j - 1] = trapzd3D_dREa(a, b, j, ph1);
        if (j >= K)
        {
            for (i = 0; i < K; i++)
            {
                h_t[i] = h[j - K + i];
                s_t[i] = s[j - K + i];
            }
            polint(h_t, s_t, 0.0, ss, dss);
            if (fabs(dss) <= EPS * fabs(ss))
                return ss;
        }
        h[j] = 0.25 * h[j - 1];
    }
    //::MessageBox(NULL, "Too many steps in routine qromb", "Message", MB_OK);
    errorMessage_magMom = "Too many steps in routine qromb";
    return 0.0;
}
double qromb3D_dREb(double a, double b, double ph1, double ph2)
{
    const int JMAX = 20, JMAXP = JMAX + 1, K = 5;
    const double EPS = 0.00001;
    double ss, dss;
    double s[JMAX], h[JMAXP], s_t[K], h_t[K];
    int i, j;

    h[0] = 1.0;
    for (j = 1; j <= JMAX; j++)
    {
        s[j - 1] = trapzd3D_dREb(a, b, j, ph1, ph2);
        if (j >= K)
        {
            for (i = 0; i < K; i++)
            {
                h_t[i] = h[j - K + i];
                s_t[i] = s[j - K + i];
            }
            polint(h_t, s_t, 0.0, ss, dss);
            if (fabs(dss) <= EPS * fabs(ss))
                return ss;
        }
        h[j] = 0.25 * h[j - 1];
    }
    //::MessageBox(NULL, "Too many steps in routine qromb", "Message", MB_OK);
    errorMessage_magMom = "Too many steps in routine qromb";
    return 0.0;
}

double qromb3D_dIMa(double a, double b, double ph1)
{
    const int JMAX = 20, JMAXP = JMAX + 1, K = 5;
    const double EPS = 0.00001;
    double ss, dss;
    double s[JMAX], h[JMAXP], s_t[K], h_t[K];
    int i, j;

    h[0] = 1.0;
    for (j = 1; j <= JMAX; j++)
    {
        s[j - 1] = trapzd3D_dIMa(a, b, j, ph1);
        if (j >= K)
        {
            for (i = 0; i < K; i++)
            {
                h_t[i] = h[j - K + i];
                s_t[i] = s[j - K + i];
            }
            polint(h_t, s_t, 0.0, ss, dss);
            if (fabs(dss) <= EPS * fabs(ss))
                return ss;
        }
        h[j] = 0.25 * h[j - 1];
    }
    //::MessageBox(NULL, "Too many steps in routine qromb", "Message", MB_OK);
    errorMessage_magMom = "Too many steps in routine qromb";
    return 0.0;
}
double qromb3D_dIMb(double a, double b, double ph1, double ph2)
{
    const int JMAX = 20, JMAXP = JMAX + 1, K = 5;
    const double EPS = 0.00001;
    double ss, dss;
    double s[JMAX], h[JMAXP], s_t[K], h_t[K];
    int i, j;

    h[0] = 1.0;
    for (j = 1; j <= JMAX; j++)
    {
        s[j - 1] = trapzd3D_dIMb(a, b, j, ph1, ph2);
        if (j >= K)
        {
            for (i = 0; i < K; i++)
            {
                h_t[i] = h[j - K + i];
                s_t[i] = s[j - K + i];
            }
            polint(h_t, s_t, 0.0, ss, dss);
            if (fabs(dss) <= EPS * fabs(ss))
                return ss;
        }
        h[j] = 0.25 * h[j - 1];
    }
    //::MessageBox(NULL, "Too many steps in routine qromb", "Message", MB_OK);
    errorMessage_magMom = "Too many steps in routine qromb";
    return 0.0;
}

// ---------- End functions for calculating integrals

// ---------- begin to calculate f_ij and their magnitude squares ---------- Do NOT change without discussion

double signalfunction3D(double p, double S1_S2, double S2_S3, double R1cube, double R2cube, double R3cube)
{
    double ph1 = p / R1cube;
    double ph2 = p / R2cube;
    double ph3 = p / R3cube;

    double ph12 = R1cube / R2cube;
    double ph23 = R2cube / R3cube;

    double REa = R1cube * qromb3D_REa(1, ph12, ph1) + qromb3D_REb(-1, 2, p, ph1, ph2);
    double IMa = R1cube * qromb3D_IMa(1, ph12, ph1) + qromb3D_IMb(-1, 2, p, ph1, ph2);

    double REb = R2cube * qromb3D_REa(1, ph23, ph2) + qromb3D_REb(-1, 2, p, ph2, ph3);
    double IMb = R2cube * qromb3D_IMa(1, ph23, ph2) + qromb3D_IMb(-1, 2, p, ph2, ph3);

    double A = (REa * REa) + (IMa * IMa);
    double B = (REb * REb) + (IMb * IMb);

    return (S1_S2)*B - (S2_S3)*A;
}
// ---------- end to calculate f_ij and their magnitude squares

// ---------- begin to set up Eq. 1 in CISSCO 2015 paper ---------- Do NOT change without discussion
double signalfunction3D_NoDiff(double g, double p, double rho, double R, double ReR)
{

    double lambda = p / (g * R * R * R);

    double S = ((4.0 * PI) / (9.0 * sqrt(3.0))) * rho * (p / g) * (qromb3D_REa_NoDiff(-1.0, -1.0 * lambda, g) + qromb3D_REb_NoDiff(-1.0 * lambda, 2.0 * lambda, lambda, g) + qromb3D_REc_NoDiff(2.0 * lambda, 2.0, g));

    return (S - ReR);
}
// ---------- end to set up Eq. 1 in CISSCO 2015 paper

double SIGN(double A, double B)
{
    if ((B > 0 && A > 0) || (B < 0 && A < 0))
        return A;
    else if ((B > 0 && A < 0) || (B < 0 && A > 0))
        return A * -1;
    else
        return 0;
}

double zbrent3D(const double x1, const double x2, const double tol, double S1_S2, double S2_S3, double R1cube, double R2cube, double R3cube)
{
    const int ITMAX = 100;
    const double EPS = 1.0e-10;
    int iter;
    // CString messagea = TEXT("");
    double a = x1, b = x2, c = x2, d, e, min1, min2;
    double fa = signalfunction3D(a, S1_S2, S2_S3, R1cube, R2cube, R3cube), fb = signalfunction3D(b, S1_S2, S2_S3, R1cube, R2cube, R3cube), fc, p, q, r, s, tol1, xm;

    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    {
        //::MessageBox(NULL, "No solution exists. Please verify center and radii.", "Message", MB_OK);
        errorMessage_magMom = "No solution exists. Please verify center and radii.";
        return 0.0;
    }
    fc = fb;
    for (iter = 0; iter < ITMAX; iter++)
    {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
        {
            c = a;
            fc = fa;
            e = d = b - a;
        }
        if (fabs(fc) < fabs(fb))
        {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        tol1 = 2.0 * EPS * fabs(b) + 0.5 * tol;
        xm = 0.5 * (c - b);
        if (fabs(xm) <= tol1 || fb == 0.0)
            return b;
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
        {
            s = fb / fa;
            if (a == c)
            {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            }
            else
            {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0)
                q = -q;
            p = fabs(p);
            min1 = 3.0 * xm * q - fabs(tol1 * q);
            min2 = fabs(e * q);
            if (2.0 * p < (min1 < min2 ? min1 : min2))
            {
                e = d;
                d = p / q;
            }
            else
            {
                d = xm;
                e = d;
            }
        }
        else
        {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        if (fabs(d) > tol1)
            b += d;
        else
            b += SIGN(tol1, xm);
        fb = signalfunction3D(b, S1_S2, S2_S3, R1cube, R2cube, R3cube);
    }
    //::MessageBox(NULL, "Maximum number of iterations exceeded in zbrent", "Message", MB_OK);
    errorMessage_magMom = "Maximum number of iterations exceeded in zbrent";
    return 0.0;
}

// ---------- begin to solve g using zbrent ---------- Do NOT change without discussion
double zbrent3D_NoDiff(const double x1, const double x2, const double tol, double mp, double rho, double R, double ReR)
{
    const int ITMAX = 100;
    const double EPS = 1.0e-10;
    int iter;
    // CString messagea = TEXT("");
    double a = x1, b = x2, c = x2, d, e, min1, min2;
    double fa = signalfunction3D_NoDiff(a, mp, rho, R, ReR), fb = signalfunction3D_NoDiff(b, mp, rho, R, ReR), fc, p, q, r, s, tol1, xm;

    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    {
        //::MessageBox(NULL, "No solution exists. Please verify center and radii.", "Message", MB_OK);
        errorMessage_magMom = "No solution exists. Please verify center and radii.";
        return 0.0;
    }
    fc = fb;
    for (iter = 0; iter < ITMAX; iter++)
    {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
        {
            c = a;
            fc = fa;
            e = d = b - a;
        }
        if (fabs(fc) < fabs(fb))
        {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        tol1 = 2.0 * EPS * fabs(b) + 0.5 * tol;
        xm = 0.5 * (c - b);
        if (fabs(xm) <= tol1 || fb == 0.0)
            return b;
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
        {
            s = fb / fa;
            if (a == c)
            {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            }
            else
            {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0)
                q = -q;
            p = fabs(p);
            min1 = 3.0 * xm * q - fabs(tol1 * q);
            min2 = fabs(e * q);
            if (2.0 * p < (min1 < min2 ? min1 : min2))
            {
                e = d;
                d = p / q;
            }
            else
            {
                d = xm;
                e = d;
            }
        }
        else
        {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        if (fabs(d) > tol1)
            b += d;
        else
            b += SIGN(tol1, xm);
        fb = signalfunction3D_NoDiff(b, mp, rho, R, ReR);
    }
    //::MessageBox(NULL, "Maximum number of iterations exceeded in zbrent", "Message", MB_OK);
    errorMessage_magMom = "Maximum number of iterations exceeded in zbrent";
    return 0.0;
}
// ---------- end to solve g using zbrent

// ---------- begin to calculate the magnetic moment ---------- Do NOT change without discussion

double FindMoment3D(double RES1, double RES2, double RES3, double IMS1, double IMS2, double IMS3, double R1, double R2, double R3)
{
    double R1cube = R1 * R1 * R1;
    double R2cube = R2 * R2 * R2;
    double R3cube = R3 * R3 * R3;

    double p1 = 0.001;
    double p2 = PI * R3cube;
    double tol = 0.0001;

    double RES2_S3 = RES2 - RES3;
    double RE1_S2 = RES1 - RES2;
    double IMS2_S3 = IMS2 - IMS3;
    double IM1_S2 = IMS1 - IMS2;
    double S2_S3 = (RES2_S3 * RES2_S3) + (IMS2_S3 * IMS2_S3);
    double S1_S2 = (RE1_S2 * RE1_S2) + (IM1_S2 * IM1_S2);

    return zbrent3D(p1, p2, tol, S1_S2, S2_S3, R1cube, R2cube, R3cube);
}
// ---------- end to calculate the magnetic moment

// ---------- begin to calculate the spin density and bkg phase ---------- Do NOT change without discussion

void CalculateSpinDensity(double RES2, double RES3, double IMS2, double IMS3, double R2, double R3, double *SpinDensity, double *BkgPhase)
{
    double p = m_MagMoment;
    double R2cube = R2 * R2 * R2;
    double R3cube = R3 * R3 * R3;

    double ph23 = R2cube / R3cube;
    double ph2 = p / R2cube;
    double ph3 = p / R3cube;

    double reS23 = RES2 - RES3;
    double imS23 = IMS2 - IMS3;
    double magS23 = sqrt((reS23 * reS23) + (imS23 * imS23));

    double ref23 = R2cube * qromb3D_REa(1, ph23, ph2) + qromb3D_REb(-1, 2, p, ph2, ph3);
    double imf23 = R2cube * qromb3D_IMa(1, ph23, ph2) + qromb3D_IMb(-1, 2, p, ph2, ph3);
    double magf23 = sqrt((ref23 * ref23) + (imf23 * imf23));

    *SpinDensity = (9 * sqrt(3.0) / (4.0 * PI)) * magS23 / magf23;

    // exp(i*BkgPhase) is proportional to (f*23)*(S23)
    double reBkgPhase = ref23 * reS23 + imf23 * imS23;
    double imBkgPhase = ref23 * imS23 - reS23 * imf23;

    *BkgPhase = atan2(imBkgPhase, reBkgPhase);
}
// ---------- end to calculate the spin density and bkg phase

// ---------- begin to calculate g and set up limits in zbrent ---------- Do NOT change without discussion
// All the _NoDiff functions below are for quantifying radius and susceptibility.
// If B0 = 7T, TE = 5 ms, and Delta Chi = 1 ppm, then g = 3.1. Thus choose p2 = 6.2
double Calculate_g(double p, double rho, double R, double ReR)
{

    double p2 = 6.2;             // 2.5*((R*R*R)/p); //upper limit of zbrent
    double p1 = p / (R * R * R); // lower limit of zbrent
    double tol = 0.0001;

    return zbrent3D_NoDiff(p1, p2, tol, p, rho, R, ReR);
}
// ---------- end to calculate g and set up limits in zbrent

double FindUncertainty3D(double RES1, double RES2, double RES3, double IMS1, double IMS2, double IMS3)
{

    double num1, num2, num, den, REf12, REf23, dREf12, dREf23,
        IMf12, IMf23, dIMf12, dIMf23, S1_S2, S2_S3, RES1_S2, RES2_S3, IMS1_S2, IMS2_S3, dS1_S2, dS2_S3, e12, e23;

    e12 = m_e12;
    e23 = m_e23;

    REf12 = (m_MagMoment / m_ROuterPhase) * qromb3D_REa(1, (m_RMiddlePhase / m_ROuterPhase), m_ROuterPhase) + qromb3D_REb(-1, 2, m_MagMoment, m_ROuterPhase, m_RMiddlePhase);
    REf23 = (m_MagMoment / m_RMiddlePhase) * qromb3D_REa(1, (m_RInnerPhase / m_RMiddlePhase), m_RMiddlePhase) + qromb3D_REb(-1, 2, m_MagMoment, m_RMiddlePhase, m_RInnerPhase);

    IMf12 = (m_MagMoment / m_ROuterPhase) * qromb3D_IMa(1, (m_RMiddlePhase / m_ROuterPhase), m_ROuterPhase) + qromb3D_IMb(-1, 2, m_MagMoment, m_ROuterPhase, m_RMiddlePhase);
    IMf23 = (m_MagMoment / m_RMiddlePhase) * qromb3D_IMa(1, (m_RInnerPhase / m_RMiddlePhase), m_RMiddlePhase) + qromb3D_IMb(-1, 2, m_MagMoment, m_RMiddlePhase, m_RInnerPhase);

    dREf12 = -2 * qromb3D_dREa(1, (m_RMiddlePhase / m_ROuterPhase), m_ROuterPhase) + qromb3D_dREb(-1, 2, m_ROuterPhase, m_RMiddlePhase);
    dREf23 = -2 * qromb3D_dREa(1, (m_RInnerPhase / m_RMiddlePhase), m_RMiddlePhase) + qromb3D_dREb(-1, 2, m_RMiddlePhase, m_RInnerPhase);

    dIMf12 = 2 * qromb3D_dIMa(1, (m_RMiddlePhase / m_ROuterPhase), m_ROuterPhase) + qromb3D_dIMb(-1, 2, m_ROuterPhase, m_RMiddlePhase);
    dIMf23 = 2 * qromb3D_dIMa(1, (m_RInnerPhase / m_RMiddlePhase), m_RMiddlePhase) + qromb3D_dIMb(-1, 2, m_RMiddlePhase, m_RInnerPhase);

    RES1_S2 = RES1 - RES2;
    RES2_S3 = RES2 - RES3;

    IMS1_S2 = IMS1 - IMS2;
    IMS2_S3 = IMS2 - IMS3;

    S1_S2 = (RES1_S2 * RES1_S2) + (IMS1_S2 * IMS1_S2);
    S2_S3 = (RES2_S3 * RES2_S3) + (IMS2_S3 * IMS2_S3);

    num1 = (2 / m_SNR) * (2 / m_SNR) * (19683 / (64 * PI * PI * PI)) * m_resx * m_resy * m_resz * ((m_R1 * m_R1 * m_R1 - m_R2 * m_R2 * m_R2) * (REf23 * REf23 + IMf23 * IMf23) * (REf23 * REf23 + IMf23 * IMf23) * (S1_S2) + (m_R2 * m_R2 * m_R2 - m_R3 * m_R3 * m_R3) * (REf12 * REf12 + IMf12 * IMf12) * (REf12 * REf12 + IMf12 * IMf12) * (S2_S3));
    num2 = (e23 * (REf12 * REf12 + IMf12 * IMf12) * S2_S3) * (e23 * (REf12 * REf12 + IMf12 * IMf12) * S2_S3) +
           (e12 * (REf23 * REf23 + IMf23 * IMf23) * S1_S2) * (e12 * (REf23 * REf23 + IMf23 * IMf23) * S1_S2);

    den = abs((2 * REf23 * dREf23 + 2 * IMf23 * dIMf23) * S1_S2 - (2 * REf12 * dREf12 + 2 * IMf12 * dIMf12) * S2_S3);

    num = (1 / m_MagMoment) * sqrt(num1 + num2);

    return num / den;
}

// This is the function used to calculate the magnetic moment in step 5
void OnBnClickedCalcmagmoment()
{
    // UpdateData(true);

    errorMessage_magMom = "";

    if ((m_R1 == 0) || (m_R2 == 0) || (m_R3 == 0))
    {
        //::MessageBox(NULL, "Cannot have values of 0 for radii.", "Message", MB_OK);
        errorMessage_magMom = "Cannot have values of 0 for radii.";
        return;
    }

    if ((m_CenterX == 0) || (m_CenterY == 0))
    {
        //::MessageBox(NULL, "Must find subpixel center first.", "Message", MB_OK);
        errorMessage_magMom = "Must find subpixel center first.";
        return;
    }

    if (SubpixelPhaseMatrix.empty())
    {
        //::MessageBox(NULL, "Must generate subpixel data first.", "Message", MB_OK);
        errorMessage_magMom = "Must generate subpixel data first.";
        return;
    }

    // ZoomedX = (int)(((m_CenterX - Xfirst + halfdisplay + 0.46) / 0.1));
    // ZoomedY = (int)(((m_CenterY - Yfirst + halfdisplay + 0.46) / 0.1));
    ZoomedX = pixelToSubpixel(m_CenterX, 0);
    ZoomedY = pixelToSubpixel(m_CenterY, 1);

    double R1, R2, R3, p, error, R_test, S_test, S_test2, ReR;
    R_test = m_R3;
    R1 = m_R1;
    R2 = m_R2;
    R3 = m_R3;

    int subpixRadius1 = round(R1 * m_SubPixels);
    int subpixRadius3 = round(R3 * m_SubPixels);
    int subpixRadius2 = round(R2 * m_SubPixels);
    int subpixRadius_test = round(R_test * m_SubPixels);
    double p0;

    // lastValueSlice = (int)(((m_CenterZ - Zfirst + halfdisplay + 0.46) / 0.1));
    lastValueSlice = pixelToSubpixel(m_CenterZ, 2);

    if (((ZoomedX + subpixRadius1) >= (2 * m_R0 * m_SubPixels)) || ((ZoomedX - subpixRadius1) <= 0))
    {
        // MessageBox("R1 is out of bounds.");
        errorMessage_magMom = "R1 is out of bounds.";
        return;
    }
    if (((ZoomedY + subpixRadius1) >= (2 * m_R0 * m_SubPixels)) || ((ZoomedY - subpixRadius1) <= 0))
    {
        // MessageBox("R1 is out of bounds.");
        errorMessage_magMom = "R1 is out of bounds.";
        return;
    }

    double RES1 = SumCircleElementsReal3D(subpixRadius1, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);
    double RES2 = SumCircleElementsReal3D(subpixRadius2, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);
    double RES3 = SumCircleElementsReal3D(subpixRadius3, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);

    double IMS1 = SumCircleElementsImag3D(subpixRadius1, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);
    double IMS2 = SumCircleElementsImag3D(subpixRadius2, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);
    double IMS3 = SumCircleElementsImag3D(subpixRadius3, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);

    //		double S_test, R_test;

    //		S_test=SumCircleElementsImag3D(subpixRadius2, (int)ZoomedX ,(int)ZoomedY, (int)lastValueSlice);
    //		S_test2=SumCircleElementsImag3D(subpixRadius3, (int)ZoomedX ,(int)ZoomedY, (int)lastValueSlice);

    // for finding epsilon
    // for (R_test=11;R_test<51;R_test++)
    //	S_test=SumCircleElementsReal3D(R_test, (int)ZoomedX ,(int)ZoomedY, (int)lastValueSlice);

    p = FindMoment3D(RES1, RES2, RES3, IMS1, IMS2, IMS3, R1, R2, R3);

    m_MagMoment = p;
    //		m_p.Format(TEXT("%3.2f"),p);
    m_p_last = p;

    if (p != 0)
    {
        m_RInnerPhase = p / (R3 * R3 * R3);
        m_RMiddlePhase = p / (R2 * R2 * R2);
        m_ROuterPhase = p / (R1 * R1 * R1);

        // m_R1CalcPhase.Format(TEXT("%3.2f"), m_ROuterPhase);
        // m_R2CalcPhase.Format(TEXT("%3.2f"), m_RMiddlePhase);
        // m_R3CalcPhase.Format(TEXT("%3.2f"), m_RInnerPhase);

        m_R1 = m_ROuterPhase;
        m_R2 = m_RMiddlePhase;
        m_R3 = m_RInnerPhase;

        CalculateSpinDensity(RES2, RES3, IMS2, IMS3, R2, R3, &rho, &BkgPhase);
        m_rho = rho;
        // m_p0.Format(TEXT("%3.2f"), m_rho);

        BackPhase = BkgPhase;
        m_BackPhase = BackPhase;
        // m_BackPhase.Format(TEXT("%3.2f"), BackPhase);

        if ((m_RChi != 0) && (m_B0 != 0) && (m_TE_first != 0))
        {
            ReR = SumCircleElementsReal3D(m_RChi, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);
            /*
            m_g = Calculate_g(p, m_rho, m_RChi, ReR);
            m_a = pow(m_p / m_g, 1.0 / 3.0);
            m_Chi = m_g / (0.08167 * m_B0 * m_TE);
            */
            m_g = Calculate_g(p, m_rho, m_RChi, ReR); // m_p_first is p at first echo time
            m_a = pow(p / m_g, 1.0 / 3.0);
            m_Chi = m_g / (0.08918167 * m_B0 * m_TE_first);
        }
        else
        {
            m_Chi = 0;
            m_a = 0;
        }

        if (m_SNR == 0)
        {
            // m_Uncertainty.Format(TEXT("Enter SNR"));
            m_Uncertainty = -1.0;
        }
        else
        {
            error = FindUncertainty3D(RES1, RES2, RES3, IMS1, IMS2, IMS3);
            error = error * 100;
            // m_Uncertainty.Format(TEXT("%3.2f"), error);
            m_Uncertainty = error;
        }
    }

    // UpdateData(false);
}

// This is the function used to get the imaginary sum within an inputted radius in step 5
void OnBnClickedImaginarysum()
{

    // UpdateData(true);

    errorMessage_sums = "";

    if (m_Ri == 0)
    {
        //::MessageBox(NULL, "Cannot have values of 0 for radii.", "Message", MB_OK);
        errorMessage_sums = "Cannot have values of 0 for radii.";
        return;
    }

    if ((m_CenterX == 0) || (m_CenterY == 0))
    {
        //::MessageBox(NULL, "Must find subpixel center first.", "Message", MB_OK);
        errorMessage_sums = "Must find subpixel center first.";
        return;
    }

    if (SubpixelPhaseMatrix.empty())
    {
        //::MessageBox(NULL, "Must generate subpixel data first.", "Message", MB_OK);
        errorMessage_sums = "Must generate subpixel data first.";
        return;
    }

    // ZoomedX = (int)(((m_CenterX - Xfirst + halfdisplay + 0.46) / 0.1));
    // ZoomedY = (int)(((m_CenterY - Yfirst + halfdisplay + 0.46) / 0.1));
    ZoomedX = pixelToSubpixel(m_CenterX, 0);
    ZoomedY = pixelToSubpixel(m_CenterY, 1);

    int subpixRadiusi = round(m_Ri * m_SubPixels);
    double Si;

    // ZoomedZ = (int)(((m_CenterZ - Zfirst + halfdisplay + 0.46) / 0.1));
    ZoomedZ = pixelToSubpixel(m_CenterZ, 2);
    Si = SumCircleElementsImag3D(subpixRadiusi, ZoomedX, ZoomedY, ZoomedZ);

    // m_Si.Format(TEXT("%3.4f"), Si);
    m_Si = Si;

    // UpdateData(false);
}

// This is the function used to get the rael sum within an inputted radius in step 5
void OnBnClickedRealsum()
{

    // UpdateData(true);

    errorMessage_sums = "";

    if (m_Ri == 0)
    {
        //::MessageBox(NULL, "Cannot have values of 0 for radii.", "Message", MB_OK);
        errorMessage_sums = "Cannot have values of 0 for radii.";
        return;
    }

    if ((m_CenterX == 0) || (m_CenterY == 0))
    {
        //::MessageBox(NULL, "Must find subpixel center first.", "Message", MB_OK);
        errorMessage_sums = "Must find subpixel center first.";
        return;
    }

    if (SubpixelPhaseMatrix.empty())
    {
        //::MessageBox(NULL, "Must generate subpixel data first.", "Message", MB_OK);
        errorMessage_sums = "Must generate subpixel data first.";
        return;
    }

    // ZoomedX = (int)(((m_CenterX - Xfirst + halfdisplay + 0.46) / 0.1));
    // ZoomedY = (int)(((m_CenterY - Yfirst + halfdisplay + 0.46) / 0.1));
    ZoomedX = (int)pixelToSubpixel(m_CenterX, 0);
    ZoomedY = (int)pixelToSubpixel(m_CenterY, 1);

    int subpixRadiusi = round(m_Ri * m_SubPixels);
    double Si2;

    // ZoomedZ = (int)(((m_CenterZ - Zfirst + halfdisplay + 0.46) / 0.1));
    ZoomedZ = (int)pixelToSubpixel(m_CenterZ, 2);

    Si2 = SumCircleElementsReal3D(subpixRadiusi, ZoomedX, ZoomedY, ZoomedZ);

    // m_Si2.Format(TEXT("%3.4f"), Si2);
    m_Si2 = Si2;

    // UpdateData(false);
}

// This is the function used to calculate g, a and chi in step 6
void OnBnClickedCalcSusceptibility()
{
    double ReR;
    ReR = SumCircleElementsReal3D(m_RChi, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);
    m_g = Calculate_g(m_p_first, m_rho, m_RChi, ReR); // m_p_first is p at first echo time
    m_a = pow(m_p_first / m_g, 1.0 / 3.0);
    m_Chi = m_g / (0.08918 * m_B0 * m_TE_first);
}

/*
These JNI functions do pretty much everything else except for pass data from Java to C++.
The functions below either let Java call the above C++ functions, or they give data from C++ back to Java
*/
/******************************************************************************************************************/
/******************************************************************************************************************/
/**********************************************    BEGIN    *******************************************************/
/**********************************************JNI FUNCTIONS*******************************************************/
/**********************************************             *******************************************************/
/******************************************************************************************************************/
/******************************************************************************************************************/

JNIEXPORT void JNICALL Java_JNIMethods_generateSubpixelArray(
    JNIEnv *env,     // Java environment
    jobject thisObj) // Reserved
{
    // generating subpixel
    OnBnClickedGenerateSubpixel();
}

JNIEXPORT jstring JNICALL Java_JNIMethods_estimateSubpixelCenter(JNIEnv *env, jobject thisObj)
{
    OnBnClickedEstimatecenter();

    const char *ermsg = errorMessage.c_str();

    return env->NewStringUTF(ermsg);
}

JNIEXPORT jstring JNICALL Java_JNIMethods_calculateMagneticMoment(JNIEnv *env, jobject thisObj)
{

    OnBnClickedCalcmagmoment();

    const char *ermsgP = errorMessage_magMom.c_str();

    return env->NewStringUTF(ermsgP);
}

// Java function to remove background phase - will be called when the button in the GUI is clicked
JNIEXPORT void JNICALL Java_JNIMethods_removeBackgroundPhase(JNIEnv *env, jobject thisObj, jdouble bPhase)
{

    removeBGPhaseAndInterpolateVoxels(bPhase);

    return;
}

JNIEXPORT jstring JNICALL Java_JNIMethods_calculateRealSum(JNIEnv *env, jobject thisObj)
{

    OnBnClickedRealsum();

    const char *erMsgS = errorMessage_sums.c_str();

    return env->NewStringUTF(erMsgS);
}

JNIEXPORT jstring JNICALL Java_JNIMethods_calculateImagSum(JNIEnv *env, jobject thisObj)
{

    OnBnClickedImaginarysum();

    const char *erMsgS = errorMessage_sums.c_str();

    return env->NewStringUTF(erMsgS);
}

JNIEXPORT void JNICALL Java_JNIMethods_estBkgAndSpinDensity(JNIEnv *env, jobject thisObj)
{

    ZoomedX = pixelToSubpixel(m_CenterX, 0);
    ZoomedY = pixelToSubpixel(m_CenterY, 1);
    lastValueSlice = pixelToSubpixel(m_CenterZ, 2);

    int subpixRadius2 = round(m_R2 * m_SubPixels);
    int subpixRadius3 = round(m_R3 * m_SubPixels);

    /*
        double RES2 = SumCircleElementsReal3D(subpixRadius2, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);
        double RES3 = SumCircleElementsReal3D(subpixRadius3, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);

        double IMS2 = SumCircleElementsImag3D(subpixRadius2, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);
        double IMS3 = SumCircleElementsImag3D(subpixRadius3, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);
    */

    double RES2 = SumCircleElementsReal3D_S5(subpixRadius2, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);
    double RES3 = SumCircleElementsReal3D_S5(subpixRadius3, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);

    double IMS2 = SumCircleElementsImag3D_S5(subpixRadius2, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);
    double IMS3 = SumCircleElementsImag3D_S5(subpixRadius3, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);

    CalculateSpinDensity(RES2, RES3, IMS2, IMS3, m_R2, m_R3, &rho, &BkgPhase);

    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_calcSusceptibility(JNIEnv *env, jobject thisObj)
{
    OnBnClickedCalcSusceptibility();
    return;
}

JNIEXPORT void JNICALL Java_JNIMethods_interpolateVoxelsSIM(JNIEnv *env, jobject thisObj, jint subpixelsize)
{
    SubpixelSimulatedRealMatrix3D.resize(subpixelsize, vector<vector<float>>(subpixelsize, vector<float>(subpixelsize, 0)));

    interpolateVoxels_SIM(subpixelsize / 10);
    return;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_SumCircleElementsReal3DSIMMED(JNIEnv *env, jobject thisObj, jint radius, jint x, jint y, jint z)
{
    return SumCircleElementsReal3D_SIM(radius, x, y, z);
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_equation10(JNIEnv *env, jobject thisObj, jdouble jp, jdouble jphi_i, jdouble jphi_j)
{
    return (jp / jphi_i) * qromb3D_REa(1, jphi_j / jphi_i, jphi_i) + qromb3D_REb(-1, 2, jp, jphi_i, jphi_j);
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_calculateUncertainty(JNIEnv *env, jobject thisObj, jdouble je12, jdouble je23)
{
    m_e12 = je12;
    m_e23 = je23;

    double R1 = m_R1;
    double R2 = m_R2;
    double R3 = m_R3;

    int subpixRadius1 = round(R1 * m_SubPixels);
    int subpixRadius2 = round(R2 * m_SubPixels);
    int subpixRadius3 = round(R3 * m_SubPixels);

    double RES1 = SumCircleElementsReal3D(subpixRadius1, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);
    double RES2 = SumCircleElementsReal3D(subpixRadius2, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);
    double RES3 = SumCircleElementsReal3D(subpixRadius3, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);

    double IMS1 = SumCircleElementsImag3D(subpixRadius1, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);
    double IMS2 = SumCircleElementsImag3D(subpixRadius2, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);
    double IMS3 = SumCircleElementsImag3D(subpixRadius3, (int)ZoomedX, (int)ZoomedY, (int)lastValueSlice);

    m_Uncertainty = FindUncertainty3D(RES1, RES2, RES3, IMS1, IMS2, IMS3);
    return m_Uncertainty;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getSubX(JNIEnv *env, jobject thisObj)
{
    return m_CenterX;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getSubY(JNIEnv *env, jobject thisObj)
{
    return m_CenterY;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getSubZ(JNIEnv *env, jobject thisObj)
{
    return m_CenterZ;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getSubXOther(JNIEnv *env, jobject thisObj)
{
    return m_CenterX3;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getSubYOther(JNIEnv *env, jobject thisObj)
{
    return m_CenterY3;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getSubZOther(JNIEnv *env, jobject thisObj)
{
    return m_CenterZ3;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getMR1Calc(JNIEnv *env, jobject thisObj)
{
    return m_ROuterPhase;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getMR2Calc(JNIEnv *env, jobject thisObj)
{
    return m_RMiddlePhase;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getMR3Calc(JNIEnv *env, jobject thisObj)
{
    return m_RInnerPhase;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getSNR(JNIEnv *env, jobject thisObj)
{
    return m_SNR;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getE12(JNIEnv *env, jobject thisObj)
{
    return m_e12;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getE23(JNIEnv *env, jobject thisObj)
{
    return m_e23;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getB0(JNIEnv *env, jobject thisObj)
{
    return m_B0;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getTE(JNIEnv *env, jobject thisObj)
{
    return m_TE_first;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getRChi(JNIEnv *env, jobject thisObj)
{
    return m_RChi;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getRho(JNIEnv *env, jobject thisObj)
{
    return m_rho;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getChi(JNIEnv *env, jobject thisObj)
{
    return m_Chi;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getA(JNIEnv *env, jobject thisObj)
{
    return m_a;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getUncertainty(JNIEnv *env, jobject thisObj)
{
    return m_Uncertainty;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getP(JNIEnv *env, jobject thisObj)
{
    return m_p_first;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getP0(JNIEnv *env, jobject thisObj)
{
    return m_p0;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getMagMoment(JNIEnv *env, jobject thisObj)
{
    return m_MagMoment;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getResX(JNIEnv *env, jobject thisObj)
{
    return m_resx;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getResY(JNIEnv *env, jobject thisObj)
{
    return m_resy;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getResZ(JNIEnv *env, jobject thisObj)
{
    return m_resz;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getRealSum(JNIEnv *env, jobject thisObj)
{
    return m_Si2;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getImagSum(JNIEnv *env, jobject thisObj)
{
    return m_Si;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getBkg(JNIEnv *env, jobject thisObj)
{
    return BkgPhase;
}

JNIEXPORT jdouble JNICALL Java_JNIMethods_getSpinDensity(JNIEnv *env, jobject thisObj)
{
    return rho;
}

/******************************************************************************************************************/
/******************************************************************************************************************/
/**********************************************     END     *******************************************************/
/**********************************************JNI FUNCTIONS*******************************************************/
/**********************************************             *******************************************************/
/******************************************************************************************************************/
/******************************************************************************************************************/