/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class JNIMethods */

#ifndef _Included_JNIMethods
#define _Included_JNIMethods
#ifdef __cplusplus
extern "C" {
#endif
/*
 * Class:     JNIMethods
 * Method:    setmVariables
 * Signature: (IDDDDDD)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_setmVariables
  (JNIEnv *, jobject, jint, jdouble, jdouble, jdouble, jdouble, jdouble, jdouble);

/*
 * Class:     JNIMethods
 * Method:    setMagMoment
 * Signature: (D)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_setMagMoment
  (JNIEnv *, jobject, jdouble);

/*
 * Class:     JNIMethods
 * Method:    setBackPhase
 * Signature: (D)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_setBackPhase
  (JNIEnv *, jobject, jdouble);

/*
 * Class:     JNIMethods
 * Method:    setRealImagNumbers
 * Signature: ([[[F[[[F)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_setRealImagNumbers
  (JNIEnv *, jobject, jobjectArray, jobjectArray);

/*
 * Class:     JNIMethods
 * Method:    setXYZ
 * Signature: (DDD)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_setXYZ
  (JNIEnv *, jobject, jdouble, jdouble, jdouble);

/*
 * Class:     JNIMethods
 * Method:    setPhaseXYMatrix
 * Signature: ([[F)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_setPhaseXYMatrix
  (JNIEnv *, jobject, jobjectArray);

/*
 * Class:     JNIMethods
 * Method:    setPhaseXZMatrix
 * Signature: ([[F)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_setPhaseXZMatrix
  (JNIEnv *, jobject, jobjectArray);

/*
 * Class:     JNIMethods
 * Method:    setMagXYMatrix
 * Signature: ([[F)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_setMagXYMatrix
  (JNIEnv *, jobject, jobjectArray);

/*
 * Class:     JNIMethods
 * Method:    setMagXZMatrix
 * Signature: ([[F)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_setMagXZMatrix
  (JNIEnv *, jobject, jobjectArray);

/*
 * Class:     JNIMethods
 * Method:    setSmallBox
 * Signature: (IIIIII)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_setSmallBox
  (JNIEnv *, jobject, jint, jint, jint, jint, jint, jint);

/*
 * Class:     JNIMethods
 * Method:    setCenterL
 * Signature: (DDD)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_setCenterL
  (JNIEnv *, jobject, jdouble, jdouble, jdouble);

/*
 * Class:     JNIMethods
 * Method:    setCenterM
 * Signature: (DDD)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_setCenterM
  (JNIEnv *, jobject, jdouble, jdouble, jdouble);

/*
 * Class:     JNIMethods
 * Method:    setCenterS
 * Signature: (DDD)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_setCenterS
  (JNIEnv *, jobject, jdouble, jdouble, jdouble);

/*
 * Class:     JNIMethods
 * Method:    setmR123
 * Signature: (DDD)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_setmR123
  (JNIEnv *, jobject, jdouble, jdouble, jdouble);

/*
 * Class:     JNIMethods
 * Method:    setR123PhaseCalc
 * Signature: (DDD)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_setR123PhaseCalc
  (JNIEnv *, jobject, jdouble, jdouble, jdouble);

/*
 * Class:     JNIMethods
 * Method:    setMagMomentVariables
 * Signature: (DDDDDD)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_setMagMomentVariables
  (JNIEnv *, jobject, jdouble, jdouble, jdouble, jdouble, jdouble, jdouble);

/*
 * Class:     JNIMethods
 * Method:    setRi
 * Signature: (D)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_setRi
  (JNIEnv *, jobject, jdouble);

/*
 * Class:     JNIMethods
 * Method:    setStep6Variables
 * Signature: (DDDD)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_setStep6Variables
  (JNIEnv *, jobject, jdouble, jdouble, jdouble, jdouble);

/*
 * Class:     JNIMethods
 * Method:    setSimulatedMatrices
 * Signature: ([[[FI)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_setSimulatedMatrices
  (JNIEnv *, jobject, jobjectArray, jint);

/*
 * Class:     JNIMethods
 * Method:    removeBackgroundPhase
 * Signature: (D)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_removeBackgroundPhase
  (JNIEnv *, jobject, jdouble);

/*
 * Class:     JNIMethods
 * Method:    generateSubpixelArray
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_JNIMethods_generateSubpixelArray
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    estimateSubpixelCenter
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_JNIMethods_estimateSubpixelCenter
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    calculateRealSum
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_JNIMethods_calculateRealSum
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    calculateImagSum
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_JNIMethods_calculateImagSum
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    estBkgAndSpinDensity
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_JNIMethods_estBkgAndSpinDensity
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    calcSusceptibility
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_JNIMethods_calcSusceptibility
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    interpolateVoxelsSIM
 * Signature: (I)V
 */
JNIEXPORT void JNICALL Java_JNIMethods_interpolateVoxelsSIM
  (JNIEnv *, jobject, jint);

/*
 * Class:     JNIMethods
 * Method:    SumCircleElementsReal3DSIMMED
 * Signature: (IIII)D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_SumCircleElementsReal3DSIMMED
  (JNIEnv *, jobject, jint, jint, jint, jint);

/*
 * Class:     JNIMethods
 * Method:    equation10
 * Signature: (DDD)D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_equation10
  (JNIEnv *, jobject, jdouble, jdouble, jdouble);

/*
 * Class:     JNIMethods
 * Method:    calculateUncertainty
 * Signature: (DD)D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_calculateUncertainty
  (JNIEnv *, jobject, jdouble, jdouble);

/*
 * Class:     JNIMethods
 * Method:    getSpinDensity
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getSpinDensity
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getBkg
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getBkg
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getSubX
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getSubX
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getSubY
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getSubY
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getSubZ
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getSubZ
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getSubXOther
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getSubXOther
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getSubYOther
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getSubYOther
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getSubZOther
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getSubZOther
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    calculateMagneticMoment
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_JNIMethods_calculateMagneticMoment
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getMR1Calc
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getMR1Calc
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getMR2Calc
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getMR2Calc
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getMR3Calc
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getMR3Calc
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getSNR
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getSNR
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getE12
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getE12
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getE23
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getE23
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getB0
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getB0
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getTE
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getTE
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getRChi
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getRChi
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getRho
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getRho
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getChi
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getChi
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getA
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getA
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getUncertainty
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getUncertainty
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getP
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getP
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getP0
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getP0
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getMagMoment
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getMagMoment
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getResX
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getResX
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getResY
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getResY
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getResZ
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getResZ
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getRealSum
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getRealSum
  (JNIEnv *, jobject);

/*
 * Class:     JNIMethods
 * Method:    getImagSum
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_JNIMethods_getImagSum
  (JNIEnv *, jobject);

#ifdef __cplusplus
}
#endif
#endif
