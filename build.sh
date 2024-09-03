# build.sh
#!/bin/bash

mkdir -p bin

echo "cleaning..."
rm -f bin/*.*
rm -f ext/ImageJ/plugins/CISSCO/*.*

echo "compiling java files[1/2]..."
javac -cp lib/ij.jar:lib/ml.jar src/java/*.java
if [ $? -ne 0 ]; then
    echo "Failed to compile java files"
    exit 1
fi

echo "compiling java files[2/2]..."
javah -cp src/java JNIMethods
if [ $? -ne 0 ]; then
    echo "Failed to compile java files"
    exit 1
fi

mv JNIMethods.h src/cpp 2>/dev/null

echo "compiling cpp files[1/2]..."
g++ -c -I"$JAVA_HOME/include" -I"$JAVA_HOME/include/linux" -m64 -std=c++11 -fPIC src/cpp/Calculate_Magnetic_Moment_3D.cpp -o Calculate_Magnetic_Moment_3D.o
if [ $? -ne 0 ]; then
    echo "Failed to compile cpp files"
    exit 1
fi

echo "compiling cpp files[2/2]..."
g++ -shared -m64 -o Calculate_Magnetic_Moment_3D_Native.so Calculate_Magnetic_Moment_3D.o
if [ $? -ne 0 ]; then
    echo "Failed to compile cpp files"
    exit 1
fi

mv src/java/*.class bin 2>/dev/null
mv Calculate_Magnetic_Moment_3D_Native.so lib 2>/dev/null
mv Calculate_Magnetic_Moment_3D.o bin 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Unable to move files to bin & lib"
    exit 1
fi

mkdir -p ext/ImageJ.app/plugins/CISSCO
cp bin/* ext/ImageJ.app/plugins/CISSCO/ 2>/dev/null
cp lib/ml.jar ext/ImageJ.app/plugins/CISSCO/ 2>/dev/null
cp lib/Calculate_Magnetic_Moment_3D_Native.so ext/ImageJ.app/plugins/CISSCO/ 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Unable to copy files to ImageJ.app/plugins/CISSCO"
    exit 1
fi

echo "done"

echo "Running ImageJ..."
open ./ext/ImageJ.app
