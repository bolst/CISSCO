REM build.cmd

@echo off

@echo cleaning...
del /q bin\*.*
del /q ext\ImageJ\plugins\CISSCO\*.*

echo compiling java files[1/2]...
javac -cp lib/ij.jar;lib/ml.jar src/java/*.java
if %ERRORLEVEL% neq 0 (
    echo Failed to compile java files
    exit /b 1
)

echo compiling java files[2/2]...
javah -cp src/java JNIMethods
if %ERRORLEVEL% neq 0 (
    echo Failed to compile java files
    exit /b 1
)

move /y JNIMethods.h src/cpp >NUL 2>NUL

echo compiling cpp files[1/2]...
g++ -c -I"%JAVA_HOME%/include" -I"%JAVA_HOME%/include/win32" -m64 -fPIC src/cpp/Calculate_Magnetic_Moment_3D.cpp -o Calculate_Magnetic_Moment_3D.o
if %ERRORLEVEL% neq 0 (
    echo Failed to compile cpp files
    exit /b 1
)

echo compiling cpp files[2/2]...
g++ -shared -m64 -o Calculate_Magnetic_Moment_3D_Native.dll Calculate_Magnetic_Moment_3D.o -Wl,-add-stdcall-alias
if %ERRORLEVEL% neq 0 (
    echo Failed to compile cpp files
    exit /b 1
)

echo done!

move /y src\java\*.class bin >NUL 2>NUL
move /y Calculate_Magnetic_Moment_3D_Native.dll lib >NUL 2>NUL
move /y Calculate_Magnetic_Moment_3D.o bin >NUL 2>NUL

copy /y bin\* ext\ImageJ\plugins\CISSCO\ >NUL 2>NUL
copy /y lib\ml.jar ext\ImageJ\plugins\CISSCO\ >NUL 2>NUL
copy /y lib\Calculate_Magnetic_Moment_3D_Native.dll ext\ImageJ\plugins\CISSCO\ >NUL 2>NUL