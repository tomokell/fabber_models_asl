# Building the ASL models #

The current version of the ASL models are designed to build with the latest version of Fabber. Although it should be possible to build them against FSL 5.0 you 
may need to edit the source slightly. So your first step should be to build fabber-core from Git. Build instructions for fabber-core are in the Wiki on its 
Git page. You will make the process straightforward if you create a directory (e.g. 'fabber') and download the Git repositories for fabber-core and 
fabber-models-cest into this directory. Then the ASL models build should be able to find the Fabber libraries without you having to install them or tell the 
build scripts where they are.

The ASL models use cmake as the build tool. CMake is designed for out-of-source builds, so you create a separate build directory and all the compiled files 
end up there. CMake is installed on many Linux distributions by default, or can easily be added. It is also readily available for OSX and Windows.

If you want to build the ASL models within an existing FSL source tree, see [Below](#fslbuild)

You need to ensure that `FSLDIR` is set to point to wherever the FSL dependencies are installed, following that the basic steps to make a build are:

    mkdir build
    cd build
    cmake ..
    make

After the `cmake` command, information on dependencies will be displayed, for example:

    -- FSL headers in /home/martinc/dev/fsl/include /home/martinc/dev/fsl/extras/include/newmat /home/martinc/dev/fsl/extras/include /home/martinc/dev/fsl/extras/include/boost
    -- Fabber headers in /home/martinc/dev/fabber_core
    -- Using Fabber libraries: /home/martinc/dev/fabber_core/Debug/libfabbercore.a /home/martinc/dev/fabber_core/Debug/libfabberexec.a
    -- Using libznz: /home/martinc/dev/fsl/lib/libznz.a
    -- Using libutils: /home/martinc/dev/fsl/lib/libutils.a 
    -- Using miscmaths: /home/martinc/dev/fsl/lib/libmiscmaths.a
    -- Using fslio: /home/martinc/dev/fsl/lib/libfslio.a
    -- Using newimage: /home/martinc/dev/fsl/lib/libnewimage.a
    -- Using niftiio: /home/martinc/dev/fsl/lib/libniftiio.a
        -- Using newmat: /home/martinc/dev/fsl/extras/lib/libnewmat.a /home/martinc/dev/fsl/extras/include/newmat
    -- Using newimage: /home/martinc/dev/fsl/lib/libnewimage.a
    -- Using prob: /home/martinc/dev/fsl/extras/lib/libprob.a
    -- Using zlib: /usr/lib/x86_64-linux-gnu/libz.so

Check that these locations seem reasonable. In particular, make sure the Fabber headers and libraries are for the latest version of Fabber which you have probably just built. If it is picking up the FSL-5.0 version you may need to edit the file CMakeLists.txt and re-run cmake.

After running `make`, if all goes well, the build should conclude with the message

    [  9%] Building CXX object CMakeFiles/fabber_asl.dir/fwdmodel_asl_rest.cc.o
    [ 18%] Building CXX object CMakeFiles/fabber_asl.dir/fwdmodel_asl_grase.cc.o
    [ 27%] Building CXX object CMakeFiles/fabber_asl.dir/fwdmodel_asl_multiphase.cc.o
    [ 36%] Building CXX object CMakeFiles/fabber_asl.dir/asl_models.cc.o
    [ 45%] Building CXX object CMakeFiles/fabber_asl.dir/fabber_client.cc.o
    [ 54%] Linking CXX executable fabber_asl
    [ 54%] Built target fabber_asl
    Scanning dependencies of target fabber_models_asl
    [ 63%] Building CXX object CMakeFiles/fabber_models_asl.dir/fwdmodel_asl_rest.cc.o
    [ 72%] Building CXX object CMakeFiles/fabber_models_asl.dir/fwdmodel_asl_grase.cc.o
    [ 81%] Building CXX object CMakeFiles/fabber_models_asl.dir/fwdmodel_asl_multiphase.cc.o
    [ 90%] Building CXX object CMakeFiles/fabber_models_asl.dir/asl_models.cc.o
    [100%] Linking CXX shared library libfabber_models_asl.so
    [100%] Built target fabber_models_asl

### It failed with something about `recompile with -fPIC`

In order to build the shared library, the FSL libraries have to contain 'Position-independent code'. If they don't you can't build the library and will need to 
use the `fabber_asl` executable instead. The only way around this is to recompile the FSL libraries with the -fPIC flag but this may be a challenge. 
You're better off just running `make fabber_asl` which will just build the executable. 

### Successful build?

Two objects are built:

1. An executable fabber_asl, which is a version of the Fabber executable with the ASL models built in
2. A shared library libfabber_models_asl.so, which can be loaded in to the generic Fabber executable using the --loadmodels option.

You may use whichever of these you prefer, however if you are intending to use the command line tool you will probably find the fabber_asl executable more convenient. The shared library is intended for when you want to embed the Fabber library in another application.

You can verify that the ASL models are available whichever method you choose:

    fabber_asl --listmodels

    asl_multiphase
    aslrest
    buxton
    linear
    poly
    trivial

or

    fabber --loadmodels=libfabber_models_asl.so --listmodels
    
    asl_multiphase
    aslrest
    buxton
    linear
    poly
    trivial
  
Note that you still need to specify the model (e.g. --model=aslrest) on the command line to tell Fabber to use one of the ASL models.

### <a name="fslbuild"></a>Building in an FSL environment

If you have built `fabber_core` within an FSL source distribution which you may prefer to build the ASL models within the same system. You will need to do the following:

1. Move the fabber_models_asl directory into the FSL source tree `mv fabber_models_asl $FSLDIR/src/`
2. `cd $FSLDIR/src/fabber_models_asl`
3. make

`fabber_asl` should now be set up to build, provided you have the following environment variables set up:

    FSLDIR
    FSLCONFDIR
    FSLMACHTYPE

Note that the FSL makefile does not attempt to build the shared library.