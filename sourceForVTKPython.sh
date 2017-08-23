#!/bin/bash

# see http://www.itk.org/Wiki/VTK/Tutorials/PythonEnvironmentSetup


# vtkBuildDir=/Users/paul/work/packages/vtk/build-5.6-py
# vtkBuildDir=/Users/paul/work/packages/vtk/build-5.10-py
# vtkBuildDir="/Users/paulaljabar/work/packages/vtk/build-5.10.1-64"
vtkBuildDir="/Users/paulaljabar/work/packages/vtk/build-6.3.0"
vtkBuildDir='/d/Users_D/paul.aljabar.MIRADA-MEDICAL/software/pa-builds/VTK-6.3.0-SRC/install'
vtkBuildDir='/d/software/pa-builds/VTK-6.3.0-SRC/install'

vtk_bin_dir="$vtkBuildDir/bin"
vtk_lib_dir="$vtkBuildDir/lib"
vtk_lib_py_dir="${vtk_lib_dir}/python2.7/site-packages"


allPaths="${vtkBuildDir}:${vtk_bin_dir}:${vtk_lib_dir}:${vtk_lib_py_dir}"
export PATH="${allPaths}:${PATH}"

# export PYTHONPATH=${vtkBuildDir}/Wrapping/Python/:${vtkBuildDir}/bin:${vtkBuildDir}/lib:${PYTHONPATH}
export PYTHONPATH="${allPaths};${PYTHONPATH}"

#export LD_LIBRARY_PATH=${vtkBuildDir}/lib
export LD_LIBRARY_PATH="${vtk_lib_dir}:${vtk_bin_dir}:${PYTHONPATH}"

#export DYLD_LIBRARY_PATH=${vtkBuildDir}/lib

# then run vtkpython

# Can call import vtk within


