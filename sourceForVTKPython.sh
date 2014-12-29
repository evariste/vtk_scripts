# see http://www.itk.org/Wiki/VTK/Tutorials/PythonEnvironmentSetup


#vtkBuildDir=/Users/paul/work/packages/vtk/build-5.6-py
#vtkBuildDir=/Users/paul/work/packages/vtk/build-5.10-py
vtkBuildDir="/Users/paulaljabar/work/packages/vtk/build-5.10.1-64"

export PATH=${vtkBuildDir}/bin:${PATH}

export PYTHONPATH=${vtkBuildDir}/Wrapping/Python/:${vtkBuildDir}/bin

export LD_LIBRARY_PATH=${vtkBuildDir}/bin

# then run vtkpython

# Can call import vtk within


