

source ~/scripts/sourceForVTKPython.sh

# image=target/mr/mr-reset-roi-resamp.nii
image=SOME_IMAGE_TO_DEFINE_DOF_ROI.NII.GZ

spacing=20

ffdcreate random.dof
ffdadd $image  random.dof random.dof -ds $spacing

sizeSigma=2
blurSigma=20

ffdrandom random.dof random.dof -sigma $sizeSigma  -blur $blurSigma

# Make structured grid
ffd2vtk random.dof random-ffd-sg.vtk

# Polydata if wanted:
vtkpython ~/scripts/vtkStructuredGridToPolydata.py random-ffd-sg.vtk random-ffd-pd.vtk
