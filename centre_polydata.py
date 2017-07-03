import sys

py_script_path = 'D:/Users_D/paul.aljabar.MIRADA-MEDICAL/scripts/python'
sys.path.append(py_script_path)

from env_for_vtk import *

import vtk

import argparse

import numpy as np

#
# def trace(frame, event, arg):
#     print "%s, %s:%d" % (event, frame.f_code.co_filename, frame.f_lineno)
#     return trace
#
#
# sys.settrace(trace)

def main():

    helpText = (" Stuff")
    parser = argparse.ArgumentParser(description=helpText)

    helpText = "Input vtk file [filename]"
    parser.add_argument('input', type=str, help=helpText)

    helpText = "Output vtk file [filename]"
    parser.add_argument('output', type=str, help=helpText)

    args = parser.parse_args()


    inputName  = args.input
    outputName = args.output


    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(inputName)
    reader.Update()


    pd = reader.GetOutput()


    xlo, xhi, ylo, yhi, zlo, zhi = pd.GetBounds()

    xcent = 0.5 * (xlo + xhi)
    ycent = 0.5 * (ylo + yhi)
    zcent = 0.5 * (zlo + zhi)

    offset = np.asarray((-xcent, -ycent, -zcent))

    n_pts = pd.GetNumberOfPoints()


    for n in range(n_pts):
        pt = np.asarray(pd.GetPoint(n))
        pd.GetPoints().SetPoint(n, pt + offset)


    pdWriter = vtk.vtkPolyDataWriter()
    pdWriter.SetFileName(outputName)
    pdWriter.SetInputData(pd)
    # pdWriter.Update()
    pdWriter.Write()

    print ('xxx')
    #
    # for p in sys.path:
    #     print p.replace('\\','/')
    # return 0


if __name__ == '__main__':
    retval = main()
    print('done')
    # sys.exit(retval)

