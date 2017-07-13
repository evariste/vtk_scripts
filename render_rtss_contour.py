import dicom
import sys
import os
import numpy as np

import argparse

os.sys.path.insert(0, 'D:/Users_D/paul.aljabar.MIRADA-MEDICAL/scripts/python')

from dicom_keys import *

import env_for_vtk
import vtk

# TODO: Write tubes around the contour in each slice, treat as polyline?

# Functions to extract roi numbers/names
getROINos = lambda x: int(x[kROINumber].value)
getROINames = lambda x: x[kROIName].value

def clean_roi_names(roiNames):
    for c in [',', ' ', '-', '=']:
        roiNames = map(lambda x: x.replace(c, '_'), roiNames)
    return roiNames




def main():

    helpText = ('Render RTSS contours for a given organ. Two outputs are generated, one showing per-slice polygons for the contour and one showing glyphs over the points')
    parser = argparse.ArgumentParser(description=helpText)

    helpText = "Input DICOM RTSS file"
    parser.add_argument('-dcm_file', type=str, help=helpText, required=True)

    helpText = "Target organ name"
    parser.add_argument('-organ', type=str, help=helpText, required=True)

    helpText = "Output vtk file [filename]"
    parser.add_argument('-output', type=str, help=helpText, required=True)

    helpText = "Scale factor for glyph output, default = 1, increase or decrease to suit"
    parser.add_argument('-glyph_scale_factor', type=float, default=1.0, help=helpText)

    helpText = "Maximum number of points for glyph output, default = 2500, increase or decrease to suit"
    parser.add_argument('-glyph_max_points', type=int, default=2500, help=helpText)

    args = parser.parse_args()

    dcm_rtss_file = args.dcm_file
    output_polygons = args.output # type: str
    target_label = args.organ

    user_scale = args.glyph_scale_factor
    max_pts = args.glyph_max_points

    if output_polygons[-4:] == '.vtk':
        output_glyphs = output_polygons.replace('.vtk', '-glyphs.vtk')
    else:
        output_glyphs = output_polygons + '-glyphs.vtk'


    # dcm_rtss_file = 'D:/Users_D/paul.aljabar.MIRADA-MEDICAL/data/Nijmegen/PROSTATE/prostate_3_anon/104836515_4259415179/20160923/1/1.2.826.0.1.3680043.2.135.736415.63184520.7.1490799097.913.96.dcm'
    # target_label = 'Prostaat'

    # dcm_rtss_file = 'D:/Users_D/paul.aljabar.MIRADA-MEDICAL/projects/aapm_challenge/predictions/stage2_full/iter_last/ct_LCTSC-Test-S2-104-heart.dcm'
    # dcm_rtss_file = 'D:/Users_D/paul.aljabar.MIRADA-MEDICAL/projects/aapm_challenge/predictions/stage2_full/iter_last/dcm_single/LCTSC-Test-S2-104.dcm'
    # target_label = 'esoph'
    # output_polygons = 'pol.vtk'
    # output_glyphs = 'pol2.vtk'

    # user_scale = 1
    # max_pts = 2500


    print('Extracting organ {:s} from file:'.format(target_label))
    print('    {:s}'.format(dcm_rtss_file))
    print('Saving result to file:')
    print('    {:s}'.format(output_polygons))
    print('    {:s}'.format(output_glyphs))



    if not os.path.exists(dcm_rtss_file):
        raise Exception('Cannot locate dicom file')

    dcm = dicom.read_file(dcm_rtss_file)

    if not( dcm[kModality].value == 'RTSTRUCT' ):
        raise Exception('DICOM file does not seem to represent RTSTRUCT data')


    structSetROIs = dcm[kStructSetROI]

    roiNumbers = map(getROINos, structSetROIs)
    roiNames = map(getROINames, structSetROIs)

    roiNames = clean_roi_names(roiNames)

    roiNumber2Name = dict(zip(roiNumbers, roiNames))




    contourSequences = dcm[kROIContSeqs]

    found = False
    for k_target, s in enumerate(contourSequences):
        roiRefNo = int(s[kContSeqRefROINo].value)
        roiRefName = roiNumber2Name[roiRefNo]
        if roiRefName == target_label:
            print('Found structure "{:s}" in dicom RTSS at index {:d}'.format(target_label, k_target))
            found = True
            break

    if not found:
        raise Exception('Cannot find structure {:s} in DICOM RTSS data'.format(target_label))


    s = contourSequences[k_target]

    assert  ( kContSeqRefROINo in s.keys() ) and ( kROIContSeq in s.keys() )

    contSeq = s[kROIContSeq]

    n_slices = len(contSeq.value)

    all_pts = []
    n_pts_total = 0

    for x in contSeq:
        n_pts = x[kNoOfContPts].value
        contData = x[kContData].value
        contData = map(float, contData)
        xs = np.asarray(contData).reshape((-1, 3))
        assert xs.shape[0] == n_pts

        n_pts_total += n_pts
        all_pts.append(xs)


    assert len(all_pts) == n_slices


    cellArr = vtk.vtkCellArray()
    cellArr.SetNumberOfCells(n_slices)
    cellArr.InitTraversal()

    pts = vtk.vtkPoints()

    pts.SetNumberOfPoints(n_pts_total)

    ix_pt = 0

    for k in range(n_slices):
        pol = vtk.vtkPolygon()
        xs = all_pts[k]
        n_pts = xs.shape[0]
        pol.GetPointIds().SetNumberOfIds(n_pts)

        for n in range(n_pts):
            pol.GetPointIds().SetId(n, ix_pt)
            x, y, z = xs[n]
            pts.InsertPoint(ix_pt, x, y, z)
            ix_pt += 1

        cellArr.InsertNextCell(pol)

    # Polydata containing polygons
    pd = vtk.vtkPolyData()
    pd.SetPoints(pts)
    pd.SetPolys(cellArr)


    # Write polygon output:
    pdWr = vtk.vtkPolyDataWriter()
    pdWr.SetInputData(pd)
    pdWr.SetFileName(output_polygons)
    pdWr.Update()
    pdWr.Write()


    # Generate glyphed version

    # Bounds
    xlo, xhi, ylo, yhi, zlo, zhi = pd.GetBounds()
    box_vol = (xhi-xlo) * (yhi - ylo) * (zhi - zlo)
    r = np.power(box_vol, 1/3.0) * 3.0 / 4.0  / np.pi

    sph_scale_factor = 25.0 / user_scale

    # The glyph
    sph = vtk.vtkSphereSource()
    sph.SetThetaResolution(6)
    sph.SetPhiResolution(6)
    sph.SetRadius(r / sph_scale_factor)
    sph.Update()

    gl = vtk.vtkGlyph3D()


    # Decimate if necessary:
    n_pts = pd.GetNumberOfPoints()

    if n_pts > max_pts:
        # The input to the glyphs filter is the result of decimating the original data
        red_fac = (n_pts - max_pts) / float(n_pts)

        tri_filt = vtk.vtkTriangleFilter()
        tri_filt.SetInputData(pd)
        tri_filt.Update()

        dec_filt = vtk.vtkDecimatePro()
        dec_filt.SetTargetReduction(red_fac)
        dec_filt.SetInputData(tri_filt.GetOutput())
        dec_filt.Update()

        gl.SetInputData(dec_filt.GetOutput())
    else:
        # Pass original data to glyph filter
        gl.SetInputData(pd)


    gl.SetSourceData(sph.GetOutput())
    gl.SetScaleModeToScaleByScalar()
    gl.SetScaleFactor(1)
    gl.Update()

    wr2 = vtk.vtkPolyDataWriter()
    wr2.SetInputData(gl.GetOutput())
    wr2.SetFileName(output_glyphs)
    wr2.Update()
    wr2.Write()


    return  0

    #
    # # SURFACE RECONSTRUCTION. DOESN'T REALLY WORK
    #
    # pd2 = vtk.vtkPolyData()
    # pd2.SetPoints(pts)
    #
    # srecon = vtk.vtkSurfaceReconstructionFilter()
    # srecon.SetSampleSpacing(0.7)
    # srecon.SetNeighborhoodSize(4)
    # srecon.SetInputData(pd2)
    # srecon.Update()
    #
    # cf = vtk.vtkContourFilter()
    # cf.SetInputData(srecon.GetOutput())
    # cf.SetValue(0, 0.0)
    # cf.Update()
    #
    # hf = vtk.vtkFillHolesFilter()
    # hf.SetInputData(cf.GetOutput())
    # hf.SetHoleSize(70)
    # hf.Update()
    #
    # nf = vtk.vtkPolyDataNormals()
    # nf.SetInputData(hf.GetOutput())
    # nf.Update()
    #
    # pdWr2 = vtk.vtkPolyDataWriter()
    # pdWr2.SetInputData(nf.GetOutput())
    # pdWr2.SetFileName('pol2.vtk')
    # pdWr2.Update()
    # pdWr2.Write()
    #



if __name__ == '__main__':
  sys.exit(main())

