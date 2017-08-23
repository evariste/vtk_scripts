import dicom
import sys
import os
import numpy as np

import argparse

import pyclipper

from PIL import Image, ImageDraw


os.sys.path.insert(0, 'D:/scripts/python')

from dicom_keys import *

# import env_for_vtk
#
# # import vtk
# from vtk import vtkPolyData, vtkCellArray, vtkPoints, vtkPolygon, vtkPolyDataWriter, vtkGlyph3D, vtkSphereSource, vtkTriangleFilter, vtkDecimatePro


# Functions to extract roi numbers/names
getROINos = lambda x: int(x[kROINumber].value)
getROINames = lambda x: x[kROIName].value

def clean_roi_names(roiNames):
    for c in [',', ' ', '-', '=']:
        roiNames = map(lambda x: x.replace(c, '_'), roiNames)
    return roiNames


############################################################################################


def get_contour_point_data(dcm_rtss_file, target_labels, verbose=False):
    """

    :param dcm_rtss_file:
    :param target_labels: List of possible variations of the name.
    :type target_labels list
    :return:
    """

    if verbose:
        print('Extracting organ from file:')
        print('    {:s}'.format(dcm_rtss_file))

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
    k_target = None

    target_labels_lower = map(lambda x: x.lower(), target_labels)

    for k_target, s in enumerate(contourSequences):
        roiRefNo = int(s[kContSeqRefROINo].value)
        roiRefName = roiNumber2Name[roiRefNo] # type: str
        if roiRefName.lower() in target_labels_lower:
            if verbose:
                print('Found structure "{:s}" in dicom RTSS at index {:d}'.format(roiRefName, k_target))
            found = True
            break

    if not found:
        print('Structure names:')
        print(target_labels)
        raise Exception('Cannot find structure in DICOM RTSS data with a name matching this list')

    s = contourSequences[k_target]

    assert  ( kContSeqRefROINo in s.keys() ) and ( kROIContSeq in s.keys() )

    contSeq = s[kROIContSeq]

    n_slices = len(contSeq.value)

    all_pts = []
    n_pts_total = 0

    for contour in contSeq:
        n_pts = contour[kNoOfContPts].value
        contData = contour[kContData].value
        contData = map(float, contData)
        xs = np.asarray(contData).reshape((-1, 3))
        assert xs.shape[0] == n_pts

        n_pts_total += n_pts
        all_pts.append(xs)


    assert len(all_pts) == n_slices

    return all_pts, n_pts_total


def draw(cnt_pts):

    N_pix = 400

    pad = 10

    im = Image.new('RGBA', (N_pix, N_pix), (255, 255, 255, 0))
    draw = ImageDraw.Draw(im)

    k = 0

    N = len(cnt_pts)

    pts = cnt_pts[N//2]

    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]


    x_min, x_max = np.min(xs), np.max(xs)

    y_min, y_max = np.min(ys), np.max(ys)

    x_range = x_max - x_min
    y_range = y_max - y_min


    p_prev = pts[0]
    for p in pts[1:]:
        x_prev = p_prev[0] - x_min
        y_prev = p_prev[1] - y_min
        x = p[0] - x_min
        y = p[1] - y_min

        x_prev = pad + int(np.round( (N_pix - 2*pad) * x_prev / x_range))
        y_prev = pad + int(np.round(  (N_pix - 2*pad) * y_prev / y_range ))
        x = pad + int(np.round( (N_pix - 2*pad) * x / x_range))
        y = pad + int(np.round(  (N_pix - 2*pad) * y / y_range ))
        draw.line((x_prev, y_prev, x, y), fill=128)

        p_prev = p



    im.show()

############################################################################################

def main():

    helpText = ('Help text')
    parser = argparse.ArgumentParser(description=helpText)

    helpText = 'DICOM RTSS file for prediction'
    parser.add_argument('-dcm_file_pred', type=str, help=helpText, required=True)

    helpText = 'DICOM RTSS file for ground truth'
    parser.add_argument('-dcm_file_gt', type=str, help=helpText, required=True)

    helpText = 'Target organ name in prediction file'
    parser.add_argument('-organ_pred', type=str, help=helpText, required=True)

    helpText = 'Target organ name in ground truth file'
    parser.add_argument('-organ_gt', type=str, help=helpText, required=True, action='append')

    helpText = 'Verbose'
    parser.add_argument('-verbose', help=helpText, action='store_true')

    helpText = 'Spacing in dimension Z (slice thickness, default = 1.0)'
    parser.add_argument('-z', help=helpText, type=float, default=1.0)
    args = parser.parse_args()

    helpText = 'Spacing in dimension X (default = 1.0)'
    parser.add_argument('-x', help=helpText, type=float, default=1.0)
    args = parser.parse_args()

    helpText = 'Spacing in dimension Y (default = 1.0)'
    parser.add_argument('-y', help=helpText, type=float, default=1.0)
    args = parser.parse_args()

    dcm_file_pred = args.dcm_file_pred
    dcm_file_gt = args.dcm_file_gt

    target_label_pred = args.organ_pred
    target_labels_gt = args.organ_gt

    verbose = args.verbose # type: bool


    sp_x = args.x
    sp_y = args.y
    sp_z = args.z



    voxel_vol = sp_x * sp_y * sp_z


    dcm_file_basename = os.path.basename(dcm_file_pred)

    cnts_pred, n_pts_pred = get_contour_point_data(dcm_file_pred, [target_label_pred], verbose=verbose)

    cnts_gt, n_pts_gt = get_contour_point_data(dcm_file_gt, target_labels_gt, verbose=verbose)

    draw(cnts_pred)


    z_gt = [x[0][2] for x in cnts_gt]
    z_pred = [x[0][2] for x in cnts_pred]

    z_all = np.sort(  np.unique(z_gt + z_pred) )

    # Store the in-slice 2-D coordinates of each slice's contour in a dictionary.

    polys_pred = dict.fromkeys(z_all)
    polys_gt = dict.fromkeys(z_all)

    for contour in cnts_pred:
        z = contour[0][2]
        pts_2d = [(np.round(x[0]), np.round(x[1])) for x in contour]
        polys_pred[z] = pts_2d

    # Slices where there is no prediction.
    for z in z_all:
        if polys_pred[z] is None:
            polys_pred[z] = []

    for contour in cnts_gt:
        z = contour[0][2]
        pts_2d = [(np.round(x[0]), np.round(x[1])) for x in contour]
        polys_gt[z] = pts_2d

    # Slices where there is no ground truth.
    for z in z_all:
        if polys_gt[z] is None:
            polys_gt[z] = []


    # Store the areas ...

    areas_pred = dict.fromkeys(z_all)
    areas_gt = dict.fromkeys(z_all)
    areas_intersection = dict.fromkeys(z_all)
    areas_union = dict.fromkeys(z_all)

    for z in z_all:

        p_pred = polys_pred[z]
        p_gt = polys_gt[z]

        if p_pred == []:
            areas_pred[z] = 0
        else:
            areas_pred[z] = pyclipper.Area(p_pred)

        if p_gt == []:
            areas_gt[z] = 0
        else:
            areas_gt[z] = pyclipper.Area(p_gt)

        if p_gt == []:
            areas_intersection[z] = 0
            areas_union[z] = areas_pred[z]
            continue

        if p_pred == []:
            areas_intersection[z] = 0
            areas_union[z] = areas_gt[z]
            continue

        # Both prediction and GT contour are present.

        pc = pyclipper.Pyclipper()

        pc.AddPath(p_pred, pyclipper.PT_SUBJECT, True)
        pc.AddPath(p_gt, pyclipper.PT_CLIP, True)

        # A set of paths defining the regions where the prediction and GT intersect.
        intersections = pc.Execute(pyclipper.CT_INTERSECTION, pyclipper.PFT_EVENODD, pyclipper.PFT_EVENODD)

        areas_intersection[z] = 0
        for region in intersections:
            areas_intersection[z] += pyclipper.Area(region)


        union = pc.Execute(pyclipper.CT_UNION, pyclipper.PFT_EVENODD, pyclipper.PFT_EVENODD)

        areas_union[z] = 0
        for region in union:
            areas_union[z] += pyclipper.Area(region)

    area_total_pred = np.sum(areas_pred.values())
    area_total_gt = np.sum(areas_gt.values())
    area_total_intersection = np.sum(areas_intersection.values())
    area_mean = 0.5 * (area_total_pred + area_total_gt)
    dice = area_total_intersection / area_mean

    if verbose:
        print('{:<15s} {:<10s} {:>7s} {:>7s} {:>6s}'.format('File', 'Organ', 'AP', 'AG', 'DSC'))

    print('{:<15s} {:<10s} {:>4.1f} {:>4.1f} {:>6.4f} '.format(dcm_file_basename,
                                                                       target_label_pred,
                                                                       area_total_pred*voxel_vol,
                                                                       area_total_gt*voxel_vol,
                                                                       dice))

    if verbose:

        print('Z: slice')
        print('AP: Area prediction')
        print('AG: Area GT')
        print('AI: Area intersection')
        print('AU: Area union')
        print('CH: Check AP+AG-AI == AU')
        print('\n{:>5s} : {:>6s} {:>6s} {:>6s}  {:>6s}    {:>6s}'.format('Z', 'AP', 'AG', 'AI', 'AU', 'CH'))
        for z in z_all:

            ap = areas_pred[z]
            ag = areas_gt[z]
            ai = areas_intersection[z]
            au = areas_union[z]

            print('{:>4.1f} : {:>6.1f} {:>6.1f} {:>6.1f}  {:>6.1f}    {:>6.1f}'.format(z, ap, ag, ai, au, ap + ag - ai))



    return  0



#####################################################################

if __name__ == '__main__':
  sys.exit(main())

