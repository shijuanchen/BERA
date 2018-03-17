#----------------------------------------------------------------------------------------------------------------------#
# This script summarize BERA functions in Point_clouds_processing_v1.0
#----------------------------------------------------------------------------------------------------------------------#
# This script creates buffer shapefiles of the centerline of seismic lines, centerline is defined by the start and end point recorded in rtk txt.

from osgeo import ogr, osr

class rtk_crd(object):
    def __init__(self, site_id, point_id, x, y, ter_h):
        self.site_id = site_id
        self.point_id = point_id
        self.x = x
        self.y = y
        self.ter_h = ter_h

def create_clip_shp(rtk_folder_path, shp_folder_path, bufferDistance):
    if not os.path.exists(rtk_folder_path):
        print "Could not find ", rtk_folder_path, "!"
        return

    if not os.path.exists(shp_folder_path):
        os.makedirs(shp_folder_path)

    rtk_file_list = glob.glob(os.path.join(rtk_folder_path, '*.txt'))
    file_num = len(rtk_file_list)
    print "There are ", file_num, "files in ", rtk_folder_path

    for k in range(0, file_num):
        fr = open(rtk_file_list[k], 'r')
        start_pt = rtk_crd('0', '0', 0, 0, 0)
        end_pt = rtk_crd('0', '0', 0, 0, 0)
        row = fr.readline().split(',')
        site_id = row[0][0] + row[0][1] + row[0][2]
        print site_id
        fr.seek(0) #newly added
        for line in fr:
            row = line.split(',')
            point_id = row[0]
            x = row[2]
            y = row[1]
            ter_h = row[3]
            if 'L000' in line:
                start_pt = rtk_crd(str(site_id), str(point_id), float(x), float(y), float(ter_h))

            if 'L150' in line:
                end_pt = rtk_crd(str(site_id), str(point_id), float(x), float(y), float(ter_h))

            else:
                end_pt = rtk_crd(str(site_id), str(point_id), float(x), float(y), float(ter_h))

        print start_pt.site_id, start_pt.point_id, start_pt.x, start_pt.y, start_pt.ter_h
        print end_pt.site_id, end_pt.point_id, end_pt.x, end_pt.y, end_pt.ter_h

        Cline = ogr.Geometry(ogr.wkbLineString)
        Cline.AddPoint(start_pt.x, start_pt.y)
        Cline.AddPoint(end_pt.x, end_pt.y)
        print Cline

        cl_wkt = str(Cline)
        cl = ogr.CreateGeometryFromWkt(cl_wkt)

        clip_bf = cl.Buffer(bufferDistance)
        clip_bf_wkt = clip_bf.ExportToWkt()
        print "%s buffered by %d is %s" % (cl.ExportToWkt(), bufferDistance, clip_bf_wkt)

        # see https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html
        # set up the shapefile driver
        driver = ogr.GetDriverByName('ESRI Shapefile')

        # create the data source, CreateDataSource(directory of output shapefile(folder path))
        layer_name = 'clip_bf_' + str(site_id)
        file_name = shp_folder_path + '\\' + layer_name + '.shp'
        print file_name
        if os.path.exists(file_name):
            driver.DeleteDataSource(file_name)

        data_source = driver.CreateDataSource(shp_folder_path)

        # create the spatial reference: NAD83 UTM Z12, EPSG Projection 26912
        srf = osr.SpatialReference()
        srf.ImportFromEPSG(26912)

        # create the layer
        layer = data_source.CreateLayer(layer_name, srf, ogr.wkbPolygon)

        # create the feature
        feature = ogr.Feature(layer.GetLayerDefn())

        # Set the feature geometry using the polygon
        feature.SetGeometry(clip_bf)
        layer.CreateFeature(feature)
        # Destroy the feature to free resources
        feature.Destroy()
        # Destroy the data source to free resources
        data_source.Destroy()
        fr.close()

#----------------------------------------------------------------------------------------------------------------------#
# This script writes a bat. file that clips las files by shp files
# Input: 1.Folder path of las. files. 2.Folder path of bat. file to save. 3. lastool path 4. output las. folder	5. version number Output:bat file
# lastool_path = r'E:\LAStools\LAStools\LAStools\bin\lasclip.exe'

def clip_las_bat(lastool_path, las_ifolder_path, shp_folder_path, las_ofolder_path, version, bat_file_path):

    if not os.path.exists(las_ifolder_path):
        print "Could not find ", las_ifolder_path, "!"
        return

    las_file_list = glob.glob(os.path.join(las_ifolder_path, '*.las'))
    shp_file_list = glob.glob(os.path.join(shp_folder_path, '*.shp'))

    file_num = len(las_file_list)
    print "There are ", file_num, "files in", las_ifolder_path

    f = open(bat_file_path, 'wb')

    f.write('rem @echo off\n')

    for k in range(0, file_num):
        old_file_name = os.path.basename(las_file_list[k]).split('.')[0]
        old_file_name = old_file_name.split('v')[0]
        las_ofile_path = las_ofolder_path + '\\' + old_file_name + 'v' + str(version) + '.las'
        line = lastool_path + ' ' + '-i' + ' ' + las_file_list[k] + ' ' + '-poly ' + shp_file_list[k] + ' -o' + ' ' + las_ofile_path + '\n'
        f.write(line)

    f.write('PAUSE')

    f.close()

#----------------------------------------------------------------------------------------------------------------------#
# The function of this script is to link ground truth of vegetation height to RTK coordinates by matching points ID.
# Link_Gr_Veg_Height v1.6 Shijuan Chen 2016/09/13
# Input data requirements: 1. veg sheets are csv files in a common folder 2. RTK points are saved in txt files in a common folder 3. veg files were arranged in the same order as RTK txt files.

import csv


class rtk_crd(object):
    def __init__(self, site_id, point_id, x, y, ter_h):
        self.site_id = site_id
        self.point_id = point_id
        self.x = x
        self.y = y
        self.ter_h = ter_h


class veg(object):
    def __init__(self, site_id, point_id, gr_veg_h):
        self.site_id = site_id
        self.point_id = point_id
        self.gr_veg_h = gr_veg_h


class rtk_veg(rtk_crd):
    def __init__(self, site_id, point_id, x, y, ter_h, gr_veg_h):
        super(rtk_veg, self).__init__(site_id, point_id, x, y, ter_h)
        self.gr_veg_h = gr_veg_h


def link_gr_veg_h(rtk_folder_path, veg_folder_path, output_folder_path):
    if not os.path.exists(veg_folder_path):
        print 'Could not find ', veg_folder_path, '!'
        if not os.path.exists(rtk_folder_path):
            print 'Could not find ', rtk_folder_path, '!'
            return
        return

    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    #read csv

    veg_file_list = glob.glob(os.path.join(veg_folder_path,'*.csv'))

    #read RTK txt
    rtk_file_list = glob.glob(os.path.join(rtk_folder_path,'*.txt'))
    site_number = len(rtk_file_list)

    print "There are", site_number, "sites."

    for index in range(0, site_number):
        with open(veg_file_list[index],'r') as fv:
            veg_reader = csv.reader(fv)
            row1 = next(veg_reader)
            v_list = list(veg_reader)
            site_id =  v_list[0][0][0] + v_list[0][0][1] + v_list[0][0][2]

            veg_pt_number = len(v_list)
            print "There are", veg_pt_number, 'veg points in file', veg_file_list[index]
            veg_list = []
            for i in range(0, veg_pt_number):
                point_id = v_list[i][0]
                gr_veg_h = v_list[i][1]
                veg_list.append(veg(site_id, point_id, gr_veg_h))
        fv.close()

        with open(rtk_file_list[index], 'r') as fr:
            rtk_line = fr.read().splitlines()
            rtk_list = []
            for line in rtk_line:
                row = line.split(',')
                site_id = row[0][0] + row[0][1] + row[0][2]
                point_id = row[0]
                x = row[2]
                y = row[1]
                ter_h = row[3]
                rtk_list.append(rtk_crd(site_id, point_id, x, y, ter_h))
            rtk_pt_number = len(rtk_list)
            print "There are", rtk_pt_number, "rtk points in file", rtk_file_list[index]
        fr.close()

        # matching by point_id
        rtk_veg_list = []
        output_file_path = output_folder_path + "\\" + site_id + '_rtk_veg.txt'
        with open(output_file_path, 'w') as fw:
            fw.write('site_id,point_id,x,y,ter_h,gr_veg_h\n')
            for i in  range(0, rtk_pt_number):
                for j in range(0, veg_pt_number):
                    if  (rtk_list[i].point_id == veg_list[j].point_id) and (rtk_list[i].site_id == veg_list[j].site_id) and \
                        (not(int(float(veg_list[j].gr_veg_h)) == 999)) and (not(int(float(veg_list[j].gr_veg_h)) == -1)):
                        gr_veg_h = float(veg_list[j].gr_veg_h) / 100.0
                        rtk_veg_list.append(rtk_veg(rtk_list[i].site_id, rtk_list[i].point_id, rtk_list[i].x, rtk_list[i].y, rtk_list[i].ter_h, gr_veg_h))
                        rtk_veg_number = len(rtk_veg_list)

                        line_list = [rtk_list[i].site_id, rtk_list[i].point_id, rtk_list[i].x, rtk_list[i].y, rtk_list[i].ter_h, str(gr_veg_h)]
                        line = ','.join(line_list)
                        fw.writelines(line)
                        fw.write('\n')
            print "There are", rtk_veg_number, "points in file", output_file_path
            fw.close()

        print "Output file was saved in", output_file_path

#----------------------------------------------------------------------------------------------------------------------#
# this script compares height of pt in las file and veg height measured in field
from liblas import file

import BERA_classes as bc

def cmp_veg_ht(las_ifolder_path, txt_ifolder_path, txt_ofolder_path):

    if not os.path.exists(las_ifolder_path):
        print "Could not find ", las_ifolder_path, "!"
        if not os.path.exists(txt_ifolder_path):
            print "Could not find ", txt_ifolder_path, "!"
            return
        return
    if not os.path.exists(txt_ofolder_path):
        os.makedirs(txt_ofolder_path)

    las_file_list = glob.glob(os.path.join(las_ifolder_path, '*.las'))
    txt_file_list = glob.glob(os.path.join(txt_ifolder_path, '*.txt'))

    las_file_num = len(las_file_list)
    txt_file_num = len(txt_file_list)
    if not (las_file_num == txt_file_num):
        print "Number of las and txt files does not match!"
        return

    for k in range(0, las_file_num):
        print las_file_num
        with open(txt_file_list[k]) as fr_txt:
            with open(las_file_list[k]) as fr_las:
               # read las
                pt_3d_las = []
                fr_las = file.File(las_file_list[k], mode='r')

                pt_list_x = []
                pt_list_y = []
                pt_list_h = []

                print len(fr_las)
                for f in fr_las:
                    pt_list_x.append(f.x)
                    pt_list_y.append(f.y)
                    pt_list_h.append(f.z)

                print "Reading las file ", las_file_list[k]

                fr_txt.readline()
                site_id = fr_txt.readline().split(',')[0]
                fr_txt.seek(0)

                output_file_name = txt_ofolder_path + '\\' + site_id + '_cmp.txt'

                with open(output_file_name, 'wb') as fw:
                    fw.write('site_id,point_id,x,y,ter_h,gr_veg_h,las_x,las_y,las_h,pt_veg_h,delta_h\n')

                    #read txt, search, and write to txt file
                    fr_txt.readline()
                    for line in fr_txt:
                        row = line.split(',')
                        txt_pt = bc.rtk_veg(str(row[0]), str(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]))

                        min_pt_x, min_pt_y, min_pt_h, min_dis = bc.search_ht_ary(txt_pt.x, txt_pt.y, pt_list_x, pt_list_y, pt_list_h)
                        pt_veg_h = min_pt_h - txt_pt.ter_h
                        delta_h = pt_veg_h - txt_pt.gr_veg_h
                        # print txt_pt.site_id, txt_pt.point_id, txt_pt.x, txt_pt.y, txt_pt.ter_h, '\n'
                        # print min_pt_x, min_pt_y, min_pt_h, min_dis, pt_veg_h, delta_h
                        # print '-----------------------------------'
                        w_line = '%s,%s,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f\n' % \
                                 (str(txt_pt.site_id), str(txt_pt.point_id),  txt_pt.x,  txt_pt.y, txt_pt.ter_h,  txt_pt.gr_veg_h\
                                , min_pt_x, min_pt_y, min_pt_h, pt_veg_h, delta_h)

                        fw.write(w_line)
                    print 'Output file saved in ', output_file_name
                fw.close()
            fr_las.close()
        fr_txt.close()

#----------------------------------------------------------------------------------------------------------------------#
# This script writes a bat. file that classifies points as ground and non-ground using lasground.exe tool.
# Input: 1.Folder path of las. files. 2.Folder path of bat. file to save. 3. lastool path 4. output las. folder	5. version number Output:bat file
# Make sure no space in the path.

def classify_las_bat(lastool_path, las_ifolder_path, las_ofolder_path, version, bat_file_path):

    if not os.path.exists(las_ifolder_path):
        print "Could not find ", las_ifolder_path, "!"
        return
    las_file_list = glob.glob(os.path.join(las_ifolder_path, '*.las'))
    file_num = len(las_file_list)
    print "There are ", file_num, "files in", las_ifolder_path

    f = open(bat_file_path, 'wb')

    f.write('rem @echo off\n')

    for k in range(0, file_num):
        old_file_name = os.path.basename(las_file_list[k]).split('.')[0]
        old_file_name = old_file_name.split('v')[0]
        las_ofile_path = las_ofolder_path + '\\' + old_file_name + 'v' + str(version) + '.las'
        line = lastool_path + ' ' + '-i' + ' ' + las_file_list[k] + ' ' + '-o' + ' ' + las_ofile_path + '\n'
        f.write(line)

    f.write('PAUSE')

    f.close()
#----------------------------------------------------------------------------------------------------------------------#
# This script writes a batch file that separates ground and veg points
# lastool_path = r'E:\LAStools\LAStools\LAStools\bin\las2las.exe'
def spt_class_las_bat(lastool_path, las_ifolder_path, las_ofolder_path_veg, las_ofolder_path_grd, veg_v, grd_v , bat_file_path):

    if not os.path.exists(las_ifolder_path):
        print "Could not find ", las_ifolder_path, "!"
        return
    las_file_list = glob.glob(os.path.join(las_ifolder_path, '*.las'))
    file_num = len(las_file_list)
    print "There are ", file_num, "files in", las_ifolder_path

    f = open(bat_file_path, 'wb')

    f.write('rem @echo off\n')

    for k in range(0, file_num):
        old_file_name = os.path.basename(las_file_list[k]).split('.')[0]
        old_file_name = old_file_name.split('v')[0]
        las_ofile_path_veg = las_ofolder_path_veg + '\\' + old_file_name + 'v' + str(veg_v) + '.las'
        line = lastool_path + ' ' + '-i' + ' ' + las_file_list[k] + ' ' + '-o' + ' ' + las_ofile_path_veg + ' -keep_class 1' +'\n'
        f.write(line)

        las_ofile_path_grd = las_ofolder_path_grd + '\\' + old_file_name + 'v' + str(grd_v) + '.las'
        line = lastool_path + ' ' + '-i' + ' ' + las_file_list[k] + ' ' + '-o' + ' ' + las_ofile_path_grd + ' -keep_class 2' +'\n'
        f.write(line)

    f.write('PAUSE')

    f.close()

#----------------------------------------------------------------------------------------------------------------------#
# This script writes a bat. file that display las files
def view_las_bat(lastool_path, las_ifolder_path, bat_file_path):

    if not os.path.exists(las_ifolder_path):
        print "Could not find ", las_ifolder_path, "!"
        return
    las_file_list = glob.glob(os.path.join(las_ifolder_path, '*.las'))
    file_num = len(las_file_list)
    print "There are ", file_num, "files in", las_ifolder_path

    f = open(bat_file_path, 'wb')

    f.write('rem @echo off\n')

    for k in range(0, file_num):
        line = lastool_path + ' ' + '-i' + ' ' + las_file_list[k] + '\n'
        f.write(line)

    f.write('PAUSE')

    f.close()

#----------------------------------------------------------------------------------------------------------------------#
# this script removes points for compare where the vegetation heigth is higher than 5 meters, _cmp file

def rmv_hpt(txt_folder_path, output_folder_path):
    if not os.path.exists(txt_folder_path):
        print "Could not find ", txt_folder_path
        return

    txt_file_list = glob.glob(os.path.join(txt_folder_path, '*.txt'))

    for file in txt_file_list:
        with open(file, 'rb') as fr:
            output_file_path = output_folder_path + '\\' + os.path.basename(file).split('.')[0] + '_new.txt'
            print output_file_path
            with open(output_file_path, 'wb') as fw:
                fw.write(fr.readline())
                site_id = fr.readline().split(',')[0]
                fr.seek(0)
                fr.readline()
                for line in fr:
                    row = line.split(',')
                    if float(row[9]) < 5.0:
                        fw.write(line)
                fw.close()
            fr.close()

#----------------------------------------------------------------------------------------------------------------------#
# this script conducts simple statistical analysis of errors.

import os
import glob
import numpy

def error_sts(txt_folder_path, output_file_path):
    if not os.path.exists(txt_folder_path):
        print "Could not find ", txt_folder_path
        return

    txt_file_list = glob.glob(os.path.join(txt_folder_path, '*.txt'))
    with open(output_file_path, 'w') as fw:
        fw.write('site_id,all_rms,L_rms,C_rms, all_nrms(%),L_nrms(%),C_nrms(%),all_error_m,L_error_m,C_error_m,all_error_std,L_error_std,C_error_std\n')
        for file in txt_file_list:
            with open(file, 'rb') as fr:
                all_error = []
                L_error = []
                C_error = []
                fr.readline()
                all_min = L_min  = C_min = 999.0
                all_max = L_max= C_max = -999.0

                for line in fr:
                    row = line.split(',')
                    site_id = str(row[0])
                    all_error.append(float(row[10]))

                    if float(row[5]) > all_max:
                        all_max = float(row[5])
                    if float(row[5]) < all_min:
                        all_min = float(row[5])

                    if ('L' in row[1]) and ('.' not in row[1]):
                        L_error.append(float(row[10]))
                        if float(row[5]) > L_max:
                            L_max = float(row[5])
                        if float(row[5]) < L_min:
                            L_min = float(row[5])

                    if ('.' in row[1]):
                        C_error.append(float(row[10]))
                        if float(row[5]) > C_max:
                            C_max = float(row[5])
                        if float(row[5]) < C_min:
                            C_min = float(row[5])

                all_rms = numpy.sqrt((numpy.square(all_error)).mean())
                L_rms = numpy.sqrt((numpy.square(L_error)).mean())
                C_rms = numpy.sqrt((numpy.square(C_error)).mean())

                all_nrms = all_rms/(all_max - all_min)
                L_nrms = L_rms/(L_max - L_min)
                C_nrms = C_rms/(C_max - C_min)

                all_mean = numpy.mean(all_error)
                L_mean = numpy.mean(L_error)
                C_mean = numpy.mean(C_error)

                all_std = numpy.std(all_error)
                L_std = numpy.std(L_error)
                C_std = numpy.std(C_error)

                w_line = '%s,%.3f,%.3f,%.3f,%.2f,%.2f,%.2f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n' % \
                         (site_id, all_rms, L_rms, C_rms, all_nrms*100, L_nrms*100, C_nrms*100, all_mean, L_mean, C_mean, all_std, L_std, C_std)
                fw.write(w_line)
            fr.close()
        fw.close()

#----------------------------------------------------------------------------------------------------------------------#
# this script draws profile comparison plots

import os
import glob
import matplotlib.pyplot as plt
import numpy

def draw_cmp_plot(txt_folder_path, output_folder_path):
    if not os.path.exists(txt_folder_path):
        print "Could not find ", txt_folder_path
        return

    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    txt_file_list = glob.glob(os.path.join(txt_folder_path, '*.txt'))

    f1 = output_folder_path + '\\all'
    f2 = output_folder_path + '\\Ltran'
    f3 = output_folder_path + '\\Ctran'
    f4 = output_folder_path + '\\6Ctran'
    f5 = output_folder_path + '\\7Ctran'
    f6 = output_folder_path + '\\9Ctran'
    if not os.path.exists(f1):
        os.makedirs(f1)
    if not os.path.exists(f2):
        os.makedirs(f2)
    if not os.path.exists(f3):
        os.makedirs(f3)
    if not os.path.exists(f4):
        os.makedirs(f4)
    if not os.path.exists(f5):
        os.makedirs(f5)
    if not os.path.exists(f6):
        os.makedirs(f6)

    for file in txt_file_list:
        all_output_file_path = f1 + '\\' + os.path.basename(file).split('.')[0] + '.png'
        L_output_file_path = f2 + '\\' + os.path.basename(file).split('.')[0] + '.png'
        C_output_file_path = f3 + '\\' + os.path.basename(file).split('.')[0] + '.png'
        C6_output_file_path = f4 + '\\' + os.path.basename(file).split('.')[0] + '.png'
        C7_output_file_path = f5 + '\\' + os.path.basename(file).split('.')[0] + '.png'
        C9_output_file_path = f6 + '\\' + os.path.basename(file).split('.')[0] + '.png'

        fw1 = open(all_output_file_path, 'wb')
        fw2 = open(L_output_file_path, 'wb')
        fw3 = open(C_output_file_path, 'wb')
        fw4 = open(all_output_file_path, 'wb')
        fw5 = open(L_output_file_path, 'wb')
        fw6 = open(C_output_file_path, 'wb')

        with open(file, 'rb') as fr:
            fr.readline()
            site_id = fr.readline().split(',')[0]
            print "Processing site", site_id
            fr.seek(0)
            fr.readline()
            all_xlist, all_y1list, all_y2list = [], [], []
            L_xlist, L_y1list, L_y2list = [], [], []
            C_xlist, C_y1list, C_y2list = [], [], []
            C6_xlist, C6_y1list, C6_y2list = [], [], []
            C7_xlist, C7_y1list, C7_y2list = [], [], []
            C9_xlist, C9_y1list, C9_y2list = [], [], []
            index = C_index = 0
            L_dis = C6_dis = C7_dis = C9_dis = 0.0
            for line in fr:
                all_xlist.append(index)
                row = line.split(',')
                indct = row[1].split('_')[1]
                y1 = float(row[5])
                y2 = float(row[9])
                all_y1list.append(y1)
                all_y2list.append(y2)
                index += 1

                if ('L' in indct) and ('.' not in indct):
                    # create long transect list
                    L_y1list.append(y1)
                    L_y2list.append(y2)
                    L_dis = float(indct.split('L')[1])
                    L_xlist.append(L_dis)

                if ('.' in indct):
                    # create cross transect list
                    C_y1list.append(y1)
                    C_y2list.append(y2)
                    C_index += 1

                if ('.' in indct) and (indct[0] == '6'):
                    # create 60-m cross transect list
                    C6_dis = float(indct[2:5])
                    C6_y1list.append(y1)
                    C6_y2list.append(y2)
                    C6_xlist.append(C6_dis)

                if ('.' in indct) and (indct[0] == '7'):
                    # create 75-m cross transect list
                    C7_dis = float(indct[2:5])
                    C7_y1list.append(y1)
                    C7_y2list.append(y2)
                    C7_xlist.append(C7_dis)

                if ('.' in indct) and (indct[0] == '9'):
                    C9_dis = float(indct[2:5])
                    C9_y1list.append(y1)
                    C9_y2list.append(y2)
                    C9_xlist.append(C9_dis)

            all_xlist = range(0, index)
            C_xlist = range(0, C_index)

            L_xlistnew, L_y1listnew, L_y2listnew = resort(L_xlist, L_y1list, L_y2list)
            C6_xlistnew, C6_y1listnew, C6_y2listnew = resort(C6_xlist, C6_y1list, C6_y2list)
            C7_xlistnew, C7_y1listnew, C7_y2listnew = resort(C7_xlist, C7_y1list, C7_y2list)
            C9_xlistnew, C9_y1listnew, C9_y2listnew = resort(C9_xlist, C9_y1list, C9_y2list)

            draw_plt(site_id, all_xlist, all_y1list, all_y2list, all_output_file_path, 0, 150, 'Point No.')
            draw_plt(site_id, C_xlist, C_y1list, C_y2list, C_output_file_path, 0, 10, 'Point No.')

            draw_plt(site_id, L_xlistnew, L_y1listnew, L_y2listnew, L_output_file_path, 0, 150)
            draw_plt(site_id, C6_xlistnew, C6_y1listnew, C6_y2listnew, C6_output_file_path, 0, 6)
            draw_plt(site_id, C7_xlistnew, C7_y1listnew, C7_y2listnew, C7_output_file_path, 0, 6)
            draw_plt(site_id, C9_xlistnew, C9_y1listnew, C9_y2listnew, C9_output_file_path, 0, 6)

            fr.close()
        fw1.close()
        fw2.close()
        fw3.close()
        fw4.close()
        fw5.close()
        fw6.close()

def draw_plt(site_id, xlist, y1list, y2list, output_file_path, xmin = 0, xmax = 30, x_label = 'distance(meters)'):
    fig = plt.figure(figsize=(xmax+10, 20))
    fig.canvas.set_window_title(site_id)

    plt.ylim(0, 6)
    plt.xlim(xmin, xmax)
    plot_gr_veg, = plt.plot(xlist, y1list, '-bo', linewidth=4.0, markersize=20)
    plot_pt_veg, = plt.plot(xlist, y2list, '-ro', linewidth=4.0,  markersize=20)

    plt.legend([plot_gr_veg, plot_pt_veg], ["Field Veg heights", "Point-Cloud Veg heights"], fontsize=40)

    plt.xlabel(x_label, fontsize=40)
    y_label = 'veg_height(meters)'
    plt.ylabel(y_label, fontsize=40)
    plt.xticks(xlist, fontsize=40)
    plt.yticks(fontsize=80)

    plt.savefig(output_file_path)
#    plt.show()
    plt.close()

def resort(xlist, y1list, y2list):
    newxlist, newy1list, newy2list = [], [], []
    newxindex = numpy.argsort(xlist)
    for i in range(0, len(xlist)):
        newxlist.append(xlist[newxindex[i]])
        newy1list.append(y1list[newxindex[i]])
        newy2list.append(y2list[newxindex[i]])
    return newxlist, newy1list, newy2list

#----------------------------------------------------------------------------------------------------------------------#
# this script conduct a batch conversion from txt to csv

import csv
import os
import glob

def txt2csv(txt_folder_path, csv_folder_path):
    if not os.path.exists(txt_folder_path):
        print "Could not find ", txt_folder_path
        return

    if not os.path.exists(csv_folder_path):
        os.makedirs(csv_folder_path)

    txt_file_list = glob.glob(os.path.join(txt_folder_path, '*.txt'))

    for txt_file_path in txt_file_list:
        with open(txt_file_path, 'rb') as fr:
            csv_file_path = csv_folder_path + '\\' + os.path.basename(txt_file_path).split('.')[0] + '.csv'
            with open(csv_file_path, 'wb') as fw:
                txt_read = csv.reader(fr, delimiter = ',')
                csv_write = csv.writer(fw)
                csv_write.writerows(txt_read)
                wr = csv.writer(fw)
            fw.close()
        fr.close()

#---------------------------------------------------------------------------------------------------------------------#
# The function of this script is to convert .xlsx files to csv. files

import xlrd
import csv
import glob
import os

def xlsx_to_csv(xls_folder_path, csv_folder_path):

    if not os.path.exists(csv_folder_path):
        os.makedirs(csv_folder_path)

    xlsx_file_list = glob.glob(os.path.join(xls_folder_path, "*.xlsx"))

    for xlsx_file_path in xlsx_file_list:

        wb = xlrd.open_workbook(xlsx_file_path)
        sh = wb.sheet_by_name('Sheet1')

        csv_file_path = csv_folder_path + '\\' + os.path.basename(xlsx_file_path).split('.')[0] + '.csv'
        print csv_file_path
        your_csv_file = open(csv_file_path, 'wb')
        wr = csv.writer(your_csv_file, quoting=csv.QUOTE_ALL)

        for rownum in xrange(sh.nrows):
            wr.writerow(sh.row_values(rownum))

        your_csv_file.close()

#---------------------------------------------------------------------------------------------------------------------#
# this script append all points, L transect points, and C-transect point seperately

import csv
import os
import glob

def txt_append(txt_folder_path, output_file_path):
    if not os.path.exists(txt_folder_path):
        print "Could not find ", txt_folder_path
        return
    txt_file_list = glob.glob(os.path.join(txt_folder_path, '*.txt'))

    with open(output_file_path, 'wb') as fw:
        fr0 = open(txt_file_list[0])
        header = fr0.readline()
        fw.writelines(header)
        for txt_file_path in txt_file_list:
            with open(txt_file_path, 'rb') as fr:
                fr.readline()
                for line in fr:
                    fw.writelines(line)
        fw.close()
    fr.close()

def txt_Lappend(txt_folder_path, output_file_path):
    if not os.path.exists(txt_folder_path):
        print "Could not find ", txt_folder_path
        return
    txt_file_list = glob.glob(os.path.join(txt_folder_path, '*.txt'))

    with open(output_file_path, 'wb') as fw:
        fr0 = open(txt_file_list[0])
        header = fr0.readline()
        fw.writelines(header)
        for txt_file_path in txt_file_list:
            with open(txt_file_path, 'rb') as fr:
                fr.readline()
                for line in fr:
                    row = line.split(',')
                    indct = row[1].split('_')[1]
                    if ('L' in indct) and ('.' not in indct):
                        fw.writelines(line)
        fw.close()
    fr.close()


def txt_Cappend(txt_folder_path, output_file_path):
    if not os.path.exists(txt_folder_path):
        print "Could not find ", txt_folder_path
        return
    txt_file_list = glob.glob(os.path.join(txt_folder_path, '*.txt'))

    with open(output_file_path, 'wb') as fw:
        fr0 = open(txt_file_list[0])
        header = fr0.readline()
        fw.writelines(header)
        for txt_file_path in txt_file_list:
            with open(txt_file_path, 'rb') as fr:
                fr.readline()
                for line in fr:
                    row = line.split(',')
                    indct = row[1].split('_')[1]
                    if ('.' in indct):
                        fw.writelines(line)
        fw.close()
    fr.close()

#---------------------------------------------------------------------------------------------------------------------#
# This script splits GCPs and other RTK points into separate files
import glob
import os

def extract_gcp(data_folder_path, output_folder_path):
    file_path_list = glob.glob(os.path.join(data_folder_path, "*.txt")) #***

    site_number = len(file_path_list)

    print 'Number of sites: ', site_number

    for k in range(0, site_number):
        with open(file_path_list[k]) as fr:
            old_file_name = os.path.basename(file_path_list[k]).split('.')[0]  #***
            GCP_folder_path = output_folder_path + "\\output_GCP"
            if not os.path.exists(GCP_folder_path):
                os.makedirs(GCP_folder_path)
            GCP_file_path = GCP_folder_path + "\\" + old_file_name + "_GCP" + ".txt"

            RTK_folder_path = output_folder_path + "\\output_RTK"
            RTK_file_path = RTK_folder_path + "\\" + old_file_name + "_RTK" + ".txt"
            if not os.path.exists(RTK_folder_path):
                os.makedirs(RTK_folder_path)

            with open(GCP_file_path, 'w') as fw1:
                with open(RTK_file_path,'w') as fw2:
                    GCP_count = 0
                    RTK_count = 0
                    for line in fr:
                        if 'GCP' in line:
                            fw1.write(line)
                            GCP_count += 1
                        else:
                            fw2.write(line)
                            RTK_count += 1
                print 'GCP coordinates of file ', old_file_name, ' were saved in ', GCP_file_path
                print 'Coordinates of other RTK points of file ', old_file_name, ' were saved in ', RTK_file_path
                print "Numbers of RTK points: ", RTK_count
                print "Numbers of GCP points: ", GCP_count
                fw2.close()
            fw1.close()
        fr.close()
#---------------------------------------------------------------------------------------------------------------------#
# The function of this script is to apply a 3D shift for RTK points after post-processing of base station.
# 3D shift v1.10  Shijuan Chen 2016.09.09
import csv
import os

def ThreeD_shift(base_data_path, rtk_folder_path, output_folder_path):

    class site_base_crd(object):
        def __init__(self, no, site_id, ex_bf, ny_bf, oh_bf, ex_af, ny_af, oh_af, dx, dy, dh):
            self.no = no
            self.site_id = site_id
            self.ex_bf = ex_bf
            self.ny_bf = ny_bf
            self.oh_bf = oh_bf
            self.ex_af = ex_af
            self.ny_af = ny_af
            self.oh_af = oh_af
            self.dx = dx
            self.dy = dy
            self.dh = dh

    site_base_list = []

    with open(base_data_path, 'rb') as fb:
        base_reader = csv.reader(fb)
        row1 = next(base_reader)

        b_list  = list(base_reader)

        site_num = base_reader.line_num - 1
        print "number of sites: ", site_num

        for i in range(0, site_num):
            site_base_list.append(site_base_crd(0,0,0,0,0,0,0,0,0,0,0))

        for i in range(0, site_num):
            site_base_list[i].no = b_list[i][0]
            site_base_list[i].site_id = b_list[i][1]
            site_base_list[i].ex_bf = float(b_list[i][2])
            site_base_list[i].ny_bf = float(b_list[i][3])
            site_base_list[i].oh_bf = float(b_list[i][4])
            site_base_list[i].ex_af = float(b_list[i][5])
            site_base_list[i].ny_af = float(b_list[i][6])
            site_base_list[i].oh_af = float(b_list[i][7])
            site_base_list[i].dx = site_base_list[i].ex_af - site_base_list[i].ex_bf
            site_base_list[i].dy = site_base_list[i].ny_af - site_base_list[i].ny_bf
            site_base_list[i].dh = site_base_list[i].oh_af - site_base_list[i].oh_bf

    fb.close()

    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)
    base_newfile_path = output_folder_path + "\\Base_info.csv"

    with open(base_newfile_path,'wb') as fb_w:
        base_writer = csv.writer(fb_w)
        base_writer.writerow(['No', 'Site_ID', 'EX_BF', 'NY_BF', 'OH_BF', 'EX_AF', 'NY_AF', 'OH_AF','DX','DY','DH'])

        for i in range(0, site_num):
            base_writer.writerow([site_base_list[i].no, site_base_list[i].site_id,
                                 site_base_list[i].ex_bf, site_base_list[i].ny_bf, site_base_list[i].oh_bf,
                                 site_base_list[i].ex_af, site_base_list[i].ny_af, site_base_list[i].oh_af,
                                 site_base_list[i].dx, site_base_list[i].dy, site_base_list[i].dh])

    print 'Base info was saved in', base_newfile_path

    fb.close()

    class rtk_crd(object):
        def __init__(self, site_id, point_id, x, y, h):
            self.site_id = site_id
            self.point_id = point_id
            self.x = x
            self.y = y
            self.h = h


    rtk_newfolder_path = output_folder_path + '\\'

    for k in range(0, site_num):
        rtk_file_path = rtk_folder_path + '\\' + site_base_list[k].site_id + '.txt'
        with open(rtk_file_path) as fr:
            rtk_point_list = []
            rtk_point_newlist = []
            rtk_point_num = 0
            for line in fr:
                rtk_row = line.split(',')
                rtk_point_list.append(rtk_crd(site_base_list[k].site_id, rtk_row[0],rtk_row[2], rtk_row[1], rtk_row[3]))
                rtk_point_num += 1
            print "Number of RTK points in site ", site_base_list[k].site_id, ' is ', rtk_point_num, '.'

            for i in range (0, rtk_point_num):
                rtk_af_x = float(rtk_point_list[i].x) + site_base_list[k].dx
                rtk_af_y = float(rtk_point_list[i].y) + site_base_list[k].dy
                rtk_af_h = float(rtk_point_list[i].h) + site_base_list[k].dh
                rtk_point_newlist.append(rtk_crd(rtk_point_list[i].site_id, rtk_point_list[i].point_id, rtk_af_x, rtk_af_y, rtk_af_h))

            rtk_newfile_path = rtk_newfolder_path + site_base_list[k].site_id + '_new.txt'

            with open(rtk_newfile_path,'w') as fw:
                for i in range(0, rtk_point_num):
                    fw_line = rtk_point_list[i].point_id + ',' + '%.4f' % rtk_point_newlist[i].y + ',' + '%.4f' % rtk_point_newlist[i].x + ',' + '%.4f' % rtk_point_newlist[i].h + '\n'
                    fw.write(fw_line)
            fw.close()
            print 'New coordinates were saved in ', rtk_newfile_path
        fr.close()

#----------------------------------------------------------------------------------------------------------------------#
# This script includes new functions in Point clouds processing v1.1
#--------------------------------------------------------------------------------------------------------------------#
# this script write a bat. file that creates las files normalized by ground height. lasheight.exe -replace_z
# Make sure the input files are classified. For more info, please see lasheight_README

import os
import glob

def CHM_las_bat(lastool_path, las_ifolder_path, las_ofolder_path, version, bat_file_path):

    if not os.path.exists(las_ifolder_path):
        print "Could not find ", las_ifolder_path, "!"
        return
    las_file_list = glob.glob(os.path.join(las_ifolder_path, '*.las'))
    file_num = len(las_file_list)
    print "There are ", file_num, "files in", las_ifolder_path

    f = open(bat_file_path, 'wb')

    f.write('rem @echo off\n')

    for k in range(0, file_num):
        old_file_name = os.path.basename(las_file_list[k]).split('.')[0]
        old_file_name = old_file_name.split('v')[0]
        las_ofile_path = las_ofolder_path + '\\' + old_file_name + 'v' + str(version) + '.las'
        line = lastool_path + ' ' + '-i' + ' ' + las_file_list[k] + ' ' + '-replace_z -o' + ' ' + las_ofile_path + '\n'
        f.write(line)

    f.write('PAUSE')

    f.close()
#--------------------------------------------------------------------------------------------------------------------#
# this script creates a circle buffer of RTK points
from osgeo import ogr, osr
import os
import glob

class rtk_crd(object):
    def __init__(self, site_id, point_id, x, y, ter_h):
        self.site_id = site_id
        self.point_id = point_id
        self.x = x
        self.y = y
        self.ter_h = ter_h

def crcl_shp(rtk_folder_path, shp_folder_path, bufferDistance):
    if not os.path.exists(rtk_folder_path):
        print "Could not find ", rtk_folder_path, "!"
        return

    if not os.path.exists(shp_folder_path):
        os.makedirs(shp_folder_path)

    rtk_file_list = glob.glob(os.path.join(rtk_folder_path, '*.txt'))
    file_num = len(rtk_file_list)
    print "There are ", file_num, "files in ", rtk_folder_path

    for k in range(0, file_num):
        fr = open(rtk_file_list[k], 'r')
        start_pt = rtk_crd('0', '0', 0, 0, 0)
        end_pt = rtk_crd('0', '0', 0, 0, 0)
        row = fr.readline().split(',')
        site_id = row[0][0] + row[0][1] + row[0][2]
        print site_id
        # set up shapefile driver
        driver = ogr.GetDriverByName('ESRI Shapefile')
        layer_name = 'clip_cbf_' + str(site_id)
        file_name = shp_folder_path + '\\' + layer_name + '.shp'
        print file_name

        if os.path.exists(file_name):
            driver.DeleteDataSource(file_name)

        data_source = driver.CreateDataSource(shp_folder_path)

        # create the spatial reference: NAD83 UTM Z12, EPSG Projection 26912
        srf = osr.SpatialReference()
        srf.ImportFromEPSG(26912)

        # create the layer
        layer = data_source.CreateLayer(layer_name, srf, ogr.wkbPolygon)

        #create a geometry collection
        #geocl = ogr.Geometry(ogr.wkbGeometryCollection)

        feature_list = []
        fr.seek(0) # newly added
        for line in fr:
            row = line.split(',')
            point_id = row[0]
            x = row[2]
            y = row[1]
            ter_h = row[3]
            rtk_point  = rtk_crd(str(site_id), str(point_id), float(x), float(y), float(ter_h))

            wkt = 'POINT (' + str(rtk_point.x) + ' ' + str(rtk_point.y) +')'
            pt = ogr.CreateGeometryFromWkt(wkt)
            poly = pt.Buffer(bufferDistance)

            # create the feature
            feature = ogr.Feature(layer.GetLayerDefn())
            feature.SetGeometry(poly)


            #Set the feature geometry using the polygon
            #layer.CreateFeature(feature)
            feature_list.append(feature)

        # merge polygons
        union_poly = ogr.Geometry(ogr.wkbPolygon)

        for feature in feature_list:
            geom = feature.GetGeometryRef()
            union_poly = union_poly.Union(geom)

        feature2 = ogr.Feature(layer.GetLayerDefn())
        feature2.SetGeometry(union_poly)
        layer.CreateFeature(feature2)

        #destroy the data source to free resources
        data_source.Destroy()

#--------------------------------------------------------------------------------------------------------------------#
# this scripts calculates height metrics of point clouds in buffers of RTK point with certain distance.
import os
import glob
from liblas import file

import BERA_classes as bc

def cal_metrics(las_ifolder_path, txt_ifolder_path, txt_ofolder_path, bufferDistance = 0.1):
    if not os.path.exists(las_ifolder_path):
        print "Could not find ", las_ifolder_path, "!"
        if not os.path.exists(txt_ifolder_path):
            print "Could not find ", txt_ifolder_path, "!"
            return
        return
    if not os.path.exists(txt_ofolder_path):
        os.makedirs(txt_ofolder_path)

    las_file_list = glob.glob(os.path.join(las_ifolder_path, '*.las'))
    txt_file_list = glob.glob(os.path.join(txt_ifolder_path, '*.txt'))

    las_file_num = len(las_file_list)
    txt_file_num = len(txt_file_list)
    if not (las_file_num == txt_file_num):
        print "Number of las and txt files does not match!"
        return

    for k in range(0, las_file_num):
        print las_file_num
        with open(txt_file_list[k]) as fr_txt:
            with open(las_file_list[k]) as fr_las:
               # read las
                pt_3d_las = []
                fr_las = file.File(las_file_list[k], mode='r')
                pt_list_x = []
                pt_list_y = []
                pt_list_h = []
                print len(fr_las)
                for f in fr_las:
                    pt_list_x.append(f.x)
                    pt_list_y.append(f.y)
                    pt_list_h.append(f.z)
                print "Reading las file ", las_file_list[k]

                # read txt files
                fr_txt.readline()
                site_id = fr_txt.readline().split(',')[0]
                fr_txt.seek(0)
                output_file_name = txt_ofolder_path + '\\' + site_id + '_metrics.txt'

                with open(output_file_name, 'w') as fw:
                    fw.write('site_id,point_id,x,y,ter_h,gr_veg_h,h_max,h95,h90,h85,h80,h75,h70,h65,h60,h55,h_ median,h_mean,nn\n')

                    #read txt, search, and write to txt file
                    fr_txt.readline()
                    for line in fr_txt:
                        row = line.split(',')
                        txt_pt = bc.rtk_veg(str(row[0]), str(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]))

                        # get height list within the buffer
                        hlist = bc.hlist_inBuffer(txt_pt.x, txt_pt.y, pt_list_x, pt_list_y, pt_list_h, bufferDistance)

                        # calculate height metrics
                        metrics_array = bc.height_metrics(hlist)

                        # get nearest-neighbour height
                        min_pt_x, min_pt_y, min_pt_h, min_dis = bc.search_ht_ary(txt_pt.x, txt_pt.y, pt_list_x,
                                                                                 pt_list_y, pt_list_h)
                        nn = min_pt_h
                        if metrics_array == -9999:
                            print line, "Error!"
                            metrics_array = []

                        metrics_array.append(nn)

                        #calculate veg height
                        for j in range(0, len(metrics_array)):
                            metrics_array[j] = metrics_array[j] - txt_pt.ter_h
                            if metrics_array[j] <0:
                                metrics_array[j] = 0

                        str_list = []
                        for m in metrics_array:
                            str_m = '%.3f' % m
                            str_list.append(str_m)
                        m_line = ','.join(str(item) for item in str_list)

                        w1_line = '%s,%s,%0.3f,%0.3f,%0.3f,%0.3f,' % \
                                 (str(txt_pt.site_id), str(txt_pt.point_id), txt_pt.x, txt_pt.y, txt_pt.ter_h,  txt_pt.gr_veg_h)
                        w_line = w1_line + m_line + '\n'
                        if metrics_array:
                            fw.write(w_line)
                    print 'Output file saved in ', output_file_name
                fw.close()
            fr_las.close()
        fr_txt.close()
#--------------------------------------------------------------------------------------------------------------------#
# this script compares veg height of UAV-based point clouds when using UAv point clouds to create terrain
import os
import glob
from liblas import file

import BERA_classes as bc

def CmpVeght_UAV_UAV(las_ifolder_path, txt_ifolder_path, txt_ofolder_path):
    if not os.path.exists(las_ifolder_path):
        print "Could not find ", las_ifolder_path, "!"
        if not os.path.exists(txt_ifolder_path):
            print "Could not find ", txt_ifolder_path, "!"
            return
        return
    if not os.path.exists(txt_ofolder_path):
        os.makedirs(txt_ofolder_path)

    las_file_list = glob.glob(os.path.join(las_ifolder_path, '*.las'))
    txt_file_list = glob.glob(os.path.join(txt_ifolder_path, '*.txt'))

    las_file_num = len(las_file_list)
    txt_file_num = len(txt_file_list)
    if not (las_file_num == txt_file_num):
        print "Number of las and txt files does not match!"
        return

    for k in range(0, las_file_num):
        print las_file_num
        with open(txt_file_list[k]) as fr_txt:
            with open(las_file_list[k]) as fr_las:
               # read las
                pt_3d_las = []
                fr_las = file.File(las_file_list[k], mode='r')

                pt_list_x = []
                pt_list_y = []
                pt_list_h = []

                print len(fr_las)
                for f in fr_las:
                    pt_list_x.append(f.x)
                    pt_list_y.append(f.y)
                    pt_list_h.append(f.z)

                print "Reading las file ", las_file_list[k]

                fr_txt.readline()
                site_id = fr_txt.readline().split(',')[0]
                fr_txt.seek(0)

                output_file_name = txt_ofolder_path + '\\' + site_id + '_cmp_UU.txt'

                with open(output_file_name, 'wb') as fw:
                    fw.write('site_id,point_id,x,y,ter_h,gr_veg_h,las_x,las_y,las_h,pt_veg_h,delta_h\n')

                    #read txt, search, and write to txt file
                    fr_txt.readline()
                    for line in fr_txt:
                        row = line.split(',')
                        txt_pt = bc.rtk_veg(str(row[0]), str(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]))

                        min_pt_x, min_pt_y, min_pt_h, min_dis = bc.search_ht_ary(txt_pt.x, txt_pt.y, pt_list_x, pt_list_y, pt_list_h)
                        pt_veg_h = min_pt_h
                        delta_h = pt_veg_h - txt_pt.gr_veg_h
                        # print txt_pt.site_id, txt_pt.point_id, txt_pt.x, txt_pt.y, txt_pt.ter_h, '\n'
                        # print min_pt_x, min_pt_y, min_pt_h, min_dis, pt_veg_h, delta_h
                        # print '-----------------------------------'
                        w_line = '%s,%s,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f\n' % \
                                 (str(txt_pt.site_id), str(txt_pt.point_id),  txt_pt.x,  txt_pt.y, txt_pt.ter_h,  txt_pt.gr_veg_h\
                                , min_pt_x, min_pt_y, min_pt_h, pt_veg_h, delta_h)

                        fw.write(w_line)
                    print 'Output file saved in ', output_file_name
                fw.close()
            fr_las.close()
        fr_txt.close()

#--------------------------------------------------------------------------------------------------------------------
# this script converts _cmp files to _smy2 files.
import os
import glob

def cmp2smy2(UAV_RTK_folder, UAV_UAV_folder, output_folder_path):
    if not os.path.exists(UAV_RTK_folder):
        print "Could not find ", UAV_RTK_folder
        if not os.path.exists(UAV_UAV_folder):
            print "Could not find ", UAV_UAV_folder
            return
        return

    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    UR_file_list = glob.glob(os.path.join(UAV_RTK_folder, '*.txt'))
    UU_file_list = glob.glob(os.path.join(UAV_UAV_folder, '*.txt'))
    if len(UR_file_list) <> len(UU_file_list):
        print "File number does not match!"
        return

    for k in range(0,len(UR_file_list)):
        with open(UR_file_list[k]) as fur:
            with open(UU_file_list[k]) as fuu:

                fur.readline()
                site_id = fur.readline().split(',')[0]
                print "Processing site", site_id
                fur.seek(0)

                fur.readline()
                fuu.readline()

                output_file_path = output_folder_path + '//' + site_id + '_smy2.txt'
                with open(output_file_path, 'w') as fw:
                    fw.write('site_id,point_id,x,y,gr_veg_h,UAV_RTK_vh,UAV_UAV_vh\n')

                    for line in fur:
                        row = line.split(',')
                        site_id = str(row[0])
                        point_id = str(row[1])
                        x = float(row[2])
                        y = float(row[3])
                        gr_veg_h = float(row[5])
                        ur_vh = float(row[9])

                        urow = fuu.readline().split(',')

                        uu_vh = float(urow[9])

                        if ur_vh < 5.0 and uu_vh < 5.0:
                            wline = site_id + ',' + point_id + ',' + str(x) + ',' + str(y) + ',' + str(gr_veg_h) + ',' + str(ur_vh) + ',' + str(uu_vh) + '\n'
                            fw.write(wline)
                    fw.close()
                fuu.close()
            fur.close()
#--------------------------------------------------------------------------------------------------------------------
# this script draws comparison plot of veg height measured from field, by UAV-UAV method, and by UAV-RTK method

import os
import glob
import matplotlib.pyplot as plt
import numpy

def plot_3p(txt_folder_path, output_folder_path):

    if not os.path.exists(txt_folder_path):
        print "Could not find ", txt_folder_path
        return

    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    txt_file_list = glob.glob(os.path.join(txt_folder_path, '*.txt'))

    fL = output_folder_path + '\\Ltran'
    f6 = output_folder_path + '\\6Ctran'
    f7 = output_folder_path + '\\7Ctran'
    f9 = output_folder_path + '\\9Ctran'

    if not os.path.exists(fL):
        os.makedirs(fL)
    if not os.path.exists(f6):
        os.makedirs(f6)
    if not os.path.exists(f7):
        os.makedirs(f7)
    if not os.path.exists(f9):
        os.makedirs(f9)

    for k in range(0, len(txt_file_list)):
        with open(txt_file_list[k], 'rb') as fr:
            fr.readline()
            site_id = fr.readline().split(',')[0]
            print "Processing site", site_id
            fr.seek(0)
            fr.readline()

            L_output_file_path = fL + '\\' + site_id + '.png'
            C6_output_file_path = f6 + '\\' + site_id + '.png'
            C7_output_file_path = f7 + '\\' + site_id + '.png'
            C9_output_file_path = f9 + '\\' + site_id + '.png'

            L_xlist, L_y1list, L_y2list, L_y3list = [], [], [], []
            C6_xlist, C6_y1list, C6_y2list, C6_y3list = [], [], [], []
            C7_xlist, C7_y1list, C7_y2list, C7_y3list = [], [], [], []
            C9_xlist, C9_y1list, C9_y2list, C9_y3list = [], [], [], []
            L_dis = C6_dis = C7_dis = C9_dis = 0.0

            for line in fr:
                row = line.split(',')
                point_id = row[1]
                indct = row[1].split('_')[1]
                y1 = float(row[4])
                y2 = float(row[5])
                y3 = float(row[6])

                if ('L' in indct) and ('.' not in indct):
                    # create long transect list
                    L_y1list.append(y1)
                    L_y2list.append(y2)
                    L_y3list.append(y3)
                    L_dis = float(indct.split('L')[1])
                    L_xlist.append(L_dis)

                if ('.' in indct) and (indct[0] == '6'):
                    # create 60-m cross transect list
                    C6_dis = float(indct[2:5])
                    C6_y1list.append(y1)
                    C6_y2list.append(y2)
                    C6_y3list.append(y3)
                    C6_xlist.append(C6_dis)

                if ('.' in indct) and (indct[0] == '7'):
                    # create 75-m cross transect list
                    C7_dis = float(indct[2:5])
                    C7_y1list.append(y1)
                    C7_y2list.append(y2)
                    C7_y3list.append(y3)
                    C7_xlist.append(C7_dis)

                if ('.' in indct) and (indct[0] == '9'):
                    C9_dis = float(indct[2:5])
                    C9_y1list.append(y1)
                    C9_y2list.append(y2)
                    C9_y3list.append(y3)
                    C9_xlist.append(C9_dis)

            L_xlistnew, L_y1listnew, L_y2listnew, L_y3listnew = resort3(L_xlist, L_y1list, L_y2list, L_y3list)
            C6_xlistnew, C6_y1listnew, C6_y2listnew, C6_y3listnew = resort3(C6_xlist, C6_y1list, C6_y2list, C6_y3list)
            C7_xlistnew, C7_y1listnew, C7_y2listnew, C7_y3listnew = resort3(C7_xlist, C7_y1list, C7_y2list, C7_y3list )
            C9_xlistnew, C9_y1listnew, C9_y2listnew, C9_y3listnew = resort3(C9_xlist, C9_y1list, C9_y2list, C9_y3list)

            draw_Lplt(site_id, L_xlistnew, L_y1listnew, L_y2listnew, L_y3listnew, L_output_file_path)
            draw_Cplt(site_id, C6_xlistnew, C6_y1listnew, C6_y2listnew, C6_y3listnew, C6_output_file_path)
            draw_Cplt(site_id, C7_xlistnew, C7_y1listnew, C7_y2listnew, C7_y3listnew, C7_output_file_path)
            draw_Cplt(site_id, C9_xlistnew, C9_y1listnew, C9_y2listnew, C9_y3listnew, C9_output_file_path)


def draw_Lplt(site_id, xlist, y1list, y2list, y3list, output_file_path):
    fig = plt.figure(figsize=(150, 20))
    fig.canvas.set_window_title(site_id)

    plt.ylim(0, 6)
    plt.xlim(0, 150)
    plot_gr_veg, = plt.plot(xlist, y1list, '-bo', linewidth=4.0, markersize=20)
    plot_pt2_veg, = plt.plot(xlist, y2list, '-ro', linewidth=4.0,  markersize=20)
    plot_pt3_veg, = plt.plot(xlist, y3list, '-go', linewidth=4.0,  markersize=20)

    plt.legend([plot_gr_veg, plot_pt2_veg, plot_pt3_veg], ["Field Veg heights", "UAV-RTK Veg heights", "UAV-UAV Veg heights"], fontsize=40)

    plt.xlabel('distance(meters)', fontsize=40)
    plt.ylabel('veg_height(meters)', fontsize=40)
    plt.xticks(xlist, fontsize=40)
    plt.yticks(fontsize=60)

    plt.savefig(output_file_path)
#    plt.show()
    plt.close()

def draw_Cplt(site_id, xlist, y1list, y2list, y3list, output_file_path):
    fig = plt.figure(figsize=(30, 20))
    fig.canvas.set_window_title(site_id)

    plt.ylim(0, 6)
    plt.xlim(0, 6)
    plot_gr_veg, = plt.plot(xlist, y1list, '-bo', linewidth=4.0, markersize=20)
    plot_pt2_veg, = plt.plot(xlist, y2list, '-ro', linewidth=4.0,  markersize=20)
    plot_pt3_veg, = plt.plot(xlist, y3list, '-go', linewidth=4.0,  markersize=20)

    plt.legend([plot_gr_veg, plot_pt2_veg, plot_pt3_veg], ["Field Veg heights", "UAV-RTK Veg heights", "UAV-UAV Veg heights"], fontsize=30)

    plt.xlabel('distance(meters)', fontsize=30)
    plt.ylabel('veg_height(meters)', fontsize=30)
    plt.xticks(xlist, fontsize=60)
    plt.yticks(fontsize=60)

    plt.savefig(output_file_path)
#    plt.show()
    plt.close()

def resort3(xlist, y1list, y2list, y3list):
    newxlist, newy1list, newy2list, newy3list = [], [], [], []
    newxindex = numpy.argsort(xlist)
    for i in range(0, len(xlist)):
        newxlist.append(xlist[newxindex[i]])
        newy1list.append(y1list[newxindex[i]])
        newy2list.append(y2list[newxindex[i]])
        newy3list.append(y3list[newxindex[i]])
    return newxlist, newy1list, newy2list, newy3list

#----------------------------------------------------------------------------------------------------------------------#
# this script summarize all new functions in Point_clouds_processing_v1.2
#------------------------------------------------------------------------------------------------------------------
# this scripts compares vegetation height of UAV-based point clouds when using lidar to create terrain when terrain is interpolated by nearest neighbour

import os
import glob
from liblas import file

import BERA_classes as bc

def CmpVeght_UAV_LDA(Ulas_ifolder_path, Llas_ifolder_path, txt_ifolder_path, txt_ofolder_path):
    if not os.path.exists(Ulas_ifolder_path):
        print "Could not find ", Ulas_ifolder_path, "!"
        if not os.path.exists(Llas_ifolder_path):
            print "Could not find ", Llas_ifolder_path, "!"
            if not os.path.exists(txt_ifolder_path):
                print "Could not find ", txt_ifolder_path, "!"
                return
            return
        return

    if not os.path.exists(txt_ofolder_path):
        os.makedirs(txt_ofolder_path)

    Ulas_file_list = glob.glob(os.path.join(Ulas_ifolder_path, '*.las'))
    Llas_file_list = glob.glob(os.path.join(Llas_ifolder_path, '*.las'))
    txt_file_list = glob.glob(os.path.join(txt_ifolder_path, '*.txt'))

    Ulas_file_num = len(Ulas_file_list)
    Llas_file_num = len(Llas_file_list)
    txt_file_num = len(txt_file_list)

    if (Ulas_file_num <> txt_file_num) or (Llas_file_num <> txt_file_num) or (Ulas_file_num <> Llas_file_num):
        print "Number of las and txt files does not match!"
        return

    for k in range(0, Ulas_file_num):
        print Ulas_file_num

        with open(txt_file_list[k]) as fr_txt:
            with open(Ulas_file_list[k]) as fr_Ulas:
               # read Ulas
                fr_Ulas = file.File(Ulas_file_list[k], mode='r')

                Upt_list_x = []
                Upt_list_y = []
                Upt_list_h = []

                print len(fr_Ulas)
                for f in fr_Ulas:
                    Upt_list_x.append(f.x)
                    Upt_list_y.append(f.y)
                    Upt_list_h.append(f.z)

                print "Reading las file ", Ulas_file_list[k]

                with open(Llas_file_list[k]) as fr_Llas:
                    # read Llas
                    fr_Llas = file.File(Llas_file_list[k], mode='r')

                    Lpt_list_x = []
                    Lpt_list_y = []
                    Lpt_list_h = []

                    print len(fr_Llas)
                    for f in fr_Llas:
                        Lpt_list_x.append(f.x)
                        Lpt_list_y.append(f.y)
                        Lpt_list_h.append(f.z)

                    print "Reading las file ", Llas_file_list[k]

                    fr_txt.readline()
                    site_id = fr_txt.readline().split(',')[0]
                    fr_txt.seek(0)

                    output_file_name = txt_ofolder_path + '\\' + site_id + '_cmp_UL.txt'

                    with open(output_file_name, 'wb') as fw:
                        fw.write('site_id,point_id,x,y,ter_h,gr_veg_h,las_x,las_y,las_h,pt_veg_h,delta_h\n')

                        #read txt, search, and write to txt file
                        fr_txt.readline()
                        for line in fr_txt:
                            row = line.split(',')
                            txt_pt = bc.rtk_veg(str(row[0]), str(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]))

                            min_Upt_x, min_Upt_y, min_Upt_h, min_Udis = bc.search_ht_ary(txt_pt.x, txt_pt.y, Upt_list_x,
                                                                                         Upt_list_y, Upt_list_h)
                            min_Lpt_x, min_Lpt_y, min_Lpt_h, min_Ldis = bc.search_ht_ary(txt_pt.x, txt_pt.y, Lpt_list_x,
                                                                                         Lpt_list_y, Lpt_list_h)
                            pt_veg_h = min_Upt_h - min_Lpt_h
                            delta_h = pt_veg_h - txt_pt.gr_veg_h
                            # print txt_pt.site_id, txt_pt.point_id, txt_pt.x, txt_pt.y, txt_pt.ter_h, '\n'
                            # print min_pt_x, min_pt_y, min_pt_h, min_dis, pt_veg_h, delta_h
                            # print '-----------------------------------'
                            w_line = '%s,%s,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f\n' % \
                                     (str(txt_pt.site_id), str(txt_pt.point_id),  txt_pt.x,  txt_pt.y, txt_pt.ter_h,  txt_pt.gr_veg_h\
                                    , min_Upt_x, min_Upt_y, min_Upt_h, pt_veg_h, delta_h)

                            fw.write(w_line)
                        print 'Output file saved in ', output_file_name
                    fw.close()
                fr_Llas.close()
            fr_Ulas.close()
        fr_txt.close()

#-------------------------------------------------------------------------------------------------------------------#
# this script converts _cmp files to _smy3 files.

def cmp2smy3(UAV_RTK_folder, UAV_UAV_folder, UAV_LDA_folder, output_folder_path):
    if not os.path.exists(UAV_RTK_folder):
        print "Could not find ", UAV_RTK_folder
        if not os.path.exists(UAV_UAV_folder):
            print "Could not find ", UAV_UAV_folder
            return
        return

    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    UR_file_list = glob.glob(os.path.join(UAV_RTK_folder, '*.txt'))
    UU_file_list = glob.glob(os.path.join(UAV_UAV_folder, '*.txt'))
    UL_file_list = glob.glob(os.path.join(UAV_LDA_folder, '*.txt'))
    if len(UR_file_list) <> len(UU_file_list) \
            or len(UU_file_list) <> len(UL_file_list) \
            or len(UR_file_list) <> len(UL_file_list):
        print "File number does not match!"
        return

    for k in range(0,len(UR_file_list)):
        with open(UR_file_list[k]) as fur:
            with open(UU_file_list[k]) as fuu:
                with open(UL_file_list[k]) as ful:
                    fur.readline()
                    site_id = fur.readline().split(',')[0]
                    print "Processing site", site_id
                    fur.seek(0)

                    fur.readline()
                    fuu.readline()
                    ful.readline()

                    output_file_path = output_folder_path + '//' + site_id + '_smy3.txt'
                    with open(output_file_path, 'w') as fw:
                        fw.write('site_id,point_id,x,y,gr_veg_h,UAV_RTK_vh,UAV_UAV_vh,UAV_LDA_vh\n')

                        for line in fur:
                            row = line.split(',')
                            site_id = str(row[0])
                            point_id = str(row[1])
                            x = float(row[2])
                            y = float(row[3])
                            gr_veg_h = float(row[5])
                            ur_vh = float(row[9])

                            urow = fuu.readline().split(',')
                            lrow = ful.readline().split(',')

                            uu_vh = float(urow[9])
                            ul_vh = float(lrow[9])

                            if ur_vh < 0:
                                ur_vh = 0
                            if uu_vh < 0:
                                uu_vh = 0
                            if ul_vh < 0:
                                ul_vh = 0

                            if ur_vh < 5.0 and uu_vh < 5.0 and ul_vh < 5.0:
                                wline = site_id + ',' + point_id + ',' + str(x) + ',' + str(y) + ',' + str(gr_veg_h) +\
                                        ',' + str(ur_vh) + ',' + str(uu_vh) + ','  + str(ul_vh) + '\n'
                                fw.write(wline)
                    fw.close()
                fuu.close()
            fur.close()

#---------------------------------------------------------------------------------------------------------------------#
# this script assign zeros to negative vegetation heights for cmp files

def dl_zero(txt_folder_path, output_folder_path):
    if not os.path.exists(txt_folder_path):
        print "Could not find ", txt_folder_path
        return

    txt_file_list = glob.glob(os.path.join(txt_folder_path, '*.txt'))

    for file in txt_file_list:
        with open(file, 'rb') as fr:
            output_file_path = output_folder_path + '\\' + os.path.basename(file).split('.')[0] + '_new.txt'
            print output_file_path
            with open(output_file_path, 'wb') as fw:
                fw.write(fr.readline())
                site_id = fr.readline().split(',')[0]
                fr.seek(0)
                fr.readline()
                for line in fr:
                    row = line.split(',')
                    if float(row[9]) < 0:
                        row[9] = 0.000
                        row[10] = '%0.3f' % (float(row[9]) - float(row[5]))
                        row[9] = '0.000'
                        newline = ','.join(row)
                        fw.write(newline)
                        fw.write('\n')
                    fw.write(line)
                fw.close()
            fr.close()

#--------------------------------------------------------------------------------------------------------------------#
# this script draws comparison plot of veg height measured from field, by UAV-UAV, UAV-RTK, UAV-LIDAR method

import os
import glob
import matplotlib.pyplot as plt
import numpy

def plot_4p(txt_folder_path, output_folder_path):

    if not os.path.exists(txt_folder_path):
        print "Could not find ", txt_folder_path
        return

    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    txt_file_list = glob.glob(os.path.join(txt_folder_path, '*.txt'))

    fL = output_folder_path + '\\Ltran'
    f6 = output_folder_path + '\\6Ctran'
    f7 = output_folder_path + '\\7Ctran'
    f9 = output_folder_path + '\\9Ctran'

    if not os.path.exists(fL):
        os.makedirs(fL)
    if not os.path.exists(f6):
        os.makedirs(f6)
    if not os.path.exists(f7):
        os.makedirs(f7)
    if not os.path.exists(f9):
        os.makedirs(f9)

    for k in range(0, len(txt_file_list)):
        with open(txt_file_list[k], 'rb') as fr:
            fr.readline()
            site_id = fr.readline().split(',')[0]
            print "Processing site", site_id
            fr.seek(0)
            fr.readline()

            L_output_file_path = fL + '\\' + site_id + '.png'
            C6_output_file_path = f6 + '\\' + site_id + '.png'
            C7_output_file_path = f7 + '\\' + site_id + '.png'
            C9_output_file_path = f9 + '\\' + site_id + '.png'

            L_xlist, L_y1list, L_y2list, L_y3list, L_y4list = [], [], [], [], []
            C6_xlist, C6_y1list, C6_y2list, C6_y3list, C6_y4list = [], [], [], [], []
            C7_xlist, C7_y1list, C7_y2list, C7_y3list, C7_y4list = [], [], [], [], []
            C9_xlist, C9_y1list, C9_y2list, C9_y3list, C9_y4list = [], [], [], [], []
            L_dis = C6_dis = C7_dis = C9_dis = 0.0

            for line in fr:
                row = line.split(',')
                point_id = row[1]
                indct = row[1].split('_')[1]
                y1 = float(row[4])
                y2 = float(row[5])
                y3 = float(row[6])
                y4 = float(row[7])

                if ('L' in indct) and ('.' not in indct):
                    # create long transect list
                    L_y1list.append(y1)
                    L_y2list.append(y2)
                    L_y3list.append(y3)
                    L_y4list.append(y4)
                    L_dis = float(indct.split('L')[1])
                    L_xlist.append(L_dis)

                if ('.' in indct) and (indct[0] == '6'):
                    # create 60-m cross transect list
                    C6_dis = float(indct[2:5])
                    C6_y1list.append(y1)
                    C6_y2list.append(y2)
                    C6_y3list.append(y3)
                    C6_y4list.append(y4)
                    C6_xlist.append(C6_dis)

                if ('.' in indct) and (indct[0] == '7'):
                    # create 75-m cross transect list
                    C7_dis = float(indct[2:5])
                    C7_y1list.append(y1)
                    C7_y2list.append(y2)
                    C7_y3list.append(y3)
                    C7_y4list.append(y4)
                    C7_xlist.append(C7_dis)

                if ('.' in indct) and (indct[0] == '9'):
                    C9_dis = float(indct[2:5])
                    C9_y1list.append(y1)
                    C9_y2list.append(y2)
                    C9_y3list.append(y3)
                    C9_y4list.append(y4)
                    C9_xlist.append(C9_dis)

            L_xlistnew, L_y1listnew, L_y2listnew, L_y3listnew, L_y4listnew \
                = resort4(L_xlist, L_y1list, L_y2list, L_y3list, L_y4list)
            C6_xlistnew, C6_y1listnew, C6_y2listnew, C6_y3listnew, C6_y4listnew \
                = resort4(C6_xlist, C6_y1list, C6_y2list, C6_y3list, C6_y4list)
            C7_xlistnew, C7_y1listnew, C7_y2listnew, C7_y3listnew, C7_y4listnew \
                = resort4(C7_xlist, C7_y1list, C7_y2list, C7_y3list, C7_y4list )
            C9_xlistnew, C9_y1listnew, C9_y2listnew, C9_y3listnew, C9_y4listnew \
                = resort4(C9_xlist, C9_y1list, C9_y2list, C9_y3list, C9_y4list)

            draw_4Lplt(site_id, L_xlistnew, L_y1listnew, L_y2listnew, L_y3listnew, L_y4listnew, L_output_file_path)
            draw_4Cplt(site_id, C6_xlistnew, C6_y1listnew, C6_y2listnew, C6_y3listnew, C6_y4listnew, C6_output_file_path)
            draw_4Cplt(site_id, C7_xlistnew, C7_y1listnew, C7_y2listnew, C7_y3listnew, C7_y4listnew, C7_output_file_path)
            draw_4Cplt(site_id, C9_xlistnew, C9_y1listnew, C9_y2listnew, C9_y3listnew, C9_y4listnew, C9_output_file_path)


def draw_4Lplt(site_id, xlist, y1list, y2list, y3list, y4list, output_file_path):
    if not xlist:
        return
    fig = plt.figure(figsize=(150, 20))
    fig.canvas.set_window_title(site_id)
    plt.ylim(0, 6)
    plt.xlim(0, 150)
    plot_gr_veg, = plt.plot(xlist, y1list, '-bs', linewidth=8.0, markersize=35)
    plot_pt2_veg, = plt.plot(xlist, y2list, '-ro', linewidth=8.0,  markersize=35)
    plot_pt3_veg, = plt.plot(xlist, y3list, '-gd', linewidth=8.0,  markersize=35)
    plot_pt4_veg, = plt.plot(xlist, y4list, '-y^', linewidth=8.0,  markersize=35)

    plt.legend([plot_gr_veg, plot_pt2_veg, plot_pt3_veg, plot_pt4_veg], ["Field", "UAV_RTK", "UAV_UAV", "UAV_LiDAR"], fontsize=120)

    plt.xlabel('distance (m)', fontsize=140)
    plt.ylabel('veg height (m)', fontsize=140)
    xticks_list = range(10, 160, 10)

    plt.xticks(xticks_list, fontsize=140)
    plt.yticks(fontsize=140)
    plt.savefig(output_file_path, bbox_inches='tight')
#    plt.show()
    plt.close()

def draw_4Cplt(site_id, xlist, y1list, y2list, y3list, y4list, output_file_path):
    if not xlist:
        return
    fig = plt.figure(figsize=(35, 25))
    fig.canvas.set_window_title(site_id)

    plt.ylim(0, 6)
    plt.xlim(0, 6)
    plot_gr_veg, = plt.plot(xlist, y1list, '-bs', linewidth=8.0, markersize=30)
    plot_pt2_veg, = plt.plot(xlist, y2list, '-ro', linewidth=8.0,  markersize=30)
    plot_pt3_veg, = plt.plot(xlist, y3list, '-gd', linewidth=8.0,  markersize=30)
    plot_pt4_veg, = plt.plot(xlist, y4list, '-y^', linewidth=8.0,  markersize=30)

    plt.legend([plot_gr_veg, plot_pt2_veg, plot_pt3_veg, plot_pt4_veg], ["Field", "UAV_RTK", "UAV_UAV", "UAV_LiDAR"], fontsize=100)

    plt.xlabel('distance (m)', fontsize=120)
    plt.ylabel('veg height (m)', fontsize=120)
    xticks_list = range(1, 7, 1)
    plt.xticks(xticks_list, fontsize=120)
    plt.yticks(fontsize=120)

    plt.savefig(output_file_path, bbox_inches='tight')
#    plt.show()
    plt.close()

def resort4(xlist, y1list, y2list, y3list, y4list):
    newxlist, newy1list, newy2list, newy3list, newy4list = [], [], [], [], []
    newxindex = numpy.argsort(xlist)
    for i in range(0, len(xlist)):
        newxlist.append(xlist[newxindex[i]])
        newy1list.append(y1list[newxindex[i]])
        newy2list.append(y2list[newxindex[i]])
        newy3list.append(y3list[newxindex[i]])
        newy4list.append(y4list[newxindex[i]])
    return newxlist, newy1list, newy2list, newy3list, newy4list

#--------------------------------------------------------------------------------------------------------------------
# this script conducts t test between GROUND, UAV_RTK and UAV_UAV
#
from scipy import stats
import glob
import os

def ttest2(txt_folder_path, output_file_path):
    if not os.path.exists(txt_folder_path):
        print "Could not find ", txt_folder_path
        return

    txt_file_list = glob.glob(os.path.join(txt_folder_path, '*.txt'))

    with open(output_file_path, 'w') as fw:
        fw.write('site_id,p_GRUR,p_GRUU,p_URUU,t_GRUR,t_GRUU,t_URUU\n')

        for file in txt_file_list:
            with open(file, 'rb') as fr:
                fr.readline()
                site_id = fr.readline().split(',')[0]
                fr.seek(0)
                fr.readline()

                gr_veg_h = []
                UAV_RTK_vh = []
                UAV_UAV_vh = []
                for line in fr:
                    row = line.split(',')
                    gr_veg_h.append(float(row[4]))
                    UAV_RTK_vh.append(float(row[5]))
                    UAV_UAV_vh.append(float(row[6]))

                t_GRUR, p_GRUR = stats.ttest_ind(gr_veg_h, UAV_RTK_vh)
                t_GRUU, p_GRUU = stats.ttest_ind(gr_veg_h, UAV_UAV_vh)
                t_URUU, p_URUU = stats.ttest_ind(UAV_UAV_vh, UAV_RTK_vh)

                wline = '%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n' % (site_id, p_GRUR,p_GRUU,p_URUU,t_GRUR,t_GRUU,t_URUU)
                fw.write(wline)
            fr.close()
        fw.close()

#--------------------------------------------------------------------------------------------------------------------
# this script conducts t test between GROUND, UAV_RTK, UAV_UAV, UAV_LIDAR

# this script conducts t test between GROUND, UAV_RTK, UAV_UAV, UAV_LIDAR

from scipy import stats
import glob
import os

def ttest3(txt_folder_path, output_file_path):
    if not os.path.exists(txt_folder_path):
        print "Could not find ", txt_folder_path
        return

    txt_file_list = glob.glob(os.path.join(txt_folder_path, '*.txt'))

    with open(output_file_path, 'w') as fw:
        fw.write('site_id,p_GRUR, p_GRUU, p_URUU, p_GRUL, p_URUL, p_UUUL,t_GRUR, t_GRUU, t_URUU, t_GRUL, t_URUL, t_UUUL\n')

        for file in txt_file_list:
            with open(file, 'rb') as fr:
                fr.readline()
                site_id = fr.readline().split(',')[0]
                fr.seek(0)
                fr.readline()

                gr_veg_h = []
                UAV_RTK_vh = []
                UAV_UAV_vh = []
                UAV_LDA_vh = []
                for line in fr:
                    row = line.split(',')
                    gr_veg_h.append(float(row[4]))
                    UAV_RTK_vh.append(float(row[5]))
                    UAV_UAV_vh.append(float(row[6]))
                    UAV_LDA_vh.append(float(row[7]))

                t_GRUR, p_GRUR = stats.ttest_ind(gr_veg_h, UAV_RTK_vh)
                t_GRUU, p_GRUU = stats.ttest_ind(gr_veg_h, UAV_UAV_vh)
                t_URUU, p_URUU = stats.ttest_ind(UAV_UAV_vh, UAV_RTK_vh)

                t_GRUL, p_GRUL = stats.ttest_ind(gr_veg_h, UAV_LDA_vh)
                t_URUL, p_URUL = stats.ttest_ind(UAV_RTK_vh, UAV_LDA_vh)
                t_UUUL, p_UUUL = stats.ttest_ind(UAV_UAV_vh, UAV_LDA_vh)

                wline = '%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n' % (site_id,
                            p_GRUR, p_GRUU, p_URUU, p_GRUL, p_URUL, p_UUUL,
                            t_GRUR, t_GRUU, t_URUU, t_GRUL, t_URUL, t_UUUL)
                fw.write(wline)
            fr.close()
        fw.close()
#--------------------------------------------------------------------------------------------------------------------
# This script converts ref. file (lat long and elevation of base station) to csv./txt. file of utm coordinates

import os
import glob
from osgeo import ogr
from osgeo import osr

def ref2utm(ref_folder_path, utm_file_path):
    ref_file_list = glob.glob(os.path.join(ref_folder_path, '*.ref'))
    with open(utm_file_path, 'w') as fw:
        fw.write('site_id,EastingX,NorthingY,ellipsoidalHeight\n')
        for file in ref_file_list:
            with open(file) as fr:
                site_id = os.path.basename(file).split('.')[0]
                fr.readline()
                lat = float(fr.readline())
                long = float(fr.readline())
                elHgt = float(fr.readline())

                # source: epsg projection 4617 - nad83(csrs)
                source = osr.SpatialReference()
                source.ImportFromEPSG(4617)

                # target: epsg projection 2956 - nad83(csrs) / utm zone 12n
                target = osr.SpatialReference()
                target.ImportFromEPSG(2956)

                transform = osr.CoordinateTransformation(source, target)
                pointWkt = 'POINT (' + str(long) + ' ' + str(lat) + ')'
                point = ogr.CreateGeometryFromWkt(pointWkt)
                point.Transform(transform)
                utmEX = point.GetX()
                utmNY = point.GetY()
                wline = '%s,%s,%s,%s\n' % (site_id, str(utmEX), str(utmNY), elHgt)
                fw.write(wline)
                print wline
    fw.close()
    fr.close()

#--------------------------------------------------------------------------------------------------------------------
# This script select lidar tiles for our study sites.

import os
import glob
from liblas import file

def select_lidar_tiles(rtk_folder_path, lidar_folder_path, output_file_path):

    if not os.path.exists(rtk_folder_path):
        print "Could not find ", rtk_folder_path, "!"
        if not os.path.exists(lidar_folder_path):
            print "Could not find ", lidar_folder_path, "!"
            return
        return

    rtk_file_list = glob.glob(os.path.join(rtk_folder_path, '*.txt'))
    lidar_file_list = glob.glob(os.path.join(lidar_folder_path, '*.las'))
    with open(output_file_path, 'w') as fw:
        fw.write('site_id,las_file_name\n')

        for lasfile in lidar_file_list:
            flr = file.File(lasfile, mode='r')
            h = flr.header
            minX = h.min[0]
            minY = h.min[1]
            maxX = h.max[0]
            maxY = h.max[1]
            flr.close() # NOTES: CLOSE THE FILE ONCE YOU DO NOT NEED IT!!!
            for rtkfile in rtk_file_list:
                with open(rtkfile) as frr:
                    line1 = frr.readline()
                    site_id = line1[0:3]
                    pointY = float(line1.split(',')[1])
                    pointX = float(line1.split(',')[2])
                    print lasfile, pointX, pointY, minX, minY, maxX, maxY
                    if (pointY >= minY) and (pointY <= maxY) and (pointX >= minX) and (pointX <= maxX):
                        wline = site_id + ',' + lasfile + '\n'
                        fw.write(wline)
                    frr.close()
        fw.close()

#--------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------
#This script summarize all functions in preprocessing
#----------------------------------------------------------------------------------------------------------------------
# This script add a constant to GCPs, because U of A measured the point under the umbrella, so we need to add the height of the umbrella
# The height of the umbrella is 30" = 0.762 meter

import os
import glob

def GCP_addHeight(rtk_folder_path, output_folder_path):
    rtk_file_list = glob.glob(os.path.join(rtk_folder_path, '*.txt'))
    for file in rtk_file_list:
        with open(file) as fr:
            output_file_path = output_folder_path + '//' + os.path.basename(file)
            with open(output_file_path, 'w') as fw:
                for line in fr:
                    if 'GCP' in line:
                        row = line.split(',')
                        point_id = row[0]
                        y = float(row[1])
                        x = float(row[2])
                        h = float(row[3])
                        note = row[4]
                        new_h = h + 0.762
                        wline = '%s,%.4f,%.4f,%.4f,%s' % (point_id, y, x, new_h, note)
                    else:
                        row = line.split(',')
                        point_id = row[0]
                        y = float(row[1])
                        x = float(row[2])
                        h = float(row[3])
                        note = row[4]
                        wline = '%s,%.4f,%.4f,%.4f,%s' % (point_id, y, x, h, note)
                    fw.write(wline)
                fw.close()
            fr.close()

#----------------------------------------------------------------------------------------------------------------------
# This script matches base station

def match_base(bin_file_path, ref_file_path, output_file_path):
    with open(output_file_path, 'w') as fw:
        fw.write('bin_file_name,site_id,delta_EX,delta_NY,delta_EH\n')
        with open(bin_file_path) as fbr:
            with open(ref_file_path) as frr:
                fbr.readline()
                frr.readline()
                for bline in fbr:
                    brow = bline.split(',')
                    bin_file_name = brow[0]
                    bin_NY = float(brow[1])
                    bin_EX = float(brow[2])
                    bin_EH = float(brow[3])
                    frr.seek(0)
                    frr.readline()
                    for rline in frr:
                        rrow = rline.split(',')
                        site_id = rrow[0]
                        ref_EX = float(rrow[1])
                        ref_NY = float(rrow[2])
                        ref_EH = float(rrow[3])
                        if abs(bin_NY - ref_NY)<5.0 and abs(bin_EX - ref_EX)<5.0 and abs(bin_EH - ref_EH)<5.0:
                            wline = '%s,%s,%.3f,%.3f,%.3f\n' % (bin_file_name, site_id,
                                    abs(bin_NY - ref_NY), abs(bin_EX - ref_EX), abs(bin_EH - ref_EH))
                            fw.write(wline)
                frr.close()
            fbr.close()
        fw.close()

# ----------------------------------------------------------------------------------------------------------------------
# This convert ref. files into csv.

import os
import glob

def ref2csv(ref_folder_path, csv_file_path):
    ref_file_list = glob.glob(os.path.join(ref_folder_path, '*.ref'))
    with open(csv_file_path, 'w') as fw:
        fw.write('latitude,longitude,ellipsoidalHeight\n')
        for file in ref_file_list:
            with open(file) as fr:
                fr.readline()
                lat = float(fr.readline())
                long = float(fr.readline())
                elHgt = float(fr.readline())
                wline = str(lat) + ',' + str(long) + ',' + str(elHgt) + '\n'
                fw.write(wline)
    fw.close()
    fr.close()

#----------------------------------------------------------------------------------------------------------------------
#  this script converts rtk data from U of A to certain format
import os
import glob

def rtkID_convert(txt_folder_path, output_folder_path):
    if not os.path.exists(output_folder_path):
        os.makedirs(txt_folder_path)

    if not os.path.exists(txt_folder_path):
        print 'Could not find', txt_folder_path

    txt_file_list = glob.glob(os.path.join(txt_folder_path, '*.txt'))

    for file in txt_file_list:
        site_id = os.path.basename(file).split('.')[0]
        print 'Processing', site_id
        with open(file, 'rb') as fr:
            output_file_path = output_folder_path + '//' + str(site_id) + '.txt'
            with open(output_file_path, 'w') as fw:
                index = 0
                for line in fr:

                    index += 1
                    if index <= 43:
                        wline = site_id + '_L' + line
                        fw.write(wline)
                    if 'GCP' in line:
                        point_id = line.split(',')[0]
                        new_id = 'GCP_' + site_id + '_' + point_id[3]
                        wline = new_id + ',' + ','.join(line.split(',')[1:])
                        fw.write(wline)

#----------------------------------------------------------------------------------------------------------------------
# This script keeps long transect points and remove cross transect ponits.

import os
import glob

def keeplongT(txt_folder_path, output_folder_path):
    if not os.path.exists(txt_folder_path):
        print "Could not find ", txt_folder_path
        return
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    txt_file_list = glob.glob(os.path.join(txt_folder_path, '*.txt'))

    for txt_file_path in txt_file_list:
        with open(txt_file_path, 'rb') as fr:
            output_file_path = output_folder_path + '//' + os.path.basename(txt_file_path)
            with open(output_file_path, 'wb') as fw:
                fr0 = open(txt_file_list[0])
                header = fr0.readline()
                fw.writelines(header)
                fr.readline()
                for line in fr:
                    row = line.split(',')
                    indct = row[1].split('_')[1]
                    if ('L' in indct) and ('.' not in indct):
                        fw.writelines(line)
        fw.close()
    fr.close()

#----------------------------------------------------------------------------------------------------------------------
# This script computes vegetation cover by height strata for each site.
import os
import glob

def VegCoverOld(cmp_folder_path, output_file_path):
    if not os.path.exists(cmp_folder_path):
        print "Could not find ", cmp_folder_path
        return
    cmp_file_list = glob.glob(os.path.join(cmp_folder_path, '*.txt'))
    with open(output_file_path, 'w') as fw:
        # variable order are: vc,vc_00_05,vc_05_10,vc_10_15,vc_15_20,vc_20_25,vc_25_30,vc_30up
        var_name_list = ['vc', 'vc_00_05', 'vc_05_10', 'vc_10_15', 'vc_15_20', 'vc_20_25', 'vc_25_30', 'vc_30up']
        wline = 'site_id,pt' + ',pt'.join(var_name_list) + ',fd' + ',fd'.join(var_name_list) \
                + ',delta_' + ',delta_'.join(var_name_list) + '\n'
        fw.write(wline)
        for file in cmp_file_list:
            with open(file) as fr:
                fr.readline()

                ptvc_list = [0, 0, 0, 0, 0, 0, 0, 0]
                fdvc_list = [0, 0, 0, 0, 0, 0, 0, 0]
                delta_list = [0, 0, 0, 0, 0, 0, 0, 0]
                index = 0
                for line in fr:
                    index += 1
                    row = line.split(',')
                    site_id = row[0]
                    pt_veg_h = float(row[9])
                    if pt_veg_h > 0.02:
                        ptvc_list[0] += 1
                    if  pt_veg_h > 0.02 and pt_veg_h <= 0.5:
                        ptvc_list[1] += 1
                    for i in range(2, 7):
                        if pt_veg_h > (i-1) * 0.5 and pt_veg_h <= i*0.5:
                            ptvc_list[i] += 1
                    if pt_veg_h > 3.0:
                        ptvc_list[7] += 1

                    fd_veg_h = float(row[5])
                    if fd_veg_h > 0.0:
                        fdvc_list[0] += 1
                    for i in range(1, 7):
                        if fd_veg_h > (i-1) * 0.5 and fd_veg_h <= i*0.5:
                            fdvc_list[i] += 1
                    if fd_veg_h > 3.0:
                        fdvc_list[7] += 1

                for k in range(0, 8):
                    ptvc_list[k] = '%.2f' % (100*float(ptvc_list[k])/float(index))
                    fdvc_list[k] = '%.2f' % (100*float(fdvc_list[k])/float(index))
                    delta_list[k] = '%.2f' % (float(ptvc_list[k]) - float(fdvc_list[k]))
                    print var_name_list[k], site_id, ptvc_list[k], fdvc_list[k], delta_list[k]
                wline = site_id + ',' + ','.join(ptvc_list) +',' + ','.join(fdvc_list) +',' + ','.join(delta_list) + '\n'
                fw.write(wline)
            fr.close()
        fw.close()

#----------------------------------------------------------------------------------------------------------------------
# This script calculates the RMSE of estimated vegetation cover.

import csv
import numpy

def VegCoverErrorsOld(vc_file_path, output_file_path):
    with open(vc_file_path, 'rb') as fr:
        vc_reader = csv.reader(fr)
        row1 = next(vc_reader)
        csv_list = list(vc_reader)
        site_num = vc_reader.line_num - 1
        print 'Number of sites:', site_num
        var_name_list = ['vc', 'vc_00_05', 'vc_05_10', 'vc_10_15', 'vc_15_20', 'vc_20_25', 'vc_25_30', 'vc_30up']

        var_rmse = []

        for k in range(0, len(var_name_list)):
            errors = []
            for i in range(0, site_num):
                errors.append(float(csv_list[i][17+k]))
            rmse = numpy.sqrt((numpy.square(errors)).mean())
            rmse = '%.2f' % (rmse)
            var_rmse.append(rmse)

        for j in range(0, len(var_name_list)):
            print var_name_list[j], var_rmse[j]

        with open(output_file_path, 'w') as fw:
            fw.write(','.join(var_name_list))
            fw.write('\n')
            fw.write(','.join(var_rmse))

#----------------------------------------------------------------------------------------------------------------------
# This script adds color info/rgb to cmp files, output are cmpRGB files.
# This script revises cmp_veg_ht function.

from liblas import file
import os
import glob
import BERA_classes as bc

def cmpRGB_veg_ht(las_ifolder_path, txt_ifolder_path, txt_ofolder_path):

    if not os.path.exists(las_ifolder_path):
        print "Could not find ", las_ifolder_path, "!"
        if not os.path.exists(txt_ifolder_path):
            print "Could not find ", txt_ifolder_path, "!"
            return
        return
    if not os.path.exists(txt_ofolder_path):
        os.makedirs(txt_ofolder_path)

    las_file_list = glob.glob(os.path.join(las_ifolder_path, '*.las'))
    txt_file_list = glob.glob(os.path.join(txt_ifolder_path, '*.txt'))

    las_file_num = len(las_file_list)
    txt_file_num = len(txt_file_list)
    if not (las_file_num == txt_file_num):
        print "Number of las and txt files does not match!"
        return

    for k in range(0, las_file_num):
        print las_file_num
        with open(txt_file_list[k]) as fr_txt:
            with open(las_file_list[k]) as fr_las:
               # read las
                pt_3d_las = []
                fr_las = file.File(las_file_list[k], mode='r')

                pt_list_x, pt_list_y, pt_list_h = [], [], []
                pt_rlist, pt_glist, pt_blist = [], [], []

                print len(fr_las)
                for f in fr_las:
                    pt_list_x.append(f.x)
                    pt_list_y.append(f.y)
                    pt_list_h.append(f.z)
                    pt_rlist.append(f.color.red)
                    pt_glist.append(f.color.green)
                    pt_blist.append(f.color.blue)

                print "Reading las file ", las_file_list[k]

                fr_txt.readline()
                site_id = fr_txt.readline().split(',')[0]
                fr_txt.seek(0)

                output_file_name = txt_ofolder_path + '\\' + site_id + '_cmpRGB.txt'

                with open(output_file_name, 'wb') as fw:
                    fw.write('site_id,point_id,x,y,ter_h,gr_veg_h,las_x,las_y,las_h,pt_veg_h,delta_h,r,g,b\n')

                    #read txt, search, and write to txt file
                    fr_txt.readline()
                    for line in fr_txt:
                        row = line.split(',')
                        txt_pt = bc.rtk_veg(str(row[0]), str(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]))

                        min_pt_x, min_pt_y, min_pt_h, min_dis, min_r, min_g, min_b = bc.search_ht_aryRGB(txt_pt.x, txt_pt.y,
                                                                pt_list_x, pt_list_y, pt_list_h, pt_rlist, pt_glist, pt_blist)
                        pt_veg_h = min_pt_h - txt_pt.ter_h
                        delta_h = pt_veg_h - txt_pt.gr_veg_h
                        # print txt_pt.site_id, txt_pt.point_id, txt_pt.x, txt_pt.y, txt_pt.ter_h, '\n'
                        # print min_pt_x, min_pt_y, min_pt_h, min_dis, pt_veg_h, delta_h
                        # print '-----------------------------------'
                        w_line = '%s,%s,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%s,%s,%s\n' % \
                                 (str(txt_pt.site_id), str(txt_pt.point_id),  txt_pt.x,  txt_pt.y, txt_pt.ter_h,  txt_pt.gr_veg_h
                                , min_pt_x, min_pt_y, min_pt_h, pt_veg_h, delta_h, str(min_r), str(min_g), str(min_b))

                        fw.write(w_line)
                    print 'Output file saved in ', output_file_name
                fw.close()
            fr_las.close()
        fr_txt.close()

#----------------------------------------------------------------------------------------------------------------------
# This script adds color info/rgb to cmp files, output are cmpRGB files.
# This script revises CmpVeght_UAV_UAV function.

import os
import glob
from liblas import file

import BERA_classes as bc

def CmpRGBVeght_UAV_UAV(las_ifolder_path, txt_ifolder_path, txt_ofolder_path):
    if not os.path.exists(las_ifolder_path):
        print "Could not find ", las_ifolder_path, "!"
        if not os.path.exists(txt_ifolder_path):
            print "Could not find ", txt_ifolder_path, "!"
            return
        return
    if not os.path.exists(txt_ofolder_path):
        os.makedirs(txt_ofolder_path)

    las_file_list = glob.glob(os.path.join(las_ifolder_path, '*.las'))
    txt_file_list = glob.glob(os.path.join(txt_ifolder_path, '*.txt'))

    las_file_num = len(las_file_list)
    txt_file_num = len(txt_file_list)
    if not (las_file_num == txt_file_num):
        print "Number of las and txt files does not match!"
        return

    for k in range(0, las_file_num):
        print las_file_num
        with open(txt_file_list[k]) as fr_txt:
            with open(las_file_list[k]) as fr_las:
               # read las
                pt_3d_las = []
                fr_las = file.File(las_file_list[k], mode='r')

                pt_list_x, pt_list_y, pt_list_h = [], [], []
                pt_rlist, pt_glist, pt_blist = [], [], []

                print len(fr_las)
                for f in fr_las:
                    pt_list_x.append(f.x)
                    pt_list_y.append(f.y)
                    pt_list_h.append(f.z)
                    pt_rlist.append(f.color.red)
                    pt_glist.append(f.color.green)
                    pt_blist.append(f.color.blue)

                print "Reading las file ", las_file_list[k]

                fr_txt.readline()
                site_id = fr_txt.readline().split(',')[0]
                fr_txt.seek(0)

                output_file_name = txt_ofolder_path + '\\' + site_id + '_cmp_UU.txt'

                with open(output_file_name, 'wb') as fw:
                    fw.write('site_id,point_id,x,y,ter_h,gr_veg_h,las_x,las_y,las_h,pt_veg_h,delta_h,r,g,b\n')

                    #read txt, search, and write to txt file
                    fr_txt.readline()
                    for line in fr_txt:
                        row = line.split(',')
                        txt_pt = bc.rtk_veg(str(row[0]), str(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]))

                        min_pt_x, min_pt_y, min_pt_h, min_dis, min_r, min_g, min_b = bc.search_ht_aryRGB(txt_pt.x, txt_pt.y,
                                                                pt_list_x, pt_list_y, pt_list_h, pt_rlist, pt_glist, pt_blist)
                        pt_veg_h = min_pt_h
                        delta_h = pt_veg_h - txt_pt.gr_veg_h
                        # print txt_pt.site_id, txt_pt.point_id, txt_pt.x, txt_pt.y, txt_pt.ter_h, '\n'
                        # print min_pt_x, min_pt_y, min_pt_h, min_dis, pt_veg_h, delta_h
                        # print '-----------------------------------'
                        w_line = '%s,%s,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%s,%s,%s\n' % \
                                 (str(txt_pt.site_id), str(txt_pt.point_id),  txt_pt.x,  txt_pt.y, txt_pt.ter_h,  txt_pt.gr_veg_h
                                , min_pt_x, min_pt_y, min_pt_h, pt_veg_h, delta_h, str(min_r), str(min_g), str(min_b))

                        fw.write(w_line)
                    print 'Output file saved in ', output_file_name
                fw.close()
            fr_las.close()
        fr_txt.close()

#----------------------------------------------------------------------------------------------------------------------
#This script add color info/rgb to cmp files
import os
import glob
from liblas import file

import BERA_classes as bc

def CmpRGBVeght_UAV_LDA(Ulas_ifolder_path, Llas_ifolder_path, txt_ifolder_path, txt_ofolder_path):
    if not os.path.exists(Ulas_ifolder_path):
        print "Could not find ", Ulas_ifolder_path, "!"
        if not os.path.exists(Llas_ifolder_path):
            print "Could not find ", Llas_ifolder_path, "!"
            if not os.path.exists(txt_ifolder_path):
                print "Could not find ", txt_ifolder_path, "!"
                return
            return
        return

    if not os.path.exists(txt_ofolder_path):
        os.makedirs(txt_ofolder_path)

    Ulas_file_list = glob.glob(os.path.join(Ulas_ifolder_path, '*.las'))
    Llas_file_list = glob.glob(os.path.join(Llas_ifolder_path, '*.las'))
    txt_file_list = glob.glob(os.path.join(txt_ifolder_path, '*.txt'))

    Ulas_file_num = len(Ulas_file_list)
    Llas_file_num = len(Llas_file_list)
    txt_file_num = len(txt_file_list)

    if (Ulas_file_num <> txt_file_num) or (Llas_file_num <> txt_file_num) or (Ulas_file_num <> Llas_file_num):
        print "Number of las and txt files does not match!"
        return

    for k in range(0, Ulas_file_num):
        print Ulas_file_num

        with open(txt_file_list[k]) as fr_txt:
            with open(Ulas_file_list[k]) as fr_Ulas:
               # read Ulas
                fr_Ulas = file.File(Ulas_file_list[k], mode='r')

                Upt_list_x, Upt_list_y, Upt_list_h = [], [], []
                Upt_rlist, Upt_glist, Upt_blist = [], [], []

                print len(fr_Ulas)
                for f in fr_Ulas:
                    Upt_list_x.append(f.x)
                    Upt_list_y.append(f.y)
                    Upt_list_h.append(f.z)
                    Upt_rlist.append(f.color.red)
                    Upt_glist.append(f.color.green)
                    Upt_blist.append(f.color.blue)

                print "Reading las file ", Ulas_file_list[k]

                with open(Llas_file_list[k]) as fr_Llas:
                    # read Llas
                    fr_Llas = file.File(Llas_file_list[k], mode='r')

                    Lpt_list_x = []
                    Lpt_list_y = []
                    Lpt_list_h = []

                    print len(fr_Llas)
                    for f in fr_Llas:
                        Lpt_list_x.append(f.x)
                        Lpt_list_y.append(f.y)
                        Lpt_list_h.append(f.z)

                    print "Reading las file ", Llas_file_list[k]

                    fr_txt.readline()
                    site_id = fr_txt.readline().split(',')[0]
                    fr_txt.seek(0)

                    output_file_name = txt_ofolder_path + '\\' + site_id + '_cmp_UL.txt'

                    with open(output_file_name, 'wb') as fw:
                        fw.write('site_id,point_id,x,y,ter_h,gr_veg_h,las_x,las_y,las_h,pt_veg_h,delta_h,r,g,b\n')

                        #read txt, search, and write to txt file
                        fr_txt.readline()
                        for line in fr_txt:
                            row = line.split(',')
                            txt_pt = bc.rtk_veg(str(row[0]), str(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]))

                            min_Upt_x, min_Upt_y, min_Upt_h, min_Udis, min_r, min_g, min_b  = bc.search_ht_aryRGB(txt_pt.x, txt_pt.y, Upt_list_x,
                                                                        Upt_list_y, Upt_list_h, Upt_rlist, Upt_glist, Upt_blist)
                            min_Lpt_x, min_Lpt_y, min_Lpt_h, min_Ldis = bc.search_ht_ary(txt_pt.x, txt_pt.y, Lpt_list_x,
                                                                                         Lpt_list_y, Lpt_list_h)
                            pt_veg_h = min_Upt_h - min_Lpt_h
                            delta_h = pt_veg_h - txt_pt.gr_veg_h
                            # print txt_pt.site_id, txt_pt.point_id, txt_pt.x, txt_pt.y, txt_pt.ter_h, '\n'
                            # print min_pt_x, min_pt_y, min_pt_h, min_dis, pt_veg_h, delta_h
                            # print '-----------------------------------'
                            w_line = '%s,%s,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%s,%s,%s\n' % \
                                     (str(txt_pt.site_id), str(txt_pt.point_id),  txt_pt.x,  txt_pt.y, txt_pt.ter_h,  txt_pt.gr_veg_h\
                                    , min_Upt_x, min_Upt_y, min_Upt_h, pt_veg_h, delta_h, str(min_r), str(min_g), str(min_b))

                            fw.write(w_line)
                        print 'Output file saved in ', output_file_name
                    fw.close()
                fr_Llas.close()
            fr_Ulas.close()
        fr_txt.close()

#----------------------------------------------------------------------------------------------------------------------

# This script calculate errors of estimated vegetation heights from smy3 files.

import glob
import os
import numpy

def sym3errors(sym3_folder_path, output_file_path):
    if not os.path.exists(sym3_folder_path):
        print "Could not find ", sym3_folder_path
        return

    sym3_file_list = glob.glob(os.path.join(sym3_folder_path, '*.txt'))
    with open(output_file_path, 'w') as fw:
        fw.write('site_id,sample_size,RMSE_UAV_RTK(m),nRMSE_UAV_RTK(%),bias_UAV_RTK(m),'
                 'RMSE_UAV_LiDAR(m),nRMSE_UAV_LiDAR(%),bias_UAV_LiDAR(m),'
                 'RMSE_UAV_UAV(m),nRMSE_UAV_UAV(%),bias_UAV_UAV(m),best,Second,Third\n')
        Alldelta_UAV_RTK, Alldelta_UAV_UAV, Alldelta_UAV_LDA= [], [], []

        site_id = 'NULL'

        for file in sym3_file_list:
            with open(file) as fr:
                fr.readline()
                delta_UAV_RTK, delta_UAV_UAV, delta_UAV_LDA = [], [], []
                for line in fr:
                    if not line:
                        continue
                    row = line.split(',')
                    site_id = row[0]
                    gr_veg_h = float(row[4])
                    UAV_RTK_vh = float(row[5])
                    UAV_UAV_vh = float(row[6])
                    UAV_LDA_vh = float(row[7])

                    delta_UAV_RTK.append((UAV_RTK_vh - gr_veg_h))
                    Alldelta_UAV_RTK.append((UAV_RTK_vh - gr_veg_h))
                    delta_UAV_UAV.append((UAV_UAV_vh - gr_veg_h))
                    Alldelta_UAV_UAV.append((UAV_UAV_vh - gr_veg_h))
                    delta_UAV_LDA.append((UAV_LDA_vh - gr_veg_h))
                    Alldelta_UAV_LDA.append((UAV_LDA_vh - gr_veg_h))

                sample_size = len(delta_UAV_RTK)
                RMSE_UAV_RTK, nRMSE_UAV_RTK, bias_UAV_RTK = stats_para(delta_UAV_RTK)
                RMSE_UAV_UAV, nRMSE_UAV_UAV, bias_UAV_UAV = stats_para(delta_UAV_UAV)
                RMSE_UAV_LDA, nRMSE_UAV_LDA, bias_UAV_LDA = stats_para(delta_UAV_LDA)
                bestMd, SecondMd, ThirdMd = MethodCmp(RMSE_UAV_RTK, RMSE_UAV_UAV, RMSE_UAV_LDA)
                wline = '%s,%d,%.2f,%d,%.2f,%.2f,%d,%.2f,%.2f,%d,%.2f,%s,%s,%s\n' % (site_id, sample_size,
                            RMSE_UAV_RTK,nRMSE_UAV_RTK,bias_UAV_RTK,
                            RMSE_UAV_LDA,nRMSE_UAV_LDA,bias_UAV_LDA,
                            RMSE_UAV_UAV,nRMSE_UAV_UAV,bias_UAV_UAV,
                            bestMd, SecondMd, ThirdMd)
                fw.write(wline)
                fr.close()
        Allsample_size = len(Alldelta_UAV_RTK)
        AllRMSE_UAV_RTK, AllnRMSE_UAV_RTK, Allbias_UAV_RTK = stats_para(Alldelta_UAV_RTK)
        AllRMSE_UAV_UAV, AllnRMSE_UAV_UAV, Allbias_UAV_UAV = stats_para(Alldelta_UAV_UAV)
        AllRMSE_UAV_LDA, AllnRMSE_UAV_LDA, Allbias_UAV_LDA = stats_para(Alldelta_UAV_LDA)
        AllbestMd, AllSecondMd, AllThirdMd = MethodCmp(AllRMSE_UAV_RTK, AllRMSE_UAV_UAV, AllRMSE_UAV_LDA)
        lastline = '%s,%d,%.2f,%d,%.2f,%.2f,%d,%.2f,%.2f,%d,%.2f,%s,%s,%s\n' % ('Overall', Allsample_size,
                       AllRMSE_UAV_RTK,AllnRMSE_UAV_RTK,Allbias_UAV_RTK,
                       AllRMSE_UAV_LDA,AllnRMSE_UAV_LDA,Allbias_UAV_LDA,
                       AllRMSE_UAV_UAV,AllnRMSE_UAV_UAV,Allbias_UAV_UAV,
                       AllbestMd, AllSecondMd, AllThirdMd)
        fw.write(lastline)

# It is better to use RMSE divided by mean of field reference value
def stats_para(delta_list):
    if delta_list:
        RMSE = numpy.sqrt((numpy.square(delta_list)).mean())
        if (max(delta_list) - min(delta_list)):
            nRMSE = 100.0 * RMSE / (max(delta_list) - min(delta_list))
        else:
            nRMSE = -999
        bias = numpy.mean(delta_list)
        return RMSE, nRMSE, bias
    if not delta_list:
        return -999, -999, -999
    
def MethodCmp(RMSE_UAV_RTK, RMSE_UAV_UAV, RMSE_UAV_LDA):
    Mdname = ['UR', 'UU', 'UL']
    Md = [RMSE_UAV_RTK, RMSE_UAV_UAV, RMSE_UAV_LDA]
    newMdlist = []
    index = numpy.argsort(Md)
    for i in index:
        newMdlist.append(Mdname[i])
    return newMdlist[0], newMdlist[1], newMdlist[2]

#----------------------------------------------------------------------------------------------------------------------
# This script stratifies points into low, medium and high

import os

def strata_all(smy3_file_path, output_folder_path):
    if not os.path.exists(smy3_file_path):
        print "Could not find ", smy3_file_path, "!"
        return
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    with open(smy3_file_path) as fr:
        header = 'site_id,point_id,x,y,gr_veg_h,UAV_RTK_vh,UAV_UAV_vh,UAV_LDA_vh\n'
        lp = output_folder_path + '//' + 'low.txt'
        fl = open(lp, 'w')
        fl.write(header)
        fm = output_folder_path + '//' + 'medium.txt'
        fm = open(fm, 'w')
        fm.write(header)
        fh = output_folder_path + '//' + 'high.txt'
        fh = open(fh,'w')
        fh.write(header)

        fr.readline()
        for line in fr:
            row = line.split(',')
            gr_veg_h = float(row[4])
            if gr_veg_h >= 0 and gr_veg_h < 0.5:
                fl.write(line)
            if gr_veg_h >= 0.5 and gr_veg_h < 2.0:
                fm.write(line)
            if gr_veg_h >= 2.0:
                fh.write(line)
        fl.close()
        fm.close()
        fh.close()
        fr.close()

#----------------------------------------------------------------------------------------------------------------------
# this script computes veg mean height for every site

import glob
import os
import numpy as np

def VegMeanHeight(sym3_folder_path, output_file_path):
    if not os.path.exists(sym3_folder_path):
        print "Could not find ", sym3_folder_path
        return

    sym3_file_list = glob.glob(os.path.join(sym3_folder_path, '*.txt'))
    with open(output_file_path, 'w') as fw:
        fw.write('site_id,gr_veg_mean,UAV_RTK_mean,UAV_LDA_mean,UAV_UAV_mean,delta_UAV_RTK,delta_UAV_UAV,delta_UAV_LDA\n')
        for file in sym3_file_list:
            with open(file) as fr:
                fr.readline()
                gr_veg_h, UAV_RTK_vh, UAV_UAV_vh, UAV_LDA_vh = [], [], [], []
                site_id = 0
                for line in fr:
                    row = line.split(',')
                    site_id = row[0]
                    gr_veg_h.append(float(row[4]))
                    UAV_RTK_vh.append(float(row[5]))
                    UAV_UAV_vh.append(float(row[6]))
                    UAV_LDA_vh.append(float(row[7]))

                gr_veg_avg = np.average(gr_veg_h)
                UAV_RTK_avg = np.average(UAV_RTK_vh)
                delta_UAV_RTK = UAV_RTK_avg - gr_veg_avg
                UAV_UAV_avg = np.average(UAV_UAV_vh)
                delta_UAV_UAV = UAV_UAV_avg - gr_veg_avg
                UAV_LDA_avg = np.average(UAV_LDA_vh)
                delta_UAV_LDA = UAV_LDA_avg - gr_veg_avg
                wline = '%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n' % (site_id, gr_veg_avg, UAV_RTK_avg,UAV_LDA_avg, UAV_UAV_avg,
                                                                     delta_UAV_RTK, delta_UAV_UAV, delta_UAV_LDA)
                fw.write(wline)

#----------------------------------------------------------------------------------------------------------------------

import numpy

def VegMeanHeightErrors(txt_file_path, output_file_path):
    with open(txt_file_path) as fr:
        fr.readline()
        delta_UAV_RTK, delta_UAV_UAV, delta_UAV_LDA = [], [], []
        for line in fr:
            row = line.split(',')
            delta_UAV_RTK.append(float(row[5]))
            delta_UAV_UAV.append(float(row[6]))
            delta_UAV_LDA.append(float(row[7]))
        RMSE_UAV_RTK, nRMSE_UAV_RTK, bias_UAV_RTK = stats_para(delta_UAV_RTK)
        RMSE_UAV_UAV, nRMSE_UAV_UAV, bias_UAV_UAV = stats_para(delta_UAV_UAV)
        RMSE_UAV_LDA, nRMSE_UAV_LDA, bias_UAV_LDA = stats_para(delta_UAV_LDA)
        fr.close()
    with open(output_file_path, 'w') as fw:
        fw.write('RMSE_UAV_RTK(m),nRMSE_UAV_RTK(%),bias_UAV_RTK(m),RMSE_UAV_LiDAR(m),nRMSE_UAV_LiDAR(%),bias_UAV_LiDAR(m),RMSE_UAV_UAV(m),nRMSE_UAV_UAV(%),bias_UAV_UAV(m)\n')
        wline = '%.2f,%d,%.2f,%.2f,%d,%.2f,%.2f,%d,%.2f\n'% (RMSE_UAV_RTK,nRMSE_UAV_RTK,bias_UAV_RTK,
                                                            RMSE_UAV_LDA, nRMSE_UAV_LDA,bias_UAV_LDA,
                                                            RMSE_UAV_UAV,nRMSE_UAV_UAV,bias_UAV_UAV,)
        fw.write(wline)
        fw.close()
#----------------------------------------------------------------------------------------------------------------------
# this script remove the points when the veg height is higher than certain height.
import os
import glob

def rmv_smy3(smy3_folder_path, output_folder_path, height_limit = 3.0):
    if not os.path.exists(smy3_folder_path):
        print "Cannot find ", smy3_folder_path, "!"
        return
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    smy3_file_list = glob.glob(os.path.join(smy3_folder_path, '*.txt'))
    for file in smy3_file_list:
        output_file = output_folder_path + '//' + os.path.basename(file)
        with open(output_file, 'w') as fw:
            with open(file) as fr:
                header = fr.readline()
                fw.write(header)
                for line in fr:
                    row = line.split(',')
                    gr_veg_h = float(row[4])
                    UAV_RTK_vh = float(row[5])
                    UAV_UAV_vh = float(row[6])
                    UAV_LDA_vh = float(row[7])
                    if gr_veg_h < height_limit and UAV_RTK_vh < height_limit and UAV_UAV_vh < height_limit and UAV_LDA_vh < height_limit:
                        fw.write(line)
                fr.close()
            fw.close()
#----------------------------------------------------------------------------------------------------------------------
# This script computes correlation coefficient

import numpy as np

def corrcof(smy3_file_path, output_file_path):
    with open(smy3_file_path) as fr:
        gr_veg_h, UAV_RTK_vh, UAV_UAV_vh, UAV_LDA_vh = [], [], [], []
        fr.readline()
        for line in fr:
            row = line.split(',')
            gr_veg_h.append(float(row[4]))
            UAV_RTK_vh.append(float(row[5]))
            UAV_UAV_vh.append(float(row[6]))
            UAV_LDA_vh.append(float(row[7]))
        fr.close()
        UAV_RTK_cf = np.corrcoef(gr_veg_h, UAV_RTK_vh)
        UAV_LDA_cf = np.corrcoef(gr_veg_h, UAV_LDA_vh)
        UAV_UAV_cf = np.corrcoef(gr_veg_h, UAV_UAV_vh)
    with open(output_file_path, 'w') as fw:
        fw.write('UAV_RTK_cf,UAV_LDA_cf,UAV_UAV_cf\n')
        wline = '%.3f,%.3f,%.3f\n' % (UAV_RTK_cf[0][1], UAV_LDA_cf[0][1], UAV_UAV_cf[0][1])
        fw.write(wline)
        fw.close()

#----------------------------------------------------------------------------------------------------------------------
# This script computes vegetation cover by height strata for each site.

import os
import glob

def VegCover(smy3_folder_path, output_folder_path):
    if not os.path.exists(smy3_folder_path):
        print 'Cannot find', smy3_folder_path, '!'
        return
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    smy3_file_list = glob.glob(os.path.join(smy3_folder_path, '*.txt'))
    UR_file = output_folder_path + '//' + 'UAV_RTK.csv'
    fur = open(UR_file, 'w')
    UL_file = output_folder_path + '//' + 'UAV_LiDAR.csv'
    ful = open(UL_file, 'w')
    UU_file = output_folder_path + '//' + 'UAV_UAV.csv'
    fuu = open(UU_file, 'w')
    var_name_list = ['vc', 'vc_00_05', 'vc_05_20', 'vc_20up']
    wline = 'site_id,pt' + ',pt'.join(var_name_list) + ',fd' + ',fd'.join(var_name_list) \
            + ',delta_' + ',delta_'.join(var_name_list) + '\n'
    fur.write(wline)
    ful.write(wline)
    fuu.write(wline)

    for file in smy3_file_list:
        with open(file) as fr:
            fr.readline()
            fdvc_list = [0, 0, 0, 0]
            UR_ptvc_list, UR_delta_list = [0, 0, 0, 0], [0, 0, 0, 0]
            UL_ptvc_list, UL_delta_list = [0, 0, 0, 0], [0, 0, 0, 0]
            UU_ptvc_list, UU_delta_list = [0, 0, 0, 0], [0, 0, 0, 0]
            index = 0
            for line in fr:
                index += 1
                row = line.split(',')
                site_id = str(row[0])
                fd_veg_h = float(row[4])
                UAV_RTK_vh = float(row[5])
                UAV_UAV_vh = float(row[6])
                UAV_LDA_vh = float(row[7])

                vc_condition(fd_veg_h, fdvc_list)
                vc_condition(UAV_RTK_vh, UR_ptvc_list)
                vc_condition(UAV_LDA_vh, UL_ptvc_list)
                vc_condition(UAV_UAV_vh, UU_ptvc_list)

            fr.close()
            for i in range(0, 4):
                fdvc_list[i] = '%d' % (100*float(fdvc_list[i])/float(index))
                UR_ptvc_list[i] = '%d' % (100*float(UR_ptvc_list[i])/float(index))
                UR_delta_list[i] = '%d' % (float(UR_ptvc_list[i]) - float(fdvc_list[i]))
                UL_ptvc_list[i] = '%d' % (100 * float(UL_ptvc_list[i]) / float(index))
                UL_delta_list[i] = '%d' % (float(UL_ptvc_list[i]) - float(fdvc_list[i]))
                UU_ptvc_list[i] = '%d' % (100 * float(UU_ptvc_list[i]) / float(index))
                UU_delta_list[i] = '%d' % (float(UU_ptvc_list[i]) - float(fdvc_list[i]))
            ur_line = site_id + ',' + ','.join(UR_ptvc_list) +',' + ','.join(fdvc_list) +',' + ','.join(UR_delta_list) + '\n'
            ul_line = site_id + ',' + ','.join(UL_ptvc_list) +',' + ','.join(fdvc_list) +',' + ','.join(UL_delta_list) + '\n'
            uu_line = site_id + ',' + ','.join(UU_ptvc_list) +',' + ','.join(fdvc_list) +',' + ','.join(UU_delta_list) + '\n'
            fur.write(ur_line)
            ful.write(ul_line)
            fuu.write(uu_line)
    fur.close()
    ful.close()
    fuu.close()

def vc_condition(veg_h, vc_list):
    if veg_h > 0.02:
        vc_list[0] += 1
    if veg_h > 0.02 and veg_h < 0.5:
        vc_list[1] += 1
    if veg_h >= 0.5 and veg_h <= 2.0:
        vc_list[2] += 1
    if veg_h > 2.0:
        vc_list[3] += 1


#----------------------------------------------------------------------------------------------------------------------
# This script calculates the RMSE of estimated vegetation cover.

import csv
import numpy

def VegCoverErrors(vc_file_path, output_file_path):
    with open(vc_file_path, 'rb') as fr:
        vc_reader = csv.reader(fr)
        row1 = next(vc_reader)
        csv_list = list(vc_reader)
        site_num = vc_reader.line_num - 1
        print 'Number of sites:', site_num
        var_name_list = ['vc', 'vc_00_05', 'vc_05_20', 'vc_20up']

        var_rmse = []

        for k in range(0, len(var_name_list)):
            errors = []
            for i in range(0, site_num):
                errors.append(float(csv_list[i][9+k]))

            rmse = numpy.sqrt((numpy.square(errors)).mean())
            rmse = '%d' % (rmse)
            var_rmse.append(rmse)

        for j in range(0, len(var_name_list)):
            print var_name_list[j], var_rmse[j]

        with open(output_file_path, 'w') as fw:
            fw.write(','.join(var_name_list))
            fw.write('\n')
            fw.write(','.join(var_rmse))

#----------------------------------------------------------------------------------------------------------------------
# this script computes statistics of field measurements of veg height. smy3 as input files.
import os
import glob
import numpy

def field_stats(smy3_folder_path, output_file_path):
    if not os.path.exists(smy3_folder_path):
        print "Cannot find", smy3_folder_path, "!"
        return
    smy3_file_list = glob.glob(os.path.join(smy3_folder_path, '*.txt'))
    with open(output_file_path, 'w') as fw:
        fw.write('site_id,mean,max,min,range,standard_deviation\n')
        for file in smy3_file_list:
            with open(file) as fr:
                fr.readline()
                fd_veg_list = []
                for line in fr:
                    row = line.split(',')
                    site_id = row[0]
                    fd_veg_h = float(row[4])
                    fd_veg_list.append(fd_veg_h)
                fd_mean = numpy.mean(fd_veg_list)
                fd_max = max(fd_veg_list)
                fd_min = min(fd_veg_list)
                fd_range = fd_max - fd_min
                fd_std = numpy.std(fd_veg_list)
                wline = '%s,%.3f,%.3f,%.3f,%.3f,%.3f\n' % (site_id, fd_mean, fd_max, fd_min, fd_range, fd_std)
                fw.write(wline)
            fr.close()
        fw.close()

#----------------------------------------------------------------------------------------------------------------------
# The following functions are used in factor_analysis_v1.0
#----------------------------------------------------------------------------------------------------------------------
# This script converts bearings from numbers to characters.
# 0: NS;   45: NE;    90: WE;    135: SE;

import os
import glob
import csv

def convert_bearing(csv_folder_path, output_folder_path):
    if not os.path.exists(csv_folder_path):
        print "Cannot find", csv_folder_path, "!"
        return
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    csv_file_list = glob.glob(os.path.join(csv_folder_path, '*.csv'))
    for file in csv_file_list:
        with open(file, 'r') as fr:
            csv_reader = csv.reader(fr)
            header = next(csv_reader)
            output_file_path = output_folder_path + '//' + os.path.basename(file)

            with open(output_file_path, 'wb') as fw:
                wheader = ','.join(header) + '\n'
                fw.write(wheader)
                csv_writer = csv.writer(fw)
                for row in csv_reader:
                    bearing = int(float(row[6]))
                    new_bearing = 'NULL'
                    if bearing == 0:
                        new_bearing = 'NS'
                    if bearing == 45:
                        new_bearing = 'NE'
                    if bearing == 90:
                        new_bearing = 'WE'
                    if bearing == 135:
                        new_bearing = 'SE'
                    row[6] = new_bearing
                    wline = ','.join(row) + '\n'
                    fw.write(wline)
                fw.close()
            fr.close()

#----------------------------------------------------------------------------------------------------------------------
# This script appends csv files that have the same header. Only keep one header
import os
import glob
import csv

def csv_append(csv_folder_path, output_file_path):
    if not os.path.exists(csv_folder_path):
        print 'Cannot find', csv_folder_path, "!"
        return
    csv_file_list = glob.glob(os.path.join(csv_folder_path,'*.csv'))
    if not csv_file_list:
        return
    with open(csv_file_list[0], 'r') as f0:
        csv_reader0 = csv.reader(f0)
        header0 = next(csv_reader0)

    with open(output_file_path, 'wb') as fw:
        line0 = ','.join(header0) + '\n'
        fw.write(line0)
        for file in csv_file_list:
            with open(file, 'r') as fr:
                csv_reader = csv.reader(fr)
                next(csv_reader)
                csv_writer = csv.writer(fw)
                csv_writer.writerows(csv_reader)

#----------------------------------------------------------------------------------------------------------------------
# This script computes bearing of the long transect of seismic line.
# Bearing = arctan(deltaX / deltaY)

import math
import os
import glob

def line_bearing(rtk_folder_path, output_file_path):
    if not os.path.exists(rtk_folder_path):
        print "Cannot find", rtk_folder_path, "!"
        return

    rtk_file_list = glob.glob(os.path.join(rtk_folder_path, '*.txt'))
    with open(output_file_path, 'w') as fw:
        fw.write('site_id,bearing,category\n')
        for file in rtk_file_list:
            startY, startX, endX, endY = 0, 0, 0, 0
            with open(file, 'r') as fr:
                row1 = fr.readline()
                site_id = row1[0:3]
                fr.seek(0)
                for line in fr:
                    row = line.split(',')
                    point_id = row[0]
                    if 'L000' in point_id:
                        startY = float(row[1])
                        startX = float(row[2])
                    if 'L150' in point_id:
                        endY = float(row[1])
                        endX = float(row[2])
                if (startX * startY * endX * endY) == 0:
                    bearing = -9999
                    category = -9999
                else:
                    ratio = (endX - startX) / (endY - startY)
                    bearing = math.atan(ratio)
                    bearing = math.degrees(bearing)
                    if bearing < 0:
                        bearing += 360
                    if bearing > 180:
                        bearing -= 180
                    if bearing >= 0 and bearing < 22.5:
                        category = 0
                    if bearing >= 22.5 and bearing < 67.5:
                        category = 45
                    if bearing >= 67.5 and bearing <= 112.5:
                        category = 90
                    if bearing > 112.5 and bearing <= 157.5:
                        category = 135
                    if bearing > 157.5 and bearing <= 180:
                        category = 0
                wline = '%s, %.2f, %.f\n' % (site_id, bearing, category)
                fw.write(wline)
#----------------------------------------------------------------------------------------------------------------------
# This script compute the percent of shadow for each site.

import os
import glob
import csv

def shadow_percent(csv_folder_path, output_file_path):
    if not os.path.exists(csv_folder_path):
        print "Cannot find", csv_folder_path, "!"
    csv_file_list = glob.glob(os.path.join(csv_folder_path, '*.csv'))

    with open(output_file_path, 'w') as fw:
        fw.write('site_id,shadow_percent\n')
        for file in csv_file_list:
            with open(file, 'r') as fr:
                shadow_index = 0
                csv_reader = csv.reader(fr)
                header = next(csv_reader)
                csv_list = list(csv_reader)
                site_id = csv_list[0][0][0] + csv_list[0][0][1] + csv_list[0][0][2]
                for row in csv_list:
                    if row[1] == 'Y':
                        shadow_index += 1
                shadow_per = float(shadow_index)/float(len(csv_list))
                wline = '%s,%.3f\n' % (site_id, shadow_per)
                fw.write(wline)
                fr.close()
        fw.close()
#----------------------------------------------------------------------------------------------------------------------
# This script link veg field height (extract from smy3 files) to factor analysis based on point id.

import glob
import os
import csv  #Note that when reading csv files, should use csv module. Do not use the way of reading txt files.

def veg_height_factor(smy3_folder_path, factor_folder_path, output_folder_path):
    if not os.path.exists(smy3_folder_path):
        print "Cannot find", smy3_folder_path, "!"
        return
    if not os.path.exists(factor_folder_path):
        print "Cannot find", factor_folder_path, "!"
        return
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    smy3_file_list = glob.glob(os.path.join(smy3_folder_path, '*.txt'))
    factor_file_list = glob.glob(os.path.join(factor_folder_path, '*.csv'))
    if len(smy3_file_list) != len(factor_file_list):
        print "Number of files does not match!"
        return
    num = len(smy3_file_list)
    print "There are", num, "sites."

    for i in range(0, num):
        with open(factor_file_list[i], 'r') as ff:
            f_reader = csv.reader(ff)
            header = next(f_reader)
            f_list = list(f_reader)
            f_num = len(f_list)
            site_id = f_list[0][0][0] + f_list[0][0][1] + f_list[0][0][2]
            output_file_path = output_folder_path + '//' + site_id + ".csv"
            with open(output_file_path, 'w') as fw:
                header_line = ','.join(header) + ',VEG_HEIGHT,UR_ERROR,UL_ERROR,UU_ERROR\n'
                fw.write(header_line)
                with open(smy3_file_list[i],'r') as fs:
                    fs.readline()
                    for sline in fs:
                        row = sline.split(',')
                        spoint_id = row[1]
                        fd_height = row[4]
                        ur_error = float(row[5]) - float(row[4])
                        ul_error = float(row[7]) - float(row[4])
                        uu_error = float(row[6]) - float(row[4])
                        add_line = ',%s,%.3f,%.3f,%.3f\n' % (fd_height, ur_error, ul_error, uu_error)
                        for k in range(0, f_num):
                            if f_list[k][0] == spoint_id:
                                wline = ','.join(f_list[k]) + add_line
                                fw.write(wline)
                    fs.close()
                fw.close()
            ff.close()

#----------------------------------------------------------------------------------------------------------------------
# This script add a column of site ID of factor_point_append file.

import csv

def add_siteID(csv_file_path, output_file_path):
    with open(csv_file_path, 'rb') as fr:
        csv_reader = csv.reader(fr)
        header = next(csv_reader)
        csv_list = list(csv_reader)
        with open(output_file_path, 'wb') as fw:
            new_header = "SITE_ID," + ','.join(header) + '\n'
            fw.write(new_header)
            for row in csv_list:
                site_id = row[0][0] + row[0][1] + row[0][2]
                new_row = site_id + ',' + ','.join(row) + '\n'
                fw.write(new_row)
            fw.close()
        fr.close()

#----------------------------------------------------------------------------------------------------------------------
# This script converts categorical data to numerical data for point level data.

import csv

def categ2num_point(csv_file_path, output_file_path):
    with open(csv_file_path, 'rb') as fr:
        csv_reader = csv.reader(fr)
        header = next(csv_reader)
        csv_list = list(csv_reader)
        with open(output_file_path, 'wb') as fw:
            new_header = 'SITE_ID,POINT_ID,SHADOW_TYPE,CL_TYPE,GCP_NUM,ADJ_TYPE,MOIST_TYPE,LINE_TYPE,VEG_HEIGHT,UR_ERROR\n'
            fw.write(new_header)
            for row in csv_list:
                is_shadow = row[2]
                shadow_type = shadow_type_rule(is_shadow)
                row[2] = str(shadow_type)

                cross_long = row[3]
                cl_type = cl_type_rule(cross_long)
                row[3] = str(cl_type)

                adj_frt = row[5]
                adj_type = adj_type_rule(adj_frt)
                row[5] = str(adj_type)

                moist = row[6][1]
                moist_type = moist_type_rule(moist)
                row[6] = str(moist_type)

                line_direct = row[7]
                line_type = line_type_rule(line_direct)
                row[7] = str(line_type)

                wline = ','.join(row) + '\n'
                fw.write(wline)
            fw.close()
        fr.close()


def shadow_type_rule(is_shadow):
    if is_shadow == 'N':
        shadow_type = 1
    elif is_shadow == 'Y':
        shadow_type = 2
    else:
        print "Shadow type error!"
        shadow_type = 0
    return shadow_type

def cl_type_rule(cross_long):
    if cross_long == 'C':
        cl_type = 1
    elif cross_long == 'L':
        cl_type = 2
    else:
        print "CL type error!"
        cl_type = 0
    return cl_type

def moist_type_rule(moist):
    if moist == 'X':
        moist_type = 1
    elif moist == 'M':
        moist_type = 2
    elif moist == 'G' or moist == 'D':
        moist_type = 3
    else:
        print 'Moist type error!'
        moist_type = 0
    return moist_type

def adj_type_rule(adj_frt):
    if adj_frt == 'OP':
        adj_type = 1
    elif adj_frt == 'MD':
        adj_type = 2
    elif adj_frt == 'CS':
        adj_type = 3
    else:
        print 'Adj type error!'
        adj_type = 0
    return  adj_type

def line_type_rule(line_direct):
    if line_direct == 'NS':
        line_type = 1
    elif line_direct == 'NE':
        line_type = 2
    elif line_direct == 'WE':
        line_type = 3
    elif line_direct == 'SE':
        line_type = 4
    else:
        print "Line type error!"
        line_type = 0
    return line_type

#----------------------------------------------------------------------------------------------------------------------
# This script converts categorical data to numerical data for site level data.

import csv

def categ2num_site(csv_file_path, output_file_path):
    with open(csv_file_path, 'rb') as fr:
        csv_reader = csv.reader(fr)
        header = next(csv_reader)
        csv_list = list(csv_reader)
        with open(output_file_path, 'wb') as fw:
            new_header = 'SITE_ID,SHADOW_PER,MEAN_HEIGHT,MOIST_TYPE,GCP_NUM,ADJ_TYPE,LINE_TYPE,UR_ERROR\n'
            fw.write(new_header)
            for row in csv_list:
                moist = row[3][1]
                moist_type = moist_type_rule(moist)
                row[3] = str(moist_type)

                adj_frt = row[5]
                adj_type = adj_type_rule(adj_frt)
                row[5] = str(adj_type)

                line_direct = row[6]
                line_type = line_type_rule(line_direct)
                row[6] = str(line_type)

                wline = ','.join(row) + '\n'
                fw.write(wline)
            fw.close()
        fr.close()

#----------------------------------------------------------------------------------------------------------------------
# This script compares field veg height and UAV_UAV estimates of veg height using a 'Cylinder' method.

import os
import glob
from liblas import file
import numpy

import BERA_classes as bc

def CmpVegHtCyl_UAV_UAV(las_ifolder_path, txt_ifolder_path, txt_ofolder_path, top_radius = 0.2, bottom_radius = 1.0):
    if not os.path.exists(las_ifolder_path):
        print "Could not find ", las_ifolder_path, "!"
        if not os.path.exists(txt_ifolder_path):
            print "Could not find ", txt_ifolder_path, "!"
            return
        return
    if not os.path.exists(txt_ofolder_path):
        os.makedirs(txt_ofolder_path)

    las_file_list = glob.glob(os.path.join(las_ifolder_path, '*.las'))
    txt_file_list = glob.glob(os.path.join(txt_ifolder_path, '*.txt'))

    las_file_num = len(las_file_list)
    txt_file_num = len(txt_file_list)
    if not (las_file_num == txt_file_num):
        print "Number of las and txt files does not match!"
        return

    for k in range(0, las_file_num):
        print las_file_num
        with open(txt_file_list[k]) as fr_txt:
            with open(las_file_list[k]) as fr_las:
               # read las
                pt_3d_las = []
                fr_las = file.File(las_file_list[k], mode='r')

                pt_list_x = []
                pt_list_y = []
                pt_list_h = []

                print len(fr_las)
                for f in fr_las:
                    pt_list_x.append(f.x)
                    pt_list_y.append(f.y)
                    pt_list_h.append(f.z)

                print "Reading las file ", las_file_list[k]

                fr_txt.readline()
                site_id = fr_txt.readline().split(',')[0]
                fr_txt.seek(0)

                output_file_name = txt_ofolder_path + '\\' + site_id + '_cmp_UU.txt'

                with open(output_file_name, 'wb') as fw:
                    fw.write('site_id,point_id,x,y,ter_h,gr_veg_h,las_x,las_y,las_h,pt_veg_h,delta_h\n')

                    #read txt, search, and write to txt file
                    fr_txt.readline()
                    for line in fr_txt:
                        row = line.split(',')
                        txt_pt = bc.rtk_veg(str(row[0]), str(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]))

                        top_hlist = bc.hlist_inBuffer(txt_pt.x, txt_pt.y, pt_list_x, pt_list_y, pt_list_h, top_radius)
                        if top_hlist:
                            top_h = numpy.percentile(top_hlist, 99)
                        else:
                            top_h = 999
                        bottom_hlist = bc.hlist_inBuffer(txt_pt.x, txt_pt.y, pt_list_x, pt_list_y, pt_list_h, bottom_radius)
                        if bottom_hlist:
                            bottom_h = min(bottom_hlist)
                        else:
                            bottom_h = -999
                        pt_veg_h = top_h - bottom_h

                        delta_h = pt_veg_h - txt_pt.gr_veg_h
                        # print txt_pt.site_id, txt_pt.point_id, txt_pt.x, txt_pt.y, txt_pt.ter_h, '\n'
                        # print min_pt_x, min_pt_y, min_pt_h, min_dis, pt_veg_h, delta_h
                        # print '-----------------------------------'
                        w_line = '%s,%s,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f\n' % \
                                 (str(txt_pt.site_id), str(txt_pt.point_id),  txt_pt.x,  txt_pt.y, txt_pt.ter_h,  txt_pt.gr_veg_h\
                                , txt_pt.x, txt_pt.y, top_h, pt_veg_h, delta_h)

                        fw.write(w_line)
                    print 'Output file saved in ', output_file_name
                fw.close()
            fr_las.close()
        fr_txt.close()

#----------------------------------------------------------------------------------------------------------------------
# This script compares field veg height and UAV_LiDAR estimates of veg height using a 'Cylinder' method.

import os
import glob
from liblas import file
import numpy

import BERA_classes as bc

def CmpVeghtCyl_UAV_LDA(Ulas_ifolder_path, Llas_ifolder_path, txt_ifolder_path, txt_ofolder_path, top_radius = 0.2, bottom_radius = 1.0):
    if not os.path.exists(Ulas_ifolder_path):
        print "Could not find ", Ulas_ifolder_path, "!"
        if not os.path.exists(Llas_ifolder_path):
            print "Could not find ", Llas_ifolder_path, "!"
            if not os.path.exists(txt_ifolder_path):
                print "Could not find ", txt_ifolder_path, "!"
                return
            return
        return

    if not os.path.exists(txt_ofolder_path):
        os.makedirs(txt_ofolder_path)

    Ulas_file_list = glob.glob(os.path.join(Ulas_ifolder_path, '*.las'))
    Llas_file_list = glob.glob(os.path.join(Llas_ifolder_path, '*.las'))
    txt_file_list = glob.glob(os.path.join(txt_ifolder_path, '*.txt'))

    Ulas_file_num = len(Ulas_file_list)
    Llas_file_num = len(Llas_file_list)
    txt_file_num = len(txt_file_list)

    if (Ulas_file_num <> txt_file_num) or (Llas_file_num <> txt_file_num) or (Ulas_file_num <> Llas_file_num):
        print "Number of las and txt files does not match!"
        return

    for k in range(0, Ulas_file_num):
        print Ulas_file_num

        with open(txt_file_list[k]) as fr_txt:
            with open(Ulas_file_list[k]) as fr_Ulas:
               # read Ulas
                fr_Ulas = file.File(Ulas_file_list[k], mode='r')

                Upt_list_x = []
                Upt_list_y = []
                Upt_list_h = []

                print len(fr_Ulas)
                for f in fr_Ulas:
                    Upt_list_x.append(f.x)
                    Upt_list_y.append(f.y)
                    Upt_list_h.append(f.z)

                print "Reading las file ", Ulas_file_list[k]

                with open(Llas_file_list[k]) as fr_Llas:
                    # read Llas
                    fr_Llas = file.File(Llas_file_list[k], mode='r')

                    Lpt_list_x = []
                    Lpt_list_y = []
                    Lpt_list_h = []

                    print len(fr_Llas)
                    for f in fr_Llas:
                        Lpt_list_x.append(f.x)
                        Lpt_list_y.append(f.y)
                        Lpt_list_h.append(f.z)

                    print "Reading las file ", Llas_file_list[k]

                    fr_txt.readline()
                    site_id = fr_txt.readline().split(',')[0]
                    fr_txt.seek(0)

                    output_file_name = txt_ofolder_path + '\\' + site_id + '_cmp_UL.txt'

                    with open(output_file_name, 'wb') as fw:
                        fw.write('site_id,point_id,x,y,ter_h,gr_veg_h,las_x,las_y,las_h,pt_veg_h,delta_h\n')

                        #read txt, search, and write to txt file
                        fr_txt.readline()
                        for line in fr_txt:

                            row = line.split(',')
                            txt_pt = bc.rtk_veg(str(row[0]), str(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]))

                            top_hlist = bc.hlist_inBuffer(txt_pt.x, txt_pt.y, Upt_list_x, Upt_list_y, Upt_list_h, top_radius)
                            if top_hlist:
                                top_h = numpy.percentile(top_hlist, 99)
                            else:
                                top_h = 999
                            bottom_hlist = bc.hlist_inBuffer(txt_pt.x, txt_pt.y, Lpt_list_x, Lpt_list_y, Lpt_list_h, bottom_radius)
                            if bottom_hlist:
                                bottom_h = min(bottom_hlist)
                            else:
                                bottom_h = -999
                            pt_veg_h = top_h - bottom_h

                            delta_h = pt_veg_h - txt_pt.gr_veg_h
                            # print txt_pt.site_id, txt_pt.point_id, txt_pt.x, txt_pt.y, txt_pt.ter_h, '\n'
                            # print min_pt_x, min_pt_y, min_pt_h, min_dis, pt_veg_h, delta_h
                            # print '-----------------------------------'
                            w_line = '%s,%s,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f\n' % \
                                     (str(txt_pt.site_id), str(txt_pt.point_id),  txt_pt.x,  txt_pt.y, txt_pt.ter_h,  txt_pt.gr_veg_h\
                                    , txt_pt.x, txt_pt.y, top_h, pt_veg_h, delta_h)

                            fw.write(w_line)
                        print 'Output file saved in ', output_file_name
                    fw.close()
                fr_Llas.close()
            fr_Ulas.close()
        fr_txt.close()

#----------------------------------------------------------------------------------------------------------------------
# This script compares field veg height and UAV_RTK estimates of veg height using a 'Cylinder' method.

import glob
import os
from liblas import file
import numpy

import BERA_classes as bc

def CmpVegHtCyl_UAV_RTK(las_ifolder_path, txt_ifolder_path, txt_ofolder_path, top_radius = 0.2):

    if not os.path.exists(las_ifolder_path):
        print "Could not find ", las_ifolder_path, "!"
        if not os.path.exists(txt_ifolder_path):
            print "Could not find ", txt_ifolder_path, "!"
            return
        return
    if not os.path.exists(txt_ofolder_path):
        os.makedirs(txt_ofolder_path)

    las_file_list = glob.glob(os.path.join(las_ifolder_path, '*.las'))
    txt_file_list = glob.glob(os.path.join(txt_ifolder_path, '*.txt'))

    las_file_num = len(las_file_list)
    txt_file_num = len(txt_file_list)
    if not (las_file_num == txt_file_num):
        print "Number of las and txt files does not match!"
        return

    for k in range(0, las_file_num):
        print las_file_num
        with open(txt_file_list[k]) as fr_txt:
            with open(las_file_list[k]) as fr_las:
               # read las
                pt_3d_las = []
                fr_las = file.File(las_file_list[k], mode='r')

                pt_list_x = []
                pt_list_y = []
                pt_list_h = []

                print len(fr_las)
                for f in fr_las:
                    pt_list_x.append(f.x)
                    pt_list_y.append(f.y)
                    pt_list_h.append(f.z)

                print "Reading las file ", las_file_list[k]

                fr_txt.readline()
                site_id = fr_txt.readline().split(',')[0]
                fr_txt.seek(0)

                output_file_name = txt_ofolder_path + '\\' + site_id + '_cmp_UR.txt'

                with open(output_file_name, 'wb') as fw:
                    fw.write('site_id,point_id,x,y,ter_h,gr_veg_h,las_x,las_y,las_h,pt_veg_h,delta_h\n')

                    #read txt, search, and write to txt file
                    fr_txt.readline()
                    for line in fr_txt:
                        row = line.split(',')
                        txt_pt = bc.rtk_veg(str(row[0]), str(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]))

                        top_hlist = bc.hlist_inBuffer(txt_pt.x, txt_pt.y, pt_list_x, pt_list_y, pt_list_h, top_radius)
                        if top_hlist:
                            top_h = numpy.percentile(top_hlist, 99)
                        else:
                            top_h = 999

                        pt_veg_h = top_h - txt_pt.ter_h
                        delta_h = pt_veg_h - txt_pt.gr_veg_h
                        # print txt_pt.site_id, txt_pt.point_id, txt_pt.x, txt_pt.y, txt_pt.ter_h, '\n'
                        # print min_pt_x, min_pt_y, min_pt_h, min_dis, pt_veg_h, delta_h
                        # print '-----------------------------------'
                        w_line = '%s,%s,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f\n' % \
                                 (str(txt_pt.site_id), str(txt_pt.point_id),  txt_pt.x,  txt_pt.y, txt_pt.ter_h,  txt_pt.gr_veg_h\
                                , txt_pt.x, txt_pt.y, top_h, pt_veg_h, delta_h)

                        fw.write(w_line)
                    print 'Output file saved in ', output_file_name
                fw.close()
            fr_las.close()
        fr_txt.close()
