# import numpy
import math

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

class rtk_veg_cmp(rtk_veg):
    def __init__(self, site_id, point_id, x, y, ter_h, gr_veg_h, pt_veg_h):
        super(rtk_veg, self).__init__(site_id, point_id, x, y, ter_h, gr_veg_h)
        self.pt_veg_h = pt_veg_h
    def cal_delta_h(self, pt_veg_h, gr_veg_h):
        # cal_delta_h calculates the differences between point cloud measurement and ground veg height
        delta_h = pt_veg_h - gr_veg_h


class pt_2d(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y

class pt_3d(pt_2d):
    def __init__(self, x, y, h):
        super(pt_3d, self).__init__(x, y)
        self.h = h

def search_ht_cl(pt_2d_a, pt_3d_list_b):
#function: Search the nearest point of pt_2d in the list of pt_3d_list, and return the nearest pt_3d
#input 1: a pt_2d, contains the x, y coordinates input 2: a pt_3d list, contains all points with x, y, h coordinates
    min_pt = pt_3d(0, 0, 0)
    min_dis = 9999
    for b in pt_3d_list_b:
        distance2 = (pt_2d_a.x - b.x)**2 + (pt_2d_a.y - b.y)**2
        if distance2 < min_dis:
            min_dis = distance2
            min_pt.x = b.x
            min_pt.y = b.y
            min_pt.h = b.h
    return min_pt, math.sqrt(min_dis)

def search_ht_ary(pt_2d_x, pt_2d_y, pt_3d_list_x, pt_3d_list_y, pt_3d_list_h):
#function: Search the nearest point of (pt_2d_x, pt_2d_y) in points where pt_3d_list_x, pt_3d_list_y, pt_3d_list_h are x, y ,h respectively
    num = len(pt_3d_list_x)
    if not ( num == len(pt_3d_list_y)):
        if not (len(pt_3d_list_y) == len(pt_3d_list_h)):
            return
        return

    min_dis = 9999
    min_pt_x, min_pt_y, min_pt_h = 0.0, 0.0, 0.0

    for i in range(0, num):
        distance2 = ( pt_2d_x- pt_3d_list_x[i])**2 + (pt_2d_y - pt_3d_list_y[i])**2
        if distance2 < min_dis:
            min_dis = distance2
            min_pt_x = pt_3d_list_x[i]
            min_pt_y = pt_3d_list_y[i]
            min_pt_h = pt_3d_list_h[i]
    return min_pt_x, min_pt_y, min_pt_h, math.sqrt(min_dis)

def search_ht_aryRGB(pt_2d_x, pt_2d_y, pt_3d_list_x, pt_3d_list_y, pt_3d_list_h, pt_3d_rlist, pt_3d_glist, pt_3d_blist):
    # function: In addition to x,y,h,distance, it also returns rgb.
    num = len(pt_3d_list_x)
    if not (num == len(pt_3d_list_y)):
        if not (len(pt_3d_list_y) == len(pt_3d_list_h)):
            return
        return

    min_dis = 9999
    min_pt_x, min_pt_y, min_pt_h, min_pt_r, min_pt_g, min_pt_b = 0.0, 0.0, 0.0, 0, 0, 0

    for i in range(0, num):
        distance2 = (pt_2d_x - pt_3d_list_x[i]) ** 2 + (pt_2d_y - pt_3d_list_y[i]) ** 2
        if distance2 < min_dis:
            min_dis = distance2
            min_pt_x = pt_3d_list_x[i]
            min_pt_y = pt_3d_list_y[i]
            min_pt_h = pt_3d_list_h[i]
            min_pt_r = pt_3d_rlist[i]
            min_pt_g = pt_3d_glist[i]
            min_pt_b = pt_3d_blist[i]
    return min_pt_x, min_pt_y, min_pt_h, math.sqrt(min_dis), min_pt_r, min_pt_g, min_pt_b
#---------------------------------------------------------------------------------------------------------------------#
# classes added for point_clouds_processing_v1.1
# this script calculates height metrics. it returns a list of metrics names and values of height metrics
import numpy as np

def height_metrics(h_list):
    if not h_list:
        return -9999
    #metrics_names = ['h_max','h95','h90','h85', 'h80', 'h75',
     #                'h70', 'h65', 'h60', 'h55', 'h_median', 'h_mean']
    h_max = np.max(h_list)
    h95 = np.percentile(h_list, 95)
    h90 = np.percentile(h_list, 90)
    h85 = np.percentile(h_list, 85)
    h80 = np.percentile(h_list, 80)
    h75 = np.percentile(h_list, 75)
    h70 = np.percentile(h_list, 70)
    h65 = np.percentile(h_list, 65)
    h60 = np.percentile(h_list, 60)
    h55 = np.percentile(h_list, 55)
    h_median = np.median(h_list)
    h_mean = np.mean(h_list)
    metrics_array = [h_max, h95, h90, h85, h80, h75, h70, h65, h60, h55, h_median, h_mean]
    return metrics_array

#----------------------------------------------------------------------------#
# this script gets height list within the buffer
def hlist_inBuffer(pt_2d_x, pt_2d_y, pt_3d_list_x, pt_3d_list_y, pt_3d_list_h, bufferDistance):
    num = len(pt_3d_list_x)
    if not (num == len(pt_3d_list_y)):
        if not (len(pt_3d_list_y) == len(pt_3d_list_h)):
            return
        return

    hlist = []
    for i in range(0, num):
        dist = (pt_2d_x - pt_3d_list_x[i]) ** 2 + (pt_2d_y - pt_3d_list_y[i]) ** 2
        if math.sqrt(dist) < bufferDistance:
            hlist.append(pt_3d_list_h[i])
    return hlist


