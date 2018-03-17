import BERA_function19 as bf19
import time

def UU_create_circle_buffer():
    bf19.crcl_shp(rtk_folder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_UAV\clip\rtk',
                  shp_folder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_UAV\clip\shp',
                  bufferDistance=3.0)

def UU_clip():
    bf19.clip_las_bat(lastool_path=r'E:\LAStools\LAStools\LAStools\bin\lasclip.exe',
                      las_ifolder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_UAV\clip\las',
                      shp_folder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_UAV\clip\shp',
                      las_ofolder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_UAV\cmp_data\las',
                      version='1.2',
                      bat_file_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_UAV\clip\clip.bat')

def UL_create_circle_buffer():
    bf19.crcl_shp(rtk_folder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_LDA\clip_lidar\rtk',
                  shp_folder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_LDA\clip_lidar\shp',
                  bufferDistance=3.0)

def UL_clip():
    bf19.clip_las_bat(lastool_path=r'E:\LAStools\LAStools\LAStools\bin\lasclip.exe',
                      las_ifolder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_LDA\clip_lidar\lidar',
                      shp_folder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_LDA\clip_lidar\shp',
                      las_ofolder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_LDA\cmp_data\lidar',
                      version='1.2',
                      bat_file_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_LDA\clip_lidar\clip_lidar.bat')

def cmp_UU(top, bottom):
    bf19.CmpVegHtCyl_UAV_UAV(las_ifolder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_UAV\cmp_data\las',
                    txt_ifolder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_UAV\cmp_data\txt',
                    txt_ofolder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_UAV\cmp_output',
                    top_radius=top, bottom_radius=bottom)

def cmp_UL(top, bottom):
    bf19.CmpVeghtCyl_UAV_LDA(Ulas_ifolder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_UAV\cmp_data\las',
                        Llas_ifolder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_LDA\cmp_data\lidar',
                        txt_ifolder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_LDA\cmp_data\txt',
                        txt_ofolder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_LDA\cmp_output',
                        top_radius=top, bottom_radius=bottom)

def cmp_UR(top):
    bf19.CmpVegHtCyl_UAV_RTK(las_ifolder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_UAV\cmp_data\las',
                    txt_ifolder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_UAV\cmp_data\txt',
                    txt_ofolder_path=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_RTK\cmp_output',
                    top_radius=top)

def cmp2smy3():
    bf19.cmp2smy3(UAV_RTK_folder=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_RTK\cmp_output',
                  UAV_UAV_folder=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_UAV\cmp_output',
                  UAV_LDA_folder=r'E:\BERA\Point_clouds_processing_v2.0\CmpVegHtCyl_UAV_LDA\cmp_output',
                  output_folder_path=r'E:\BERA\Point_clouds_processing_v2.0\Cmp2smy3\smy3_output')

def remove_smy3():
    bf19.rmv_smy3(smy3_folder_path=r'E:\BERA\Point_clouds_processing_v2.0\Cmp2smy3\smy3_output',
                  output_folder_path=r'E:\BERA\Point_clouds_processing_v2.0\rmv_smy3',
                  height_limit = 3.0)

def Smy3errors(top, bottom):
    sym3_folder = r'E:\BERA\Point_clouds_processing_v2.0\rmv_smy3'
    output_file = r'E:\BERA\Point_clouds_processing_v2.0\Smy3errors\errors_' + str(top) + '_' + str(bottom) + '.csv'
    bf19.sym3errors(sym3_folder,output_file)


def plot():
    bf19.plot_4p(txt_folder_path=r'E:\BERA\Point_clouds_processing_v2.0\rmv_smy3',
                 output_folder_path=r'E:\BERA\Point_clouds_processing_v2.0\Plots')

def meanHeight():
    bf19.VegMeanHeight(sym3_folder_path=r'E:\BERA\Point_clouds_processing_v2.0\rmv_smy3',
                       output_file_path=r'E:\BERA\Point_clouds_processing_v2.0\MeanHeight\VegMeanHeight.csv')

def meanHeightError(top, bottom):
    txt_file = r'E:\BERA\Point_clouds_processing_v2.0\MeanHeight\VegMeanHeight.csv'
    output_file = r'E:\BERA\Point_clouds_processing_v2.0\MeanHeightErrors\MeanHeightErrors_' + str(top) + '_' + str(bottom) + '.csv'
    bf19.VegMeanHeightErrors(txt_file, output_file)

def VegCover():
    bf19.VegCover(smy3_folder_path=r'E:\BERA\Point_clouds_processing_v2.0\rmv_smy3',
                  output_folder_path=r'E:\BERA\Point_clouds_processing_v2.0\VegCover')

def VegCoverErrors():
    bf19.VegCoverErrors(vc_file_path=r'E:\BERA\Point_clouds_processing_v2.0\VegCover\UAV_UAV.csv',
                        output_file_path=r'E:\BERA\Point_clouds_processing_v2.0\VegCoverErrors\VegCoverErrorsUU.csv')

def main():
    top = 0.2
    bottom_list = [0.25, 0.5, 0.75, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3]
    for bottom in bottom_list:
        cmp_UU(top, bottom)
        cmp_UL(top, bottom)
        cmp_UR(top)
        cmp2smy3()
        remove_smy3()
        Smy3errors(top, bottom)
        meanHeight()
        meanHeightError(top, bottom)

start_time = time.time()
main()
print time.time() - start_time







