import numpy as np
from math import *
import scipy.io as scio
import cv2
"""
# Class ReadDemPar
# 功能:用于读取并存储DEM头文件,文件格式见说明书
# author:csuhuacan@163.com,2022.1.14
"""
class ReadDemPar:
        def __init__(self,DEM_par):
            with open(DEM_par) as f:
                for line in f:
                    if 'title' in line:
                        s = line.split()
                        self.title = s[1]

                    if 'DEM_projection' in line:
                        s = line.split()
                        self.DEM_projection = s[1]

                    if 'data_format' in line:
                        s = line.split()
                        self.data_format = s[1]

                    if 'DEM_hgt_offset' in line:
                        s = line.split()
                        self.DEM_hgt_offset = float(s[1])

                    if 'DEM_scale' in line:
                        s = line.split()
                        self.DEM_scale = float(s[1])

                    if 'width' in line:
                        s = line.split()
                        self.width = int(s[1])

                    if 'nlines' in line:
                        s = line.split()
                        self.nlines = int(s[1])

                    if 'corner_lat' in line:
                        s = line.split()
                        self.corner_lat = float(s[1])

                    if 'corner_lon' in line:
                        s = line.split()
                        self.corner_lon = float(s[1])

                    if 'post_lat' in line:
                        s = line.split()
                        self.post_lat = float(s[1])

                    if 'post_lon' in line:
                        s = line.split()
                        self.post_lon = float(s[1])

                    if 'ellipsoid_name' in line:
                        s = line.split()
                        self.ellipsoid_name = s[1]
                    
                    if 'ellipsoid_ra' in line:
                        s = line.split()
                        self.ellipsoid_ra = float(s[1])
                    
                    if 'ellipsoid_reciprocal_flattening' in line:
                        s = line.split()
                        self.ellipsoid_reciprocal_flattening = float(s[1])

"""
# class ReadSlcPar
# 功能:用于读取并存储slc SAR数据的头文件,文件格式见说明书
# author:csuhuacan@163.com,2022.1.14
"""
class ReadSlcPar:
        def __init__(self,SLC_par):
            with open(SLC_par) as f:
                i = 0
                for line in f:
                    if 'title' in line:
                        s = line.split()
                        self.title = s[1]

                    if 'sensor' in line:
                        s = line.split()
                        self.sensor = s[1]

                    if 'date' in line:
                        s = line.split()
                        self.date = s[1:5]

                    if 'start_time' in line:
                        s = line.split()
                        self.start_time = float(s[1])

                    if 'center_time' in line:
                        s = line.split()
                        self.center_time = float(s[1])

                    if 'end_time' in line:
                        s = line.split()
                        self.end_time = float(s[1])

                    if 'azimuth_line_time' in line:
                        s = line.split()
                        self.azimuth_line_time = float(s[1])

                    if 'line_header_size' in line:
                        s = line.split()
                        self.line_header_size = int(s[1])

                    if 'range_samples' in line:
                        s = line.split()
                        self.range_samples = int(s[1])

                    if 'azimuth_lines' in line:
                        s = line.split()
                        self.azimuth_lines = int(s[1])

                    if 'range_looks' in line:
                        s = line.split()
                        self.range_looks = int(s[1])

                    if 'azimuth_looks' in line:
                        s = line.split()
                        self.azimuth_looks = int(s[1])   

                    if 'image_format' in line:
                        s = line.split()
                        self.image_format = s[1]   

                    if 'image_geometry' in line:
                        s = line.split()
                        self.image_geometry = s[1]  

                    if 'range_scale_factor' in line:
                        s = line.split()
                        self.range_scale_factor = float(s[1])
                    
                    if 'azimuth_scale_factor' in line:
                        s = line.split()
                        self.azimuth_scale_factor = float(s[1])

                    if 'center_latitude' in line:
                        s = line.split()
                        self.center_latitude = float(s[1])

                    if 'center_longitude' in line:
                        s = line.split()
                        self.center_longitude = float(s[1])

                    if 'heading' in line:
                        s = line.split()
                        self.heading = float(s[1])


                    if 'range_pixel_spacing' in line:
                        s = line.split()
                        self.range_pixel_spacing = float(s[1])

                    if 'azimuth_pixel_spacing' in line:
                        s = line.split()
                        self.azimuth_pixel_spacing = float(s[1])

                    if 'near_range_slc' in line:
                        s = line.split()
                        self.near_range_slc = float(s[1])

                    if 'center_range_slc' in line:
                        s = line.split()
                        self.center_range_slc = float(s[1])

                    if 'far_range_slc' in line:
                        s = line.split()
                        self.far_range_slc = float(s[1])

                    if 'first_slant_range_polynomial' in line:
                        s = line.split()
                        dpo = np.zeros((6),dtype='float64')
                        dpo[0] = float(s[1])
                        dpo[1] = float(s[2])
                        dpo[2] = float(s[3])
                        dpo[3] = float(s[4])
                        dpo[4] = float(s[5])
                        dpo[5] = float(s[6])
                        self.fsrp = dpo
                    
                    if 'center_slant_range_polynomial' in line:
                        s = line.split()
                        dpo = np.zeros((6),dtype='float64')
                        dpo[0] = float(s[1])
                        dpo[1] = float(s[2])
                        dpo[2] = float(s[3])
                        dpo[3] = float(s[4])
                        dpo[4] = float(s[5])
                        dpo[5] = float(s[6])
                        self.csrp = dpo
                    
                    if 'last_slant_range_polynomial' in line:
                        s = line.split()
                        dpo = np.zeros((6),dtype='float64')
                        dpo[0] = float(s[1])
                        dpo[1] = float(s[2])
                        dpo[2] = float(s[3])
                        dpo[3] = float(s[4])
                        dpo[4] = float(s[5])
                        dpo[5] = float(s[6])
                        self.lsrp = dpo

                    if 'incidence_angle' in line:
                        s = line.split()
                        self.incidence_angle = float(s[1])

                    if 'azimuth_deskew' in line:
                        s = line.split()
                        self.azimuth_deskew = s[1]
                    
                    if 'azimuth_angle' in line:
                        s = line.split()
                        self.azimuth_angle = float(s[1])
                    
                    if 'radar_frequency' in line:
                        s = line.split()
                        self.radar_frequency = float(s[1])

                    if 'adc_sampling_rate' in line:
                        s = line.split()
                        self.adc_sampling_rate = float(s[1])

                    if 'chirp_bandwidth' in line:
                        s = line.split()
                        self.chirp_bandwidth = float(s[1])

                    if 'prf' in line:
                        s = line.split()
                        self.prf = float(s[1])

                    if 'azimuth_proc_bandwidth' in line:
                        s = line.split()
                        self.azimuth_proc_bandwidth = float(s[1])


                    if 'doppler_polynomial' in line:
                        s = line.split()
                        dpo = np.zeros((4),dtype='float64')
                        dpo[0] = float(s[1])
                        dpo[1] = float(s[2])
                        dpo[2] = float(s[3])
                        dpo[3] = float(s[4])
                        self.dp = dpo

                    if 'doppler_poly_dot' in line:
                        s = line.split()
                        dpo = np.zeros((4),dtype='float64')
                        dpo[0] = float(s[1])
                        dpo[1] = float(s[2])
                        dpo[2] = float(s[3])
                        dpo[3] = float(s[4])
                        self.dpd = dpo

                    if 'doppler_poly_ddot' in line:
                        s = line.split()
                        dpo = np.zeros((4),dtype='float64')
                        dpo[0] = float(s[1])
                        dpo[1] = float(s[2])
                        dpo[2] = float(s[3])
                        dpo[3] = float(s[4])
                        self.dpdd = dpo

                    if 'receiver_gain' in line:
                        s = line.split()
                        self.receiver_gain = float(s[1])

                    if 'calibration_gain' in line:
                        s = line.split()
                        self.calibration_gain = float(s[1])
                    
                    if 'sar_to_earth_center' in line:
                        s = line.split()
                        self.sar_to_earth_center = float(s[1])

                    if 'earth_radius_below_sensor' in line:
                        s = line.split()
                        self.earth_radius_below_sensor = float(s[1])

                    if 'earth_semi_major_axis' in line:
                        s = line.split()
                        self.earth_semi_major_axis = float(s[1])

                    if 'earth_semi_minor_axis' in line:
                        s = line.split()
                        self.earth_semi_minor_axis = float(s[1])

                    if 'number_of_state_vectors' in line:
                        s = line.split()
                        self.number_of_state_vectors = int(s[1])
                        # 创建轨道数组和速度数组
                        position = np.zeros((int(s[1]),3),dtype='float64')
                        velocity = np.zeros((int(s[1]),3),dtype='float64')
                    
                    if 'time_of_first_state_vector' in line:
                        s = line.split()
                        self.time_of_first_state_vector = float(s[1])

                    if 'state_vector_interval' in line:
                        s = line.split()
                        self.state_vector_interval = float(s[1])

                    # 读取轨道数据                    
                    if 'state_vector_position' in line:
                        s = line.split()
                        position[i,0] = s[1]
                        position[i,1] = s[2]
                        position[i,2] = s[3]
                    
                    if 'state_vector_velocity' in line:
                        s = line.split()
                        velocity[i,0] = s[1]
                        velocity[i,1] = s[2]
                        velocity[i,2] = s[3]
                        i = i+1
            self.position = position
            self.velocity = velocity

"""
# class ReadOffPar
# 功能:用于读取并存储off_par数据的头文件,文件格式见说明书
# author:csuhuacan@163.com,2022.1.14
"""
class ReadOffPar:
        def __init__(self,OFF_par):
            with open(OFF_par) as f:
                i = 0
                for line in f:
                    if 'title' in line:
                        s = line.split()
                        self.title = s[1]

                    if 'initial_range_offset' in line:
                        s = line.split()
                        self.initr = float(s[1])

                    if 'initial_azimuth_offset' in line:
                        s = line.split()
                        self.initaz = float(s[1])

                    if 'slc1_starting_range_pixel' in line:
                        s = line.split()
                        self.nofstr = float(s[1])

                    if 'number_of_slc_range_pixels' in line:
                        s = line.split()
                        self.npr = float(s[1])

                    if 'offset_estimation_starting_range' in line:
                        s = line.split()
                        self.rstr = float(s[1])

                    if 'offset_estimation_ending_range' in line:
                        s = line.split()
                        self.rend = float(s[1])
                    if 'offset_estimation_range_samples' in line:
                        s = line.split()
                        self.nr = float(s[1])

                    if 'offset_estimation_range_spacing' in line:
                        s = line.split()
                        self.rsp = float(s[1])

                    if 'offset_estimation_starting_azimuth' in line:
                        s = line.split()
                        self.azstr = float(s[1])

                    if 'offset_estimation_ending_azimuth' in line:
                        s = line.split()
                        self.azend = float(s[1])

                    if 'offset_estimation_azimuth_samples' in line:
                        s = line.split()
                        self.naz = float(s[1])

                    if 'offset_estimation_azimuth_spacing' in line:
                        s = line.split()
                        self.azsp = float(s[1])

                    if 'offset_estimation_window_width' in line:
                        s = line.split()
                        self.rwin = float(s[1])

                    if 'offset_estimation_window_height' in line:
                        s = line.split()
                        self.azwin = float(s[1])

                    if 'offset_estimation_threshhold' in line:
                        s = line.split()
                        self.thres = float(s[1])

                    if 'range_offset_polynomial' in line:
                        s = line.split()
                        dpo = np.zeros((6),dtype='float64')
                        dpo[0] = float(s[1])
                        dpo[1] = float(s[2])
                        dpo[2] = float(s[3])
                        dpo[3] = float(s[4])
                        dpo[4] = float(s[5])
                        dpo[5] = float(s[6])
                        self.rop = dpo
                    
                    if 'azimuth_offset_polynomial' in line:
                        s = line.split()
                        dpo = np.zeros((6),dtype='float64')
                        dpo[0] = float(s[1])
                        dpo[1] = float(s[2])
                        dpo[2] = float(s[3])
                        dpo[3] = float(s[4])
                        dpo[4] = float(s[5])
                        dpo[5] = float(s[6])
                        self.aop = dpo

                    if 'slc1_starting_azimuth_line' in line:
                        s = line.split()
                        self.lbegin = float(s[1])
                    
                    if 'interferogram_azimuth_lines' in line:
                        s = line.split()
                        self.nls = float(s[1])

                    if 'interferogram_width' in line:
                        s = line.split()
                        self.nr_int = float(s[1])

                    if 'first_nonzero_range_pixel' in line:
                        s = line.split()
                        self.nrb = float(s[1])

                    if 'number_of_nonzero_range_pixels' in line:
                        s = line.split()
                        self.nrps = float(s[1])
                    
                    if 'interferogram_range_looks' in line:
                        s = line.split()
                        self.rlks = float(s[1])

                    if 'interferogram_azimuth_looks' in line:
                        s = line.split()
                        self.azlks = float(s[1])

                    if 'interferogram_range_pixel_spacing' in line:
                        s = line.split()
                        self.rps_int = float(s[1])

                    if 'interferogram_azimuth_pixel_spacing' in line:
                        s = line.split()
                        self.azps_int = float(s[1])
                    
                    if 'resampled_range_pixel_spacing' in line:
                        s = line.split()
                        self.rps_res = float(s[1])
                    
                    if 'resampled_azimuth_pixel_spacing' in line:
                        s = line.split()
                        self.azps_res = float(s[1])

                    if 'resampled_starting_ground_range' in line:
                        s = line.split()
                        self.grg_start = float(s[1])

                    if 'resampled_pixels_per_line' in line:
                        s = line.split()
                        self.ngrg = float(s[1])

                    if 'resampled_number_of_lines' in line:
                        s = line.split()
                        self.ngaz = float(s[1])       

"""
# class WGS84
# 功能:WGS84空间直角坐标X,Y,Z
# author:csuhuacan@163.com,2022.1.14    
"""     
class WGS84:
        def __init__(self):
            self.x = ''
            self.y = ''
            self.z = ''

""" 
#class Ctcn
#功能:TCN坐标系
#author:csuhuacan@163.com,2022.1.14   
# """       
class Ctcn:
        def __init__(self):
            self.t = ''
            self.c = ''
            self.n = ''

""" 
#class Vec Vec2
#功能:空间直角坐标X,Y,Z
#author:csuhuacan@163.com:2022.1.14    
"""      
class Vec:
        def __init__(self,a):
            self.x = np.zeros((a),dtype = 'float64')
            self.y = np.zeros((a),dtype = 'float64')
            self.z = np.zeros((a),dtype = 'float64')

class Vec2:
        def __init__(self,a):
            self.x = np.zeros((2,a),dtype = 'float64')
            self.y = np.zeros((2,a),dtype = 'float64')
            self.z = np.zeros((2,a),dtype = 'float64')          

""" 
#class LLH
#功能:大地经纬度坐标L,L,H
#author:csuhuacan@163.com,2022.1.14  
"""        
class LLH:
        def __init__(self):
            self.lon = ''
            self.lat = ''
            self.alt = ''

def freadbk(filename,datatype,row,flag):
    """
    Function:
      读取二进制文件
    Author:
      csuhuacan@163.com,2022.11.09 
    Parameters:
      filename - 二进制文件路径
      datatype - 输入的数据类型,gamma对应的为:'>f',python中fwrite输入输出对应的是:'<f'
      row - 数据的行号
      flag - 0:实数数据,1:复数数据
    Returns: 
      二维矩阵 
    Raises:
      KeyError - raises an exception 
    """
    pass

    # 读取实数文件
    if flag == 0 :  
        data = np.fromfile(filename, dtype=datatype)
        count = data.shape[0]
        b=(count % row == 0)
        if b!=1:
            print('输入的行号有误')
            return
        else:
            output = data.reshape(row, data.shape[0] // row)
    
    # 读取复数文件
    elif flag == 1 :
        data = np.fromfile(filename, dtype=datatype)
        count = data.shape[0]
        b=(count % row == 0)
        if b!=1:
            print('输入的行号有误')
            return    
        else:
            a = np.array(data[0:count:2])
            real = a.reshape(row, a.shape[0] // row)

            b = np.array(data[1:count:2])
            imag = b.reshape(row, a.shape[0] // row)
            output = real + 1j*imag

    return output  

def p_interp(ds,flag,t):
    """
    Function:
      利用三次多项式对卫星轨道位置进行拟合
    Author:
      csuhuacan@163.com,2022.1.14 
    Parameters:
      ds - slc.par文件参数
      flag - 0:计算整个影像的每个像素的位置,1:计算任意时间点的位置,需要输入时间t
    Returns: 
      flag=1 时xyz返回矩阵, flag=0 时返回xyz单点坐标
    Raises:
      KeyError - raises an exception 
    """
    pass
    if flag == 0:#计算整个影像的每个像素的位置
        p_x = ds.position[:,0]
        p_y = ds.position[:,1]
        p_z = ds.position[:,2]

        state_time = np.zeros((ds.number_of_state_vectors),dtype='float64')
        sar_t = np.arange(ds.start_time,ds.end_time,ds.azimuth_line_time)
        position1 = np.zeros((ds.azimuth_lines,ds.range_samples),dtype='float64')
        position2 = np.zeros((ds.azimuth_lines,ds.range_samples),dtype='float64')
        position3 = np.zeros((ds.azimuth_lines,ds.range_samples),dtype='float64')

        for i in range(ds.number_of_state_vectors):
            state_time[i] = ds.time_of_first_state_vector + i*ds.state_vector_interval

        xa=np.polyfit(state_time,p_x,3) # 用3次多项式拟合x，y数组
        xb=np.poly1d(xa)                # 拟合完之后用这个函数来生成多项式对象

        ya=np.polyfit(state_time,p_y,3) # 用3次多项式拟合x，y数组
        yb=np.poly1d(ya)                # 拟合完之后用这个函数来生成多项式对象

        za=np.polyfit(state_time,p_z,3) # 用3次多项式拟合x，y数组
        zb=np.poly1d(za)                # 拟合完之后用这个函数来生成多项式对象

        position1[:,0] = xb(sar_t)
        position2[:,0] = yb(sar_t)
        position3[:,0] = zb(sar_t)

        for i in range(ds.azimuth_lines):
            position1[i,:] = position1[i,0] 
            position2[i,:] = position2[i,0] 
            position3[i,:] = position3[i,0] 

        a = WGS84()
        a.x = position1
        a.y = position2
        a.z = position3
        return a
    else: #计算时刻t的卫星位置
        p_x = ds.position[:,0]
        p_y = ds.position[:,1]
        p_z = ds.position[:,2]
        state_time = np.zeros((ds.number_of_state_vectors),dtype='float64')
        for i in range(ds.number_of_state_vectors):
            state_time[i] = ds.time_of_first_state_vector + i*ds.state_vector_interval

        xa=np.polyfit(state_time,p_x,3)#用2次多项式拟合x，y数组
        xb=np.poly1d(xa)#拟合完之后用这个函数来生成多项式对象

        ya=np.polyfit(state_time,p_y,3)#用2次多项式拟合x，y数组
        yb=np.poly1d(ya)#拟合完之后用这个函数来生成多项式对象

        za=np.polyfit(state_time,p_z,3)#用2次多项式拟合x，y数组
        zb=np.poly1d(za)#拟合完之后用这个函数来生成多项式对象
        p = WGS84()
        p.x = xb(t)
        p.y = yb(t)
        p.z = zb(t)
        return p
        
def v_interp(ds,flag,t):
    """
    Function:
      利用三次多项式对卫星轨道速度进行拟合
    Author:
      csuhuacan@163.com,2022.1.14 
    Parameters:
      ds - slc.par文件参数
      flag - 0:计算整个影像的每个像素的速度,1:计算任意时间点的速度,需要输入时间t
    Returns: 
      flag=1 时Vxyz返回矩阵, flag=0 时返回Vxyz单点坐标
    Raises:
      KeyError - raises an exception 
    """
    pass
    if flag == 0:#计算整个影像的每个像素的速度
        v_x = ds.velocity[:,0]
        v_y = ds.velocity[:,1]
        v_z = ds.velocity[:,2]

        state_time = np.zeros((ds.number_of_state_vectors),dtype='float64')
        sar_t = np.arange(ds.start_time,ds.end_time,ds.azimuth_line_time)
        velocity = np.zeros((ds.azimuth_lines,3),dtype='float64')

        for i in range(ds.number_of_state_vectors):
            state_time[i] = ds.time_of_first_state_vector + i*ds.state_vector_interval

        xa=np.polyfit(state_time,v_x,3)#用2次多项式拟合x，y数组
        xb=np.poly1d(xa)#拟合完之后用这个函数来生成多项式对象

        ya=np.polyfit(state_time,v_y,3)#用2次多项式拟合x，y数组
        yb=np.poly1d(ya)#拟合完之后用这个函数来生成多项式对象

        za=np.polyfit(state_time,v_z,3)#用2次多项式拟合x，y数组
        zb=np.poly1d(za)#拟合完之后用这个函数来生成多项式对象

        velocity[:,0] = xb(sar_t)
        velocity[:,1] = yb(sar_t)
        velocity[:,2] = zb(sar_t)

        return velocity
    else:
        v_x = ds.velocity[:,0]
        v_y = ds.velocity[:,1]
        v_z = ds.velocity[:,2]
        state_time = np.zeros((ds.number_of_state_vectors),dtype='float64')
        for i in range(ds.number_of_state_vectors):
            state_time[i] = ds.time_of_first_state_vector + i*ds.state_vector_interval

        xa=np.polyfit(state_time,v_x,3) # 用2次多项式拟合x，y数组
        xb=np.poly1d(xa)                # 拟合完之后用这个函数来生成多项式对象

        ya=np.polyfit(state_time,v_y,3) # 用2次多项式拟合x，y数组
        yb=np.poly1d(ya)                # 拟合完之后用这个函数来生成多项式对象

        za=np.polyfit(state_time,v_z,3) # 用2次多项式拟合x，y数组
        zb=np.poly1d(za)                # 拟合完之后用这个函数来生成多项式对象

        v = WGS84()
        v.x = xb(t)
        v.y = yb(t)
        v.z = zb(t)
        return v
 
    
def R_interp(ds):
    """
    Function:
      计算每个像素的斜距
    Author:
      csuhuacan@163.com,2022.1.14 
    Parameters:
      ds - slc.par文件参数
    Returns: 
      二维矩阵:每个像素的斜距
    Raises:
      KeyError - raises an exception 
    """
    pass
    dR = ds.range_pixel_spacing
    R = np.zeros((ds.azimuth_lines,ds.range_samples),dtype='float64')
    R[:,0] = ds.near_range_slc
    for i in range(1,ds.range_samples):
        R[:,i] = R[:,0] + i * dR
    return R

def xyz_ll(x,y,z):
    """
    Function:
      空间直角坐标转大地经纬度坐标xyz--->L，B
    Author:
      csuhuacan@163.com,2022.1.14 
    Parameters:
      x,y,z:空间直角坐标的位置
    Returns: 
      包含经纬度和海拔的结构体
    Raises:
      KeyError - raises an exception 
    """
    pass
    a = 6378137
    e2 = 0.0066943799013
    b = a * sqrt(1.0-e2)
    p81 = 57.29577951308232

    e1 = e2 * a*a / (b*b)

    p = sqrt(x*x + y*y)
    theta = atan((z * a) / (p * b))
    st = sin(theta)
    ct = cos(theta)
    phi = atan((z+e1*b*st*st*st)/(p-e2*a*ct*ct*ct))
    lam = atan2(y,x)

    L = LLH()
    L.lat = phi*p81
    L.lon = lam*p81
    nu = a/sqrt(1.0 - e2*sin(phi)*sin(phi))
    L.alt = p/cos(phi) - nu
    return L

def ll_xyz(lat,lon,H):
    """
    Function:
      大地经纬度坐标转空间直角坐标L,B,H--->xyz
    Author:
      csuhuacan@163.com,2022.1.14 
    Parameters:
      lat,lon,H:经纬度和海拔
    Returns: 
      x,y,z:空间直角坐标的位置
    Raises:
      KeyError - raises an exception 
    """
    pass
    a = 6378137
    b = 6356752.3142
    e2 = 0.0066943799013
    p18 = 0.017453292519943295
    lat = lat * p18
    lon = lon * p18
    N = a/np.sqrt(1-e2*np.sin(lat)*np.sin(lat))

    W = WGS84()
    W.x = (N + H)*np.cos(lat)*np.cos(lon)
    W.y = (N + H)*np.cos(lat)*np.sin(lon)
    W.z = (N*(1-e2)+H)*np.sin(lat)
    return W

def ij_ll(dp):
    """
    Function:
      计算DEM每个像素的经纬度
    Author:
      csuhuacan@163.com,2022.1.14 
    Parameters:
      dp:DEM头文件
    Returns: 
      二维矩阵:DEM每个像素的经纬度
    Raises:
      KeyError - raises an exception 
    """
    pass
    lat = np.zeros((dp.nlines,dp.width),dtype='float64')
    lon = np.zeros((dp.nlines,dp.width),dtype='float64')
    for i in range(dp.nlines):
        lat[i,:] = dp.corner_lat + dp.post_lat * i
    for j in range(dp.width):
        lon[:,j] = dp.corner_lon + dp.post_lon * j

    return lat,lon 

def ll_ij(dp,lat,lon):
    """
    Function:
      经纬度对应的dem的行列号
    Author:
      csuhuacan@163.com,2022.1.14 
    Parameters:
      dp:DEM头文件,lat,lon:经纬度
    Returns: 
      行列号,vflag=0:超出了当前DEM的边界,vflag=1:正常使用
    Raises:
      KeyError - raises an exception 
    """
    pass
    irow = ((lat - dp.corner_lat) / dp.post_lat)
    jcol = ((lon - dp.corner_lon) / dp.post_lon)

    if((irow < 0) or (irow > dp.nlines - 1) or (jcol < 0) or (jcol > dp.width - 1) ):
        vflag = 0
    else:
        vflag = 1

    return irow,jcol,vflag 

def cross(a, b):
    c = WGS84()
    c.x = a.y * b.z - a.z * b.y
    c.y = a.z * b.x - a.x * b.z
    c.z = a.x * b.y - a.y * b.x
    return c

def unit(a):
    m = sqrt(a.x * a.x + a.y * a.y + a.z * a.z)
    n = WGS84()
    n.x = a.x / m
    n.y = a.y / m
    n.z = a.z / m
    return n.x,n.y,n.z

def DEM_read(dp, demr4, line, col):
    """
    Function:
      按照输入的经纬度和行列号准确的读取DEM
    Author:
      csuhuacan@163.com,2022.1.14 
    Parameters:
      dp:DEM头文件,
      demr4: DEM,
      line、col:行列号
    Returns: 
      DEM
    Raises:
      KeyError - raises an exception 
    """
    pass
    if ((line >= dp.nlines) or (col >= dp.width) or (line < 0) or (col < 0)):
        hgt = 0.0
        return hgt
    hgt = demr4[line][col] * dp.DEM_scale - dp.DEM_hgt_offset
    return hgt

def DEM_interp(dp, dem, north, east):
    """
    Function:
      对DEM进行插值,插值得出成像平面
    Author:
      csuhuacan@163.com,2022.1.14 
    Parameters:
      dp:DEM头文件,
      demr: DEM,
      north, east:经纬度
    Returns: 
      当前经纬度对应的DEM
    Raises:
      KeyError - raises an exception 
    """
    pass

    s4 = np.zeros((4),dtype='float64')
    A = 1.0
    B = 1.0 + (A - 1.0) / 2.0
    
    n_lin = (north - dp.corner_lat) / dp.post_lat
    e_col = (east - dp.corner_lon) / dp.post_lon

    n_start = floor(n_lin)
    e_start = floor(e_col)
    n_frac = n_lin - n_start
    e_frac = e_col - e_start
    for i1 in range(-1,3):
        d0 = DEM_read(dp, dem, n_start + i1, e_start - 1)
        if d0 == 0.0:
            hgt = 0.0

        d1 = DEM_read(dp, dem, n_start + i1, e_start)
        if d1 == 0.0:
            hgt = 0.0

        d2 = DEM_read(dp, dem, n_start + i1, e_start + 1)
        if d2 == 0.0:
            hgt = 0.0

        d3 = DEM_read(dp, dem, n_start + i1, e_start + 2)
        if d3 == 0.0:
            hgt = 0.0
 
        Abb = (d2 - d0) / 2.0
        Acc = d2 - d1 - Abb
        Bbb = (d3 - d1) / 2.0
        Bcc = d3 - d2 - Bbb
        Aweight = max(0.0, B - A * e_frac)
        Aweight = min(1.0, Aweight)
        Bweight = 1.0 - Aweight
        s4[i1+1] = Aweight * (d1 + Abb * e_frac + Acc * e_frac * e_frac) + Bweight * (d2 + Bbb * (e_frac - 1.0) + Bcc * (e_frac - 1.0) * (e_frac - 1.0))

    Abb = (s4[2] - s4[0]) / 2.
    Acc = s4[2] - s4[1] - Abb
    Bbb = (s4[3] - s4[1]) / 2.
    Bcc = s4[3] - s4[2] - Bbb
    Aweight = max(0.0, B - A * n_frac)
    Aweight = min(1.0, Aweight)
    Bweight = 1. - Aweight
    hgt = Aweight * (s4[1] + Abb * n_frac + Acc * n_frac * n_frac) + Bweight * (s4[2] + Bbb * (n_frac - 1.0) + Bcc * (n_frac - 1.0) * (n_frac - 1.0)) 
    return hgt

def dopc(ds, r):
    dr = r - ds.center_range_slc
    dr = r - ds.center_range_slc
    fd = ds.dp[0] + dr * (ds.dp[1] + dr * (ds.dp[2] + dr * ds.dp[3]))
    return fd

def vec_mat( a, b, c):
    z = np.zeros((3,3),dtype='float64')
    z[0][0] = a.x
    z[0][1] = a.y
    z[0][2] = a.z
    z[1][0] = b.x
    z[1][1] = b.y
    z[1][2] = b.z
    z[2][0] = c.x
    z[2][1] = c.y
    z[2][2] = c.z
    return z

def c_tcn(pc, vc):
    t = WGS84()
    c = WGS84()
    n = WGS84()
    tmp= WGS84()

    tmp.x = (-pc.x)
    tmp.y = (-pc.y)
    tmp.z = (-pc.z)
    n.x,n.y,n.z = unit(tmp)

    tmp = cross(n, vc)
    c.x,c.y,c.z = unit(tmp)

    t = cross(c, n)
    tcn = vec_mat(t, c, n)
    return tcn

def matvm(a, b):
    c = Ctcn()
    c.t = a[0][0] * b.x + a[0][1] * b.y + a[0][2] * b.z
    c.c = a[1][0] * b.x + a[1][1] * b.y + a[1][2] * b.z
    c.n = a[2][0] * b.x + a[2][1] * b.y + a[2][2] * b.z
    return c

def smul(sm, a):
    b = WGS84()
    b.x = a.x * sm
    b.y = a.y * sm
    b.z = a.z * sm
    return  b.x,b.y,b.z

def norm(a,i,j):
    mag = sqrt(a.x[i,j] * a.x[i,j] + a.y[i,j] * a.y[i,j] + a.z[i,j] * a.z[i,j])
    return mag

def sub(a,i,j,b,m,n):
    c = WGS84()
    if i==0 and j==0:
        c.x = a.x - b.x[m,n]
        c.y = a.y - b.y[m,n]
        c.z = a.z - b.z[m,n]      
    else:
        c.x = a.x[i,j] - b.x[m,n]
        c.y = a.y[i,j] - b.y[m,n]
        c.z = a.z[i,j] - b.z[m,n]
    
    return c
def sind(x):
    y = np.sin(x*pi/180)
    return y

def bilinear_interpolate(data, x, y):
    """Function to perform bilinear interpolation on the input array data, at
    the image coordinates given by input arguments x and y.
    
    Arguments
        data (array): 2D array containing raster data to interpolate.
        x (array): the X coordinate values at which to interpolate (in array
            indices, starting at zero).  Note that X refers to the second
            dimension of data (e.g., the columns).
        y (array): the Y coordinate values at which to interpolate (in array
            indices, starting at zero).  Note that Y refers to the first
            dimension of data (e.g., the rows).
            
    Returns:
        intdata (array): The 2D interpolated array, with same dimensions as
            x and y.

    """
    x = np.asarray(x)
    y = np.asarray(y)

    # Get lower and upper bounds for each pixel.
    x0 = np.floor(x).astype(int)
    x1 = x0 + 1
    y0 = np.floor(y).astype(int)
    y1 = y0 + 1

    # Clip the image coordinates to the size of the input data.
    x0 = np.clip(x0, 0, data.shape[1]-1);
    x1 = np.clip(x1, 0, data.shape[1]-1);
    y0 = np.clip(y0, 0, data.shape[0]-1);
    y1 = np.clip(y1, 0, data.shape[0]-1);

    data_ll = data[ y0, x0 ] # lower left corner image values
    data_ul = data[ y1, x0 ] # upper left corner image values
    data_lr = data[ y0, x1 ] # lower right corner image values
    data_ur = data[ y1, x1 ] # upper right corner image values

    w_ll = (x1-x) * (y1-y) # weight for lower left value
    w_ul = (x1-x) * (y-y0) # weight for upper left value
    w_lr = (x-x0) * (y1-y) # weight for lower right value
    w_ur = (x-x0) * (y-y0) # weight for upper right value
    
    # Where the x or y coordinates are outside of the image boundaries, set one
    # of the weights to nan, so that these values are nan in the output array.
    w_ll[np.less(x,0)] = np.nan
    w_ll[np.greater(x,data.shape[1]-1)] = np.nan
    w_ll[np.less(y,0)] = np.nan
    w_ll[np.greater(y,data.shape[0]-1)] = np.nan
    
    intdata = w_ll*data_ll + w_ul*data_ul + w_lr*data_lr + w_ur*data_ur

    return intdata

def write_dem_par(dp, str_path):
    """
    Function:
      写DEM头文件
    Author:
      csuhuacan@163.com,2022.11.10
    Parameters:
      dp:包含DEM头文件的参数,是前面所定义的类,str_path:保存路径文件名称
    Returns: 
      1:写入成功,0:写入失败 
    Raises:
      KeyError - raises an exception 
    """
    pass
    try:
        with open(str_path,"w") as f:
            f.write('DEM paramter file,  CSU-PPSM\r\n')
            f.write('title:  InSAR DEM\r\n')
            f.write('DEM_projection:     %s\n'%(dp.DEM_projection))
            f.write('data_format:        %s\n'%(dp.data_format))
            f.write('DEM_hgt_offset:          %.5f\n'%(dp.DEM_hgt_offset))
            f.write('DEM_scale:               %.5f\n'%(dp.DEM_scale))
            f.write('width:                 %d\n'%(dp.width))
            f.write('nlines:                %d\n'%(dp.nlines))
            f.write('corner_lat:      %.8f    decimal degrees\n'%(dp.corner_lat))
            f.write('corner_lon:      %.8f    decimal degrees\n'%(dp.corner_lon))
            f.write('post_lat:    -%.12f    decimal degrees\n'%(dp.post_lat))
            f.write('post_lon:    %.12f    decimal degrees\n'%(dp.post_lon))

            f.write('ellipsoid_name: WGS 84\n')
            f.write('ellipsoid_ra:        %.3f   m\n'%(dp.ellipsoid_ra))
            f.write('ellipsoid_reciprocal_flattening:   %.6f\n'%(dp.ellipsoid_reciprocal_flattening))

            f.write('datum_name: WGS 1984\n')
            f.write('datum_shift_dx:              0.000   m\n')
            f.write('datum_shift_dy:              0.000   m\n')
            f.write('datum_shift_dz:              0.000   m\n')
            f.write('datum_scale_m:         0.00000e+00\n')
            f.write('datum_rotation_alpha:  0.00000e+00   arc-sec\n')
            f.write('datum_rotation_beta:   0.00000e+00   arc-sec\n')
            f.write('datum_rotation_gamma:  0.00000e+00   arc-sec\n')
            f.write('datum_country_list Global Definition, WGS84, World\n')
            f.close()
            return 1
    except:
        return 0
  
def rascc(cc, dtype, nlines,start_cc = 0, end_cc = -1,cmin = 0.1,cmax = 0.9, scale = 1,rasf = None):
    """
    Function:
      输出相干性,cc图,    特别注明:输出路径中不能包含中文字符,且必须以*.bmp,*.jpg,*.png结尾
    Author:
      csuhuacan@163.com,2022.11.10
    Parameters:
      cc:相干性文件
      dtype:数据类型,gamma对应的为:'>f',python中fwrite输入输出对应的是:'<f'
      start_cc:输出图片开始的行数,默认为第一行:0
      end_cc:输出图片结束的行数,默认为最后行:-1
      cmin:输出图片的最小值,默认为0.1
      cmax:输出图片的最大值,默认为0.9
      scale:对所输出的图片就行缩放,默认为1,大于1进行放大,小于1进行缩小
      rasf:输出图片的路径文件
    Returns: 
      1:写入成功,0:写入失败 
    Raises:
      KeyError - raises an exception 
    """
    pass
    try:
        
        fltcc = freadbk(cc,dtype,nlines,0) # 读取数据
        temp1 = fltcc[start_cc:end_cc,:] # 获取需要显示的部分
        temp2 = (temp1-cmin)/(cmax-cmin) # 对数据进行归化
        temp2[temp2 < cmin] = cmin
        temp2[temp2 > cmax] = cmax

        if scale != 1: # 对数据进行缩放
            temp3 = (temp2 * 255).astype(np.uint8)  
            colormap1 = cv2.applyColorMap(np.array(temp3), cv2.COLORMAP_JET) 
            width = int(temp2.shape[1]*scale)
            height = int(temp2.shape[0]*scale)
            colormap = cv2.resize(colormap1,(width,height))
            
            cv2.imwrite(rasf, colormap)
        else:
            image_arry = (temp2 * 255).astype(np.uint8)  
            colormap = cv2.applyColorMap(np.array(image_arry), cv2.COLORMAP_JET) 
            cv2.imwrite(rasf, colormap)
        
        return 1           
    except:
        return 0 

def raspwr(pwr, dtype, nlines,start_nlines = 0, end_nlines = -1,scale = 1,rasf = None):
    """
    Function:
      输出强度图,    特别注明:输出路径中不能包含中文字符,且必须以*.bmp,*.jpg,*.png结尾
    Author:
      csuhuacan@163.com,2022.11.10
    Parameters:
      pwr:强度文件
      dtype:数据类型,gamma对应的为:'>f',python中fwrite输入输出对应的是:'<f'
      start_nlines:输出图片开始的行数,默认为第一行:0
      end_nlines:输出图片结束的行数,默认为最后行:-1
      scale:对所输出的图片就行缩放,默认为1,大于1进行放大,小于1进行缩小
      rasf:输出图片的路径文件
    Returns: 
      1:写入成功,0:写入失败 
    Raises:
      KeyError - raises an exception 
    """
    pass
    try:
        
        mli = freadbk(pwr,dtype,nlines,0)      # 读取数据
        temp1 = mli[start_nlines:end_nlines,:] # 获取需要显示的部分
        temp1 = 10*np.log10(temp1)               # 计算分贝值
        
        h1 = np.percentile(temp1, 0.01)
        h99 = np.percentile(temp1, 99.99)
        temp2 = (temp1-h1)/(h99-h1) # 对数据进行归化

        if scale != 1: # 对数据进行缩放
            image_arry = (temp2 * 255).astype(np.uint8)  
            width = int(temp2.shape[1]*scale)
            height = int(temp2.shape[0]*scale)
            colormap = cv2.resize(image_arry,(width,height))
            
            cv2.imwrite(rasf, colormap)
        else:
            image_arry = (temp2 * 255).astype(np.uint8)  
            cv2.imwrite(rasf, image_arry)
        
        return 1           
    except:
        return 0 
    
def rasmph(cpx, dtype, nlines,start_nlines = 0, end_nlines = -1,scale = 1,rasf = None):
    """
    Function:
      显示复数数据图,主要以相位的形式就行显示,    特别注明:输出路径中不能包含中文字符,且必须以*.bmp,*.jpg,*.png结尾
    Author:
      csuhuacan@163.com,2022.11.10
    Parameters:
      cpx:复数数据文件
      dtype:数据类型,gamma对应的为:'>f',python中fwrite输入输出对应的是:'<f'
      start_nlines:输出图片开始的行数,默认为第一行:0
      end_nlines:输出图片结束的行数,默认为最后行:-1
      scale:对所输出的图片就行缩放,默认为1,大于1进行放大,小于1进行缩小
      rasf:输出图片的路径文件
    Returns: 
      1:写入成功,0:写入失败 
    Raises:
      KeyError - raises an exception 
    """
    pass
    try:
        
        cpx = freadbk(cpx,dtype,nlines,1)      # 读取数据
        temp1 = cpx[start_nlines:end_nlines,:] # 获取需要显示的部分
        temp1 = np.angle(temp1)                # 计算相位
        
        h1 = np.min(temp1)
        h99 = np.max(temp1)
        temp2 = (temp1-h1)/(h99-h1) # 对数据进行归化

        if scale != 1: # 对数据进行缩放
            temp3 = (temp2 * 255).astype(np.uint8)  
            colormap1 = cv2.applyColorMap(np.array(temp3), cv2.COLORMAP_JET) 
            width = int(temp2.shape[1]*scale)
            height = int(temp2.shape[0]*scale)
            colormap = cv2.resize(colormap1,(width,height))
            
            cv2.imwrite(rasf, colormap)
        else:
            image_arry = (temp2 * 255).astype(np.uint8)  
            colormap = cv2.applyColorMap(np.array(image_arry), cv2.COLORMAP_JET) 
            cv2.imwrite(rasf, colormap)
        
        return 1           
    except:
        return 0                 
                        
def rasrmg(unw, dtype, nlines,start_nlines = 0, end_nlines = -1,ph_scale = 6.28318530718,scale = 1,rasf = None,show_flag = 0):
    """
    Function:
      显示相位解缠图,    特别注明:输出路径中不能包含中文字符,且必须以*.bmp,*.jpg,*.png结尾
    Author:
      csuhuacan@163.com,2022.11.10
    Parameters:
      unw:解缠后的相位文件
      dtype:数据类型,gamma对应的为:'>f',python中fwrite输入输出对应的是:'<f'
      start_nlines:输出图片开始的行数,默认为第一行:0
      end_nlines:输出图片结束的行数,默认为最后行:-1
      ph_scale:显示的周期,默认为2π
      scale:对所输出的图片就行缩放,默认为1,大于1进行放大,小于1进行缩小
      rasf:输出图片的路径文件
      show_flag:是否按照周期显示,0:按照ph_scale为周期进行显示,1不按周期,按绝对值显示
    Returns: 
      1:写入成功,0:写入失败 
    Raises:
      KeyError - raises an exception 
    """
    pass
    try:
        
        unw = freadbk(unw,dtype,nlines,0)      # 读取数据
        temp1 = unw[start_nlines:end_nlines,:] # 获取需要显示的部分
        if show_flag ==0 :
            temp1 = np.mod(temp1,ph_scale)                # 计算相位
        else:
            temp1 = temp1
        
        h1 = np.min(temp1)
        h99 = np.max(temp1)
        temp2 = (temp1-h1)/(h99-h1) # 对数据进行归化

        if scale != 1: # 对数据进行缩放
            temp3 = (temp2 * 255).astype(np.uint8)  
            colormap1 = cv2.applyColorMap(np.array(temp3), cv2.COLORMAP_JET) 
            width = int(temp2.shape[1]*scale)
            height = int(temp2.shape[0]*scale)
            colormap = cv2.resize(colormap1,(width,height))
            
            cv2.imwrite(rasf, colormap)
        else:
            image_arry = (temp2 * 255).astype(np.uint8)  
            colormap = cv2.applyColorMap(np.array(image_arry), cv2.COLORMAP_JET) 
            cv2.imwrite(rasf, colormap)
        
        return 1           
    except:
        return 0    

def rashgt(hgt, dtype, nlines,start_nlines = 0, end_nlines = -1,ph_scale = 160,scale = 1,rasf = None,show_flag = 1):
    """
    Function:
      显示高程数据图,主要以相位的形式就行显示,    特别注明:输出路径中不能包含中文字符,且必须以*.bmp,*.jpg,*.png结尾
    Author:
      csuhuacan@163.com,2022.11.10
    Parameters:
      hgt:解缠后的相位文件
      dtype:数据类型,gamma对应的为:'>f',python中fwrite输入输出对应的是:'<f'
      start_nlines:输出图片开始的行数,默认为第一行:0
      end_nlines:输出图片结束的行数,默认为最后行:-1
      ph_scale:显示的周期,默认为160m
      scale:对所输出的图片就行缩放,默认为1,大于1进行放大,小于1进行缩小
      rasf:输出图片的路径文件
      show_flag:是否按照周期显示,0:按照ph_scale为周期进行显示,1不按周期,按绝对值显示
    Returns: 
      1:写入成功,0:写入失败 
    Raises:
      KeyError - raises an exception 
    """
    pass
    try:
        
        hgt = freadbk(hgt,dtype,nlines,0)      # 读取数据
        temp1 = hgt[start_nlines:end_nlines,:] # 获取需要显示的部分
        if show_flag ==0 :                     # 按周期显示
            temp1 = np.mod(temp1,ph_scale)                # 计算相位
            h1 = np.percentile(temp1, 0.01)
            h99 = np.percentile(temp1, 99.99)
            temp2 = (temp1-h1)/(h99-h1) # 对数据进行归化

            if scale != 1: # 对数据进行缩放
                temp3 = (temp2 * 255).astype(np.uint8)  
                colormap1 = cv2.applyColorMap(np.array(temp3), cv2.COLORMAP_JET) 
                width = int(temp2.shape[1]*scale)
                height = int(temp2.shape[0]*scale)
                colormap = cv2.resize(colormap1,(width,height))
                
                cv2.imwrite(rasf, colormap)
            else:
                image_arry = (temp2 * 255).astype(np.uint8)  
                colormap = cv2.applyColorMap(np.array(image_arry), cv2.COLORMAP_JET) 
                cv2.imwrite(rasf, colormap)
        else:                                  # 不按周期显示,输出绘图高程图
            temp1 = temp1
            h1 = np.percentile(temp1, 0.01)
            h99 = np.percentile(temp1, 99.99)
            temp2 = (temp1-h1)/(h99-h1) # 对数据进行归化
            if scale != 1: # 对数据进行缩放
                image_arry = (temp2 * 255).astype(np.uint8)  
                width = int(temp2.shape[1]*scale)
                height = int(temp2.shape[0]*scale)
                colormap = cv2.resize(image_arry,(width,height))
                
                cv2.imwrite(rasf, colormap)
            else:
                image_arry = (temp2 * 255).astype(np.uint8)  
                cv2.imwrite(rasf, image_arry)
        
        
        return 1           
    except:
        return 0    


        