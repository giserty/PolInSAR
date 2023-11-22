import matplotlib.pyplot as plt
# from numpy.lib.financial import pv 
import sys
import numpy as np
from math import *
plt.rcParams['font.sans-serif']=['SimHei'] #显示中文标签
plt.rcParams['axes.unicode_minus']=False   #这两行需要手动设置
from scipy.ndimage import interpolation
##绘图显示
font1 = {'family' : 'Times New Roman',
'weight' : 'normal',
'size' : 20,
}

# 计算r_2
def computeCorrelation(x,y):
    xBar = np.mean(x)
    yBar = np.mean(y)
    SSR = 0.0
    varX = 0.0
    varY = 0.0
    for i in range(0,len(x)):
        diffXXbar = x[i] - xBar
        difYYbar = y[i] - yBar
        SSR += (diffXXbar * difYYbar)
        varX += diffXXbar**2
        varY += difYYbar**2
    SST = sqrt(varX * varY)
    return (SSR/SST)


# def plot_imagesc(im,title,xlab,ylab):
#     fig = plt.figure()
#     #plt.pcolor(im,cmap = 'jet')
#     plt.pcolor(im,cmap = 'jet')
#     plt.colorbar()
#     plt.jet
#     plt.clim(0,1)

#     plt.colorbar()
#     plt.jet() 
#     plt.xticks(fontsize=10);plt.yticks(fontsize=10)
#     plt.xlabel(xlab,font1) #为子图设置横轴标题
#     plt.ylabel(ylab,font1) #为子图设置纵轴标题
#     plt.title(title, size=12)
#     return plt

def plot_imagesc(im,title,xlab,ylab):
    fig = plt.figure()
    plt.pcolor(im,cmap = 'jet')
    plt.colorbar()
    plt.jet
    plt.clim(0,1) 
    plt.xticks(fontsize=10);plt.yticks(fontsize=10)
    plt.xlabel(xlab,font1) #为子图设置横轴标题
    plt.ylabel(ylab,font1) #为子图设置纵轴标题
    plt.title(title, size=12)
    return plt


import matplotlib as mpl
def plot_imshow(im,title,xlab,ylab):
    fig = plt.figure()
    plt.imshow(im)
    plt.colorbar()
    plt.jet() 
    plt.xticks(fontsize=10);plt.yticks(fontsize=10)
    plt.xlabel(xlab,font1) #为子图设置横轴标题
    plt.ylabel(ylab,font1) #为子图设置纵轴标题
    plt.title(title, size=12)
    return plt



import time
#进度条函数
def view_bar(num, total):
    rate = float(num) / total
    rate_num = int(rate * 100)+1
    r = '\r[%s%s]%d%%' % ("#"*rate_num, " "*(100-rate_num), rate_num, )
    sys.stdout.write(r)
    sys.stdout.flush()


def plot_density(XX,YY,title,xlab,ylab):
    XX[YY==0] = 0
    YY[XX==0] = 0
    t = np.argwhere(YY != 0)
    size = XX.shape
    if sum(size) == size[0]:
        x = XX[t[:,0]]
        y = YY[t[:,0]]
    else:
        x = XX[t[:,0],t[:,1]]
        y = YY[t[:,0],t[:,1]]
  
    size=x.shape
    
    #数据长度
    length = size[0]
    c = np.zeros((length,1))
    max_x = np.max(x);   min_x = np.min(x)   #搜索边界点
    max_y = np.max(y);   min_y = np.min(y)   #搜索边界点

    NLevel = 101   #划分等级150份;分的等级越多，颜色变化的渐变越精细，即颜色条分辨率越高
    color_Map = np.zeros((NLevel+1,NLevel+1))
    step_x = (max_x-min_x)/(NLevel-1);    # x轴步长
    step_y = (max_y-min_y)/(NLevel-1);    # y轴步长 
    for j in range(length):
        color_Map_x = int((x[j]-min_x)/step_x)
        color_Map_y = int((y[j]-min_y)/step_y)
        color_Map[color_Map_x,color_Map_y] = color_Map[color_Map_x,color_Map_y]+1;

    for j in range(length): 
        color_Map_x = int((x[j]-min_x)/step_x)
        color_Map_y = int((y[j]-min_y)/step_y)
        c[j,0]=color_Map[color_Map_x,color_Map_y]

    plt.scatter(x, y, 3, c, alpha=0.4, label='类别A')
    plt.colorbar()
    plt.jet()
    plt.xticks(fontsize=10);plt.yticks(fontsize=10)
    plt.xlabel(xlab,font1) #为子图设置横轴标题
    plt.ylabel(ylab,font1) #为子图设置纵轴标题
    plt.title(title, size=12)
    return plt


from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
def plot_density_ver(x,y,title,xlab,ylab):
    # XX[YY==0] = 0
    # YY[XX==0] = 0
    # t = np.argwhere(YY != 0)
    # size = XX.shape
    # if sum(size) == size[0]:
    #     x = XX[t[:,0]]
    #     y = YY[t[:,0]]
    # else:
    #     x = XX[t[:,0],t[:,1]]
    #     y = YY[t[:,0],t[:,1]]
  
    size=x.shape
    #数据长度
    length = size[0]
    c = np.zeros((length,1))
    max_x = np.max(x);   min_x = np.min(x)   #搜索边界点
    max_y = np.max(y);   min_y = np.min(y)   #搜索边界点

    NLevel = 101   #划分等级150份;分的等级越多，颜色变化的渐变越精细，即颜色条分辨率越高
    color_Map = np.zeros((NLevel+1,NLevel+1))
    step_x = (max_x-min_x)/(NLevel-1);    # x轴步长
    step_y = (max_y-min_y)/(NLevel-1);    # y轴步长 
    for j in range(length):
        color_Map_x = int((x[j]-min_x)/step_x)
        color_Map_y = int((y[j]-min_y)/step_y)
        color_Map[color_Map_x,color_Map_y] = color_Map[color_Map_x,color_Map_y]+1;

    for j in range(length): 
        color_Map_x = int((x[j]-min_x)/step_x)
        color_Map_y = int((y[j]-min_y)/step_y)
        c[j,0]=color_Map[color_Map_x,color_Map_y]

    #决定系数
   
    R2 = computeCorrelation(x,y)
    #均方误差、均方根误差

    mse = mean_squared_error(x,y)
    rmse = np.sqrt(mse)
    x1 = np.linspace(-100,100,40)
    y1 = np.linspace(-100,100,40)
     

    plt.scatter(x, y, 3, c, alpha=0.4, label='类别A')
    plt.colorbar()
    plt.jet()
    plt.xticks(fontsize=10);plt.yticks(fontsize=10)
    plt.xlabel(xlab,font1) #为子图设置横轴标题
    plt.ylabel(ylab,font1) #为子图设置纵轴标题
    plt.title(title, size=12)
    plt.plot(x1,y1,'b-')  
       
    return plt,R2,rmse

#有效样地统计
def Plot_statistics(winh,winl,sg,lidar,bhs):
    
    yd_num = 0
    size = sg.shape
    h = size[0]
    l = size[1]
    yd_sg = np.zeros((int((h*l)/(winh*winl)),1))
    yd_lidar = np.zeros((int((h*l)/(winh*winl)),1))
    for i in range(0, h, winh):
        for j in range(0, l, winl):
            num = 0
            sgnum = 0
            lidarnum = 0
            if i+winh < h:
                h1 = i + winh
            else: 
                h1 = h - 1
            
            if j+winl < l:
                l1 = j + winl
            else:
                l1 = l - 1
            
            for k1 in range(i, h1):
                for k2 in range(j, l1):
                    if (sg[k1,k2] > 0) and (lidar[k1,k2] > 0):
                        num = num+1
                        sgnum = sgnum + sg[k1,k2]
                        lidarnum = lidarnum + lidar[k1,k2]
            if (num/winh/winl) >= bhs:
                yd_sg[yd_num] = sgnum/num
                yd_lidar[yd_num] = lidarnum/num
                yd_num = yd_num + 1

    t1 = np.argwhere(yd_sg != 0)
    y_sg = yd_sg[t1[:,0],t1[:,1]]
    y_lidar = yd_lidar[t1[:,0],t1[:,1]]

    return y_lidar,y_sg,yd_num

def density_ver(XX,YY,title,xlab,ylab):
    XX[YY==0] = 0
    YY[XX==0] = 0
    t = np.argwhere(YY != 0)
    size = XX.shape
    if sum(size) == size[0]:
        x = XX[t[:,0]]
        y = YY[t[:,0]]
    else:
        x = XX[t[:,0],t[:,1]]
        y = YY[t[:,0],t[:,1]]

    #均方误差、均方根误差

    mse = mean_squared_error(x,y)
    rmse = np.sqrt(mse)
    return rmse


#功能：将数组写入二进制文件
