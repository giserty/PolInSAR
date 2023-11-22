%% 六维非线性迭代算法
clear;clc
lines = 11099;
cols = 1510;
% 变量初始化 目标规划 精度评价会用到
m = lines;
n = cols;
font_size = 10;

%% 加载数据
yv = "PDH";
y1 = "HH";
y2 = "HV";
y3 = "VV";
s = 'y';  %  只能填单个字符 n or y
[YV,Y1,Y2,Y3,theta,kz] = load_data(lines,yv,y1,y2,y3,s);

%% DEM差分法
% h_dem = (angle(Y2)-angle(Y1))./kz;
% fwritebk(h_dem,'../Dem/h_dem_s.bin','float32');

%% 变量初始化
% h = ones(lines,cols) .* 0.001;
% h_dem = struct2array(load('hv.mat'));  % DEM差分法生成的树高
% h = struct2array(load('./data/H23'));  % 三阶段算法生成的树高

% fai,hh,sig
% fai = struct2array(load('./data/fai2.mat'));  % 三阶段算法生成的地表相位
% fai(isnan(fai)) = 0;
fai = 0;
hh = h_dem;  % struct2array(load('./data/H23'));
hh(isnan(hh)) = 0;
hh(hh>40) = 40;
hh(hh<0) = 0;
sig = 0.1;
% 执行生成初值的函数
[h,sigma,fai0,m1,m2,m3] = initial_value(hh,sig,fai,m,n,YV,Y1,Y2,Y3);
% 上下界
lb = [0,0.0001,-pi,0.01,0.01,0.01]; 
ub = [40,0.23,pi,100,100,100];

%% 目标规划
% clear H;clear SIGMA;clear FAI0;clear M1;clear M2;clear M3
tic
% 执行目标规划函数
[H,SIGMA,FAI0,M1,M2,M3] = func(h,sigma,fai0,m1,m2,m3,lb,ub,theta,kz,Y1,Y2,Y3,m,n);
% % 保存结果变量
save('./result/H_dem_s.mat','H');
save('./result/SIGMA_dem_s','SIGMA');
save('./result/FAI0_dem_s','FAI0');
save('./result/M1_dem_s','M1');
save('./result/M2_dem_s','M2');
save('./result/M3_dem_s','M3');
toc

%% 绘图
img = H;
clims = [0 30];
% 执行绘图函数
draw(img,font_size,clims);

%% 精度评价
% 导入精度评价数据
r_truth = struct2array(load('./data/CHM_Lidar.mat'));
r_hv = freadbk('E:/ty/PolInSAR/AfriSAR/BioSAR/result/bin/h_six_three.bin',lines);
r_hv(r_hv>30) = 30;
r_hv(r_hv<0) = 0;
% 多视窗口，行列号在第一节定义
p = 48;
q = 40;
% 执行掩膜和精度评价函数
[A,B,~] = nn_cktjf(p,q,r_hv,r_truth,0.8); % 0.8是有效像素比
% 计算参数
RMSE = sqrt((sum((A-B).^2))./(length(A)-1));
RMSE = roundn(RMSE,-2); % 保留2位小数
R = corrcoef(A, B);
R2 = roundn(R(1,2),-2); % 保留2位小数
% 绘图
plot_density(A,B,RMSE,R2);
