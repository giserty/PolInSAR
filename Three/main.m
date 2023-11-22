%% 三阶段算法
clear;clc
lines = 11099;
cols = 1510;
m = lines;
n = cols;
font_size = 10;

%% 加载数据
work_path = 'F:/童谣_本科毕设/PolInSAR/FHI/data/0201_0207_FER_MLK_GSS/T6/cmplx_coh_';
PDH = freadbk([work_path,'PDHigh.bin'],lines,'cpxfloat32');
PDL = freadbk([work_path,'PDLow.bin'],lines,'cpxfloat32');
Opt1 = freadbk([work_path,'Opt1.bin'],lines,'cpxfloat32');
Opt2 = freadbk([work_path,'Opt2.bin'],lines,'cpxfloat32');
Opt3 = freadbk([work_path,'Opt3.bin'],lines,'cpxfloat32');
HH = freadbk([work_path,'HH.bin'],lines,'cpxfloat32');
HV = freadbk([work_path,'HV.bin'],lines,'cpxfloat32');
VV = freadbk([work_path,'VV.bin'],lines,'cpxfloat32');
HHpVV = freadbk([work_path,'HHpVV.bin'],lines,'cpxfloat32');
HHmVV = freadbk([work_path,'HHmVV.bin'],lines,'cpxfloat32');
Inc = struct2array(load('./data/IncAngle0201.mat'));
Local_Inc = struct2array(load('./data/Local_IncAngle0201.mat'));
theta = Local_Inc;
Rngslope = Inc - Local_Inc;
Kz = struct2array(load('./data/Kz.mat'));
kz = squeeze(Kz(:,:,3));
kz = kz.*sin(Inc)./sin(Local_Inc);
kz = abs(kz);

%% 直线拟合和地表相位
n_pol = 2;
[fai0,p] = line_fai(PDH,PDL,Opt1,Opt2,Opt3,HH,HV,VV,HHpVV,HHmVV,n_pol,m,n);

%% 查找表
fai0 = struct2array(load('./result/fai2.mat'));
fai0(isnan(fai0)) = 0;
dethv = 0.5;
detg = 0.02;
detmu = 0.09;
yv = HV;
% 执行LUT
tic
theta = Inc;
[h_min,sigma_min] = LUT(fai0,dethv,detg,yv,m,n,theta,kz);
% sig = 0.5;
% [h_min,sigma_min] = LUT_m(fai0,dethv,detmu,yv,m,n,theta,kz,sig);
% h_min = h_min ./ cos(Rngslope);
% fwritebk(h_min,'./result/three_mu_s/h_three_mu_s_5.bin');
toc

%% 精度评价
r_truth = struct2array(load('./data/CHM_Lidar.mat'));
r_hv = freadbk('./result/three_mu/h_three_mu_3.bin',lines);  %%0.1 3.44; 0.2 3.4; 0.3 3.58; 0.4 3.85    S:0.1 3.96; 0.2 3.52; 0.3 3.48; 0.4 3.57

r_hv(r_hv>40) = 40;
r_hv(r_hv<0) = 0;
p = 48;
q = 40;
% 掩膜
[A,B,~] = nn_cktjf(p,q,r_hv,r_truth,0.8); % 0.8是有效像素比
% 精度评价
[RMSE,R2] = evalue(A,B);
% 绘图
plot_density(A,B,RMSE,R2);


%% tiff to bin
h_three_mu_s = struct2array(load('./result/h_three_mu_s.mat'));
fwritebk(r_hv,'E:/ty/PolInSAR/AfriSAR/BioSAR/result/bin/h_three_mu_s.bin');

