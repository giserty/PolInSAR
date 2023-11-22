%% ��ά�����Ե����㷨
clear;clc
lines = 11099;
cols = 1510;
% ������ʼ�� Ŀ��滮 �������ۻ��õ�
m = lines;
n = cols;
font_size = 10;

%% ��������
yv = "PDH";
y1 = "HH";
y2 = "HV";
y3 = "VV";
s = 'y';  %  ֻ������ַ� n or y
[YV,Y1,Y2,Y3,theta,kz] = load_data(lines,yv,y1,y2,y3,s);

%% DEM��ַ�
% h_dem = (angle(Y2)-angle(Y1))./kz;
% fwritebk(h_dem,'../Dem/h_dem_s.bin','float32');

%% ������ʼ��
% h = ones(lines,cols) .* 0.001;
% h_dem = struct2array(load('hv.mat'));  % DEM��ַ����ɵ�����
% h = struct2array(load('./data/H23'));  % ���׶��㷨���ɵ�����

% fai,hh,sig
% fai = struct2array(load('./data/fai2.mat'));  % ���׶��㷨���ɵĵر���λ
% fai(isnan(fai)) = 0;
fai = 0;
hh = h_dem;  % struct2array(load('./data/H23'));
hh(isnan(hh)) = 0;
hh(hh>40) = 40;
hh(hh<0) = 0;
sig = 0.1;
% ִ�����ɳ�ֵ�ĺ���
[h,sigma,fai0,m1,m2,m3] = initial_value(hh,sig,fai,m,n,YV,Y1,Y2,Y3);
% ���½�
lb = [0,0.0001,-pi,0.01,0.01,0.01]; 
ub = [40,0.23,pi,100,100,100];

%% Ŀ��滮
% clear H;clear SIGMA;clear FAI0;clear M1;clear M2;clear M3
tic
% ִ��Ŀ��滮����
[H,SIGMA,FAI0,M1,M2,M3] = func(h,sigma,fai0,m1,m2,m3,lb,ub,theta,kz,Y1,Y2,Y3,m,n);
% % ����������
save('./result/H_dem_s.mat','H');
save('./result/SIGMA_dem_s','SIGMA');
save('./result/FAI0_dem_s','FAI0');
save('./result/M1_dem_s','M1');
save('./result/M2_dem_s','M2');
save('./result/M3_dem_s','M3');
toc

%% ��ͼ
img = H;
clims = [0 30];
% ִ�л�ͼ����
draw(img,font_size,clims);

%% ��������
% ���뾫����������
r_truth = struct2array(load('./data/CHM_Lidar.mat'));
r_hv = freadbk('E:/ty/PolInSAR/AfriSAR/BioSAR/result/bin/h_six_three.bin',lines);
r_hv(r_hv>30) = 30;
r_hv(r_hv<0) = 0;
% ���Ӵ��ڣ����к��ڵ�һ�ڶ���
p = 48;
q = 40;
% ִ����Ĥ�;������ۺ���
[A,B,~] = nn_cktjf(p,q,r_hv,r_truth,0.8); % 0.8����Ч���ر�
% �������
RMSE = sqrt((sum((A-B).^2))./(length(A)-1));
RMSE = roundn(RMSE,-2); % ����2λС��
R = corrcoef(A, B);
R2 = roundn(R(1,2),-2); % ����2λС��
% ��ͼ
plot_density(A,B,RMSE,R2);
