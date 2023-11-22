
r_truth = struct2array(load('./data/CHM_Lidar.mat'));
r_hv = struct2array(load('./result/H5.mat'));
r_hv = struct2array(load('./result/H6.mat'));
r_hv = h_six_three_s;
r_hv = freadbk('E:/ty/PolInSAR/AfriSAR/BioSAR/result/bin/h_three_mu_s.bin',lines);

r_hv(r_hv>40) = 40;
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


plot_density_RMSE(A,B);