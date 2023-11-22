
function [RMSE,R2] = evalue(A,B)
RMSE = sqrt((sum((A-B).^2))./(length(A)-1));
RMSE = roundn(RMSE,-2); % 保留2位小数
R = corrcoef(A, B);
R2 = roundn(R(1,2),-2); % 保留2位小数
end

%% 分区间精度评价
% % 多视窗口，行列号在第一节定义
% p = 1;
% q = 1;
% r_truth = struct2array(load('./data/CHM_Lidar.mat'));
% hv_prod = freadbk('E:/ty/PolInSAR/AfriSAR/BioSAR/result/pauli/prod/hv_prod.bin',lines);
% hv_ecc = freadbk('E:/ty/PolInSAR/AfriSAR/BioSAR/result/pauli/ecc/hh_hhpvv/hv_ecc.bin',lines);
% hv_var = freadbk('E:/ty/PolInSAR/AfriSAR/BioSAR/result/pauli/var/hv_var.bin',lines);
% 
% % 执行掩膜和精度评价函数
% lb = 30;
% ub = 40;
% r_truth(r_truth<=lb) = 0;
% r_truth(r_truth>ub) = 0;
% r_truth(r_truth==0) = nan;
% hv_prod(isnan(r_truth)) = nan;
% hv_ecc(isnan(r_truth)) = nan;
% hv_var(isnan(r_truth)) = nan;
% 
% r_hv = hv_prod;
% 
% [A,B,~] = nn_cktjf(p,q,r_hv,r_truth,0.8); % 0.8是有效像素比
% % 计算参数
% RMSE = sqrt((sum((A-B).^2))./(length(A)-1));
% RMSE = roundn(RMSE,-2); % 保留2位小数
% R = corrcoef(A, B);
% R2 = roundn(R(1,2),-2); % 保留2位小数
% 
% num_p = sum(sum(~isnan(A)));
