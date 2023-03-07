%% 三阶段算法
% 输入数据
%%相干性：Opt1, Opt2, Opt3, PDH, PDL
%%RVoG参数：Inc, Local_Inc, Kz
%%第二步计算的地表相位：fai2 or fai5
%%第三步计算的树高：H1, H2
%%CHM_LiDAR
clear;clc
lines = 11099;
cols = 1510;
font_size = 10;

%% 加载数据
Opt1 = struct2array(load('Opt1.mat'));
Opt2 = struct2array(load('Opt2.mat'));
Opt3 = struct2array(load('Opt3.mat'));
PDH = struct2array(load('PDH.mat'));
PDL = struct2array(load('PDL.mat'));
Inc = struct2array(load('IncAngle0201.mat'));
Local_Inc = struct2array(load('Local_IncAngle0201.mat'));
Rngslope = Inc - Local_Inc;
theta = Local_Inc;
Kz = struct2array(load('Kz.mat'));
kz = squeeze(Kz(:,:,3));
kz = kz.*sin(Inc)./sin(Local_Inc);
kz = abs(kz);

%% 相位幅度计算
Opt1_pha = angle(Opt1); Opt1_mod = abs(Opt1); 
Opt1_re = real(Opt1); Opt1_im = imag(Opt1);
Opt2_pha = angle(Opt2); Opt2_mod = abs(Opt2); 
Opt2_re = real(Opt2); Opt2_im = imag(Opt2);
Opt3_pha = angle(Opt3); Opt3_mod = abs(Opt3); 
Opt3_re = real(Opt3); Opt3_im = imag(Opt3);
PDH_pha = angle(PDH); PDH_mod = abs(PDH); 
PDH_re = real(PDH); PDH_im = imag(PDH);
PDL_pha = angle(PDL); PDL_mod = abs(PDL); 
PDL_re = real(PDL); PDL_im = imag(PDL);

%% 直线拟合和地表相位
tic
for i = 1 : 1ines
    for j = 1 : cols
        % 直线拟合
        x = [Opt1_re(i,j) Opt2_re(i,j) Opt3_re(i,j) PDH_re(i,j) PDL_re(i,j)];
        y = [Opt1_im(i,j) Opt2_im(i,j) Opt3_im(i,j) PDH_im(i,j) PDL_im(i,j)];
        f = polyfit(x,y,1);
        k = f(1);
        b = f(2);
        % 计算交点
        aa = 1+k*k;
        bb = 2*k*b;
        cc = b*b-1;
        xx = roots([aa bb cc]);
        yy = xx.*k + b;
        zz = xx + yy.*1i;
        pha = angle(zz);  % 某像素2个相位
        % 筛选地表相位
        dist1_pdh = abs(zz(1)-PDH(i,j));
        dist1_pdl = abs(zz(1)-PDL(i,j));
        if dist1_pdh > dist1_pdl
            fai0(i,j) = pha(1);
            p(i,j) = zz(1);
        else
            fai0(i,j) = pha(2);
            p(i,j) = zz(2);
        end
    end
end
toc

%% 查找表
% 5通道地表相位
fai0 = struct2array(load('fai5.mat'));

tic
sigma = (0.001:0.05:0.5)'.*0.4;
size_sigma = max(size(sigma));
h = 0.01:0.2:30.01;
size_h = max(size(h));
sigma = sigma*ones(1,size_h);
h = ones(size_sigma,1)*h;

parfor i = 1 : lines
    for j = 1 : cols
        try
            p = 2 * sigma / cos(theta(i,j));
            p1 = p + kz(i,j)*1i;
            y = exp(fai0(i,j)*1i) .* (p./p1) .* (exp(p1.*h)-1)./(exp(p.*h)-1);
            if abs(PDH_pha(i,j)-fai0(i,j)) > abs(PDL_pha(i,j)-fai0(i,j))
                F = abs(PDH(i,j) - y);
            else
                F = abs(PDL(i,j) - y);
            end
            F_min = min(F(:));
            [row, col] = find(F == F_min);
            sigma_min(i,j) = 0.001+(row-1)*0.05;
            h_min(i,j) = 0.01+(col-1)*0.2;
        catch
            sigma_min(i,j) = NaN;
            h_min(i,j) = NaN;
        end
    end
    if (mod(i,1000) == 0)
        disp(['已执行完',num2str(i),'行']);
    end
end

% 计算树高
h_min = h_min ./ cos(Rngslope);
disp('LUT算法已完成');
toc

%% 绘图
% 极坐标图
polarplot(Opt1(i,j),'ro','linewidth',1)
hold on
polarplot(Opt2(i,j),'bd','linewidth',1)
polarplot(Opt3(i,j),'yp','linewidth',1)
polarplot(PDH(i,j),'gs','linewidth',1)
polarplot(PDL(i,j),'k*','linewidth',1)
% 绘制直线
X = -pi : 0.01 : pi;
Y = b ./ (sin(X)-k.*cos(X));
polarplot(X,Y,'k-')
legend('Opt1','Opt2','Opt3','PDH','PDL');
rlim([0 1])
% 绘制地表相位
figure,imagesc(fai0);
colormap jet;
colorbar;
set(gca,'FontSize',font_size);
pbaspect([1 2.5 1]);

%% 精度评价
% 导入数据
CHM_LiDAR = struct2array(load('CHM_Lidar.mat'));
% 未进行坡度校正的树高
H1 = struct2array(load('h_0301_1.mat'));
% 进行了坡度校正的树高
H2 = struct2array(load('h_0305_2.mat'));
% r_hv = H2-4;
r_truth = CHM_LiDAR;
% 裁剪
mask = ones(lines,cols);
mask(abs(r_hv)<0.3) = 0;
mask(r_truth==0) = 0;
A = mask.*r_truth;
B = mask.*r_hv;
A = multilook(A,40,40);
B = multilook(B,40,40);
A(A==0) = [];
B(B==0) = [];
% 计算参数
RMSE = sqrt((sum((A(:)-B(:)).^2))./(length(A)-1));
R = corrcoef(A, B);
R2 = R(1,2)^2;
% 绘图
plot_density(A,B,RMSE,R2);
