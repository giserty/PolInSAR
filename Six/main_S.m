%% 六维非线性迭代算法
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
Rngslope  = Inc - Local_Inc;
theta = Local_Inc;
Kz = struct2array(load('Kz.mat'));
kz = squeeze(Kz(:,:,3));
kz = kz.*sin(Inc)./sin(Local_Inc);
kz = abs(kz);

%% 变量初始化
% h = ones(lines,cols) .* 0.0001;
% h = struct2array(load('hv.mat'));  % DEM差分法生成的树高
h = struct2array(load('h_0305_2'));  % 三阶段算法生成的树高(S-RVoG)
h(isnan(h)) = 0;
sigma = ones(lines,cols) .* 0.0001;
fai0 = struct2array(load('fai0.mat'));  % 三阶段算法生成的地表相位
fai0(isnan(fai0)) = 0;
m1 = abs((exp(1i.*fai0).*PDH-Opt1)./(Opt1-exp(1i.*fai0)));
m2 = abs((exp(1i.*fai0).*PDH-Opt2)./(Opt2-exp(1i.*fai0)));
m3 = abs((exp(1i.*fai0).*PDH-Opt3)./(Opt3-exp(1i.*fai0)));
% 上下界
lb = [0,0.0001,-pi,0.01,0.01,0.01]; 
ub = [40,0.23,pi,100,100,100];

%% 目标规划
% clear H;clear SIGMA;clear FAI0;clear M1;clear M2;clear M3
for i = 1 : lines
    for j = 1 : cols
        x0 = [h(i,j) sigma(i,j) fai0(i,j) m1(i,j) m2(i,j) m3(i,j)];
        y(i,j) = Obj_fun(x0,theta(i,j),kz(i,j),Opt1(i,j),Opt2(i,j),Opt3(i,j));
    end
end
clear i; clear j;clear x0;

tic
% M = ceil(lines/2);
% N = ceil(cols/2);
parfor i = 1 : lines
    for j = 1 : cols
        % 该像素各变量初值
        x0 = [h(i,j) sigma(i,j) fai0(i,j) m1(i,j) m2(i,j) m3(i,j)];
        % 初始函数值为空
        if isnan(y(i,j))
            H(i,j) = 0;
            SIGMA(i,j) = 0.0001;
            FAI0(i,j) = 0;
            M1(i,j) = 0.01;
            M2(i,j) = 0.01;
            M3(i,j) = 0.01;
        else
            % 该像素目标函数
            f = @(x)Obj_fun(x,theta(i,j),kz(i,j),Opt1(i,j),Opt2(i,j),Opt3(i,j));
            options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
            X = lsqnonlin(f,x0,lb,ub);
            % 赋值
            H(i,j) = X(1);
            SIGMA(i,j) = X(2);
            FAI0(i,j) = X(3);
            M1(i,j) = X(4);
            M2(i,j) = X(5);
            M3(i,j) = X(6);
        end
    end
    disp(i);
end

H = H2./cos(Rngslope);
H(H>40) = 40;
H(H<0) = 0;
save('H22.mat','H');
save('SIGMA2','SIGMA');
save('FAI02','FAI0');
save('M12','M1');
save('M22','M2');
save('M32','M3');

toc




%%
H2 = struct2array(load('H2.mat'));
% H2 = H2./Rngslope;
clims = [0 30];
figure,imagesc(H,clims);
colormap jet;
colorbar;
set(gca,'FontSize',font_size);

pbaspect([1 2.5 1]);

%%
CHM_LiDAR = struct2array(load('CHM_Lidar.mat'));
% H2 = struct2array(load('h_0305_2.mat'));
r_hv = H;
r_truth = CHM_LiDAR;
% 裁剪
mask = ones(lines,cols);
mask(abs(r_hv)<0.3)=0;
% mask=multilook(mask,40,40);
mask(r_truth==0)=0;
A = mask.*r_truth;
B = mask.*r_hv;
A = multilook(A,40,40);
B = multilook(B,40,40);
A(isnan(A)) = 0;
B(isnan(B)) = 0;
A(A==0) = [];
B(B==0) = [];
% 计算参数
RMSE = sqrt((sum((A(:)-B(:)).^2))./(length(A)-1));
R = corrcoef(A, B);
R2 = R(1,2)^2;
% 绘图
plot_density(A(1:100),B(1:100),RMSE,R2);
