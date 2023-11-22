%% SA模拟退火
tic
clear; clc

lines = 11099;
cols = 1510;
Opt1 = struct2array(load('Opt1.mat'));
PDLow = struct2array(load('PDLow.mat'));
hv = struct2array(load('hv.mat'));   % 森林高度
pha0 = angle(PDLow);   % 地表相位
sigma = ones(lines,cols) .* 0.1;   % 消光系数
m1 = ones(lines,cols) .* 0.1;  % 地体幅度比1
% m2 = 0.2;  % 地体幅度比2
% m3 = 0.3;  % 地体幅度比3

%% 参数初始化
narvs = 4;  % 变量个数
T0 = 1000;   % 初始温度
T = T0; % 迭代中温度会发生改变，第一次迭代时温度就是T0
maxgen = 100;  % 最大迭代次数
Lk = 30;  % 每个温度下的迭代次数
alfa = 0.95;  % 温度衰减系数
theta = sqrt(2)/2;  % 入射角
kz = 1;  % 垂直波数
% sigma = 0.1;  % 消光系数
% hv = 10;  % 树高
% pha = 0.5;  % 相位
% m1 = 0.1;  % 地体幅度比
sigma_lb = zeros(lines,cols);
sigma_ub = ones(lines,cols);
h_lb = ones(lines,cols) .* 10;
h_ub = ones(lines,cols) .* 30;
pha_lb = ones(lines,cols) .* 0.5;
pha_ub = ones(lines,cols) .* 5;
m_lb = ones(lines,cols) .* 0.1;
m_ub = ones(lines,cols);
x_lb = [sigma_lb h_lb pha_lb m_lb]; % x的下界
x_ub = [sigma_ub h_ub pha_ub m_ub]; % x的上界

%%  生成初始解
% x0 = zeros(1,narvs);  % 初始0矩阵
% for i = 1: narvs
%     x0(i) = x_lb(i) + (x_ub(i)-x_lb(i))*rand(1);  % 随机初始化粒子所在的位置在定义域内
% end
x0 = [hv,pha0,sigma,m1];
y0 = abs(Opt1 - Obj_fun3(x0));  % 计算当前解的函数值

%% 定义一些保存中间过程的量，方便输出结果
min_y = y0;     % 初始化找到的最佳的解对应的函数值为y0

%% 模拟退火过程
for iter = 1 : maxgen  % 外循环, 我这里采用的是指定最大迭代次数
    for i = 1 : Lk  %  内循环，在每个温度下开始迭代
        y = randn(lines,cols*4);  % 生成m行n列的N(0,1)随机数
        z = y ./ sqrt(sum(sum(y.^2))); % 根据新解的产生规则计算z
        x_new = x0 + z.*T; % 根据新解的产生规则计算x_new的值
        % 如果这个新解的位置超出了定义域，就对其进行调整
%         for j = 1: lines*cols*4
%             if x_new(j) < x_lb(j)
%                 r = rand(1);
%                 x_new(j) = r.*x_lb(j)+(1-r).*x0(j);
%             elseif x_new(j) > x_ub(j)
%                 r = rand(1);
%                 x_new(j) = r.*x_ub(j)+(1-r).*x0(j);
%             end
%         end
        x1 = x_new;    % 将调整后的x_new赋值给新解x1
        y1 = abs(Opt1 - Obj_fun3(x1));  % 计算新解的函数值
        if y1 < y0    % 如果新解函数值小于当前解的函数值
            x0 = x1; % 更新当前解为新解
            y0 = y1;
        else
            p = exp(-(y0 - y1)/T); % 根据Metropolis准则计算一个概率
            if rand(1) < p   % 生成一个随机数和这个概率比较，如果该随机数小于这个概率
                x0 = x1; % 更新当前解为新解
                y0 = y1;
            end
        end
        % 判断是否要更新找到的最佳的解
        if y0 < min_y  % 如果当前解更好，则对其进行更新
            min_y = y0;  % 更新最小的y
            best_x = x0;  % 更新找到的最好的x
        end
    end
    T = alfa*T;   % 温度下降
end

%%
toc
disp(['运行时间: ',num2str(toc)]);
