function [H,SIGMA,FAI0,M1,M2,M3] = func(h,sigma,fai0,m1,m2,m3,lb,ub,theta,kz,Y1,Y2,Y3,m,n)
% 功能：目标规划函数，使用非线性最小二乘拟合
% 输入参数：h,sigma,fai0,m1,m2,m3,lb,ub,M,N
% 输出参数：矩阵 H,SIGMA,FAI0,M1,M2,M3

% 计算初始函数值
for i = 1 : m
    for j = 1 : n
        x0 = [h(i,j) sigma(i,j) fai0(i,j) m1(i,j) m2(i,j) m3(i,j)];
        y(i,j) = Obj_fun(x0,theta(i,j),kz(i,j),Y1(i,j),Y2(i,j),Y3(i,j));
    end
end
clear i; clear j;clear x0;

% 循环目标规划
parfor i = 1 : m
    for j = 1 : n
        % 初始值
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
            f = @(x)Obj_fun(x,theta(i,j),kz(i,j),Y1(i,j),Y2(i,j),Y3(i,j));
            options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
            X = lsqnonlin(f,x0,lb,ub,options);
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