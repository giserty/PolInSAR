function [ground_phase]=traditional_pdh_pdl_line(ypdh,ypdl,flag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Name: traditional_pdh_pdl_line.m
%   Author:  wuchuanjun & songtianyi
%   Description: 两点直线拟合确定地表相位
%   Input: (1)复相干系数: ypdh,ypdl；(2)flag:if kz > 0,flag = 1;kz < 0,flag = -1
%   Output:(1)ground_phase: 地表相位。
%   Date：2019/08/11
%   Version: V1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%线性模型：y=bx+a
%获取矩阵大小,设置线性模型x,y大小
[row,col]=size(ypdh);
%设置地表相位ground_phase矩阵
ground_phase=zeros(row,col);
% parfor_progress(row);
parfor m=1:row
   
    for n=1:col
        
%% 第一步，相干直线确定地表相位
        %获取PDH,PDL的复相干信息，组成x,y
       x1=real(ypdh(m,n));x2=real(ypdl(m,n));
       y1=imag(ypdh(m,n));y2=imag(ypdl(m,n));
       b=(y2-y1)/(x2-x1);
       a=y1-b*x1;
%%  第二步，确定地表相位
        %计算拟合直线y=bx+a与单位圆x^2+y^2=1的两个交点
        
        delta=4*a^2*b^2-4*(1+b^2)*(a^2-1);
        fai1_x=(-2*a*b+sqrt(delta))/(2*(1+b^2));
        fai1_y=b*fai1_x+a;
        fai1=fai1_x+fai1_y*1i;
        fai2_x=(-2*a*b-sqrt(delta))/(2*(1+b^2));
        fai2_y=b*fai2_x+a;
        fai2=fai2_x+fai2_y*1i;
        
%%%%%%%%%%%% KZ正负号判断，更加稳健%%%%%%%%%      
        if flag==1
            if angle(fai1*conj(fai2))>0
                ground_phase(m,n)=angle(fai1);
            else
                ground_phase(m,n)=angle(fai2);
            end
        else
            if angle(fai1*conj(fai2))>0
                ground_phase(m,n)=angle(fai2);
            else
                ground_phase(m,n)=angle(fai1);
            end
        end
%%%%%参照论文-2015年Forest Height Estimation by Means of Pol-InSAR Data
%%%%%Inversion: The Role of the Vertical Wavenumber-IEEE.TGRS.%%%%

%         %将交点与PDH（体散射占优，代表植被相位）和PDL（二面角散射占优，代表地表相位）比较，距PDH较远、距PDL较近的点为地表相位点
%         dist1_ypdh=abs(fai1-ypdh(m,n));
%         dist1_ypdl=abs(fai1-ypdl(m,n));
%         if dist1_ypdh > dist1_ypdl 
%             ground_phase(m,n)=angle(fai1);
%         else
%             ground_phase(m,n)=angle(fai2);
%         end       
    end
%   parfor_progress;  
end
% parfor_progress(0);



