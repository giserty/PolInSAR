function [ground_phase]=traditional_pdh_pdl_line(ypdh,ypdl,flag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Name: traditional_pdh_pdl_line.m
%   Author:  wuchuanjun & songtianyi
%   Description: ����ֱ�����ȷ���ر���λ
%   Input: (1)�����ϵ��: ypdh,ypdl��(2)flag:if kz > 0,flag = 1;kz < 0,flag = -1
%   Output:(1)ground_phase: �ر���λ��
%   Date��2019/08/11
%   Version: V1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%����ģ�ͣ�y=bx+a
%��ȡ�����С,��������ģ��x,y��С
[row,col]=size(ypdh);
%���õر���λground_phase����
ground_phase=zeros(row,col);
% parfor_progress(row);
parfor m=1:row
   
    for n=1:col
        
%% ��һ�������ֱ��ȷ���ر���λ
        %��ȡPDH,PDL�ĸ������Ϣ�����x,y
       x1=real(ypdh(m,n));x2=real(ypdl(m,n));
       y1=imag(ypdh(m,n));y2=imag(ypdl(m,n));
       b=(y2-y1)/(x2-x1);
       a=y1-b*x1;
%%  �ڶ�����ȷ���ر���λ
        %�������ֱ��y=bx+a�뵥λԲx^2+y^2=1����������
        
        delta=4*a^2*b^2-4*(1+b^2)*(a^2-1);
        fai1_x=(-2*a*b+sqrt(delta))/(2*(1+b^2));
        fai1_y=b*fai1_x+a;
        fai1=fai1_x+fai1_y*1i;
        fai2_x=(-2*a*b-sqrt(delta))/(2*(1+b^2));
        fai2_y=b*fai2_x+a;
        fai2=fai2_x+fai2_y*1i;
        
%%%%%%%%%%%% KZ�������жϣ������Ƚ�%%%%%%%%%      
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
%%%%%��������-2015��Forest Height Estimation by Means of Pol-InSAR Data
%%%%%Inversion: The Role of the Vertical Wavenumber-IEEE.TGRS.%%%%

%         %��������PDH����ɢ��ռ�ţ�����ֲ����λ����PDL�������ɢ��ռ�ţ�����ر���λ���Ƚϣ���PDH��Զ����PDL�Ͻ��ĵ�Ϊ�ر���λ��
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



