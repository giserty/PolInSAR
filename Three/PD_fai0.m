function [fai0,p] = PD_fai0(ypdh,ypdl,m,n)
% ����PDHigh��PDLow��ر���λ
[row,col]=size(ypdh);    % ��ȡ�����С,��������ģ��x,y��С
fai0=zeros(row,col);     % ���������ر���λfai0����
p=zeros(row,col);

% RVOGģ�����׶η�����ֲ����
parfor i=1:m
    for j=1:n
        x1=real(ypdh(i,j));x2=real(ypdl(i,j));
        y1=imag(ypdh(i,j));y2=imag(ypdl(i,j));
        b=(y2-y1)./(x2-x1);
        a=y1-b*x1;
        r=1; 
        delta=4*a.^2.*b.^2-4.*(1+b.^2).*(a.^2-r.^2);
        fai1_x=(-2*a*b+sqrt(delta))/(2*(1+b^2));
        fai1_y=b*fai1_x+a;
        fai1=fai1_x+fai1_y*1i;
        fai2_x=(-2*a*b-sqrt(delta))/(2*(1+b^2));
        fai2_y=b*fai2_x+a;
        fai2=fai2_x+fai2_y*1i;
        
        %���������ж��ĸ�����
        %%��������HV/PDH����ɢ��ռ�ţ�����ֲ����λ����HH-VV/PDL�������ɢ��ռ�ţ�����ر���λ���Ƚϣ�����ĵ����жϵر���λ��
        dist1_pdh=abs(fai1-ypdh(i,j));
        dist1_pdl=abs(fai1-ypdl(i,j));
        if dist1_pdh>dist1_pdl
            fai0(i,j)=angle(fai1);
            p(i,j)=fai1;
        else
            fai0(i,j)=angle(fai2);
            p(i,j)=fai2;
        end
    end
end
