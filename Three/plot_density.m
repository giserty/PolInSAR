%  �����ܶȵ�ɢ��ͼ
% x�����꣬y������
function plot_density(x,y,RMSE,r2,accuracy)
if nargin==1   %����1���������
    if size(x,1)
        x=x';    %ת��Ϊ������
    end
    y=x;
    x=(1:length(y))';
else
    if nargin==2   %��2���������
        if size(x,1)
            x=x';   %ת��Ϊ������
        end
        if size(y,1)
            y=y';       %ת��Ϊ������
        end
    end
end
N=length(x);     %���ݳ���
c=zeros(N,1);
max_x=max(x);
min_x=min(x);    %�����߽��
max_y=max(y);
min_y=min(y);    %�����߽��
NLevel=300;      %���ֵȼ�150��;�ֵĵȼ�Խ�࣬��ɫ�仯�Ľ���Խ��ϸ������ɫ���ֱ���Խ��
color_Map=zeros(NLevel+1);
step_x=(max_x-min_x)/(NLevel-1);  % x�Ჽ��
step_y=(max_y-min_y)/(NLevel-1);    % y�Ჽ�� 
for j=1:N
    color_Map_x=int32((x(j)-min_x)/step_x)+1;
    color_Map_y=int32((y(j)-min_y)/step_y)+1;
    color_Map(color_Map_x,color_Map_y)=color_Map(color_Map_x,color_Map_y)+1;
end
for j=1:N
    color_Map_x=int32((x(j)-min_x)/step_x)+1;
    color_Map_y=int32((y(j)-min_y)/step_y)+1;
    c(j)=color_Map(color_Map_x,color_Map_y);
end
figure
scatter(x,y,5,'filled')
hold on  % ����Ϊ���õ�Ĵ�С
% image(linspace(min_x,max_x,NLevel+1),linspace(min_y,max_y,NLevel+1),color_Map);
max_h=max(max(max(x)),max(max(y)));
max_plot=max_h+5;
a=[0:max_plot];  %��x=y��ֱ��
b=a;
plot(a,b,'k--')  %��x=y��ֱ��
%�������ֱ��
a_fit=polyfit(x,y,1);
plot(a,a_fit(1).*a+a_fit(2),'m--')
% plot(a,a_fit(1).*a.^2+a_fit(2).*a+a_fit(3),'b--');
% plot(a,a_fit(1).*a.^4+a_fit(2).*a.^3+a_fit(3).*a.^2+a_fit(4).*a+a_fit(5),'b--');
xlabel('Lidar Forest Height  (m)','Color','k','fontsize',14,'FontName','΢���ź�'); % ���x����
ylabel('Estemated Forest Height  (m)','Color','k','fontsize',14,'FontName','΢���ź�'); % ���y����
% text(corner*2/3+3,3,['mean-Lidar-h=',num2str(mean_h)],'FontSize',9);  % �����ʵ���ߵ�ƽ��ֵ
text(1,33,[' RMSE=',num2str(RMSE),' m'],'FontSize',12);  % ���RMSEֵ
text(1,30,[' R  =',num2str(r2)],'FontSize',12);
text(1.5,31,['    2'],'FontSize',8);
a_fit=roundn(a_fit,-2);
% text(21,5,['---:y=',num2str(a_fit(1)),'x+',num2str(a_fit(2))],'FontSize',12);  
% colorbar
caxis([0,3])
axis([0,35,0,35]);
% set(gca,'FontSize',16);
% colormap jet
% colorbar;   %��ʾ��ɫ��
end





