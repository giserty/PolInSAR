%  画带密度的散点图
% x横坐标，y纵坐标
function plot_density(x,y,RMSE,r2,accuracy)
if nargin==1   %仅有1个输入变量
    if size(x,1)
        x=x';    %转置为列向量
    end
    y=x;
    x=(1:length(y))';
else
    if nargin==2   %有2个输入变量
        if size(x,1)
            x=x';   %转置为列向量
        end
        if size(y,1)
            y=y';       %转置为列向量
        end
    end
end
N=length(x);     %数据长度
c=zeros(N,1);
max_x=max(x);
min_x=min(x);    %搜索边界点
max_y=max(y);
min_y=min(y);    %搜索边界点
NLevel=300;      %划分等级150份;分的等级越多，颜色变化的渐变越精细，即颜色条分辨率越高
color_Map=zeros(NLevel+1);
step_x=(max_x-min_x)/(NLevel-1);  % x轴步长
step_y=(max_y-min_y)/(NLevel-1);    % y轴步长 
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
hold on  % 数字为设置点的大小
% image(linspace(min_x,max_x,NLevel+1),linspace(min_y,max_y,NLevel+1),color_Map);
max_h=max(max(max(x)),max(max(y)));
max_plot=max_h+5;
a=[0:max_plot];  %加x=y的直线
b=a;
plot(a,b,'k--')  %加x=y的直线
%绘制拟合直线
a_fit=polyfit(x,y,1);
plot(a,a_fit(1).*a+a_fit(2),'m--')
% plot(a,a_fit(1).*a.^2+a_fit(2).*a+a_fit(3),'b--');
% plot(a,a_fit(1).*a.^4+a_fit(2).*a.^3+a_fit(3).*a.^2+a_fit(4).*a+a_fit(5),'b--');
xlabel('Lidar Forest Height  (m)','Color','k','fontsize',14,'FontName','微软雅黑'); % 添加x轴标记
ylabel('Estemated Forest Height  (m)','Color','k','fontsize',14,'FontName','微软雅黑'); % 添加y轴标记
% text(corner*2/3+3,3,['mean-Lidar-h=',num2str(mean_h)],'FontSize',9);  % 添加真实树高的平均值
text(1,33,[' RMSE=',num2str(RMSE),' m'],'FontSize',12);  % 添加RMSE值
text(1,30,[' R  =',num2str(r2)],'FontSize',12);
text(1.5,31,['    2'],'FontSize',8);
a_fit=roundn(a_fit,-2);
% text(21,5,['---:y=',num2str(a_fit(1)),'x+',num2str(a_fit(2))],'FontSize',12);  
% colorbar
caxis([0,3])
axis([0,35,0,35]);
% set(gca,'FontSize',16);
% colormap jet
% colorbar;   %显示颜色条
end





