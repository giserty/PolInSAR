% 画带密度的散点图
% x横坐标，y纵坐标
function plot_density_RMSE(true_h,TDX_h)

zero_value1=find(true_h==0);
true_h(zero_value1)=[];
TDX_h(zero_value1)=[];

if nargin==1   %仅有1个输入变量
    if size(true_h,1)
        true_h=true_h';    %转置为列向量
    end
    TDX_h=true_h;
    true_h=(1:length(TDX_h))';
else if nargin==2   %有2个输入变量
        if size(true_h,1)
            true_h=true_h';   %转置为列向量
        end
        if size(TDX_h,1)
            TDX_h=TDX_h';       %转置为列向量
        end
    end
end
N=length(true_h);   %数据长度
c=zeros(N,1);
max_x=max(true_h);min_x=min(true_h);   %搜索边界点
max_y=max(TDX_h);min_y=min(TDX_h);   %搜索边界点
NLevel=151;   %划分等级150份;分的等级越多，颜色变化的渐变越精细，即颜色条分辨率越高
color_Map=zeros(NLevel+1);
step_x=(max_x-min_x)/(NLevel-1);    % x轴步长
step_y=(max_y-min_y)/(NLevel-1);    % y轴步长 
for j=1:N
    color_Map_x=int32((true_h(j)-min_x)/step_x)+1;
    color_Map_y=int32((TDX_h(j)-min_y)/step_y)+1;
    color_Map(color_Map_x,color_Map_y)=color_Map(color_Map_x,color_Map_y)+1;
end
for j=1:N
    color_Map_x=int32((true_h(j)-min_x)/step_x)+1;
    color_Map_y=int32((TDX_h(j)-min_y)/step_y)+1;
    c(j)=color_Map(color_Map_x,color_Map_y);
end
figure;
scatter(true_h,TDX_h,6,c,'filled');hold on; % 数字是设置点的大小为10
box on; % 打开上侧和右侧的坐标轴

set(gca,'FontSize',14);
colormap(jet);   %查阅colormap函数改变颜色变化趋势
colorbar;   %显示颜色条

max_h=max(true_h);
diff=true_h-TDX_h;
p=find(abs(diff)>max_h);
diff(p)=[];
true_h(p)=[];
TDX_h(p)=[];

% 计算RMSE的值，空值不计入统计内
[m,n]=size(TDX_h);
num=0;
temp=0;
for i=1:m
    for j=1:n
        if diff(i,j)~=0
        temp=diff(i,j)^2+temp;
        num=num+1;
        end
    end
end
RMSE=sqrt(temp/num);
RMSE=roundn(RMSE,-2); % 保留2位小数

cor=corrcoef(true_h,TDX_h);
cor(1,2)=roundn(cor(1,2),-2); % 保留2位小数

max1=max(true_h);
max2=max(TDX_h);
max3=max(max1,max2);
max4=fix(max3/5)+1;
corner=max4*5;
plot([0,corner],[0,corner],'k--','LineWidth',2);  %画一条y=x的直线
box on; 

effct_h95plot=true_h(true_h>0);
plot_count=length(effct_h95plot);

xlabel('LiDAR Forest Height(m)','Color','k','fontsize',18,'FontName','Arial');    % 添加x轴标记
ylabel('Estimated  Forest  Height(m)','Color','k','fontsize',18,'FontName','Arial');   % 添加y轴标记
text(corner*3/6,corner*5/6+3.6,['R^2=',num2str(cor(1,2))],'FontSize',20,'FontName','Arial');     % 添加相关系数值
text(corner*3/6,corner*5/6+1.8,['RMSE=',num2str(RMSE),' m'],'FontSize',20','FontName','Arial');        % 添加RMSE值
text(corner*3/6,corner*5/6,[num2str(plot_count),' Plots'],'FontSize',20,'FontName','Arial');     % 添加相关系数值

pbaspect([1 1 1]);
set(gca,'xtick',0:5:corner,'ytick',0:5:corner,'FontSize',17);

end
