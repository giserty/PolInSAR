% �����ܶȵ�ɢ��ͼ
% x�����꣬y������
function plot_density_RMSE(true_h,TDX_h)

zero_value1=find(true_h==0);
true_h(zero_value1)=[];
TDX_h(zero_value1)=[];

if nargin==1   %����1���������
    if size(true_h,1)
        true_h=true_h';    %ת��Ϊ������
    end
    TDX_h=true_h;
    true_h=(1:length(TDX_h))';
else if nargin==2   %��2���������
        if size(true_h,1)
            true_h=true_h';   %ת��Ϊ������
        end
        if size(TDX_h,1)
            TDX_h=TDX_h';       %ת��Ϊ������
        end
    end
end
N=length(true_h);   %���ݳ���
c=zeros(N,1);
max_x=max(true_h);min_x=min(true_h);   %�����߽��
max_y=max(TDX_h);min_y=min(TDX_h);   %�����߽��
NLevel=151;   %���ֵȼ�150��;�ֵĵȼ�Խ�࣬��ɫ�仯�Ľ���Խ��ϸ������ɫ���ֱ���Խ��
color_Map=zeros(NLevel+1);
step_x=(max_x-min_x)/(NLevel-1);    % x�Ჽ��
step_y=(max_y-min_y)/(NLevel-1);    % y�Ჽ�� 
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
scatter(true_h,TDX_h,6,c,'filled');hold on; % ���������õ�Ĵ�СΪ10
box on; % ���ϲ���Ҳ��������

set(gca,'FontSize',14);
colormap(jet);   %����colormap�����ı���ɫ�仯����
colorbar;   %��ʾ��ɫ��

max_h=max(true_h);
diff=true_h-TDX_h;
p=find(abs(diff)>max_h);
diff(p)=[];
true_h(p)=[];
TDX_h(p)=[];

% ����RMSE��ֵ����ֵ������ͳ����
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
RMSE=roundn(RMSE,-2); % ����2λС��

cor=corrcoef(true_h,TDX_h);
cor(1,2)=roundn(cor(1,2),-2); % ����2λС��

max1=max(true_h);
max2=max(TDX_h);
max3=max(max1,max2);
max4=fix(max3/5)+1;
corner=max4*5;
plot([0,corner],[0,corner],'k--','LineWidth',2);  %��һ��y=x��ֱ��
box on; 

effct_h95plot=true_h(true_h>0);
plot_count=length(effct_h95plot);

xlabel('LiDAR Forest Height(m)','Color','k','fontsize',18,'FontName','Arial');    % ���x����
ylabel('Estimated  Forest  Height(m)','Color','k','fontsize',18,'FontName','Arial');   % ���y����
text(corner*3/6,corner*5/6+3.6,['R^2=',num2str(cor(1,2))],'FontSize',20,'FontName','Arial');     % ������ϵ��ֵ
text(corner*3/6,corner*5/6+1.8,['RMSE=',num2str(RMSE),' m'],'FontSize',20','FontName','Arial');        % ���RMSEֵ
text(corner*3/6,corner*5/6,[num2str(plot_count),' Plots'],'FontSize',20,'FontName','Arial');     % ������ϵ��ֵ

pbaspect([1 1 1]);
set(gca,'xtick',0:5:corner,'ytick',0:5:corner,'FontSize',17);

end
