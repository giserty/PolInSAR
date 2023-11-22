function [h,gv]=SingleBaseline_fixuwAlgorith(ypdh,fai0,kz,thet,uw0)
%thet，为入射角文件
kz=abs(kz);
[row,col]=size(ypdh);
h=zeros(row,col);%树高
gv=zeros(row,col);%消光系数
uw=ones(row,col).*uw0;
%%
%构建二维查找表
%gg为消光系数，hvhv为树高
dethv=0.5;
hv=[0:dethv:30];
size_hv=max(size(hv));
detg=0.02;
g=[0:detg:0.5]'*0.23;
size_g=max(size(g));
hvhv=ones(size_g,1)*hv;
gg=g*ones(1,size_hv);

%%
parfor m=1:row
    for n=1:col
        %消光系数关系建模相干性理论值，地体幅度比关系观测想干性的值
        %disp([num2str(m) '_' num2str(n)]);
        %当地体幅度比为0时(u(w)=0),即只有植被能量贡献，只存在体散射
        %传统三阶段假设HV/PDH通道u(w)=0，可用yhv=e^(i*fai0)*y_v计算纯体去相干系数，可引入pdh替代hv
        y_v=(ypdh(m,n)*(uw(m,n)+1))/exp(fai0(m,n)*1i)-uw(m,n); %1 i代表复数i  
        %第三步，二维查表法解算植被参数解算
        %前半部分观测相干性，后半部分理论相干性，yvyv代表相干性差值
        yvyv=abs(y_v-2*gg.*(exp(2*gg.*hvhv./cos(thet(m,n))+kz(m,n)*hvhv*1i)-1)./((2*gg+kz(m,n)*cos(thet(m,n))*1i).*(exp(2*gg.*hvhv/cos(thet(m,n)))-1)));
        [min_m,min_n]=find(yvyv==min(min(yvyv))); %找到yvyv中最小值的索引   %列代表树高，行代表消光系数
        if isempty(min_n)==1
            h(m,n)=NaN;
            gv(m,n)=NaN;
        else
            h(m,n)=0+(min_n(1)-1)*dethv;
            gv(m,n)=(0+(min_m(1)-1)*detg);
        end
    end
end
end