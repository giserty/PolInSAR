function [h_min,mu_min] = LUT_m(fai0,dethv,detmu,ywv,m,n,theta,kz,sig)
% 功能：
% 输入：ground 地表相干,dethv h步长,detg sigma步长,yv 体散射相干,m,n,theta,kz,sig
% 输出：h_min 树高,mu_min 地体幅度比

sigma = sig * 0.23;
mu = (0.01:detmu:1)';
size_m = max(size(mu));
h = 0.01:dethv:30.01;
size_h = max(size(h));
mu = mu*ones(1,size_h);
h = ones(size_m,1)*h;

parfor i = 1 : m
    for j = 1 : n
        try
            y_v = (ywv(i,j).*(mu+1))./exp(fai0(i,j)*1i)-mu; %1 i代表复数i  
            %第三步，二维查表法解算植被参数解算
            %前半部分观测相干性，后半部分理论相干性，yvyv代表相干性差值
            F = abs(y_v-2*sigma.*(exp(2*sigma.*h./cos(theta(i,j))+kz(i,j)*h*1i)-1)./((2*sigma+kz(i,j)*cos(theta(i,j))*1i).*(exp(2*sigma.*h/cos(theta(i,j)))-1)));
            F_min = min(F(:));
            [row, col] = find(F == F_min);
            h_min(i,j) = 0.01+(col-1)*dethv;
            mu_min(i,j) = 0.01+(row-1)*detmu;
        catch
            mu_min(i,j) = NaN;
            h_min(i,j) = NaN;
        end
    end
end

end