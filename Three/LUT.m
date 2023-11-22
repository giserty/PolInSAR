function [h_min,sigma_min] = LUT(fai0,dethv,detg,yv,m,n,theta,kz)
% 功能：
% 输入：
% 输出：
% global theta
% global kz
sigma = (0.001:detg:0.5)'.*0.23;
size_sigma = max(size(sigma));
h = 0.01:dethv:30.01;
size_h = max(size(h));
sigma = sigma*ones(1,size_h);
h = ones(size_sigma,1)*h;

parfor i = 1 : m
    for j = 1 : n
        try
            p = 2 * sigma / cos(theta(i,j));
            p1 = p + kz(i,j)*1i;
            y = exp(fai0(i,j)*1i) .* (p./p1) .* (exp(p1.*h)-1)./(exp(p.*h)-1);
            F = abs(yv(i,j) - y);
            F_min = min(F(:));
            [row, col] = find(F == F_min);
            h_min(i,j) = 0.01+(col-1)*dethv;
            sigma_min(i,j) = 0.001+(row-1)*detg;
        catch
            sigma_min(i,j) = NaN;
            h_min(i,j) = NaN;
        end
    end
end

end