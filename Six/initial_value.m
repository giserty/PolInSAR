function [h,sigma,fai0,m1,m2,m3] = initial_value(hh,sig,fai,m,n,YV,Y1,Y2,Y3)
% 生成初值的函数
% 输入参数：hv,sig,fai,PDH,Y1,Y2,Y3，个数可以是6/3/2
% 输出参数：h,sigma,fai0,m1,m2,m3

% 如果输入参数不够
if isnumeric(hh)
   h = ones(m,n) .* hh; 
end

if isnumeric(sig)
   sigma = ones(m,n) .* sig; 
end

if isnumeric(fai)
   fai0 = ones(m,n) .* fai; 
end

if nargin == 9
    m1 = abs((exp(1i.*fai).*YV-Y1)./(Y1-exp(1i.*fai)));
    m2 = abs((exp(1i.*fai).*YV-Y2)./(Y2-exp(1i.*fai)));
    m3 = abs((exp(1i.*fai).*YV-Y3)./(Y3-exp(1i.*fai)));
else
    m1 = ones(m,n) .* 0.01;
    m2 = ones(m,n) .* 0.01;
    m3 = ones(m,n) .* 0.01;
end

end