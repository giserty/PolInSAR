function [fai0,p] = PDMCD_fai0(PDH,PDL,Opt1,Opt2,Opt3,m,n)
% 相位幅度计算
Opt1_pha = angle(Opt1); Opt1_mod = abs(Opt1); 
Opt1_re = real(Opt1); Opt1_im = imag(Opt1);
Opt2_pha = angle(Opt2); Opt2_mod = abs(Opt2); 
Opt2_re = real(Opt2); Opt2_im = imag(Opt2);
Opt3_pha = angle(Opt3); Opt3_mod = abs(Opt3); 
Opt3_re = real(Opt3); Opt3_im = imag(Opt3);
PDH_pha = angle(PDH); PDH_mod = abs(PDH); 
PDH_re = real(PDH); PDH_im = imag(PDH);
PDL_pha = angle(PDL); PDL_mod = abs(PDL); 
PDL_re = real(PDL); PDL_im = imag(PDL);

for i = 1 : m
    for j = 1 : n
        % 直线拟合
        x = [Opt1_re(i,j) Opt2_re(i,j) Opt3_re(i,j) PDH_re(i,j) PDL_re(i,j)];
        y = [Opt1_im(i,j) Opt2_im(i,j) Opt3_im(i,j) PDH_im(i,j) PDL_im(i,j)];
        f = polyfit(x,y,1);
        k = f(1);
        b = f(2);
        % 计算交点
        aa = 1+k*k;
        bb = 2*k*b;
        cc = b*b-1;
        xx = roots([aa bb cc]);
        yy = xx.*k + b;
        zz = xx + yy.*1i;
        pha = angle(zz);  % 某像素2个相位
        % 筛选地表相位
        dist1_pdh = abs(zz(1)-PDH(i,j));
        dist1_pdl = abs(zz(1)-PDL(i,j));
        if dist1_pdh > dist1_pdl
            fai0(i,j) = pha(1);
            p(i,j) = zz(1);
        else
            fai0(i,j) = pha(2);
            p(i,j) = zz(2);
        end
    end
end

end