function f = Obj_fun(x,Inc,kz,Opt1,Opt2,Opt3)
    % h,sigma,fai0,m1,m2,m3
    % 纯体去相干系数
    yv = 2.*x(2).*(exp(2.*x(2).*x(1)./cos(Inc)+1i.*kz.*x(1))-1) ./ ...
        (2.*x(2)+1i.*kz.*cos(Inc)) ./ (exp(2.*x(2).*x(1)./cos(Inc))-1);
    % 复相干系数
    y1 = exp(1i.*x(3)) .* (yv+x(4)) ./ (1+x(4));
    y2 = exp(1i.*x(3)) .* (yv+x(5)) ./ (1+x(5));
    y3 = exp(1i.*x(3)) .* (yv+x(6)) ./ (1+x(6));
%     f = sqrt(abs(y1-Opt1)).^2 + (abs(y2-Opt2)).^2 + (abs(y3-Opt3)).^2);
    f = abs(y1-Opt1) + abs(y2-Opt2) + abs(y3-Opt3);
end