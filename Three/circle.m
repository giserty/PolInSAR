function circle(i,j)

% 极坐标图
polarplot(Opt1(i,j),'ro','linewidth',1)
hold on
polarplot(Opt2(i,j),'bd','linewidth',1)
polarplot(Opt3(i,j),'yp','linewidth',1)
polarplot(PDH(i,j),'gs','linewidth',1)
polarplot(PDL(i,j),'k*','linewidth',1)
% 绘制直线
X = -pi : 0.01 : pi;
Y = b ./ (sin(X)-k.*cos(X));
polarplot(X,Y,'k-')
legend('Opt1','Opt2','Opt3','PDH','PDL');
rlim([0 1])

end