

clear;
clc;
lines=11099;

addpath 'F:\4.树高\BioSAR_2008L新生测试数据\0201_0207_FER_BOX\T6';
pdh=freadbk('cmplx_coh_PDHigh.bin',lines,'cpxfloat32');
pd1=freadbk('cmplx_coh_PDLow.bin',lines,'cpxfloat32');
opt1=freadbk('cmplx_coh_Opt1.bin',lines,'cpxfloat32');
opt2=freadbk('cmplx_coh_Opt2.bin',lines,'cpxfloat32');
opt3=freadbk('cmplx_coh_Opt3.bin',lines,'cpxfloat32');

HH=freadbk('cmplx_coh_HH.bin',lines,'cpxfloat32');
HV=freadbk('cmplx_coh_HV.bin',lines,'cpxfloat32');
VV=freadbk('cmplx_coh_VV.bin',lines,'cpxfloat32');
HHmVV=freadbk('cmplx_coh_HHmVV.bin',lines,'cpxfloat32');
HHpVV=freadbk('cmplx_coh_HHpVV.bin',lines,'cpxfloat32');

load('CHM_Lidar.mat');
load('Kz.mat');
kz3=Kz(:,:,3);
figure,imagesc(CHM_sar),title('LiDAR Height');colormap jet;colorbar;caxis([0,30]);axis off;set(gca,'FontSize',24);
figure,imagesc(kz3),title('kz3');colormap jet;colorbar;caxis([0,30]);axis off;set(gca,'FontSize',24);

%% 低矮（6m）：a=1136;b=726; midlle（15m）：a=2850;b=1101; high(30m): a=1508;b=504;
a=907;
b=801;
yv1=pdh(a,b);yv2=pd1(a,b);
yv3=opt1(a,b);yv4=opt2(a,b);yv5=opt3(a,b);
yv6=HH(a,b);yv7=HV(a,b);yv8=VV(a,b);
yv9=HHmVV(a,b);yv10=HHpVV(a,b);

figure;
polarplot(angle(yv1),abs(yv1),'k.','MarkerSize',24,'linewidth',3);hold on;
polarplot(angle(yv2),abs(yv2),'m.','MarkerSize',24,'linewidth',3);hold on;
% polarplot(angle(yv3),abs(yv3),'rp','MarkerSize',10,'linewidth',3);hold on; % pentagram,五角星
% polarplot(angle(yv4),abs(yv4),'gp','MarkerSize',10,'linewidth',3);hold on;
% polarplot(angle(yv5),abs(yv5),'bp','MarkerSize',10,'linewidth',3);hold on;
polarplot(angle(yv6),abs(yv6),'rs','MarkerSize',12,'linewidth',3);hold on; % square，正方形
polarplot(angle(yv7),abs(yv7),'gs','MarkerSize',12,'linewidth',3);hold on;
polarplot(angle(yv8),abs(yv8),'bs','MarkerSize',12,'linewidth',3);hold on;
polarplot(angle(yv9),abs(yv9),'k+','MarkerSize',8,'linewidth',3);hold on;
polarplot(angle(yv10),abs(yv10),'m+','MarkerSize',8,'linewidth',3);hold on;
set(gca,'FontSize',24);rlim([0,1]);
legend1 = legend('PD High','PD Low','HH','HV','VV','HH+HV','HH-HV');
set(legend1,'Position',[0.867530682823755 0.174509806921279 0.112419698172004 0.210457510340447]);






