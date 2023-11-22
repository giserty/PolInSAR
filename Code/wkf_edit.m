%% 导入数据
clear
clc
filename='D:\BioSAR2008L';
lines=11099;
ypdh=freadbk(strcat(filename,'\0201_0207_FER_MLK_BOX\T6\cmplx_coh_PDHigh.bin'),lines,'cpxfloat32');
ypdl=freadbk(strcat(filename,'\0201_0207_FER_MLK_BOX\T6\cmplx_coh_PDLow.bin'),lines,'cpxfloat32');
cmplx_coh_HV=freadbk(strcat(filename,'\0201_0207_FER_MLK_BOX\T6\cmplx_coh_HV.bin'),lines,'cpxfloat32');
load('D:\BioSAR2008L\data\Kz');kz=squeeze(Kz(:,:,3));

%% 计算地表相位
flag=1; 
if mean(kz(:))<0  
    flag=-1;
end

[ground_phase]=traditional_pdh_pdl_line(ypdh,ypdl,flag);
%[ground_phase,y_v]=MCD_PD_TLS_linear_fitting_ACDkz(ypdh,ypdl,flag);

figure,imagesc(ground_phase),title('ground phase');colormap jet;colorbar;caxis([-pi,pi]);
save('ground_phase.mat','ground_phase');

%% 查找表
load('D:\BioSAR2008L\三阶段\ground_phase');
load('D:\BioSAR2008L\data\IncAngle0201');
load('D:\BioSAR2008L\data\Local_IncAngle0201');
Rngslope=IncAngle0201-Local_IncAngle0201; % 坡度角=名义入射角-局部入射角
kz=kz.*sin(IncAngle0201)./sin(Local_IncAngle0201);
kz=abs(kz);

% [h,gv,y_v]=classic_SB_LUT_calculate_height(ypdh,ypdl,ground_phase,kz,IncAngle0201);
% [h,gv,y_v]=classic_SB_LUT_calculate_height(ypdh,ypdl,ground_phase,kz,Local_IncAngle0201);
% [h,gv,y_v]=Copy_of_classic_SB_LUT_calculate_height(ypdh,ypdl,ground_phase,kz,IncAngle0201);
% [h,gv,y_v]=Copy_of_classic_SB_LUT_calculate_height(ypdh,ypdl,ground_phase,kz,Local_IncAngle0201);
[h,gv]=SingleBaseline_fixuwAlgorith(ypdh,ground_phase,kz,Local_IncAngle0201,0);

% figure();
% t=fspecial('average',3);
% h_filter = filter2(t,h);imagesc(h_filter);

h1=h./cos((Rngslope));
h1(h1>30)=30;
h1(h1<0)=0;
save('h1.mat','h1'); 







