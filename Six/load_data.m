function [PDH,Y1,Y2,Y3,theta,kz] = load_data(lines,yv,y1,y2,y3,s)
% 功能：目标规划函数，使用非线性最小二乘拟合
% 输入参数：h,sigma,fai0,m1,m2,m3,lb,ub,M,N
% 输出参数：矩阵 H,SIGMA,FAI0,M1,M2,M3

% 复相干
if yv == "PDH"
    PDH = freadbk('E:/ty/PolInSAR/FHI/data/0201_0207_FER_MLK_GSS/T6/cmplx_coh_PDHigh.bin',lines,'cpxfloat32');
end
if y1 == "HH"
    Y1 = freadbk('E:/ty/PolInSAR/FHI/data/0201_0207_FER_MLK_GSS/T6/cmplx_coh_HH.bin',lines,'cpxfloat32');
elseif y1 == "Opt1"
    Y1 = freadbk('E:/ty/PolInSAR/FHI/data/0201_0207_FER_MLK_GSS/T6/cmplx_coh_Opt1.bin',lines,'cpxfloat32');
end
if y2 == "HV"
    Y2 = freadbk('E:/ty/PolInSAR/FHI/data/0201_0207_FER_MLK_GSS/T6/cmplx_coh_HV.bin',lines,'cpxfloat32');
elseif y2 == "Opt2"
    Y2 = freadbk('E:/ty/PolInSAR/FHI/data/0201_0207_FER_MLK_GSS/T6/cmplx_coh_Opt2.bin',lines,'cpxfloat32');
end
if y3 == "VV"
    Y3 = freadbk('E:/ty/PolInSAR/FHI/data/0201_0207_FER_MLK_GSS/T6/cmplx_coh_VV.bin',lines,'cpxfloat32');
elseif y3 == "Opt3"
    Y3 = freadbk('E:/ty/PolInSAR/FHI/data/0201_0207_FER_MLK_GSS/T6/cmplx_coh_Opt3.bin',lines,'cpxfloat32');
end

% 其他参数
Inc = struct2array(load('./data/IncAngle0201.mat'));
Kz = struct2array(load('./data/Kz.mat'));
kz = squeeze(Kz(:,:,3));

% 若进行坡度校正
if s == 'y'
    Local_Inc = struct2array(load('./data/Local_IncAngle0201.mat'));
    Rngslope = Inc - Local_Inc;
    theta = Local_Inc;
    kz = kz.*sin(Inc)./sin(Local_Inc);
    kz = abs(kz);
else
    theta = Inc;
    kz = abs(kz);
end

end