lines = 11099;
work_path = 'E:/ty/PolInSAR/FHI/data/0201_0207_FER_MLK_GSS/T6/cmplx_coh_';
PDH = freadbk([work_path,'PDHigh.bin'],lines,'cpxfloat32');
PDL = freadbk([work_path,'PDLow.bin'],lines,'cpxfloat32');
HH = freadbk([work_path,'HH.bin'],lines,'cpxfloat32');
HV = freadbk([work_path,'HV.bin'],lines,'cpxfloat32');
VV = freadbk([work_path,'VV.bin'],lines,'cpxfloat32');
HHpVV = freadbk([work_path,'HHpVV.bin'],lines,'cpxfloat32');
HHmVV = freadbk([work_path,'HHmVV.bin'],lines,'cpxfloat32');



img = H;
figure,imagesc(img);
colormap jet;
colorbar;
set(gca,'FontSize',font_size);
pbaspect([1 2 1]);