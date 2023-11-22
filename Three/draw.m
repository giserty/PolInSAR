
img = h_min;

img(img>35)=35;
img(img<0)=0;

figure,imagesc(img);
colormap jet;
colorbar;
set(gca,'FontSize',font_size);
pbaspect([1 2 1]);
