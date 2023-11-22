function [mtlk_A2]=multilook(A,p,q)
% ���ܣ�������ľ�����ж���
% ���������AΪ�����ӵľ���p��q�Ƕ��ӱ���
% ���������mtlk_A2Ϊ���Ӻ���¾���

[m,n]=size(A);
pp=fix(m/p);
qq=fix(n/q);
end_p=m-pp*p;
end_q=n-qq*q;

mtlk_A1=zeros(pp,n);
for i=1:p
    mtlk_A1=mtlk_A1+A(i:p:end-end_p,:)/p;
end

mtlk_A2=zeros(pp,qq);
for j=1:q
    mtlk_A2=mtlk_A2+mtlk_A1(:,j:q:end-end_q)/q;
end

% % 2:2����
% fltcc2=fltcc2(1:2:end-1,:)/2+fltcc2(2:2:end-1,:)/2;
% fltcc2=fltcc2(:,1:2:end-1)/2+fltcc2(:,2:2:end-1)/2;
% % 4:2����
% kz=kz(:,1:2:1500)/2+kz(:,2:2:1500)/2;
% kz=kz(1:4:end-1,:)/4+kz(2:4:end-1,:)/4+kz(3:4:end-1,:)/4+kz(4:4:end-1,:)/4;

% figure,imagesc(A),title('A');colormap jet;colorbar;
% figure,imagesc(mtlk_A2),title('multilook-A');colormap jet;colorbar;
end