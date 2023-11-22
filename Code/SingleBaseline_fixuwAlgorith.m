function [h,gv]=SingleBaseline_fixuwAlgorith(ypdh,fai0,kz,thet,uw0)
%thet��Ϊ������ļ�
kz=abs(kz);
[row,col]=size(ypdh);
h=zeros(row,col);%����
gv=zeros(row,col);%����ϵ��
uw=ones(row,col).*uw0;
%%
%������ά���ұ�
%ggΪ����ϵ����hvhvΪ����
dethv=0.5;
hv=[0:dethv:30];
size_hv=max(size(hv));
detg=0.02;
g=[0:detg:0.5]'*0.23;
size_g=max(size(g));
hvhv=ones(size_g,1)*hv;
gg=g*ones(1,size_hv);

%%
parfor m=1:row
    for n=1:col
        %����ϵ����ϵ��ģ���������ֵ��������ȱȹ�ϵ�۲�����Ե�ֵ
        %disp([num2str(m) '_' num2str(n)]);
        %��������ȱ�Ϊ0ʱ(u(w)=0),��ֻ��ֲ���������ף�ֻ������ɢ��
        %��ͳ���׶μ���HV/PDHͨ��u(w)=0������yhv=e^(i*fai0)*y_v���㴿��ȥ���ϵ����������pdh���hv
        y_v=(ypdh(m,n)*(uw(m,n)+1))/exp(fai0(m,n)*1i)-uw(m,n); %1 i������i  
        %����������ά�������ֲ����������
        %ǰ�벿�ֹ۲�����ԣ���벿����������ԣ�yvyv��������Բ�ֵ
        yvyv=abs(y_v-2*gg.*(exp(2*gg.*hvhv./cos(thet(m,n))+kz(m,n)*hvhv*1i)-1)./((2*gg+kz(m,n)*cos(thet(m,n))*1i).*(exp(2*gg.*hvhv/cos(thet(m,n)))-1)));
        [min_m,min_n]=find(yvyv==min(min(yvyv))); %�ҵ�yvyv����Сֵ������   %�д������ߣ��д�������ϵ��
        if isempty(min_n)==1
            h(m,n)=NaN;
            gv(m,n)=NaN;
        else
            h(m,n)=0+(min_n(1)-1)*dethv;
            gv(m,n)=(0+(min_m(1)-1)*detg);
        end
    end
end
end