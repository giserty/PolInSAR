%% SAģ���˻�
tic
clear; clc

lines = 11099;
cols = 1510;
Opt1 = struct2array(load('Opt1.mat'));
PDLow = struct2array(load('PDLow.mat'));
hv = struct2array(load('hv.mat'));   % ɭ�ָ߶�
pha0 = angle(PDLow);   % �ر���λ
sigma = ones(lines,cols) .* 0.1;   % ����ϵ��
m1 = ones(lines,cols) .* 0.1;  % ������ȱ�1
% m2 = 0.2;  % ������ȱ�2
% m3 = 0.3;  % ������ȱ�3

%% ������ʼ��
narvs = 4;  % ��������
T0 = 1000;   % ��ʼ�¶�
T = T0; % �������¶Ȼᷢ���ı䣬��һ�ε���ʱ�¶Ⱦ���T0
maxgen = 100;  % ����������
Lk = 30;  % ÿ���¶��µĵ�������
alfa = 0.95;  % �¶�˥��ϵ��
theta = sqrt(2)/2;  % �����
kz = 1;  % ��ֱ����
% sigma = 0.1;  % ����ϵ��
% hv = 10;  % ����
% pha = 0.5;  % ��λ
% m1 = 0.1;  % ������ȱ�
sigma_lb = zeros(lines,cols);
sigma_ub = ones(lines,cols);
h_lb = ones(lines,cols) .* 10;
h_ub = ones(lines,cols) .* 30;
pha_lb = ones(lines,cols) .* 0.5;
pha_ub = ones(lines,cols) .* 5;
m_lb = ones(lines,cols) .* 0.1;
m_ub = ones(lines,cols);
x_lb = [sigma_lb h_lb pha_lb m_lb]; % x���½�
x_ub = [sigma_ub h_ub pha_ub m_ub]; % x���Ͻ�

%%  ���ɳ�ʼ��
% x0 = zeros(1,narvs);  % ��ʼ0����
% for i = 1: narvs
%     x0(i) = x_lb(i) + (x_ub(i)-x_lb(i))*rand(1);  % �����ʼ���������ڵ�λ���ڶ�������
% end
x0 = [hv,pha0,sigma,m1];
y0 = abs(Opt1 - Obj_fun3(x0));  % ���㵱ǰ��ĺ���ֵ

%% ����һЩ�����м���̵���������������
min_y = y0;     % ��ʼ���ҵ�����ѵĽ��Ӧ�ĺ���ֵΪy0

%% ģ���˻����
for iter = 1 : maxgen  % ��ѭ��, ��������õ���ָ������������
    for i = 1 : Lk  %  ��ѭ������ÿ���¶��¿�ʼ����
        y = randn(lines,cols*4);  % ����m��n�е�N(0,1)�����
        z = y ./ sqrt(sum(sum(y.^2))); % �����½�Ĳ����������z
        x_new = x0 + z.*T; % �����½�Ĳ����������x_new��ֵ
        % �������½��λ�ó����˶����򣬾Ͷ�����е���
%         for j = 1: lines*cols*4
%             if x_new(j) < x_lb(j)
%                 r = rand(1);
%                 x_new(j) = r.*x_lb(j)+(1-r).*x0(j);
%             elseif x_new(j) > x_ub(j)
%                 r = rand(1);
%                 x_new(j) = r.*x_ub(j)+(1-r).*x0(j);
%             end
%         end
        x1 = x_new;    % ���������x_new��ֵ���½�x1
        y1 = abs(Opt1 - Obj_fun3(x1));  % �����½�ĺ���ֵ
        if y1 < y0    % ����½⺯��ֵС�ڵ�ǰ��ĺ���ֵ
            x0 = x1; % ���µ�ǰ��Ϊ�½�
            y0 = y1;
        else
            p = exp(-(y0 - y1)/T); % ����Metropolis׼�����һ������
            if rand(1) < p   % ����һ���������������ʱȽϣ�����������С���������
                x0 = x1; % ���µ�ǰ��Ϊ�½�
                y0 = y1;
            end
        end
        % �ж��Ƿ�Ҫ�����ҵ�����ѵĽ�
        if y0 < min_y  % �����ǰ����ã��������и���
            min_y = y0;  % ������С��y
            best_x = x0;  % �����ҵ�����õ�x
        end
    end
    T = alfa*T;   % �¶��½�
end

%%
toc
disp(['����ʱ��: ',num2str(toc)]);
