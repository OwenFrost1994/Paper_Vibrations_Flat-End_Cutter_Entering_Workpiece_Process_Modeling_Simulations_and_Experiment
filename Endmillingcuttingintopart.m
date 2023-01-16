clear;
clc;
close all;
tic
%%����ֱ�߽��빤�����̵�������ģ��
%%�������빤���������ž�������ı仯������������֮�仯�����ﲻ����Ҫ���㵶�ߵ���ת�Ƕȣ���Ҫ���㵶�ߵ�λ��
%%ģ������Ǵӵ��߿�ʼ�ƶ��������빤��һ���ľ��룬���ľ��빤��ǰ��S1��mm����֮������S2��mm��
%%�����ߵ���Ϊ�������֣�����s1������s2���ȶ�ϳ��s3
S1=7;%%��λ��mm
S2=1;%%��λ��mm

%%���߲���
D=12;%%���߰뾶
N=2;%%���߳���
B=pi/6;%%����������
Cp=2*pi/N;%%�ݼ��
dl=0.02;

%%���߸նȲ���%%
kx=1.74e7;%%x����ն�
ky=1.74e7;%%y����ն�
cx=116.4931;%%x��������
cy=116.4931;%%y��������
mx=0.2198;%%x��������
my=0.2198;%%y��������

%%���ϲ���%%����AL7075%%
% Ktc=951.751;%%���������ϵ��
% Kte=14.0371;%%�����п���ϵ��
% Krc=608.561;%%���������ϵ��
% Kre=16.5002;%%�����п���ϵ��
% Kac=288.478;%%���������ϵ��
% Kae=1.25118;%%�����п���ϵ��
%%���ϲ���%%����GH909%%
% Ktc=3250.0332;%%���������ϵ��
% Kte=90.8541;%%�����п���ϵ��
% Krc=1200.5036;%%���������ϵ��
% Kre=115.1116;%%�����п���ϵ��
% Kac=251.55;%%���������ϵ��
% Kae=52.9445;%%�����п���ϵ��
%%���ϲ���%%����45��%%
% Ktc=1807.5725;%%���������ϵ��
% Kte=36.8380;%%�����п���ϵ��
% Krc=551.0207;%%���������ϵ��
% Kre=26.2109;%%�����п���ϵ��
% Kac=261.6069;%%���������ϵ��
% Kae=-29.9157;%%�����п���ϵ��
Ktc=2129.9;%%���������ϵ��
Kte=28.6524;%%�����п���ϵ��
Krc=946.7411;%%���������ϵ��
Kre=18.1387;%%�����п���ϵ��
Kac=261.6069;%%���������ϵ��
Kae=-29.9157;%%�����п���ϵ��

%%�ӹ�����%%
%%ϳ����ʽ:˳ϳ%%
Cm=1;%%ϳ����ʽ��˳ϳΪ1����ϳΪ0
S=1000;%%����ת��
f=240;%%�����ٶ�
fs=10000;%%����Ƶ��
ap=1.5;%%���������λmm��
ae=3;%%���������λmm��

%%������������%%
R=D/2;%%���߰뾶
kb=(2*tan(B))/D;%%k�¼���
fe=f/(N*S);%%feed every tooth
w=2*pi*S/60;%%���߽��ٶ�
T=2*pi/w;%%��������
Nc=floor(60*fs/S);%%һ�������ڵĲ��������
%%��ȫ����֮��������г��Ƕ�
if Cm==1%%˳ϳ
    Cst=pi-acos((R-ae)/R);%%�����
    Cex=pi;%%�г���
else%%��ϳ
    Cst=0;%%�����
    Cex=acos((R-ae)/R);%%�г���
end
Cs=0;%%��ʼ�Ƕ�
Dt=T/Nc;%%ʱ�䲽��
DC=Dt*w;%%�Ƕ�����
Ca=Cs;%%��ʼ�Ƕ�

%%ģ�����ĸ����׶�s1��s2��s3
if Cm ==1
    if ae<=R
        s1=S1-R*sin(Cst);
    else
        s1=S1-R;
    end
else
    if ae<=R
        s1=S1-R*sin(Cex);
    else
        s1=S1-R;
    end 
end

s2=S1;

s3=S1+S2;

%%ģ���ʱ�򳤶�
Ts=60*(S1+S2)/f;%%��λ��s
%%�ܵĲ��������
Ns=Ts*fs;

%%���ִ洢��Ԫ%%
Fx=0;
Fy=0;
Fz=0;%%���������ۼӱ���
F=zeros(3,Ns);%%�洢���������������
Xx=zeros(1,Ns);%%������x�������λ��
Xy=zeros(1,Ns);%%������y�������λ��
Vx=zeros(1,Ns);%%������x��������ٶ�
Vy=zeros(1,Ns);%%������y��������ٶ�
Dx=zeros(1,Ns);%%������x�����λ�Ʊ仯����fv�ĵ�һ������
Dy=zeros(1,Ns);%%������y�����λ�Ʊ仯����fv�ĵڶ�������
FA=zeros(1,Ns);
apx=zeros(1,floor(Ns/Nc)*N);
apy=zeros(1,floor(Ns/Nc)*N);
Cstd=zeros(1,Ns);%%��¼��̬�����г��Ƕ�
Cexd=zeros(1,Ns);
%%��������ģ����Ĵ�ѭ��
%%�����ͬʱҪ����㵶�ߵ�ʵʱ�������г��Ƕ�
for i=1:1:Ns
    s=i*Dt*f/60;
    if s<=s1
        Cstd(i)=0;%%��̬�����
        Cexd(i)=0;%%��̬�г���
    else if s<=s2
            if ae<=R
                Cstd(i)=Cst;%%��ʱ��������ǿ϶���
                Cexd(i)=pi-asin((R*sin(Cst)-s+s1)/R);
            else
                Cexd(i)=pi-asin((R-s+s1)/R);
                if s<=s1+R-R*sin(Cst)
                    Cstd(i)=asin((R-s+s1)/R);
                else
                    Cstd(i)=Cst;
                end
            end
        else
            Cstd(i)=Cst;%%���һ���ǲ�ϳ
            Cexd(i)=Cex;
        end
    end
    
    Ca=Ca+DC;%%΢Ԫ�Ƕȵ��Ӽ��㵶�ߵ�ת����
    if Ca>=2*pi%%���ǵ��߶����ת���ڣ��ۼӵĵ��߽Ƕȳ���һ�ܾͼ�ȥһ��2��
        Ca=Ca-2*pi;
    else
    end
    
    %%���㸽�ӽ�����
    if i<=(Nc/N)%%����ǵ�һ���������е�ʱ����м�ϱ���û��ǰһ�������г��Ĳ��Ʊ��棬��ʱ�൱�ڵ������񶯵����
        Dx(i)=1000*(Xx(i)-0);%%ע�⵽����ѧ�����������е�λΪ�����Ƶ�λ������õ�λ�Ƶ�λΪm����ÿ�ݽ������ĵ�λ��mm
        Dy(i)=1000*(Xy(i)-0);
    else
        q=floor(i*N/Nc);
        if (i-q*Nc/N)==0
            q=q-1;
        else
        end
        for n=1:1:q
            apx(n)=Xx(i-n*Nc/N);
            apy(n)=Xy(i-n*Nc/N);
        end
        Dx(i)=1000*(Xx(i)-min(apx(1:q)));
        Dy(i)=1000*(Xy(i)-min(apy(1:q)));
    end
    
        %���Ӽ��������ݵ���ѭ��
        for j=1:1:N
            C=Ca-(j-1)*Cp;%%���Ƕ�ݴ��ڵĳݼ���ͺ�
            if C<0
                C=C+2*pi;%%������ݽǶ�С������ת�����ĽǶ�
            else
            end
            %%΢Ԫ���Ӽ���һ�������ϵ�������
            for m=1:1:ap/dl
                Cd=C-(m-1)*kb*dl;
                fa=fe*sin(Cd)+Dx(i)*sin(Cd)+Dy(i)*cos(Cd);%
                if Cd<Cexd(i)&&Cd>Cstd(i)&&fa>=0
                    fx=(-cos(Cd))*(Ktc*fa+Kte)*dl+(-sin(Cd))*(Krc*fa+Kre)*dl;
                    fy=( sin(Cd))*(Ktc*fa+Kte)*dl+(-cos(Cd))*(Krc*fa+Kre)*dl;
                    fz=(Kac*fa+Kae)*dl;
                else
                    fx=0;fy=0;fz=0;
                end
                Fx=Fx+fx;Fy=Fy+fy;Fz=Fz+fz;
            end
        end
    F(1,i)=Fx;F(2,i)=Fy;F(3,i)=Fz;%%�ھ����д洢������
    Fx=0;Fy=0;Fz=0;%%���������ۼӱ�������
        %Runge-Kutta�����㵶�ߵ���λ�ƺ��ٶ�
    if  i==1 %%�ж��Ƿ��ǳ�ʼ�㣨��һ���㣩
        K1=0;L1=(-cx*0-kx*0+0)/mx;
        K2=0+Dt*L1/2;L2=(-cx*(0+Dt*L1/2)-kx*(0+Dt*K1/2)+(0+F(1,i))/2)/mx;
        K3=0+Dt*L2/2;L3=(-cx*(0+Dt*L2/2)-kx*(0+Dt*K2/2)+(0+F(1,i))/2)/mx;
        K4=0+Dt*L3;L4=(-cx*(0+Dt*L3)-kx*(0+Dt*K3)+F(1,i))/mx;
        Xx(i)=0+Dt*(K1+2*K2+2*K3+K4)/6;
        Vx(i)=0+Dt*(L1+2*L2+2*L3+L4)/6;
        K1=0;L1=(-cy*0-ky*0+0)/my;
        K2=0+Dt*L1/2;L2=(-cy*(0+Dt*L1/2)-ky*(0+Dt*K1/2)+(0+F(2,i))/2)/my;
        K3=0+Dt*L2/2;L3=(-cy*(0+Dt*L2/2)-ky*(0+Dt*K2/2)+(0+F(2,i))/2)/my;
        K4=0+Dt*L3;L4=(-cy*(0+Dt*L3)-ky*(0+Dt*K3)+F(2,i))/my;
        Xy(i)=0+Dt*(K1+2*K2+2*K3+K4)/6;
        Vy(i)=0+Dt*(L1+2*L2+2*L3+L4)/6;
    else %%���ǳ�ʼ��
        K1=Vx(i-1);L1=(-cx*Vx(i-1)-kx*Xx(i-1)+F(1,i-1))/mx;
        K2=Vx(i-1)+Dt*L1/2;L2=(-cx*(Vx(i-1)+Dt*L1/2)-kx*(Xx(i-1)+Dt*K1/2)+(F(1,i-1)+F(1,i))/2)/mx;
        K3=Vx(i-1)+Dt*L2/2;L3=(-cx*(Vx(i-1)+Dt*L2/2)-kx*(Xx(i-1)+Dt*K2/2)+(F(1,i-1)+F(1,i))/2)/mx;
        K4=Vx(i-1)+Dt*L3;L4=(-cx*(Vx(i-1)+Dt*L3)-kx*(Xx(i-1)+Dt*K3)+F(1,i))/mx;
        Xx(i)=Xx(i-1)+Dt*(K1+2*K2+2*K3+K4)/6;
        Vx(i)=Vx(i-1)+Dt*(L1+2*L2+2*L3+L4)/6;
        K1=Vy(i-1);L1=(-cy*Vy(i-1)-ky*Xy(i-1)+F(2,i-1))/my;
        K2=Vy(i-1)+Dt*L1/2;L2=(-cy*(Vy(i-1)+Dt*L1/2)-ky*(Xy(i-1)+Dt*K1/2)+(F(2,i-1)+F(2,i))/2)/my;
        K3=Vy(i-1)+Dt*L2/2;L3=(-cy*(Vy(i-1)+Dt*L2/2)-ky*(Xy(i-1)+Dt*K2/2)+(F(2,i-1)+F(2,i))/2)/my;
        K4=Vy(i-1)+Dt*L3;L4=(-cy*(Vy(i-1)+Dt*L3)-ky*(Xy(i-1)+Dt*K3)+F(2,i))/my;
        Xy(i)=Xy(i-1)+Dt*(K1+2*K2+2*K3+K4)/6;
        Vy(i)=Vy(i-1)+Dt*(L1+2*L2+2*L3+L4)/6;
    end

end
figure(1)
plot(Dt:Dt:Ns*Dt,180*Cstd/pi,'b-','Markersize',7,'Markerface','white','linewidth',3.0);%%�����õ��ǽǶ���
hold on;
plot(Dt:Dt:Ns*Dt,180*Cexd/pi,'r-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('Angle(��)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Dynamic immerse angle and exit angle');
legend('immerse angle','exit angle')

A1=xlsread('EE8.xlsx',1);
A1=downsample(A1,2);
[n,m]=size(F);
A=zeros(m,1);
A=A1(:,2);
A1(:,2)=A1(:,3);
A1(:,3)=A;%%xy���򽻻�
F1=A1(1:m,1:3);
Xx1=zeros(1,m);
Xy1=zeros(1,m);
Vx1=zeros(1,m);
Vy1=zeros(1,m);
% A1(:,2)=A1(:,2)-(sum(A1(1:100,2))+sum(A1(m-99:m,2)))/200-2;%%��������Ҫ������Ư
% A1(:,3)=A1(:,3)-(sum(A1(1:100,3))+sum(A1(m-99:m,3)))/200-2;
for i=1:1:m
        %Runge-Kutta�����㵶�ߵ���λ�ƺ��ٶ�
    if  i==1 %%�ж��Ƿ��ǳ�ʼ�㣨��һ���㣩
        K1=0;L1=(-cx*0-kx*0+0)/mx;
        K2=0+Dt*L1/2;L2=(-cx*(0+Dt*L1/2)-kx*(0+Dt*K1/2)+(0+F1(i,2))/2)/mx;
        K3=0+Dt*L2/2;L3=(-cx*(0+Dt*L2/2)-kx*(0+Dt*K2/2)+(0+F1(i,2))/2)/mx;
        K4=0+Dt*L3;L4=(-cx*(0+Dt*L3)-kx*(0+Dt*K3)+F1(i,2))/mx;
        Xx1(i)=0+Dt*(K1+2*K2+2*K3+K4)/6;
        Vx1(i)=0+Dt*(L1+2*L2+2*L3+L4)/6;
        K1=0;L1=(-cy*0-ky*0+0)/my;
        K2=0+Dt*L1/2;L2=(-cy*(0+Dt*L1/2)-ky*(0+Dt*K1/2)+(0+F1(i,3))/2)/my;
        K3=0+Dt*L2/2;L3=(-cy*(0+Dt*L2/2)-ky*(0+Dt*K2/2)+(0+F1(i,3))/2)/my;
        K4=0+Dt*L3;L4=(-cy*(0+Dt*L3)-ky*(0+Dt*K3)+F1(i,3))/my;
        Xy1(i)=0+Dt*(K1+2*K2+2*K3+K4)/6;
        Vy1(i)=0+Dt*(L1+2*L2+2*L3+L4)/6;
    else %%���ǳ�ʼ��
        K1=Vx1(i-1);L1=(-cx*Vx1(i-1)-kx*Xx1(i-1)+F1(i-1,2))/mx;
        K2=Vx1(i-1)+Dt*L1/2;L2=(-cx*(Vx1(i-1)+Dt*L1/2)-kx*(Xx1(i-1)+Dt*K1/2)+(F1(i-1,2)+F1(i,2))/2)/mx;
        K3=Vx1(i-1)+Dt*L2/2;L3=(-cx*(Vx1(i-1)+Dt*L2/2)-kx*(Xx1(i-1)+Dt*K2/2)+(F1(i-1,2)+F1(i,2))/2)/mx;
        K4=Vx1(i-1)+Dt*L3;L4=(-cx*(Vx1(i-1)+Dt*L3)-kx*(Xx1(i-1)+Dt*K3)+F1(i,2))/mx;
        Xx1(i)=Xx1(i-1)+Dt*(K1+2*K2+2*K3+K4)/6;
        Vx1(i)=Vx1(i-1)+Dt*(L1+2*L2+2*L3+L4)/6;
        K1=Vy1(i-1);L1=(-cy*Vy1(i-1)-ky*Xy1(i-1)+F1(i-1,2))/my;
        K2=Vy1(i-1)+Dt*L1/2;L2=(-cy*(Vy1(i-1)+Dt*L1/2)-ky*(Xy1(i-1)+Dt*K1/2)+(F1(i-1,3)+F1(i,3))/2)/my;
        K3=Vy1(i-1)+Dt*L2/2;L3=(-cy*(Vy1(i-1)+Dt*L2/2)-ky*(Xy1(i-1)+Dt*K2/2)+(F1(i-1,3)+F1(i,3))/2)/my;
        K4=Vy1(i-1)+Dt*L3;L4=(-cy*(Vy1(i-1)+Dt*L3)-ky*(Xy1(i-1)+Dt*K3)+F1(i,3))/my;
        Xy1(i)=Xy1(i-1)+Dt*(K1+2*K2+2*K3+K4)/6;
        Vy1(i)=Vy1(i-1)+Dt*(L1+2*L2+2*L3+L4)/6;
    end

end
% [n,m]=size(F);
figure(2)
subplot(2,2,1);
plot(Dt:Dt:m*Dt,F(1,1:m),'b-','Markersize',7,'Markerface','white','linewidth',2.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('F_x(N)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Simulated Forces in X direction');
subplot(2,2,3);
plot(Dt:Dt:m*Dt,F(2,1:m),'r-','Markersize',7,'Markerface','white','linewidth',2.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('F_y(N)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Simulated Forces in Y direction');

subplot(2,2,2);
plot(Dt:Dt:m*Dt,F1(1:m,2),'b-','Markersize',7,'Markerface','white','linewidth',2.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('F_x(N)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Measured Forces in X direction');
subplot(2,2,4);
plot(Dt:Dt:m*Dt,F1(1:m,3),'r-','Markersize',7,'Markerface','white','linewidth',2.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('F_y(N)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Measured Forces in Y direction');
% figure(2)
% subplot(3,1,1);
% plot(Dt:Dt:Ns*Dt,F(1,:),'b-','Markersize',7,'Markerface','white','linewidth',3.0);
% hold on;
% grid on;
% xlabel('time(s)')
% ylabel('Force(N)')
% set(gca, 'FontSize', 20)
% set(get(gca,'XLabel'),'Fontsize',20)
% set(get(gca,'YLabel'),'Fontsize',20)
% title('Simulated Forces in X direction');
% subplot(3,1,2);
% plot(Dt:Dt:Ns*Dt,F(2,:),'r-','Markersize',7,'Markerface','white','linewidth',3.0);
% hold on;
% grid on;
% xlabel('time(s)')
% ylabel('Force(N)')
% set(gca, 'FontSize', 20)
% set(get(gca,'XLabel'),'Fontsize',20)
% set(get(gca,'YLabel'),'Fontsize',20)
% title('Simulated Forces in Y direction');
% subplot(3,1,3);
% plot(Dt:Dt:Ns*Dt,F(3,:),'g-','Markersize',7,'Markerface','white','linewidth',3.0);
% grid on;
% xlabel('time(s)')
% ylabel('Force(N)')
% set(gca, 'FontSize', 20)
% set(get(gca,'XLabel'),'Fontsize',20)
% set(get(gca,'YLabel'),'Fontsize',20)
% title('Simulated Forces in Z direction');

% n1=40000;n2=45000;
% figure(3)
% subplot(3,1,1);
% plot(n1*Dt:Dt:n2*Dt,F(1,n1:n2),'b-','Markersize',7,'Markerface','white','linewidth',3.0);
% hold on;
% grid on;
% xlabel('time(s)')
% ylabel('Force(N)')
% set(gca, 'FontSize', 20)
% set(get(gca,'XLabel'),'Fontsize',20)
% set(get(gca,'YLabel'),'Fontsize',20)
% title('Simulated Forces in X direction');
% subplot(3,1,2);
% plot(n1*Dt:Dt:n2*Dt,F(2,n1:n2),'r-','Markersize',7,'Markerface','white','linewidth',3.0);
% hold on;
% grid on;
% xlabel('time(s)')
% ylabel('Force(N)')
% set(gca, 'FontSize', 20)
% set(get(gca,'XLabel'),'Fontsize',20)
% set(get(gca,'YLabel'),'Fontsize',20)
% title('Simulated Forces in Y direction');
% subplot(3,1,3);
% plot(n1*Dt:Dt:n2*Dt,F(3,n1:n2),'g-','Markersize',7,'Markerface','white','linewidth',3.0);
% grid on;
% xlabel('time(s)')
% ylabel('Force(N)')
% set(gca, 'FontSize', 20)
% set(get(gca,'XLabel'),'Fontsize',20)
% set(get(gca,'YLabel'),'Fontsize',20)
% title('Simulated Forces in Z direction');

%%���Ƶ�����x�����y�����λ�ƺ��ٶ�ͼ��
figure(4)
subplot(2,2,1)
plot(Dt:Dt:m*Dt,1000*Xx(1:m),'b-','Markersize',7,'Markerface','white','linewidth',2.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('x(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Predicted displacement in x')

figure(4)
subplot(2,2,3)
plot(Dt:Dt:m*Dt,1000*Xy(1:m),'r-','Markersize',7,'Markerface','white','linewidth',2.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('y(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Predicted displacement in y')

figure(4)
subplot(2,2,2)
plot(Dt:Dt:m*Dt,1000*Xx1(1:m),'b-','Markersize',7,'Markerface','white','linewidth',2.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('x(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Measured displacement in x')

figure(4)
subplot(2,2,4)
plot(Dt:Dt:m*Dt,1000*Xy1(1:m),'r-','Markersize',7,'Markerface','white','linewidth',2.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('y(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Measured displacement in y')
% figure(5)
% subplot(2,1,1)
% plot(Dt:Dt:Ns*Dt,1000*Vx(:),'r-');
% hold on;
% grid on;
% xlabel('time(s)')
% ylabel('v_x(mm/s)')
% set(gca, 'FontSize', 20)
% set(get(gca,'XLabel'),'Fontsize',20)
% set(get(gca,'YLabel'),'Fontsize',20)
% title('velocity at x direction')
% 
% figure(5)
% subplot(2,1,2)
% plot(Dt:Dt:Ns*Dt,1000*Vy(:),'r-');
% hold on;
% grid on;
% xlabel('time(s)')
% ylabel('v_y(mm/s)')
% set(gca, 'FontSize', 20)
% set(get(gca,'XLabel'),'Fontsize',20)
% set(get(gca,'YLabel'),'Fontsize',20)
% title('velocity at y direction')

%ʱ��FFT
% N1=40000;N2=70000;
% L=length(F(1,N1:N2));
% Nx=2^nextpow2(L);
% x2=fft(F(1,N1:N2),Nx);
% ma=abs(x2)/Nx*2;
% mag=fftshift(ma);
% fx=linspace(-fs/2,fs/2,Nx);
% figure(6)
% subplot(3,1,1)
% stem(fx,mag,'r','LineWidth',2);
% xlabel('Frequency(Hz)');
% ylabel('Amplitude(N)');
% xlim([0 1000]);
% legend('FFT in x');
% set(gca, 'FontSize', 20)
% set(get(gca,'XLabel'),'Fontsize',20)
% set(get(gca,'YLabel'),'Fontsize',20)
% grid on;
% 
% L=length(F(2,N1:N2));
% Nx=2^nextpow2(L);
% x2=fft(F(2,N1:N2),Nx);
% ma=abs(x2)/Nx*2;
% mag=fftshift(ma);
% fx=linspace(-fs/2,fs/2,Nx);
% figure(6)
% subplot(3,1,2)
% stem(fx,mag,'b','LineWidth',2);
% xlabel('Frequency(Hz)');
% ylabel('Amplitude(N)');
% xlim([0 1000]);
% legend('FFT in y');
% set(gca, 'FontSize', 20)
% set(get(gca,'XLabel'),'Fontsize',20)
% set(get(gca,'YLabel'),'Fontsize',20)
% grid on;
% 
% L=length(F(3,N1:N2));
% Nx=2^nextpow2(L);
% x2=fft(F(3,N1:N2),Nx);
% ma=abs(x2)/Nx*2;
% mag=fftshift(ma);
% fx=linspace(-fs/2,fs/2,Nx);
% figure(6)
% subplot(3,1,3)
% stem(fx,mag,'k','LineWidth',2);
% xlabel('Frequency(Hz)');
% ylabel('Amplitude(N)');
% xlim([0 1000]);
% legend('FFT in z');
% set(gca, 'FontSize', 20)
% set(get(gca,'XLabel'),'Fontsize',20)
% set(get(gca,'YLabel'),'Fontsize',20)
% grid on;

% figure(8)
% plot(Dt:Dt:Ns*Dt,Dx(:),'r-');
% hold on;
% figure(9)
% plot(Dt:Dt:Ns*Dt,Dy(:),'r-');
% hold on;
% grid on;

% n1=20000;n2=24000;
% n3=20070;n4=24070;
% figure(10)
% plot(n1*Dt:Dt:n2*Dt,F(1,n1:n2),'b--','Markersize',7,'Markerface','white','linewidth',2.0);
% hold on;
% plot(n1*Dt:Dt:n2*Dt,F1(n3:n4,2),'b-','Markersize',7,'Markerface','white','linewidth',2.0);
% hold on;
% grid on;
% xlabel('time(s)')
% ylabel('F_x(N)')
% set(gca, 'FontSize', 20)
% set(get(gca,'XLabel'),'Fontsize',20)
% set(get(gca,'YLabel'),'Fontsize',20)
% legend('Simulated F_x','Measured F_x');
% % title('Simulated Forces in X direction');
% 
% figure(11)
% plot(n1*Dt:Dt:n2*Dt,F(2,n1:n2),'r--','Markersize',7,'Markerface','white','linewidth',2.0);
% hold on;
% plot(n1*Dt:Dt:n2*Dt,F1(n3:n4,3),'r-','Markersize',7,'Markerface','white','linewidth',2.0);
% hold on;
% grid on;
% xlabel('time(s)')
% ylabel('F_y(N)')
% set(gca, 'FontSize', 20)
% set(get(gca,'XLabel'),'Fontsize',20)
% set(get(gca,'YLabel'),'Fontsize',20)
% legend('Simulated F_y','Measured F_y');
% % title('Simulated Forces in X direction');
toc