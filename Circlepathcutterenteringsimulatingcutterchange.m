tic
%%Բ�������ӵ����̣�����Բ�������߲�ϳ��ģ���м任������
clear;
clc;
close all;
%%ʹ����ϳ����������Ԥ�⣬��ƫ��%%
D=12;%%���߰뾶
N=2;%%���߳���
B=pi/6;%%����������
Cp=2*pi/N;%%�ݼ��
dl=0.02;%%����΢Ԫ����

%%���߸նȲ���%%
kx=1.74e7;%%x����ն�
ky=1.74e7;%%y����ն�
cx=116.4931;%%x��������
cy=116.4931;%%y��������
mx=0.2198;%%x��������
my=0.2198;%%y��������

%%���ϲ���
%%����GH909%%
% Ktc=3250.0332;%%���������ϵ��
% Kte=90.8541;%%�����п���ϵ��
% Krc=1200.5036;%%���������ϵ��
% Kre=115.1116;%%�����п���ϵ��
% Kac=251.55;%%���������ϵ��
% Kae=52.9445;%%�����п���ϵ��
%%���ϲ���%%����45��%%
% Ktc=3225.0273;%%���������ϵ��
% Kte=82.7731;%%�����п���ϵ��
% Krc=598.3905;%%���������ϵ��
% Kre=60.1713;%%�����п���ϵ��
% Kac=-859.6977;%%���������ϵ��
% Kae=-13.3381;%%�����п���ϵ��
Ktc=2129.9;%%���������ϵ��
Kte=28.6524;%%�����п���ϵ��
Krc=946.7411;%%���������ϵ��
Kre=18.1387;%%�����п���ϵ��
Kac=261.6069;%%���������ϵ��
Kae=-29.9157;%%�����п���ϵ��

%%�ӹ�����%%
%%ϳ����ʽ:˳ϳ%%
Cm=0;%%ϳ����ʽ��˳ϳΪ1����ϳΪ0
S=1000;%%����ת��
f=80;%%�����ٶ�
fs=10000;%%����Ƶ��
ap=2;%%���������λmm��
ae=1;%%���������λmm��
s=3;%%����3mm��ϳ

%%������������%%
R=D/2;%%���߰뾶
R1=5*R;%%������ת�뾶
kb=(2*tan(B))/D;%%k�¼���
fe=f/(N*S);%%feed every tooth
w=2*pi*S/60;%%���߽��ٶ�
w1=f/(R1*60);
T=2*pi/w;%%��������
Nc=floor(60*fs/S);%%һ�������ڵĲ��������
if Cm==1%%˳ϳ
    Cst=pi-acos((R-ae)/R);%%�����
    Cex=pi;%%�г���
else%%��ϳ
    Cst=0;%%�����
    Cex=acos((R-ae)/R);%%�г���
end
Cs=0;%%��ʼ�Ƕ�
Dt=1/fs;%%ʱ�䲽��
DC=Dt*w;%%�Ƕ�����
DC1=Dt*w1;%%��ת�����Ƕ�����
Ca=Cs;%%��ʼ�Ƕ�

%%ģ�����ĸ����׶�s1��s2��s3
C1=asin((R1-ae)/R1);
C2=asin((R1+R-ae)/(R1+R));

%%�����׶η������ܸ���
N1=floor(C1/DC1);
N2=floor(pi/(2*DC1));
N3=floor(60*fs*s/f);
NS=N2+N3;

%%���ִ洢��Ԫ%%
Fx=0;
Fy=0;
Fz=0;%%���������ۼӱ���
F=zeros(3,NS);%%�洢���������������
Xx=zeros(1,NS);%%������x�������λ��
Xy=zeros(1,NS);%%������y�������λ��
Vx=zeros(1,NS);%%������x��������ٶ�
Vy=zeros(1,NS);%%������y��������ٶ�
Dx=zeros(1,NS);%%������x�����λ�Ʊ仯����fv�ĵ�һ������
Dy=zeros(1,NS);%%������y�����λ�Ʊ仯����fv�ĵڶ�������
yy1=zeros(1,NS);
yy2=zeros(1,NS);
l=zeros(1,NS);
FA=zeros(1,NS);
apx=zeros(1,floor(NS/Nc)*N);
apy=zeros(1,floor(NS/Nc)*N);
Cstd=zeros(1,NS);%%��¼��̬�����г��Ƕ�
Cexd=zeros(1,NS);

%%��������ֽ׶����
%%Բ�����߽׶�
for i=1:1:N1
    Ca=Ca+DC;
    if Ca>2*pi
        Ca=Ca-2*pi;
    else
    end
    F(:,i)=[0;0;0];
end

%%Բ������׶�
Cc=C1;
for i=(N1+1):1:N2
    Cc=Cc+DC1;
    %%�������ĵ�λ��
    Ocx=R1*(-sin(Cc));
    Ocy=R1*(-cos(Cc));
    yy1(i)=Ocy+sqrt(R^2-(R+R1+Ocx-ae)^2);
    yy2(i)=Ocy-sqrt(R^2-(R+R1+Ocx-ae)^2);
    l(i)=R+R1+Ocx-ae;
    Yc=[-sin(Cc);-cos(Cc)];
    OcB=[-(R+R1)+ae,yy1(i)]-[Ocx,Ocy];
    OcA=[-(R+R1)+ae,yy2(i)]-[Ocx,Ocy];
    if Cc<=C2
        Cstd(i)=acos((OcA*Yc)/sqrt(OcA(1)^2+OcA(2)^2));
        Cexd(i)=acos((OcB*Yc)/sqrt(OcB(1)^2+OcB(2)^2));
    else
        Cstd(i)=Cst;
        Cexd(i)=acos((OcB*Yc)/sqrt(OcB(1)^2+OcB(2)^2));
    end
    
    Ca=Ca+DC-DC1;%%΢Ԫ�Ƕȵ��Ӽ��㵶�ߵ�ת���ǣ����ڵ��߹켣��ת�Ĵ��ڵ��߽ǶȵĻ�׼y��
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
    
    %%���ϼ����ھ�̬����ϵ����ɣ�ת������������ϵ��
    DD=[-cos(Cc),sin(Cc);-sin(Cc),-cos(Cc)]*[Dx(i);Dy(i)];
        %���Ӽ��������ݵ���ѭ��
        for j=1:1:N
            C=Ca-(j-1)*Cp;
            if C<0
                C=C+2*pi;%%������ݽǶ�С������ת�����ĽǶ�
            else
            end
            %%΢Ԫ���Ӽ���һ�������ϵ�������
            for m=1:1:ap/dl
                Cd=C-(m-1)*kb*dl;
                fa=fe*sin(Cd)+DD(1)*sin(Cd)+DD(2)*cos(Cd);%ʵ�ʽ�����
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
        FF=[-cos(Cc),-sin(Cc),0;sin(Cc),-cos(Cc),0;0,0,1]*[Fx;Fy;Fz];
    F(1,i)=FF(1);F(2,i)=FF(2);F(3,i)=FF(3);%%�ھ����д洢������
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

%%ֱ���ߵ��׶�
for i=(N2+1):1:NS
    %%��ʱ����Ҫ���㶯̬�����г��Ƕ�
    Cstd(i)=Cst;
    Cexd(i)=Cex;
    
    Ca=Ca+DC;%%΢Ԫ�Ƕȵ��Ӽ��㵶�ߵ�ת���ǣ����ڵ��߹켣ֱ�߲���Ҫ��ʲô
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
    
    %%���ϼ����ھ�̬����ϵ����ɣ�ת������������ϵ��
    DD=[0,1;-1,0]*[Dx(i);Dy(i)];
        %���Ӽ��������ݵ���ѭ��
        for j=1:1:N
            C=Ca-(j-1)*Cp;
            if C<0
                C=C+2*pi;%%������ݽǶ�С������ת�����ĽǶ�
            else
            end
            %%΢Ԫ���Ӽ���һ�������ϵ�������
            for m=1:1:ap/dl
                Cd=C-(m-1)*kb*dl;
                fa=fe*sin(Cd)+DD(1)*sin(Cd)+DD(2)*cos(Cd);%ʵ�ʽ�����
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
        FF=[0,-1,0;1,0,0;0,0,1]*[Fx;Fy;Fz];
    F(1,i)=FF(1);F(2,i)=FF(2);F(3,i)=FF(3);%%�ھ����д洢������
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
plot(Dt:Dt:NS*Dt,180*Cstd/pi,'b-','Markersize',7,'Markerface','white','linewidth',3.0);%%�����õ��ǽǶ���
hold on;
plot(Dt:Dt:NS*Dt,180*Cexd/pi,'r-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('Angle(��)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Dynamic immerse angle and exit angle');
legend('immerse angle','exit angle')

A1=xlsread('CE1.xlsx',1);
A1=downsample(A1,2);
[m,n]=size(A1);
A1(:,2)=-A1(:,2);%%xy������
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

n1=307000;
figure(2)
subplot(2,2,1);
plot(Dt:Dt:m*Dt,F(1,n1+1:n1+m),'b-','Markersize',7,'Markerface','white','linewidth',2.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('F_x(N)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Simulated Forces in X direction');
subplot(2,2,3);
plot(Dt:Dt:m*Dt,F(2,n1+1:n1+m),'r-','Markersize',7,'Markerface','white','linewidth',2.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('F_y(N)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Simulated Forces in Y direction');
% subplot(3,1,3);
% plot(Dt:Dt:NS*Dt,F(3,:),'g-','Markersize',7,'Markerface','white','linewidth',3.0);
% grid on;
% xlabel('time(s)')
% ylabel('Force(N)')
% set(gca, 'FontSize', 20)
% set(get(gca,'XLabel'),'Fontsize',20)
% set(get(gca,'YLabel'),'Fontsize',20)
% title('Simulated Forces in Z direction');
subplot(2,2,2);
plot(Dt:Dt:m*Dt,F1(:,2),'b-','Markersize',7,'Markerface','white','linewidth',2.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('F_x(N)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Measured Forces in X direction');
subplot(2,2,4);
plot(Dt:Dt:m*Dt,F1(:,3),'r-','Markersize',7,'Markerface','white','linewidth',2.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('F_y(N)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Measured Forces in Y direction');

figure(3)
subplot(2,2,1)
plot(Dt:Dt:m*Dt,1000*Xx(n1+1:n1+m),'b-','Markersize',7,'Markerface','white','linewidth',2.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('x(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Predicted displacement in x')
subplot(2,2,3)
plot(Dt:Dt:m*Dt,1000*Xy(n1+1:n1+m),'r-','Markersize',7,'Markerface','white','linewidth',2.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('y(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Predicted displacement in y')
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

figure(4)
subplot(2,1,1)
plot(Dt:Dt:NS*Dt,1000*Vx(:),'r-');
hold on;
grid on;
xlabel('time(s)')
ylabel('v_x(mm/s)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('velocity at x direction')
subplot(2,1,2)
plot(Dt:Dt:NS*Dt,1000*Vy(:),'r-');
hold on;
grid on;
xlabel('time(s)')
ylabel('v_y(mm/s)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('velocity at y direction')

% clear s;
% s=zeros(1,Ns);
% for i=1:1:Ns
%     if i==1
%         s(i)=0+Dt*f(i)/60;
%     else
%         s(i)=s(i-1)+Dt*f(i)/60;
%     end
% end
% figure(7)
% plot(s(:),1000*Xy(:),'b-');
% hold on;
% grid on;
% xlabel('distance(mm)')
% ylabel('y(mm)')
% set(gca, 'FontSize', 20)
% set(get(gca,'XLabel'),'Fontsize',20)
% set(get(gca,'YLabel'),'Fontsize',20)
% title('displacement at y direction')
toc