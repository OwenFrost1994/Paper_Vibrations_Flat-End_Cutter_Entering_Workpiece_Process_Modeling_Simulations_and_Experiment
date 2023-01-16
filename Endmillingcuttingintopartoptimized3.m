clear;
clc;
close all;
%%刀具直线进入工件过程的切削力模拟
%%刀具切入工件过程随着径向切深的变化刀具切入点会随之变化，这里不仅仅要计算刀具的旋转角度，还要计算刀具的位置
%%模拟过程是从刀具开始移动，到切入工件一定的距离，刀心距离工件前沿S1（mm），之后再走S2（mm）
%%刀具走刀分为三个部分：空走s1；切入s2；稳定铣削s3
%%变齿间角刀具加工
S1=5.5;%%单位：mm
S2=1;%%单位：mm

%%刀具参数
D=10;%%刀具半径
N=4;%%刀具齿数
B=pi/6;%%刀具螺旋角
% Cp=[96*pi/180,92*pi/180,88*pi/180,84*pi/180];%%变化的齿间角
Cp=[94*pi/180,86*pi/180,94*pi/180,86*pi/180];
dl=0.02;

%%刀具刚度参数%%
kx=23030727.4827;%%x方向刚度
ky=23030727.4827;%%y方向刚度
cx=214.1390;%%x方向阻尼
cy=214.1390;%%y方向阻尼
mx=0.7037;%%x方向质量
my=0.7037;%%y方向质量

%%材料参数%%材料AL7075%%
% Ktc=951.751;%%切向剪切力系数
% Kte=14.0371;%%切向刃口力系数
% Krc=608.561;%%径向剪切力系数
% Kre=16.5002;%%径向刃口力系数
% Kac=288.478;%%轴向剪切力系数
% Kae=1.25118;%%轴向刃口力系数
%%材料参数%%材料GH909%%
Ktc=3250.0332;%%切向剪切力系数
Kte=90.8541;%%切向刃口力系数
Krc=1200.5036;%%径向剪切力系数
Kre=115.1116;%%径向刃口力系数
Kac=251.55;%%轴向剪切力系数
Kae=52.9445;%%轴向刃口力系数

%%加工参数%%
%%铣削方式:顺铣%%
Cm=1;%%铣削方式，顺铣为1，逆铣为0
S=1000;%%主轴转速
f1=200;%%进给速度
f2=200;
fs=10000;%%采样频率
ap=0.5;%%轴向切深（单位mm）
ae=1;%%径向切深（单位mm）

%%基本参数计算%%
R=D/2;%%刀具半径
kb=(2*tan(B))/D;%%kβ计算
w=2*pi*S/60;%%刀具角速度
T=2*pi/w;%%刀具周期
Nc=floor(60*fs/S);%%一个周期内的采样点个数
%%完全切入之后的切入切出角度
if Cm==1%%顺铣
    Cst=pi-acos((R-ae)/R);%%切入角
    Cex=pi;%%切出角
else%%逆铣
    Cst=0;%%切入角
    Cex=acos((R-ae)/R);%%切出角
end
Cs=0;%%开始角度
Dt=T/Nc;%%时间步长
DC=Dt*w;%%角度增量
Ca=Cs;%%初始角度

%%模拟计算的各个阶段s1，s2，s3
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

%%计算模拟时域的总仿真点数
Ns=0;
s=0;
for i=1:1:10000000
    if s<=s1
        f(i)=f1;
    else if s>=s2
            f(i)=f2;
        else
            f(i)=((f2-f1)/(s2-s1))*(s-s1)+f1;
        end
    end
    s=s+Dt*f(i)/60;
    Ns=Ns+1;
    if s>s3
        break;
    end
end

figure(1);
subplot(1,2,1)
plot(Dt:Dt:Ns*Dt,f,'k-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('Feed speed(mm/min)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Dynamic feed following time');
axis([0 Ns*Dt 0 300]);
axis square;
subplot(1,2,2)
plot([0;s1;s2;s3],[f1;f1;f2;f2],'y-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('Displacement(mm)')
ylabel('Feed speed(mm/min)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Dynamic feed following distance');
axis([0 s3 0 300]);
axis square;



%%各种存储单元%%
Fx=0;
Fy=0;
Fz=0;%%刀齿受力累加变量
F=zeros(3,Ns);%%存储三个方向的切削力
Xx=zeros(1,Ns);%%刀具在x方向的振动位移
Xy=zeros(1,Ns);%%刀具在y方向的振动位移
Vx=zeros(1,Ns);%%刀具在x方向的振动速度
Vy=zeros(1,Ns);%%刀具在y方向的振动速度
Dx=zeros(1,Ns);%%刀具在x方向的位移变化量，fv的第一个参数
Dy=zeros(1,Ns);%%刀具在y方向的位移变化量，fv的第二个参数
FA=zeros(1,Ns);
apx=zeros(1,floor(Ns/Nc)*N);
apy=zeros(1,floor(Ns/Nc)*N);
Cstd=zeros(1,Ns);%%记录动态切入切出角度
Cexd=zeros(1,Ns);
%%计算整个模拟域的大循环
%%计算的同时要求计算刀具的实时的切入切出角度
s=0;
for i=1:1:Ns
    s=s+Dt*f(i)/60;
    if s<=s1
        Cstd(i)=0;%%动态切入角
        Cexd(i)=0;%%动态切出角
    else if s<=s2
            if ae<=R
                Cstd(i)=Cst;%%此时的切入角是肯定的
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
            Cstd(i)=Cst;%%最后一段是槽铣
            Cexd(i)=Cex;
        end
    end
    
    Ca=Ca+DC;%%微元角度叠加计算刀具的转动角
    if Ca>=2*pi%%考虑刀具多个旋转周期，累加的刀具角度超过一周就减去一个2π
        Ca=Ca-2*pi;
    else
    end
    
    %%计算附加进给量
    if i<=(Nc/N)%%如果是第一个刀齿在切的时候，切屑上表面没有前一个刀齿切出的波纹表面，此时相当于刀具无振动的情况
        Dx(i)=1000*(Xx(i)-0);%%注意到动力学方程求解过程中单位为国际制单位，求解获得的位移单位为m，而每齿进给量的单位是mm
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
    
        %叠加计算多个刀齿的内循环
        for j=1:1:N
            if j==1
                C=Ca;%%考虑多齿存在的齿间角滞后
            else
                C=Ca-sum(Cp(1:j-1));
            end
            if C<0
                C=C+2*pi;%%如果刀齿角度小于零则转成正的角度
            else
            end
            %%微元叠加计算一个刀刃上的切削力
            for m=1:1:ap/dl
                Cd=C-(m-1)*kb*dl;
                fa=f(i)*sin(Cd)/(N*S)+Dx(i)*sin(Cd)+Dy(i)*cos(Cd);%
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
    F(1,i)=Fx;F(2,i)=Fy;F(3,i)=Fz;%%在矩阵中存储切削力
    Fx=0;Fy=0;Fz=0;%%刀具受力累加变量归零
        %Runge-Kutta法计算刀具的振动位移和速度
    if  i==1 %%判断是否是初始点（第一个点）
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
    else %%不是初始点
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
figure(2)
plot(Dt:Dt:Ns*Dt,180*Cstd/pi,'b-','Markersize',7,'Markerface','white','linewidth',3.0);%%这里用的是角度制
hold on;
plot(Dt:Dt:Ns*Dt,180*Cexd/pi,'r-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('Angle(。)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Dynamic immerse angle and exit angle');
legend('immerse angle','exit angle')

figure(3)
subplot(3,1,1);
plot(Dt:Dt:Ns*Dt,F(1,:),'b-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('Force(N)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Simulated Forces in X direction');
subplot(3,1,2);
plot(Dt:Dt:Ns*Dt,F(2,:),'r-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('Force(N)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Simulated Forces in Y direction');
subplot(3,1,3);
plot(Dt:Dt:Ns*Dt,F(3,:),'g-','Markersize',7,'Markerface','white','linewidth',3.0);
grid on;
xlabel('time(s)')
ylabel('Force(N)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Simulated Forces in Z direction');

% n1=22000;n2=50000;
% figure(4)
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

%%绘制刀具在x方向和y方向的位移和速度图像
figure(5)
subplot(2,1,1)
plot(Dt:Dt:Ns*Dt,1000*Xx(:),'b-');
hold on;
grid on;
xlabel('time(s)')
ylabel('x(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('displacement at x direction')
subplot(2,1,2)
plot(Dt:Dt:Ns*Dt,1000*Xy(:),'b-');
hold on;
grid on;
xlabel('time(s)')
ylabel('y(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('displacement at y direction')

figure(6)
subplot(2,1,1)
plot(Dt:Dt:Ns*Dt,1000*Vx(:),'r-');
hold on;
grid on;
xlabel('time(s)')
ylabel('v_x(mm/s)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('velocity at x direction')
figure(6)
subplot(2,1,2)
plot(Dt:Dt:Ns*Dt,1000*Vy(:),'r-');
hold on;
grid on;
xlabel('time(s)')
ylabel('v_y(mm/s)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('velocity at y direction')

clear s;
s=zeros(1,Ns);
for i=1:1:Ns
    if i==1
        s(i)=0+Dt*f(i)/60;
    else
        s(i)=s(i-1)+Dt*f(i)/60;
    end
end
figure(7)
plot(s(:),1000*Xy(:),'b-');
hold on;
grid on;
xlabel('distance(mm)')
ylabel('y(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('displacement at y direction')

%时域FFT
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