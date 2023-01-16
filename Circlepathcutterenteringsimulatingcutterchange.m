tic
%%圆弧进给接刀过程，先走圆弧，后走侧铣，模拟中间换刀过程
clear;
clc;
close all;
%%使用立铣刀的切削力预测，无偏心%%
D=12;%%刀具半径
N=2;%%刀具齿数
B=pi/6;%%刀具螺旋角
Cp=2*pi/N;%%齿间角
dl=0.02;%%刀具微元长度

%%刀具刚度参数%%
kx=1.74e7;%%x方向刚度
ky=1.74e7;%%y方向刚度
cx=116.4931;%%x方向阻尼
cy=116.4931;%%y方向阻尼
mx=0.2198;%%x方向质量
my=0.2198;%%y方向质量

%%材料参数
%%材料GH909%%
% Ktc=3250.0332;%%切向剪切力系数
% Kte=90.8541;%%切向刃口力系数
% Krc=1200.5036;%%径向剪切力系数
% Kre=115.1116;%%径向刃口力系数
% Kac=251.55;%%轴向剪切力系数
% Kae=52.9445;%%轴向刃口力系数
%%材料参数%%材料45钢%%
% Ktc=3225.0273;%%切向剪切力系数
% Kte=82.7731;%%切向刃口力系数
% Krc=598.3905;%%径向剪切力系数
% Kre=60.1713;%%径向刃口力系数
% Kac=-859.6977;%%轴向剪切力系数
% Kae=-13.3381;%%轴向刃口力系数
Ktc=2129.9;%%切向剪切力系数
Kte=28.6524;%%切向刃口力系数
Krc=946.7411;%%径向剪切力系数
Kre=18.1387;%%径向刃口力系数
Kac=261.6069;%%轴向剪切力系数
Kae=-29.9157;%%轴向刃口力系数

%%加工参数%%
%%铣削方式:顺铣%%
Cm=0;%%铣削方式，顺铣为1，逆铣为0
S=1000;%%主轴转速
f=80;%%进给速度
fs=10000;%%采样频率
ap=2;%%轴向切深（单位mm）
ae=1;%%径向切深（单位mm）
s=3;%%再走3mm侧铣

%%基本参数计算%%
R=D/2;%%刀具半径
R1=5*R;%%进刀旋转半径
kb=(2*tan(B))/D;%%kβ计算
fe=f/(N*S);%%feed every tooth
w=2*pi*S/60;%%刀具角速度
w1=f/(R1*60);
T=2*pi/w;%%刀具周期
Nc=floor(60*fs/S);%%一个周期内的采样点个数
if Cm==1%%顺铣
    Cst=pi-acos((R-ae)/R);%%切入角
    Cex=pi;%%切出角
else%%逆铣
    Cst=0;%%切入角
    Cex=acos((R-ae)/R);%%切出角
end
Cs=0;%%开始角度
Dt=1/fs;%%时间步长
DC=Dt*w;%%角度增量
DC1=Dt*w1;%%旋转进给角度增量
Ca=Cs;%%初始角度

%%模拟计算的各个阶段s1，s2，s3
C1=asin((R1-ae)/R1);
C2=asin((R1+R-ae)/(R1+R));

%%各个阶段仿真点的总个数
N1=floor(C1/DC1);
N2=floor(pi/(2*DC1));
N3=floor(60*fs*s/f);
NS=N2+N3;

%%各种存储单元%%
Fx=0;
Fy=0;
Fz=0;%%刀齿受力累加变量
F=zeros(3,NS);%%存储三个方向的切削力
Xx=zeros(1,NS);%%刀具在x方向的振动位移
Xy=zeros(1,NS);%%刀具在y方向的振动位移
Vx=zeros(1,NS);%%刀具在x方向的振动速度
Vy=zeros(1,NS);%%刀具在y方向的振动速度
Dx=zeros(1,NS);%%刀具在x方向的位移变化量，fv的第一个参数
Dy=zeros(1,NS);%%刀具在y方向的位移变化量，fv的第二个参数
yy1=zeros(1,NS);
yy2=zeros(1,NS);
l=zeros(1,NS);
FA=zeros(1,NS);
apx=zeros(1,floor(NS/Nc)*N);
apy=zeros(1,floor(NS/Nc)*N);
Cstd=zeros(1,NS);%%记录动态切入切出角度
Cexd=zeros(1,NS);

%%整个计算分阶段完成
%%圆弧空走阶段
for i=1:1:N1
    Ca=Ca+DC;
    if Ca>2*pi
        Ca=Ca-2*pi;
    else
    end
    F(:,i)=[0;0;0];
end

%%圆弧进入阶段
Cc=C1;
for i=(N1+1):1:N2
    Cc=Cc+DC1;
    %%刀具中心的位置
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
    
    Ca=Ca+DC-DC1;%%微元角度叠加计算刀具的转动角，由于刀具轨迹旋转的存在刀具角度的基准y轴
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
    
    %%以上计算在静态坐标系内完成，转换到刀具坐标系内
    DD=[-cos(Cc),sin(Cc);-sin(Cc),-cos(Cc)]*[Dx(i);Dy(i)];
        %叠加计算多个刀齿的内循环
        for j=1:1:N
            C=Ca-(j-1)*Cp;
            if C<0
                C=C+2*pi;%%如果刀齿角度小于零则转成正的角度
            else
            end
            %%微元叠加计算一个刀刃上的切削力
            for m=1:1:ap/dl
                Cd=C-(m-1)*kb*dl;
                fa=fe*sin(Cd)+DD(1)*sin(Cd)+DD(2)*cos(Cd);%实际进给量
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
    F(1,i)=FF(1);F(2,i)=FF(2);F(3,i)=FF(3);%%在矩阵中存储切削力
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

%%直线走刀阶段
for i=(N2+1):1:NS
    %%此时不需要计算动态切入切出角度
    Cstd(i)=Cst;
    Cexd(i)=Cex;
    
    Ca=Ca+DC;%%微元角度叠加计算刀具的转动角，由于刀具轨迹直线不需要减什么
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
    
    %%以上计算在静态坐标系内完成，转换到刀具坐标系内
    DD=[0,1;-1,0]*[Dx(i);Dy(i)];
        %叠加计算多个刀齿的内循环
        for j=1:1:N
            C=Ca-(j-1)*Cp;
            if C<0
                C=C+2*pi;%%如果刀齿角度小于零则转成正的角度
            else
            end
            %%微元叠加计算一个刀刃上的切削力
            for m=1:1:ap/dl
                Cd=C-(m-1)*kb*dl;
                fa=fe*sin(Cd)+DD(1)*sin(Cd)+DD(2)*cos(Cd);%实际进给量
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
    F(1,i)=FF(1);F(2,i)=FF(2);F(3,i)=FF(3);%%在矩阵中存储切削力
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

figure(1)
plot(Dt:Dt:NS*Dt,180*Cstd/pi,'b-','Markersize',7,'Markerface','white','linewidth',3.0);%%这里用的是角度制
hold on;
plot(Dt:Dt:NS*Dt,180*Cexd/pi,'r-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('Angle(。)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Dynamic immerse angle and exit angle');
legend('immerse angle','exit angle')

A1=xlsread('CE1.xlsx',1);
A1=downsample(A1,2);
[m,n]=size(A1);
A1(:,2)=-A1(:,2);%%xy方向反向
F1=A1(1:m,1:3);
Xx1=zeros(1,m);
Xy1=zeros(1,m);
Vx1=zeros(1,m);
Vy1=zeros(1,m);
% A1(:,2)=A1(:,2)-(sum(A1(1:100,2))+sum(A1(m-99:m,2)))/200-2;%%几乎不需要消除零漂
% A1(:,3)=A1(:,3)-(sum(A1(1:100,3))+sum(A1(m-99:m,3)))/200-2;
for i=1:1:m
        %Runge-Kutta法计算刀具的振动位移和速度
    if  i==1 %%判断是否是初始点（第一个点）
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
    else %%不是初始点
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