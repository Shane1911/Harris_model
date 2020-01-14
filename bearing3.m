clear all
close all
clc
global alfa A Zn F K_oj K_ij f D omega_i m J dm 
%系统输入:内圈受载：N N*m
Fa=15000;
F=Fa;
f=[0.523 0.523];%沟道曲率
dm=125;%单位：mm
D=22.23;%单位：mm
alfa=40;%初始接触角/°度
Zn=16;%滚子个数
den=7800*1e-9;%钢球密度 kg/mm^3
m=pi*den*D^3/6;%钢球质量 kg
J=m*(D*1e-3)^2/10;%钢球的惯性矩 kg*m^2
K_ij =1.3617e+06;%单位 N*mm^1.5
K_oj =1.3221e+06;%单位 N*mm^1.5
A=(f(1)+f(2)-1)*D;%曲率中心距离/mm
 %初值设定:位移/m 角位移弧度
deta_a=0.057;
X_1j=0.2571;X_2j=0.4521;deta_oj=0.0118;deta_ij=0.0096;
Y1=[X_1j X_2j deta_oj deta_ij];%单位/mm
Wx=2410;Wy=0;Wz=1600;Wm=500;
%迭代初值设定
X1=[X_1j X_2j deta_oj deta_ij Wx Wy Wz Wm deta_a];
output=[];
% for kk=1:2
omega_i=10200*pi/30;%rad/s
% omega_i=(3000+200*kk)*pi/30;%rad/s
options =optimoptions('fsolve','Algorithm','Levenberg-Marquardt');
[Y,fval,exitflag]=fsolve(@fun2,X1,options);
if exitflag==1
    disp('结果收敛');
%     output(kk,:)=Y;
else
    disp('所设初值无效');
end
X1=Y
% end
% xlswrite('date1.xlsx',Y,1,'B7')


    