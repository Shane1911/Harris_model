clear all
close all
clc
global alfa A Zn F K_oj K_ij f D omega_i m J dm 
%ϵͳ����:��Ȧ���أ�N N*m
Fa=15000;
F=Fa;
f=[0.523 0.523];%��������
dm=125;%��λ��mm
D=22.23;%��λ��mm
alfa=40;%��ʼ�Ӵ���/���
Zn=16;%���Ӹ���
den=7800*1e-9;%�����ܶ� kg/mm^3
m=pi*den*D^3/6;%�������� kg
J=m*(D*1e-3)^2/10;%����Ĺ��Ծ� kg*m^2
K_ij =1.3617e+06;%��λ N*mm^1.5
K_oj =1.3221e+06;%��λ N*mm^1.5
A=(f(1)+f(2)-1)*D;%�������ľ���/mm
 %��ֵ�趨:λ��/m ��λ�ƻ���
deta_a=0.057;
X_1j=0.2571;X_2j=0.4521;deta_oj=0.0118;deta_ij=0.0096;
Y1=[X_1j X_2j deta_oj deta_ij];%��λ/mm
Wx=2410;Wy=0;Wz=1600;Wm=500;
%������ֵ�趨
X1=[X_1j X_2j deta_oj deta_ij Wx Wy Wz Wm deta_a];
output=[];
% for kk=1:2
omega_i=10200*pi/30;%rad/s
% omega_i=(3000+200*kk)*pi/30;%rad/s
options =optimoptions('fsolve','Algorithm','Levenberg-Marquardt');
[Y,fval,exitflag]=fsolve(@fun2,X1,options);
if exitflag==1
    disp('�������');
%     output(kk,:)=Y;
else
    disp('�����ֵ��Ч');
end
X1=Y
% end
% xlswrite('date1.xlsx',Y,1,'B7')


    