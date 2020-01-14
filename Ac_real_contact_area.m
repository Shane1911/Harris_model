%椭圆接触参数
clear all
close all
clc
f=[0.523 0.523];%沟道曲率
dm=125;%单位：mm
D=22.23;%单位：mm
v=[0.3 0.3];
alfa0=40;
Q=3534;%单位：N 假设量
% E12b=1e5*[2.06 2.06];%单位：N/mm^2
% E=2/((1-v(1)^2)/E12b(1)+(1-v(2)^2)/E12b(2))%单位：N/mm^2
gama=D*cosd(alfa0)/dm;
Ry=[];%单位：mm
Ry(1)=(2/D-2/D*(gama/(1+gama)))^(-1);
Ry(2)=(2/D+2/D*(gama/(1-gama)))^(-1)
Rx=[];%单位：mm
Rx(1)=(2/D-1/f(1)/D)^(-1);
Rx(2)=(2/D-1/f(2)/D)^(-1)
sum_cur=[];
sum_cur(1)=1/D*(4-1/f(1)-2*gama/(1-gama));
sum_cur(2)=1/D*(4-1/f(2)+2*gama/(1-gama));
Ep=1.0003+0.5968*(Ry./Rx);
k=1.0339*(Rx./Ry).^0.636;
F=1.5277+0.6023*log(Rx./Ry);
a=0.0236*(2*k.^2.*Ep/pi).^(1/3).*(Q./sum_cur).^(1/3)%单位：mm
b=0.0236*(2*Ep./(pi*k)).^(1/3).*(Q./sum_cur).^(1/3)%单位：mm
xigma=2.*F/pi.*(pi./(2*k.^2.*Ep)).^(1/3);
% Ep12=2.15e5.*sum_cur.^(-0.5).*(xigma).^(-1.5);
% E12=(Ep12(1)^(-2/3)+Ep12(2)^(-2/3))^(-1.5)
E12b=1e5*[2.06 2.06];%单位：N/mm^2
E=2/((1-v(1)^2)/E12b(1)+(1-v(2)^2)/E12b(2))%单位：N/mm^2
A_o=pi.*a.*b%单位：mm^2
%%GW模型
 hA_s=[2 2];%假设值
 m0=[0.0625];m2=[0.0018];m4=[0.000104];%单位：1e-6*mm                         %%待确定参数
 Dsum=m4/(6*pi*m2*sqrt(3));
 alfaA=m0*m4/m2/m2;
 R_Ac=3/8*sqrt(pi/m4);
 S_s=sqrt((1-0.8968/alfaA)*m0);
 d_Ss=(hA_s-4/sqrt(pi*alfaA))./sqrt(1-0.8968/alfaA)
 %计算正态分布值
 F1=[];
 for jj=1:2
 y=@(x)(x-d_Ss(jj)).^(1.5).*normpdf(x);
 F1(jj)=quad(y,d_Ss(jj),10*d_Ss(jj));
 end
 Q_Ao=4/3*E*(R_Ac)^0.5*S_s^1.5*Dsum.*F1%单位:Mpa
 Qa=A_o.*Q_Ao
 %计算点接触的总面积,利用赫兹弹性理论计算
%     Ac=[];a12=1e-9*[5.0752 5.4365];K12=[K_oj K_ij];
%     for kk=1:2 %椭圆接触的面积
%     a=(a12(kk))^(1/3)*(Q12(kk))^(1/3);
%     b=a/K12(kk);
%     As=pi*a*b;%单位：m^2
%     Ac(kk)=As*pi*R_Ac*sqrt(S_s)*Dsum*F1(kk)
%     end
%     Qa=[];
%     Qa(1)=0.25*E*Ac(1)*deta_theta* I_lamta(1)/pi/pi%单位：N
%     Qa(2)=0.25*E*Ac(2)*deta_theta* I_lamta(2)/pi/pi%单位：N