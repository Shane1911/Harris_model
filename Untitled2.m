%纯滚动速度估算
clear all
close all
clc
dm=56.6;%单位：mm
D=13.5;%单位：mm
alfa=25;%初始接触角/°度
wi=2000*pi/30;
Wm=0.5*wi*(1-(D*cosd(alfa)/dm))
WR=0.5*dm/D*wi*(1-(D*cosd(alfa)/dm)^2)
aj=atan(sind(alfa)/(cosd(alfa)+D/dm))
wy=WR*sin(aj)
wz=WR*cos(aj)
wx=sqrt(WR^2-wy^2-wz^2)
Wm_Wi=Wm/wi