function [F] = fun3( x,y,W_j,a,b,n );%m
global  Q_1j Q_2j alfa_i alfa_o dm D hc omega_i 
T=30;
Po=3*Q_1j/(2*pi*a(1)*b(1))*sqrt(1-x.^2-y.^2);
Pi=3*Q_2j/(2*pi*a(2)*b(2))*sqrt(1-x.^2-y.^2);
yita_o=0.033.*exp(Po.*1.28e-8);
yita_i=0.033.*exp(Pi.*1.28e-8);
%ÇóËÙ¶È²î
wso=W_j(1)*sin(alfa_o)-W_j(3)*cos(alfa_o);
uoy=-0.5*D*1e-3*(W_j(1)*cos(alfa_o)+W_j(1)*sin(alfa_o))+W_j(4)*1e-3*(0.5*dm+0.5*D*cos(alfa_o))-wso.*a(1).*x;
uox=0.5*D*1e-3*W_j(2)+wso.*b(1).*y;
wsi=W_j(1)*sin(alfa_i)-W_j(3)*cos(alfa_i);
uiy=-0.5*D*1e-3*(W_j(1)*cos(alfa_i)+W_j(3)*sin(alfa_i))+(omega_i-W_j(4))*1e-3*(0.5*dm-0.5*D*cos(alfa_i))+wsi*a(2).*x;
uix=0.5*D*1e-3*W_j(2)-wsi*b(2).*y;
switch n
    case 1
        F=a(1)*b(1)*yita_o.*uoy/hc(1);
    case 2
        F=a(1)*b(1)*yita_o.*uox/hc(1);
    case 3
        F=a(2)*b(2)*yita_i.*uiy/hc(2);
   otherwise
        F=a(2)*b(2)*yita_i.*uix/hc(2);
end
