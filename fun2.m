function [ Y ] = fun2( X );
global  alfa A Zn F K_oj K_ij f D m J dm omega_i Q_1j Q_2j hc  alfa_i alfa_o 
    Y=zeros(9,1);
    A_1j=A*sind(alfa)+X(9);
    A_2j=A*cosd(alfa);
    cos_1j=X(2)/((f(1)-0.5)*D+X(3));
    sin_1j=X(1)/((f(1)-0.5)*D+X(3));
    cos_2j=(A_2j-X(2))/((f(2)-0.5)*D+X(4));
    sin_2j=(A_1j-X(1))/((f(2)-0.5)*D+X(4));
    %计算中心油膜厚度
    Q_1j=K_oj*X(3)^1.5;%单位/N
    Q_2j=K_ij*X(4)^1.5;%单位/N
    F_cj=0.5*dm*(1e-3)*m*X(8)^2;
    M_gy=J*X(7)*X(8);
    M_gz=J*X(6)*X(8);
    %求hc
    Q_12=[Q_1j Q_2j];
    gama=[];
    gama(1)=D*cos_1j/dm;
    gama(2)=D*cos_2j/dm;
    u=[];
    u(1)=0.5*dm*1e-3*((1+gama(1))*X(8)-gama(1)*(X(5)*cos_1j+X(7)*sin_1j));
    u(2)=0.5*dm*1e-3*((1-gama(2))*(omega_i-X(8))-gama(2)*(X(5)*cos_2j-X(7)*sin_2j));
    Rx=[];
    Rx(1)=0.5*D*1e-3*(1+D*cosd(alfa)/dm);
    Rx(2)=0.5*D*1e-3*(1-D*cosd(alfa)/dm);
    E1=2.08e11;
    E=E1/(1-0.3^2);
    yita_0=0.033;
    alfa_0=1.28e-8;
    U=yita_0.*u./(2*E.*Rx);
    G=E*alfa_0;
    W_l=Q_12./(E.*Rx.^2);
    k=[6.9520 8.2763];
    hc=2.69*Rx.*U.^(0.67).*G.^(0.53).*W_l.^(-0.067).*(1-0.61*exp(-0.73.*k));%单位：m
    %求摩擦力
    %求椭圆参数
    deta1=[1.0301 1.0230];
    R1=1e-3*[12.0282 9.2494];%m
    a=(6*k.^2.*deta1.*Q_12.*R1/(pi*E)).^(1/3);%m
    b=a./k;
    alfa_o=acos(cos_1j);
    alfa_i=acos(cos_2j);
    W_j=X(5:8);
    ymin1=@(x)-sqrt(1-x.^2);
    ymax1=@(x)sqrt(1-x.^2);
    F_deta1=integral2(@(x,y)fun3(x,y,W_j,a,b,1),-1,1,ymin1,ymax1);
    F_yita1=integral2(@(x,y)fun3(x,y,W_j,a,b,2),-1,1,ymin1,ymax1);
    F_deta2=integral2(@(x,y)fun3(x,y,W_j,a,b,3),-1,1,ymin1,ymax1);
    F_yita2=integral2(@(x,y)fun3(x,y,W_j,a,b,4),-1,1,ymin1,ymax1);
    F_v=0.05*pi*860*(D*1e-3)^2*(dm*1e-3*X(8))^1.95/32/9.8;
    Y(1,1)=-Q_2j*sin_2j+Q_1j*sin_1j+F_yita1*cos_1j-F_yita2*cos_2j;
    Y(2,1)=Q_2j*cos_2j-Q_1j*cos_1j+F_yita1*sin_1j-F_yita2*sin_2j+F_cj;
    Y(3,1)=-F_deta1-F_deta2+F_v;
    Y(4,1)=0.5*D*1e-3*(F_deta1*cos_1j-F_deta2*cos_2j);
    Y(5,1)=0.5*D*1e-3*(F_yita1+F_yita2)-M_gy;
    Y(6,1)=0.5*D*1e-3*(F_deta1*sin_1j-F_deta2*sin_2j)-M_gz;
    Y(7,1)=Q_2j*sin_2j+F_yita2*cos_2j-F/Zn;
    Y(8,1)=(A_1j-X(1))^2+(A_2j-X(2))^2-((f(2)-0.5)*D+X(4))^2;
    Y(9,1)=X(1)^2+X(2)^2-((f(1)-0.5)*D+X(3))^2;
end
