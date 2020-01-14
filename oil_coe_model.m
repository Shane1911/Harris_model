clear all
close all
clc
   mu2=[];
    E=2.26e11;
    TK=[25 75 125];
    for kk=1:3;
      TT=TK(kk);
      jj=1;
    for SR=0:0.02:.2;
    Ry=10;u=30;yita_0=0.033;beta_l=0.0215;K_l=0.0966;
    RF=Ry*1e-3;%单位：m
    Q12=29;
    W=Q12./(E.*RF.^2);
    UF=yita_0.*u./(E.*RF);
    T=TT*sqrt(K_l.*yita_0.*beta_l)/(E.*RF);
    WF=7.956015e12.*T^2+5553.09.*T+1.758455e-6;
    A1=-4.793526e-8*W^(0.0068361.*abs(WF/W-1))*UF^(-0.4047492)*T^(-0.1833848);
    B=8.37449e-15*W^(0.01409715.*abs(WF/W-1))*UF^(-0.5868325)*T^(-0.8173636);
    C=1.180823e-4*W^(0.006321.*abs(WF/W-1))*UF^(-0.20617)*T^(-0.3631740);
    DF=4.793526e-8*W^(0.0068361.*abs(WF/W-1))*UF^(-0.4047492)*T^(-0.1833848);
    mu2(jj)=(A1+B*SR)*exp(-C*SR)+DF;
    jj=jj+1;
    end
    X=0:0.02:.2;
    plot(X,mu2,'-*')
    hold on
    end
    title('润滑油4109，拖动系数随滚滑比的变化(P=20N;V=30 m/s;R取值10mm)')
    grid on
    legend('T=25 C^o','T=75 C^o','T=125 C^o')
%     axis([0 140 0 0.03])