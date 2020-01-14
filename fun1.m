function [ Y ] = fun1( X );
global K_ij K_oj Zn f A alfa R_i K_nc K_c Cp yita_0 D Dc dm omega_i m J F H_dc
%�������ת�ٶȱ仯��
    W_oj=[];
    for j=1:Zn
    W_oj(j)=X(8*j);
    end
    deta_Woj=[];
    deta_Woj(1)=0.5*W_oj(1)*(W_oj(2)-W_oj(16))*Zn/2/pi;
    deta_Woj(16)=0.5*W_oj(16)*(W_oj(1)-W_oj(15))*Zn/2/pi;
    for n=2:15
    deta_Woj(n)=0.5*W_oj(n)*(W_oj(n+1)-W_oj(n-1))*Zn/2/pi;
    end
    deta_Wxj=[];deta_Wyj=[];deta_Wzj=[];
    deta_Wxj(1)=0.5*W_oj(1)*(X(8*2-3)-X(8*16-3))*Zn/2/pi;
    deta_Wxj(16)=0.5*W_oj(16)*(X(8*1-3)-X(8*15-3))*Zn/2/pi;
    deta_Wyj(1)=0.5*W_oj(1)*(X(8*2-2)-X(8*16-2))*Zn/2/pi;
    deta_Wyj(16)=0.5*W_oj(16)*(X(8*1-2)-X(8*15-2))*Zn/2/pi;
    deta_Wzj(1)=0.5*W_oj(1)*(X(8*2-1)-X(8*16-1))*Zn/2/pi;
    deta_Wzj(16)=0.5*W_oj(16)*(X(8*1-1)-X(8*15-1))*Zn/2/pi;
    for mm=2:15
    deta_Wxj(mm)=0.5*W_oj(mm)*(X(8*(mm+1)-3)-X(8*(mm-1)-3))*Zn/2/pi;
    deta_Wyj(mm)=0.5*W_oj(mm)*(X(8*(mm+1)-2)-X(8*(mm-1)-2))*Zn/2/pi;
    deta_Wzj(mm)=0.5*W_oj(mm)*(X(8*(mm+1)-1)-X(8*(mm-1)-1))*Zn/2/pi;
    end
    Wc=sum(W_oj)/Zn;
    %���ּܶ��������ķ���������
    Zcj=[];Zcj(1)=0.08;sum_z=0;
    for nn=2:Zn
        fan=(j-1)*2*pi/Zn;
        sum_z=sum_z+0.5*(W_oj(nn-1)+W_oj(nn))/Wc-1;
        Zcj(nn)=pi*dm/Zn*sum_z-0.1-X(135)*sin(fan)+X(134)*cos(fan);
    end
    Y=[];sumF_cz=0;sumF_cy=0;sumM_x=0;sum_Fx=0;sum_Fy=0;sum_Fz=0;sum_My=0;sum_Mz=0;
    for j=1
    %����λ�ý�
    fan=(j-1)*2*pi/Zn;
    %��������
    A_1j=A*sind(alfa)+X(129)+R_i*X(132)*sin(fan)+R_i*X(133)*cos(fan);
    A_2j=A*cosd(alfa)+X(130)*cos(fan)+X(131)*sin(fan);
    %����Q1(2)
    Q1=K_oj*X(8*j-5)^1.5;
    Q2=K_ij*X(8*j-4)^1.5;
    Q12=[Q1 Q2];
    %��������Ӵ������Ǻ���
    cos_1j=X(8*j-7)/((f(1)-0.5)*D+X(8*j-5));
    sin_1j=X(8*j-6)/((f(1)-0.5)*D+X(8*j-5));
    cos_2j=(A_2j-X(8*j-7))/((f(2)-0.5)*D+X(8*j-4));
    sin_2j=(A_1j-X(8*j-6))/((f(2)-0.5)*D+X(8*j-4));
    %���������֮��ķ�����
    Z_cj=Zcj(j);
    if Z_cj<=Cp
    Q_cj=K_c*Z_cj;%��λ��N
    else
    Q_cj=K_c*Cp+K_nc*(Z_cj-Cp)^1.5;%��λ��N
    end
    %����������Ĥ���
    W2=omega_i;W1=0;%Ȧ��ת�٣�rad/s
    Wxj=X(8*j-3);Wyj=X(8*j-2);Wzj=X(8*j-1);
    Woj=W_oj(j);%����ת�ٶȣ�rad/s
    E=2.26e11;%��λ��Pa
    %R-���������Ч���ʰ뾶
    gama_12=[];
    gama_12(1)=D*cos_1j/dm;
    gama_12(2)=D*cos_2j/dm;
    %���������ϵĵ�Ч�뾶
    Rx=[];%��λ��mm
    Rx(1)=252;
    Rx(2)=252;
    Ry=[];%��λ��mm
    Ry(1)=12.43;
    Ry(2)=9.79;
    %�������k
    k=1.0339*(Rx./Ry).^0.636;
    %����������Ĥ���
    %��������ٶ�u
    u=[];%��������ЧӦ:m/s
    u(1)=0.25*D*1e-3*((1/gama_12(1)+1)*abs(Woj)*cos_1j+Wxj*cos_1j+Wyj*sin_1j);
    u(2)=0.25*D*1e-3*((1/gama_12(2)-1)*abs(W2-Woj)*cos_2j+Wxj*cos_2j+Wyj*sin_2j);
    %���������ٲ���
    %���Ͳ���-4109
    alfa_l=1.28e-8;%Pa^-1 ճѹϵ��
    beta_l=0.0215;%C^-1 ճ��ϵ��
    K_l=0.0953;%N/S*C ����ϵ��
    U=yita_0.*u./(E.*Ry*1e-3);
    G=alfa_l*E;
    W=Q12/E./((Ry*1e-3).^2);
    Hc=2.69*G^(0.53).*(U.^(0.67)).*(W.^(-0.067)).*(1-0.61.*exp(-0.73.*k));
    hc=Hc.*Ry*1e-3;%��Ĥ���ĺ��δ����ֵ:mm
    %������ϵ��
    QT=yita_0*u.*beta_l./0.13;
    CTp=2.564./(2.564+QT.^0.548);
    hc=0.8*hc.*CTp;
    %%��Ĥ�϶���
    %ʵ�ʽӴ����
    s=2;%��λ��mm
    %��������Ĥ����
    lamta=hc./s;
    mu1=0.1;%��Ħ��ϵ��                                    %%��ȷ����
    %�����Ħ����ѹ��Qa
    I_lamta=[];
    for ii=1:2
    if lamta(ii)>3||lamta(ii)<=0.4
    I_lamta(ii)=0;
    elseif lamta(ii)<=2&&lamta(ii)>=0.4
    I_lamta(ii)=2.31*exp(-1.84*lamta(ii))+0.1175*(lamta(ii)-0.4)^0.6*(2-lamta(ii))^2;
    else
    I_lamta(ii)=17*exp(-2.84*lamta(ii))+1.44e-4*(lamta(ii)-2)^1.1*(4-lamta(ii))^7.8;
    end
    end
    deta_theta=tand(2);%�ֲڷ�б��/��
    %%����ʵ�ʽӴ����,����GWģ��
    %����ֲڶȵķ��β���
    m0=[0.0625];m2=[0.0018];m4=[0.000104];%��λ��1e-6*mm                         %%��ȷ������
    Dsum=m4/(6*pi*m2*sqrt(3));
    alfaA=m0*m4/m2/m2;
    hA_s=hc/sqrt(m0);
    R_Ac=3/8*sqrt(pi/m4);
    S_s=(1-0.8968/alfaA)*m0;
    d_Ss=(hA_s-4/sqrt(pi*alfaA))./sqrt(1-0.8968/alfaA);
    %������̬�ֲ�ֵ
    F1=[];
    for jj=1:2
    y=@(x)(x-d_Ss(jj)).*normpdf(x);
    F1(jj)=quad(y,d_Ss(jj),10*d_Ss(jj));
    end
    %�����Ӵ��������,���ú��ȵ������ۼ���
    Ac=[];a12=1e-9*[5.0752 5.4365];K12=[K_oj K_ij];
    for kk=1:2 %��Բ�Ӵ������
    a=(a12(kk))^(1/3)*(Q12(kk))^(1/3);
    b=a/K12(kk);
    As=pi*a*b*1e-6;%��λ��m^2
    Ac(kk)=As*pi*R_Ac*sqrt(S_s)*Dsum*F1(kk);
    end
    Qa=[];
    Qa(1)=0.25*E*Ac(1)*deta_theta* I_lamta(1)/pi/pi;%��λ��N
    Qa(2)=0.25*E*Ac(2)*deta_theta* I_lamta(2)/pi/pi;%��λ��N
    Ta=mu1*Qa;%��λ��N
    %���㲻��ȫ��Ħ����
    Qf=Q12-Qa;
    %���������϶�ϵ��
    mu2=[];SR=[];
    SR(1)=abs(0.5*D*1e-3*((1/gama_12(1)+1)*(Woj)*cos_1j-Wxj*cos_1j-Wyj*sin_1j))/u(1);
    SR(2)=abs(0.5*D*1e-3*((1/gama_12(2)-1)*(W2-Woj)*cos_2j-Wxj*cos_2j-Wyj*sin_2j))/u(2);
    RF=Ry*1e-3;%��λ��m
    WF=Q12./(E.*RF.^2);
    UF=yita_0.*u./(E.*RF);
    T=25*sqrt(K_l*beta_l*yita_0)./(E.*RF);
    for nn=1:2
    A1=-4.793562e-8*WF(nn)^(0.0068361*abs((7.956e12*T(nn)^2+5553.09*T(nn)+1.758455e-6)/WF(nn)-1))*UF(nn)^-0.4047492*T(nn)^-0.1833848;
    B=8.37449e-15*WF(nn)^(0.01409715*abs((7.956e12*T(nn)^2+5553.09*T(nn)+1.758455e-6)/WF(nn)-1))*UF(nn)^-0.5868325*T(nn)^-0.8173636;
    C=1.180823e-4*WF(nn)^(0.0061321*abs((7.9566e12*T(nn)^2+5553.09*T(nn)+1.758455e-6)/WF(nn)-1))*UF(nn)^-0.20617*T(nn)^-0.3431740;
    DF=4.793562e-8*WF(nn)^(0.0068361*abs((7.9566e12*T(nn)^2+5553.09*T(nn)+1.758455e-6)/WF(nn)-1))*UF(nn)^-0.4047492*T(nn)^-0.1833848;
    mu2(nn)=(A1+B*SR(nn))*exp(-C*SR(nn))+DF;
    end
    Tf2=Qf.*mu2;%��λ��N
    T_delta=Ta+Tf2%��λ��N
    mu3=[];SR1=[];
    SR1(1)=2;
    SR1(2)=2;
    u1=0.25*D*1e-3*Wzj;
    RF=Rx*1e-3;%��λ��m
    WF=Q12./(E.*RF.^2);
    UF=yita_0.*u1./(E.*RF);
    T=25*sqrt(K_l*beta_l*yita_0)./(E.*RF);
    for nn=1:2
    A1=-4.793562e-8*WF(nn)^(0.0068361*abs((7.956e12*T(nn)^2+5553.09*T(nn)+1.758455e-6)/WF(nn)-1))*UF(nn)^-0.4047492*T(nn)^-0.1833848;
    B=8.37449e-15*WF(nn)^(0.01409715*abs((7.956e12*T(nn)^2+5553.09*T(nn)+1.758455e-6)/WF(nn)-1))*UF(nn)^-0.5868325*T(nn)^-0.8173636;
    C=1.180823e-4*WF(nn)^(0.0061321*abs((7.9566e12*T(nn)^2+5553.09*T(nn)+1.758455e-6)/WF(nn)-1))*UF(nn)^-0.20617*T(nn)^-0.3431740;
    DF=4.793562e-8*WF(nn)^(0.0068361*abs((7.9566e12*T(nn)^2+5553.09*T(nn)+1.758455e-6)/WF(nn)-1))*UF(nn)^-0.4047492*T(nn)^-0.1833848;
    mu3(nn)=(A1+B*SR1(nn))*exp(-C*SR1(nn))+DF;
    end
    Tf3=Qf.*mu3;%��λ��N
    T_yita=Ta+Tf3%��λ��N  
    
    %%���嶯ѹĦ����
    %ball-raceway
%     k_s12=Ry./Rx;
%     K_cur=[];
%     K_cur(1)=(1+gama_12(1))*(2*f(1)-1)/2/f(1);
%     K_cur(2)=(1+gama_12(2))*(2*f(2)-1)/2/f(2);
%     r_s12=((Ry.*1e-3)./K_cur).^(0.5).*sqrt((3+2.*K_cur).*(4/3*hc*1e-3)-2.*(2+K_cur).*hc*1e-3);
%     u_deta=u;
%     u_yita=[];
%     u_yita(1)=0.5*D*1e-3*Wzj;
%     u_yita(2)=0.5*D*1e-3*Wzj;
%     theta_FR=atan((3+2.*k_s12)./(k_s12.^0.5.*(3+2./k_s12)).*(u_yita./u_deta));
%     t_FR=r_s12.*sqrt(cos(theta_FR).^2+(sin(theta_FR).^2)./k_s12).*sqrt(2e-6*hc.*Rx);
%     if t_FR>5
%     F_R12=28.59.*log(t_FR)-10.10;
%     else
%     F_R12=36.57.*log(t_FR)-22.85;
%     end
%     C_FR=yita_0.*u_deta.*((Rx.*Ry*1e-6).^0.5).*((3+2.*k_s12).^(-2)+(u_yita./u_deta).^2.*(3+2./k_s12).^(-2)./k_s12).^0.5;
    FR_yita=[0 0];%0.5.*C_FR.*F_R12.*sin(theta_FR).*sqrt(k_s12);
    FR_deta=[0 0];%0.5.*C_FR.*F_R12.*cos(theta_FR);
%     %�������嶯ѹ������ˮƽ����
    FH_deta=[0 0];%2.*C_FR.*F_R12.*Rx.*cos(theta_FR)/D;
    FH_yita=[0 0];%2.*C_FR.*F_R12.*Ry.*sqrt(k_s12).*sin(theta_FR)/D;
%     %ball-cage
%     u_spnj=0.5*D*1e-3*Wyj;%��λ��m/s
%     u_spej=0.5*D*1e-3*Wxj;%��λ��m/s
%     u_pej=0.25*D*1e-3*Wxj;%��λ��m/s
%     u_pnj=0.25*D*1e-3*Wyj;%��λ��m/s
%     R_pe=0.5*1e-3*D;%��λ��m
%     R_pn=1e-3/(2/D-1/(0.5*D+Cp));%��λ��m
%     kp=R_pe/R_pn;
%     h_opj=abs(Cp-Z_cj)*1e-3;%��λ��m
%     theta_pj=atan(u_pnj*(3+2*kp)/(kp^0.5*u_pej*(3+2/kp)));
%     Rou_1pj=0.25*1e-3*(H_dc)/sqrt(2*R_pe*h_opj*((cos(theta_pj))^2+(sin(theta_pj))^2/kp));
%     P_Rj=34.74*log(Rou_1pj)-27.6;
%     P_Sj=0.26*P_Rj+10.9;
%     C_opj=yita_0*u_pej*sqrt(R_pe*R_pn*((3+2*kp)^-2+u_pnj^2*(3+2/kp)^(-2)/kp/(u_pej^2)));
    P_Rej=0;%0.5*C_opj*P_Rj*cos(theta_pj);
    P_Rnj=0;%0.5*C_opj*P_Rj*sqrt(R_pe/R_pn)*sin(theta_pj);
    P_Sej=0;%P_Sj*yita_0*u_spej*sqrt(R_pe*R_pn);
    P_Snj=0;%P_Sj*yita_0*u_spnj*sqrt(R_pe*R_pn);
    %%���������
    F_nj=m*Wc^2*0.5*dm*1e-3;
    F_taoj=-m*0.5*dm*1e-3*deta_Woj(j);
    G_yj=J*Wc*Wzj;
    G_zj=J*Wc*Wyj;
    C_v=0.45;Rou_ef=0.015*800;                                             %%��ȷ������
    F_Dj=Rou_ef*pi*C_v*(D*1e-3)^2*(dm*1e-3*Wc)^1.95/(32*9.8);%��λ:N
    %����ƽ�ⷽ����
    Y(8*j-7,1)=Q2*sin_2j-Q1*sin_1j-FR_yita(2)*cos_2j+FR_yita(1)*cos_1j+FH_yita(2)*cos_2j-FH_yita(1)*cos_1j+P_Sej+P_Rej+T_yita(2)*cos_2j-T_yita(1)*cos_1j;
    Y(8*j-6,1)=Q2*cos_2j-Q1*cos_1j+FR_yita(2)*sin_2j-FR_yita(1)*sin_1j-FH_yita(2)*sin_2j+FH_yita(1)*sin_1j+F_nj+P_Snj+P_Rnj;
    Y(8*j-5,1)=T_delta(2)-T_delta(1)+FR_deta(1)-FR_deta(2)-FH_deta(1)+FH_deta(2)+Q_cj-F_Dj-F_taoj;
    Y(8*j-4,1)=(T_delta(1)-FR_deta(1))*0.5*D*1e-3*cos_1j-(T_delta(2)-FR_deta(2))*0.5*1e-3*D*cos_2j-(P_Snj+P_Rnj)*0.5*1e-3*D+J*deta_Wxj(j);
    Y(8*j-3,1)=-(T_delta(1)-FR_deta(1))*0.5*D*1e-3*sin_1j-(T_delta(2)-FR_deta(2))*0.5*1e-3*D*sin_2j+G_yj-(P_Sej+P_Rej)*0.5*1e-3*D-J*deta_Wyj(j);
    Y(8*j-2,1)=(T_delta(1)-FR_deta(1))*0.5*D*1e-3+(T_delta(2)-FR_deta(2))*0.5*1e-3*D-G_zj-J*deta_Wzj(j);
    Y(8*j-1,1)=(A_1j-X(8*j-6))^2+(A_2j-X(8*j-7))^2-((f(2)-0.5)*D+X(8*j-4))^2;
    Y(8*j,1)=X(8*j-7)^2+X(8*j-6)^2-((f(1)-0.5)*D+X(8*j-5))^2;
%     sumF_cz=0;sumF_cy=0;sumM_x=0;
%     sumF_cz=sumF_cz-Q_cj*sin(fan)-(P_Snj+P_Rnj)*cos(fan);
%     sumF_cy=sumF_cy-Q_cj*cos(fan)+(P_Snj+P_Rnj)*sin(fan);
%     sumM_x=sumM_x+Q_cj*dm*0.5-0.5*(P_Snj+P_Rnj)*Dc;
%     %�����Ȧ���������
%     sum_Fx=sum_Fx+Q2*sin_2j+T_yita(2)*cos_2j-FR_yita(2)*cos_2j;
%     sum_Fz=sum_Fz+(Q2*cos_2j-T_yita(2)*sin_2j+FR_yita(2)*sin_2j)*cos(fan);
%     sum_Fy=sum_Fy+(Q2*cos_2j-T_yita(2)*sin_2j+FR_yita(2)*sin_2j)*sin(fan);
%     sum_Mz=sum_Mz+((Q2*sin_2j+T_yita(2)*cos_2j-FR_yita(2)*cos_2j)*R_i-f(2)*D*T_yita(2)*cos_2j+f(2)*D*FR_yita(2)*cos_2j)*sin(fan);
%     sum_My=sum_My+((Q2*sin_2j+T_yita(2)*cos_2j-FR_yita(2)*cos_2j)*R_i-f(2)*D*T_yita(2)*cos_2j+f(2)*D*FR_yita(2)*cos_2j)*cos(fan); 
% end
%     %���ּ���������
%    %%���ּܴ�ȷ����
%    R_c1=134.52*1e-3;C_c1=0.14*1e-3;L_c=28.2*1e-3;%��λ��m;                 %%��ȷ��
%    deta_zc=X(135);deta_yc=X(134);
%    fan_c=X(136);                     
%    delta_c=sqrt(deta_zc^2+deta_yc^2)/C_c1;
%    u_c=R_c1*(Wc+omega_i);
%    V_c=R_c1*(omega_i-Wc);
%    FC_z=yita_0*u_c*L_c^3*delta_c^2/(C_c1^2*(1-delta_c^2)^2);
%    FC_y=pi*yita_0*u_c*L_c^3*delta_c^2/(4*C_c1^2*(1-delta_c^2)^(3/2));
%    MC_x=2*pi*yita_0*V_c*L_c*R_c1/(C_c1*sqrt(1-delta_c^2));
%    MC=[1 0 0;0 cos(fan_c) -sin(fan_c);0 sin(fan_c) cos(fan_c)];
%    FC=MC*[MC_x FC_z FC_y]';
%    Y(129,1)=sumF_cz+FC(2)-0.0469*9.8;
%    Y(130,1)=sumF_cy+FC(3);
%    Y(131,1)=sumM_x+FC(1);
%    Y(132,1)=F(1)-sum_Fx;
%    Y(133,1)=F(3)-sum_Fz;
%    Y(134,1)=F(2)-sum_Fy;
%    Y(135,1)=F(4)-sum_My;
%    Y(136,1)=F(5)-sum_Mz;
end

