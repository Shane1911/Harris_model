clear all
clc
ho=0.57e-4;k=0.038;U=0.1173e-5;
f=@(x)(5.5/(ho^1.5)-12*x/(ho+0.5*x^2)^2-(3+2*k)/U);
alfa1=fzero(f,0.001)
r=alfa1*12.6
