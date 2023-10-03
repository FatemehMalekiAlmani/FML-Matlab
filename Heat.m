clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%% Heat(2) with newman condition %%%%%%%%%%%%%%%%%%%%%
T=2;
F=1/T;
w=2*pi*F;
t=-T/2:.00001*T:T/2;
s=0;
for n=1:1:100
    a(n)=2*(2*cos(n*pi/2)-cos (n*pi/4)-cos(3*n*pi/4))/((n*pi)^2);
    s=s+a(n)*cos(n*w*t);
end
s=s+(1.16);
plot(t,s,'k')
grid on