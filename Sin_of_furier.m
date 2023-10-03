clear all
close all
clc
T=1;
f=1/T;
w=2*pi*f;
t=0:.01*T:T;
x=sin(w*t);
plot(t,x);
grid on
hold on
x1=sin(2*w*t);
plot(t,x1,'r')
hold on
x2=sin(3*w*t);
plot(t,x2,'k')
hold on
x3=sin(10*w*t);
plot(t,x3,'g')


