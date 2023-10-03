clear all
close all
clc
%%%%%%%%%%%%% D'alembert(7)%%%%%%%%%%%%%%
L=1;
T=2*L;
x=-1.5*T:.01*T:1.5*T;
s=0;
w=2*pi/T;
time=5;
Nt=500;
Dt=time/Nt;
c=1;
s1=0;
for m=1:Nt
    s=0;
    s1=0;
    t(m)=(m-1)*Dt;
for n=1:2:40
s=s+(8*sin(n*pi/2)/((n*pi)^2))*sin(n*w*(x+c*t(m)));
end
for n=1:2:40
s1=s1+(8*sin(n*pi/2)/((n*pi)^2))*sin(n*w*(x-c*t(m)));
end
plot(x,s,x,s1,'r',x,(s+s1)/2,'k')
axis([x(1) x(end) -1.1 1.1])

grid on
pause(.2)
end