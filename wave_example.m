clear all
close all
clc
%%%%%%%%%%%%% Homogeneouse Wave equation Example(4)%%%%%%%%%%%%%%
l=1;
t=2*l;
w=2*pi*t;
x=0:.01*t:4*pi;
time=5;
Nt=500;
Dt=time/Nt;
c=1;
for m=1:Nt
    t(m)=(m-1)*Dt;
    F=cos(x-4*t(m));
    F1=cos(x+4*t(m));
    F2=F+F1;
    plot(x,F,'y',x,F1,'r',x,F2,'b')
    grid on
    axis([x(1) x(end) -2.1 2.1])
    pause(.05)
end