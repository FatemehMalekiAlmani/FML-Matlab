clear all
close all
clc
%%%%%%%%%%%%% standing Wave(5)%%%%%%%%%%%%%%
L=1;
x=0:0.01*L:L;
time=50;
Nt=500;
Dt=time/Nt;
c=1;
for m=1:Nt
    t(m)=(m-1)*Dt;
    s=0;
    s1=0;
    s2=0;
    for n=1:2:100
        s=s+(8*sin(n*pi/2)/(((n*pi)^2))*sin(n*pi*x/L)*cos(n*pi*t(m)*c/L));
        s1=s1+(4*sin(n*pi/2)/(((n*pi)^2))*sin(n*pi*(x+c*t(m))/L));% rolling wave in opposit direction
        s2=s2+(4*sin(n*pi/2)/(((n*pi)^2))*sin(n*pi*(x-c*t(m))/L));% rolling wave in straight direction
    end
    plot(x,s1,'b',x,s2,'r',x,s,'k')
    grid on
    axis([x(1) x(end) -1.1 1.1])
pause(.2)
end