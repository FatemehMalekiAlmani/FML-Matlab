clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%% Heat equation(1) with Drichlet condition %%%%%%%%%%%%%%%%%%%%%
L=1;
x=0:0.01*L:L;
time=50;
Nt=100;
Dt=time/Nt;
c=.05;
for m=1:Nt
    t(m)=(m-1)*Dt;
    s=0;
    for n=1:2:100
        s=s+4*L*sin(n*pi/2)/((n*pi^2))*sin(n*pi*x/L)*exp(-(n*pi*c/L)^2*t(m));
    end
plot(x,s)
grid on
axis ([0 L 0 1])
pause(.2)
end