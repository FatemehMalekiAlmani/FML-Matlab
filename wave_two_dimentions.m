clear all
close all
clc
%%%%%%%%%%%%% Wave two dimentions(6)%%%%%%%%%%%%%%
a=4;
x=0:.01*a:a;
y=x;
w=pi/a;
[X Y]=meshgrid(x,y);
time=5;
Nt=500;
Dt=time/Nt;
c=1;
for k=1:Nt
    t(k)=(k-1)*Dt;
    s=0;
    for n=1:10
    for m=1:10
        Gama=c*sqrt((n*pi/a)^2+(m*pi/a)^2);
        s=s+(8*(2*sin(n*pi/2)-sin(n*pi/4)-sin(3*n*pi/4))/((n*pi)^2))*sin(n*w*X).*...
            (8*(2*sin(m*pi/2)-sin(m*pi/4)-sin(3*m*pi/4))/((m*pi)^2)).*sin(m*w*x).*cos(Gama*t(k));
    end
    end
    surf(X,Y,s)
    pause(.02)
    axis([0 a -.3 1.1])
    grid on
    xlabel('X')
    ylabel('Y')
end