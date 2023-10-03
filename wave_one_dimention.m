 clear all
 close all
 clc
 %%%%%%%%%%%%% Homogeneouse Wave equation in one dimention(3)%%%%%%%%%%%%%%
 l=1;
 t=2*l;
 w=2*pi/t;
 x=0:.01*t:l;
 time=5;
 Nt=500;
 Dt=time/Nt;
 c=1;
 for m=1:Nt
     t(m)=(m-1)*Dt;
     s=0;
     for n=1:2:100
         s=s+(8*sin(n*pi/2)/((n*pi)^2))*sin(n*pi*x/l)*cos(n*pi*t(m)*c/l);
     end
     plot(x,s)
     grid on
     axis([x(1) x(end) -1.1 1.1])
     pause(.1)
 end