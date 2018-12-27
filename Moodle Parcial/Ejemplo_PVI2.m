%Ejemplo PVI2: Movimiento del péndulo simple
%
u0=[pi/2;0];
N=1000;
%método de Euler
A=[0];
b=[1];
c=[0];
[U,t]=RKexplicito(@fPVI2,10,0,N,u0,b,c,A);
plot(t,U(1,:))
hold on 
v=pi*ones(1,N+1)/2;
plot(t,v)
plot(t,-v)
%método de Heun
A=[0, 0,0;1/3,0,0;0,2/3,0];
b=[1/4,0,3/4];
c=[0;1/3;2/3];
[U,t]=RKexplicito(@fPVI2,10,0,N,u0,b,c,A);
figure(2)
plot(t,U(1,:))
hold on
plot(t,v)
%método de Kutta
A=[0, 0,0,0;1/2,0,0,0;0,1/2,0,0;0,0,1,0];
b=[1/6,1/3,1/3,1/6];
c=[0;1/2;1/2;1];
[U,t]=RKexplicito(@fPVI2,10,0,N,u0,b,c,A);
figure(3)
plot(t,U(1,:))
hold on
plot(t,v)