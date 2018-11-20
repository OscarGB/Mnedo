%Runge Kutta Implicitos
%
%ejemplo de euler
%aplicado al ejercicio 1

%Para graficos en ubuntu (si no me peta) (SOLO EN OCTAVE)
% graphics_toolkit("gnuplot");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Definicion del método RK Euler
A=[0];
b=[1];
c=[0];
%valor inicial
u0=[pi/2; 0]; % angulo, velocidad
t0 = 0;
tf = 10;
%Método Numérico
[U,t]=RKexplicito(@f2,tf,t0,100000,u0,b,c,A);
%Pintamos resultados del MN
figure(1);
hold on;
plot(t,U(1,:),'r');
title("Método de Euler");
plot(t, U(2,:),'b');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Definicion del método RK Heun
A=[0,0,0;1/3,0,0;0,2/3,0];
b=[1/4,0,3/4];
c=[0;1/3;2/3];
%valor inicial
u0=[pi/2; 0]; % angulo, velocidad
t0 = 0;
tf = 10;
%Método Numérico
[U,t]=RKexplicito(@f2,tf,t0,100000,u0,b,c,A);
%Pintamos resultados del MN
figure(2);
hold on;
plot(t,U(1,:),'r');
title("Método de Heun");
plot(t, U(2,:),'b');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Definicion del método RK Kutta
A=[0,0,0,0; 1/2,0,0,0;0,1/2,0,0;0,0,1,0];
b=[1/6,1/3,1/3,1/6];
c=[0;1/2;1/2;1];
%valor inicial
u0=[pi/2; 0]; % angulo, velocidad
t0 = 0;
tf = 10;
%Método Numérico
[U,t]=RKexplicito(@f2,tf,t0,100000,u0,b,c,A);
%Pintamos resultados del MN
figure(3);
hold on;
plot(t,U(1,:),'r');
title("Método de Kutta");
plot(t, U(2,:),'b');

hold off;

