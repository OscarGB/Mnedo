%%%%%%%%%%%%%%%%%%%%%%%%%%
%RK- explicito Figura (1)
%%%%%%%%%%%%%%%%%%%%%%%%%%

%orden 1 Euler

A=0; b=1; c=0;

figure(1);ezplot(@(x,y)Fest_RK(x,y,A,b,c),[-5,2,-3,3])
hold all
%orden 2: Método de Runge (1895) de orden 2
A=[0,0; 1/2,0];
c=[0; 1/2];
b=[0,1];

figure(1);ezplot(@(x,y)Fest_RK(x,y,A,b,c),[-5,2,-3,3])
hold all
%orden3: Método de Heun de orden 3

A=[0,0,0; 1/3,0,0; 0,2/3,0];
c=[0; 1/3;2/3];
b=[1/4,0,3/4];

figure(1);ezplot(@(x,y)Fest_RK(x,y,A,b,c),[-5,2,-3,3])
hold all
%orden 4: Método de Kutta de orden 4

A=[0,0,0,0; 1/2,0,0,0; 0,1/2,0,0; 0,0,1,0];
c=[0; 1/2;1/2; 1];
b=[1/6,1/3,1/3,1/6];

figure(1);ezplot(@(x,y)Fest_RK(x,y,A,b,c),[-5,2,-3,3])
hold all

% %orden 4: Método de los 3/8 de orden 4
 
 A=[0,0,0,0; 1/3,0,0,0; -1/3,1,0,0; 1,-1,1,0];
 c=[0; 1/3;2/3; 1];
 b=[1/8,3/8,3/8,1/8];

figure(1); ezplot(@(x,y)Fest_RK(x,y,A,b,c),[-5,2,-3,3])
% hold all

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Euler y implicitos Figura (2)
%%%%%%%%%%%%%%%%%%%%%%%%%%

%orden 1 Euler
figure(2);ezplot(@(x,y)Fest_Euler(x,y),[-5,2,-4,4])
hold all

%orden 1. Euler implicito

A=1; b=1; c=1;

figure(2);ezplot(@(x,y)Fest_RK(x,y,A,b,c),[-5,2,-4,4]); hold all

%orden 2. Trapecio
figure(2);ezplot(@(x,y)Fest_trapecio(x,y),[-5,2,-4,4]); hold all

%orden 2. Metodo de Gauss: Regla implícita del punto medio
A=1/2; b=1; c=1/2;

figure(2);ezplot(@(x,y)Fest_RK(x,y,A,b,c),[-5,2,-4,4]); hold all

%Lobatto IIIA. Regla Trapezoidal
A=[0,0; 1/2,1/2];
b=[1/2,1/2];
c=[0;1];
figure(2);ezplot(@(x,y)Fest_RK(x,y,A,b,c),[-5,2,-4,4]); 