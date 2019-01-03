graphics_toolkit("gnuplot");

%MATRICES PARA RK
ARunge = [0,0;1/2,0];
BRunge = [0;1/2];
CRunge = [0,1];

AHeun2 = [0,0;1,0];
BHeun2 = [0;1];
CHeun2 = [1/2,1/2];

AHeun3 = [0,0,0;1/3,0,0;0,2/3,0];
BHeun3 = [0;1/3;2/3];
CHeun3 = [1/4,0,3/4];

AKutta = [0,0,0,0;1/2,0,0,0;0,1/2,0,0;0,0,1,0];
BKutta = [0;1/2;1/2;1];
CKutta = [1/6,1/3,1/3,1/6];

A38 = [0,0,0,0;1/3,0,0,0;-1/3,1,0,0;1,-1,1,0];
B38 = [0;1/3;2/3;1];
C38 = [1/8,3/8,3/8,1/8];

u0 = [2,3]; % Valor inicial

%EULER
N=[20, 200, 2000, 20000, 200000];
e1=zeros(1,5);
e2=e1;
for n=1:5
  [U,t]=euler(@ej1, N(n), 0, 10, u0);
  uexac = feval(@exacta, t);   
  e1(n)=max(abs(U(:,N(n)+1)-uexac(:,N(n)+1)));
  e2(n)=max(max(abs(U-uexac)));
end
eeuler = e2;

figure(1);
clf;
hold on;
figure(1);
title("Método de Euler");
loglog(N,e1)
loglog(N,e2)
legend('final', 'máximo');
hold off;

%Runge
N=[20, 200, 2000, 20000, 200000];
e1=zeros(1,5);
e2=e1;
for n=1:5
  [U,t]=RKexplicito(@ej1, N(n), 0, 10, u0, ARunge, BRunge, CRunge);
  uexac = feval(@exacta, t);   
  e1(n)=max(abs(U(:,N(n)+1)-uexac(:,N(n)+1)));
  e2(n)=max(max(abs(U-uexac)));
end
erunge = e2;

figure(2);
clf;
hold on;
figure(2);
title("Método de Runge");
loglog(N,e1)
loglog(N,e2)
legend('final', 'máximo');
hold off;

%Heun2
N=[20, 200, 2000, 20000, 200000];
e1=zeros(1,5);
e2=e1;
for n=1:5
  [U,t]=RKexplicito(@ej1, N(n), 0, 10, u0, AHeun2, BHeun2, CHeun2);
  uexac = feval(@exacta, t);   
  e1(n)=max(abs(U(:,N(n)+1)-uexac(:,N(n)+1)));
  e2(n)=max(max(abs(U-uexac)));
end
eheun2 = e2;

figure(3);
clf;
hold on;
figure(3);
title("Método de Heun2");
loglog(N,e1)
loglog(N,e2)
legend('final', 'máximo');
hold off;

%Heun3
N=[20, 200, 2000, 20000, 200000];
e1=zeros(1,5);
e2=e1;
for n=1:5
  [U,t]=RKexplicito(@ej1, N(n), 0, 10, u0, AHeun3, BHeun3, CHeun3);
  uexac = feval(@exacta, t);   
  e1(n)=max(abs(U(:,N(n)+1)-uexac(:,N(n)+1)));
  e2(n)=max(max(abs(U-uexac)));
end
eheun3 = e2;

figure(4);
clf;
hold on;
figure(4);
title("Método de Heun3");
loglog(N,e1)
loglog(N,e2)
legend('final', 'máximo');
hold off;


%Kutta
N=[20, 200, 2000, 20000, 200000];
e1=zeros(1,5);
e2=e1;
for n=1:5
  [U,t]=RKexplicito(@ej1, N(n), 0, 10, u0, AKutta, BKutta, CKutta);
  uexac = feval(@exacta, t);   
  e1(n)=max(abs(U(:,N(n)+1)-uexac(:,N(n)+1)));
  e2(n)=max(max(abs(U-uexac)));
end
ekutta = e2;

figure(5);
clf;
hold on;
figure(5);
title("Método de Kutta");
loglog(N,e1)
loglog(N,e2)
legend('final', 'máximo');
hold off;


%3/8
N=[20, 200, 2000, 20000, 200000];
e1=zeros(1,5);
e2=e1;
for n=1:5
  [U,t]=RKexplicito(@ej1, N(n), 0, 10, u0, A38, B38, C38);
  uexac = feval(@exacta, t);   
  e1(n)=max(abs(U(:,N(n)+1)-uexac(:,N(n)+1)));
  e2(n)=max(max(abs(U-uexac)));
end
e38 = e2;

figure(6);
clf;
hold on;
figure(6);
title("Método de 3/8");
loglog(N,e1)
loglog(N,e2)
legend('final', 'máximo');
hold off;

figure(7);
clf;
hold on;
figure(7);
title("Comparaciones");
loglog(N,eeuler)
loglog(N,erunge)
loglog(N,eheun2)
loglog(N,eheun3)
loglog(N,ekutta)
loglog(N,e38)
legend('Euler', 'Runge', 'Heun2', 'Heun3', 'Kutta', '3/8');
hold off;


clear;
