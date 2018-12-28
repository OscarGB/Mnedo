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

%EULER
u0 = [2,3]; % Valor inicial
[u1,t] = euler(@ej1, 100, 0, 10, u0);
[u2,t] = euler(@ej1, 1000, 0, 10, u0);
[u3,t] = euler(@ej1, 10000, 0, 10, u0);

%EXACTA
uexac = feval(@exacta, t);

figure(1);
clf;
hold on;
figure(1);
title("Método de Euler");
plot(uexac(1,:),uexac(2,:));
plot(u1(1,:),u1(2,:));
plot(u2(1,:),u2(2,:));
plot(u3(1,:),u3(2,:));
legend('exacta', '100', '1000', '10000');
hold off;

%RK Runge
u0 = [2,3]; % Valor inicial
[u1,t] = RKexplicito(@ej1, 100, 0, 10, u0, ARunge, BRunge, CRunge);
[u2,t] = RKexplicito(@ej1, 1000, 0, 10, u0, ARunge, BRunge, CRunge);
[u3,t] = RKexplicito(@ej1, 10000, 0, 10, u0, ARunge, BRunge, CRunge);

figure(2);
clf;
hold on;
figure(2);
title("Método de Runge");
plot(uexac(1,:),uexac(2,:));
plot(u1(1,:),u1(2,:));
plot(u2(1,:),u2(2,:));
plot(u3(1,:),u3(2,:));
legend('exacta', '100', '1000', '10000');
hold off;

%RK Heun2
u0 = [2,3]; % Valor inicial
[u1,t] = RKexplicito(@ej1, 100, 0, 10, u0, AHeun2, BHeun2, CHeun2);
[u2,t] = RKexplicito(@ej1, 1000, 0, 10, u0, AHeun2, BHeun2, CHeun2);
[u3,t] = RKexplicito(@ej1, 10000, 0, 10, u0, AHeun2, BHeun2, CHeun2);

figure(3);
clf;
hold on;
figure(3);
title("Método de Heun2");
plot(uexac(1,:),uexac(2,:));
plot(u1(1,:),u1(2,:));
plot(u2(1,:),u2(2,:));
plot(u3(1,:),u3(2,:));
legend('exacta', '100', '1000', '10000');
hold off;

%RK Heun3
u0 = [2,3]; % Valor inicial
[u1,t] = RKexplicito(@ej1, 100, 0, 10, u0, AHeun3, BHeun3, CHeun3);
[u2,t] = RKexplicito(@ej1, 1000, 0, 10, u0, AHeun3, BHeun3, CHeun3);
[u3,t] = RKexplicito(@ej1, 10000, 0, 10, u0, AHeun3, BHeun3, CHeun3);

figure(4);
clf;
hold on;
figure(4);
title("Método de Heun3");
plot(uexac(1,:),uexac(2,:));
plot(u1(1,:),u1(2,:));
plot(u2(1,:),u2(2,:));
plot(u3(1,:),u3(2,:));
legend('exacta', '100', '1000', '10000');
hold off;

%RK Kutta
u0 = [2,3]; % Valor inicial
[u1,t] = RKexplicito(@ej1, 100, 0, 10, u0, AKutta, BKutta, CKutta);
[u2,t] = RKexplicito(@ej1, 1000, 0, 10, u0, AKutta, BKutta, CKutta);
[u3,t] = RKexplicito(@ej1, 10000, 0, 10, u0, AKutta, BKutta, CKutta);

figure(5);
clf;
hold on;
figure(5);
title("Método de Kutta");
plot(uexac(1,:),uexac(2,:));
plot(u1(1,:),u1(2,:));
plot(u2(1,:),u2(2,:));
plot(u3(1,:),u3(2,:));
legend('exacta', '100', '1000', '10000');
hold off;

%RK 38
u0 = [2,3]; % Valor inicial
[u1,t] = RKexplicito(@ej1, 100, 0, 10, u0, A38, B38, C38);
[u2,t] = RKexplicito(@ej1, 1000, 0, 10, u0, A38, B38, C38);
[u3,t] = RKexplicito(@ej1, 10000, 0, 10, u0, A38, B38, C38);

figure(6);
clf;
hold on;
figure(6);
title("Método de 3/8");
plot(uexac(1,:),uexac(2,:));
plot(u1(1,:),u1(2,:));
plot(u2(1,:),u2(2,:));
plot(u3(1,:),u3(2,:));
legend('exacta', '100', '1000', '10000');
hold off;

clear;