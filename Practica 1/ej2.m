graphics_toolkit("gnuplot");

%EULER
u0 = [2,3]; % Valor inicial
[u,t] = euler(@ej1, 100, 0, 10, u0);

%EXACTA
uexac = feval(@exacta, t);

figure(1);
clf;
hold on;
figure(1);
plot(u(1,:),u(2,:),'g');
figure(1);
plot(uexac(1,:),u(2,:),'r');
hold off;
