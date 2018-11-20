%Runge Kutta Implicitos
%
%ejemplo de euler
%aplicado al ejercicio 1

%Para graficos en ubuntu (si no me peta) (SOLO EN OCTAVE)
% graphics_toolkit("gnuplot");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Definicion del método RK GAUSS
A=[1/2];
b=[1];
c=[1/2];
%valor inicial
u0=[2;3];
t0 = 0;
tf = 10;
%Método Numérico
[U,t, it]=RKIfuncional(@f1,t0,tf,10000,u0,b,c,A, 100, 0.001);
%Pintamos resultados del MN
hold on;
figure(1)
plot(U(1,:),U(2,:));
title("Método de Gauss");
%Sol exacta
sol=f1_exac(t);
%La pintamos en la misma grafica
figure(1)
plot(sol(1,:),sol(2,:));

h = [];
diferencia = [];

% con N < 10000 no funciona (peta por alguna razon) 
for N=5000:5000:100000
  h((N-5000)/5000+1)=(tf-t0)/N;
  [U,t, it]=RKIfuncional(@f1,t0,tf,N,u0,b,c,A, 100, 0.001);
  sol = f1_exac(t);
  diferencia(1,(N-5000)/5000+1)=norm((U-sol),1);
end

pendiente(1)=mean(diff(log(diferencia(1,:)))./diff(log(h)));
figure(2); 
loglog(h, diferencia); 
title(["Método de Gauss: la pendiente es = ", num2str(pendiente)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Definicion del método RK Lobatto IIIA
A=[0,0;1/2,1/2];
c=[0;1];
b=[1/2,1/2];
%valor inicial
u0=[2;3];
t0 = 0;
tf = 10;
%Método Numérico
[U,t, it]=RKIfuncional(@f1,t0,tf,10000,u0,b,c,A, 100, 0.001);
%Pintamos resultados del MN
figure(3)
plot(U(1,:),U(2,:));
title("Método de Lobatto");
%Sol exacta
sol=f1_exac(t);
%La pintamos en la misma grafica
figure(3);
plot(sol(1,:),sol(2,:));

h = [];
diferencia = [];

% con N < 10000 no funciona (peta por alguna razon) 
for N=5000:5000:10000
  h((N-5000)/5000+1)=(tf-t0)/N;
  [U,t, it]=RKIfuncional(@f1,t0,tf,N,u0,b,c,A, 100, 0.001);
  sol = f1_exac(t);
  diferencia(1,(N-5000)/5000+1)=norm((U-sol),1);
end

pendiente(1)=mean(diff(log(diferencia(1,:)))./diff(log(h)));
figure(4);
loglog(h, diferencia);
title(["Método de Lobatto: la pendiente es = ", num2str(pendiente)]);

hold off;

