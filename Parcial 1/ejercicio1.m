%Runge Kutta Implicitos
%
%ejemplo de euler
%aplicado al ejercicio 1

%Para graficos en ubuntu (si no me peta) (MAGIA)
graphics_toolkit("gnuplot");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Definicion del método RK
A=[1/2];
b=[1];
c=[1/2];
%valor inicial
u0=[2;3];
t0 = 0
tf = 10
%Método Numérico
[U,t, it]=RKIfuncional(@f1,t0,tf,1000,u0,b,c,A, 1000, 0.1);
%Pintamos resultados del MN
figure(1)
plot(U(1,:),U(2,:));
hold on;
%Sol exacta
sol=f1_exac(t);
%La pintamos en la misma grafica
figure(1)
plot(sol(1,:),sol(2,:));

h = [];
diferencia = [];

for N=1000:5000:10000
  h((N-1000)/100+1)=(tf-t0)/N;
  [U,t, it]=RKIfuncional(@f1,t0,tf,N,u0,b,c,A, 100, 0.001);
  sol = f1_exac(t);
  diferencia(1,(N-1000)/100+1)=norm((U-sol),1);
end

pendiente(1)=mean(diff(log(diferencia(1,:)))./diff(log(h)));
figure(2); 
loglog(h, diferencia); 
title(["La pendiente es = ", num2str(pendiente)]);
hold off;


