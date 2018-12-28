function [u,t]=RKexplicito(f,N,t0,T,u0,A,b,c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta función resuelve el problema de valor inicial
% u'=f(t,u)
% u(t0)=u0
% utilizandodo un método Runge-Kutta explícito
%
% [u,t]=RKexplicito(f,N,t0,T,u0,A,b,c)
%
% Variables de Entrada:
%
% f: vector columna. función que rige el sistema de EDO,
% tiene dos argumentos f(t,u) donde t es escalar
% y u vector columna.
% N: número de pasos en los que dividimos el intervalo de
% integración
% t0: tiempo inicial
% T: tiempo final
% u0: vector columna, dato inicial
% b,c,A: coeficientes del tablero de BUTCHER.
% A: matriz cuadrada de orden s
% c: vector columna de orden s
% b: vector fila de orden s
%
% Variables de Salida:
%
% u: matriz de length(u0) x length(t) que contiene la solución
% t: vector de tiempos
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  h = (T-t0)/N;
  u(:,1) = u0;
  t = [t0:h:T];
  s = length(c);
  k = zeros(length(u0),s);
  for n=1:N
    k(:,1) = feval(f, t(n)+c(1)*h, u(:,n));
    for i=2:s
      sum = zeros(length(u0),1);
      for j=1:i-1
        sum = sum + A(i,j)*k(:,j);
      end
      k(:,i) = feval(f, t(n)+c(i)*h, u(:,n) + h*sum);
    end
    sum = zeros(length(u0),1);
    for i=1:s
      sum = sum + b(i)*k(:,i);
    end
    u(:,n+1) = u(:,n) + h*sum;
  end  
end
