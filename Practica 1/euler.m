function [u,t]=euler(f,N,t0,T,u0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta función resuelve el problema de valor inicial
% u'=f(t,u)
% u(t0)=u0
% utilizando el método de Euler
% u(n+1) = u(n) + hf(t(n),u(n))
% [u,t]=euler(N,t0,T,u0)
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
% u0: vector columna. Dato inicial
%
% Variables de Salida:
% u: matriz de length(u0) x N que contiene la solución
% t: vector de tiempos
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  h = (T-t0)/N;
  u(:,1) = u0;
  t = [t0:h:T];
  for n=1:N
    u(:,n+1) = u(:,n) + h*feval(f,t(n),u(:,n));
  end  
end