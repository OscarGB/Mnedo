function [U,t]=euler(f,N,t0,T,u0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta funci�n resuelve el problema de valor inicial
% u�=f(t,u)
% u(t0)=u0
% en [t0,T] utilizando el m�todo de Euler
%
% [u,t]=euler(f,N,t0,T,u0)
%
% Variables de Entrada:
%
% f: define un vector columna de dimensi�n u. funci�n que rige el sistema de EDO,
%    tiene dos argumentos f(t,u) donde t es escalar
%    y u vector columna.
% N: n�mero de pasos en los que dividimos el intervalo de
%    integraci�n
% t0: tiempo inicial
% T: tiempo final
% u0: vector columna. Dato inicial
%
% Variables de Salida:
% U: matriz de length(u0) x (N+1) que contiene la soluci�n
% t: vector fila de tiempos de N+1 componentes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% defino el paso
h=(T-t0)/N;
% defino u y t
t=[t0:h:T]; 
U=zeros(length(u0),length(t));  
% identifico el dato inicial
U(:,1)=u0;
% aplicamos el algoritmo de Euler
for n=1:N
    U(:,n+1)=U(:,n)+h*f(t(n),U(:,n));
end
end
















