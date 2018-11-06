function [u,t]=RKexplicito(f,tf,t0,N,u0,b,c,A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Esta función resuelve el problema de valor inicial
%       u'=f(t,u)
%       u(t0)=u0
% utilizandodo un método Runge-Kutta explicito
%
%           [u,t]=RKexplicito(f,tf,t0,N,u0,b,c,A)
%
%  Variables de Entrada:
%
%       f: vector columna. función que rige el sistema de EDO, 
%          tiene dos argumentos f(t,u) donde t es escalar
%          y u vector columna.   
%       N: número de pasos en los que dividimos el intervalo de
%          integración
%       t0: tiempo inicial 
%       tf: tiempo final
%       u0: vector columna. Dato inicial
%       b,c,A: coeficientes del tablero de BUTCHER.       
%
%  Variables de Salida:
%
%       u: matriz de length(u0) x length(t) que contiene la solución
%       t: vector de tiempos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=(tf-t0)/N;
t=t0:h:tf;
u=zeros(length(u0),N+1);
u(:,1)=u0; %En la columna 1 guardamos la solución inicial 
s=size(b,2);
k=zeros(length(u0),s);

for n=1:N
    k(:,1)=f(t(n)+c(1)*h,u(:,n));%primera fila de A todo ceros
    for i=2:s
        k(:,i)=f(t(n)+c(i)*h,u(:,n)+h*k(:,1:i-1)*transpose(A(i,1:i-1)));
    end
    u(:,n+1)=u(:,n)+h*k*transpose(b);
end

end