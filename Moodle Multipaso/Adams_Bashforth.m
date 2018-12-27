function [u,t]=Adams_Bashforth(f,tf,t0,N,u0,beta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Esta función resuelve el problema de valor inicial
%       u'=f(t,u)
%       u(t0)=u0
% utilizandodo el método 
%
%          [u,t]=Adams_Bashforth(f,tf,t0,N,u0,beta)
%
%  Variables de Entrada:
%
%       f: vector columna. función que rige el sistema de EDO, 
%          tiene dos argumentos f(t,u) donde t es escalar
%          y u vector columna.   
%       N : número de pasos 
%       t0: tiempo inicial 
%       tf: tiempo final
%       u0: vector columna. Dato inicial
%       beta: coeficientes del polinomio de interpolación       
%       
%  Variables de Salida:
%
%       u: matriz de length(u0) x length(t) que contiene la solución
%       t: vector de tiempos
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=length(beta)
h=(tf-t0)/N;
%
%En la columna 1 guardamos la solución inicial 
%llamamos a un método de un paso para calcular, partiendo de u0,
%u1,u2,..uk mediante el método de Kutta
A=[0,0,0,0; 1/2,0,0,0; 0,1/2,0,0; 0,0,1,0];
c=[0; 1/2;1/2; 1];
b=[1/6,1/3,1/3,1/6];
[u,t]=RKexplicito(f,(k-1)*h+t0,t0,k-1,u0,b,c,A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=t0:h:tf;
F=zeros(length(u0),k);
for i=1:k
    F(:,i)=f(t(i),u(:,i)); 
end 
%
%calculo de la solucion
%
for n=1:N-k+1
      u(:,n+k)=u(:,n+k-1)+h*F*transpose(beta);
      F(:,1:k-1)=F(:,2:k);
      F(:,k)=f(t(n+k),u(:,n+k)); 
end    
end