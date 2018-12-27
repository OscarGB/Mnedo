function [u,t]=BDF_pc(f,tf,t0,N,u0,alpha,beta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Esta función resuelve el problema de valor inicial
%       u'=f(t,u)
%       u(t0)=u0
% utilizandodo el método Adams-Bahfort para predecir y BDF para corregir 
%
%          [u,t]=BDF_pc(f,tf,t0,N,u0,alpha,beta)
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
%       alpha: coeficientes de interpolación BDF
%       beta: coeficientes de interpolación Adams-Bashfort      
%       
%  Variables de Salida:
%
%       u: matriz de length(u0) x length(t) que contiene la solución
%       t: vector de tiempos
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=length(alpha)-1; %para saber cuantos datos necesito para iniciar multipaso
h=(tf-t0)/N;
%En la columna 1 guardamos la solución inicial 
%llamamos a un método de un paso para calcular, partiendo de u0,
%u1,u2,..uk
A=[0,0,0,0; 1/2,0,0,0; 0,1/2,0,0; 0,0,1,0];
c=[0; 1/2;1/2; 1];
b=[1/6,1/3,1/3,1/6];
[u,t]=RKexplicito(f,(k-1)*h,t0,k-1,u0,b,c,A);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ha=h/alpha(k+1);
alpha=alpha(1:k)/alpha(k+1);
t=t0:h:tf;
F=zeros(length(u0),k+1);
for i=1:k
    F(:,i)=f(t(i),u(:,i)); 
end 
%
%calculo de la solucion
%
for n=1:N-k+1
      u(:,n+k)=u(:,n+k-1)+h*F(:,1:k)*transpose(beta); %predecir
      F(:,1+k)=f(t(n+k),u(:,n+k));%evaluar
      u(:,n+k)=-u(:,n:n+k-1)*transpose(alpha)+ha*F(:,1+k); %corregir
      F(:,1+k)=f(t(n+k),u(:,n+k));%evaluar
      F(:,1:k)=F(:,2:k+1);
end    
    
 

end