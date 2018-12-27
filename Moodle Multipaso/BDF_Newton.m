function [u,t,it]=BDF_Newton(f,df,tf,t0,N,u0,alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Esta función resuelve el problema de valor inicial
%       u'=f(t,u)
%       u(t0)=u0
% utilizandodo el método 
%
%          [u,t]=BDF_Newton(f,df,tf,t0,N,u0,alpha)
%
%  Variables de Entrada:
%
%       f: vector columna. función que rige el sistema de EDO, 
%          tiene dos argumentos f(t,u) donde t es escalar
%          y u vector columna.
%       df: matriz cuadrada. jacobiano de la función que rige las EDOs, 
%          tiene dos argumentos df(t,u) donde t es escalar
%          y u vector columna.  

%       N : número de pasos 
%       t0: tiempo inicial 
%       tf: tiempo final
%       u0: vector columna. Dato inicial
%       alpha: coeficientes  de interpolación de BDF     
%       
%  Variables de Salida:
%
%       u: matriz de length(u0) x length(t) que contiene la solución
%       t: vector de tiempos
%       it: vector de iteraciones de Newton para cada paso
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
t=t0:h:tf; 
d=length(u0);
it=zeros(N-k+1,1);
I=eye(d);

tol=1d-16 %tolerancia

%
ha=h/alpha(k+1);
alpha=alpha(1:k)/alpha(k+1);


for n=1:N-k+1
    
   %Para resolver el sistema implicito utilizamos Newton-
   %inicio con el valor antiguo
   U=u(:,n+k-1);
   
   A=u(:,n:n+k-1)*transpose(alpha(1:k));
   norma=tol+1;
   
    while (norma>tol && it(n)< 100) 
        it(n)=it(n)+1;
        %auxU=U; %la guardamos para calcular luego la norma
        F=U-ha*f(t(n+k),U)+A;
        DF=I-ha*df(t(n+k),U);
        %Resolvemos el sistema calculando Delta
        E=DF\F;
        U=U-E;
        norma=norm(E)/norm(U);
    end    
    u(:,n+k)=U;
end
end