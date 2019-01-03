function [u,t,it]=RKIfuncional(f,t0,T,N,u0,b,c,A,itmax,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta función resuelve el problema de valor inicial
%       u'=f(t,u)
%       u(t0)=u0
% utilizandodo un método Runge-Kutta ímplicito con la iteración 
% funcional (IF)
%
%           [u,t,it]=RKIfuncional(f,t0,T,N,u0,b,c,A,itmax,tol)
%
%  Variables de Entrada:
%
%       f: vector columna. función que rige el sistema de EDO, 
%          tiene dos argumentos f(t,u) donde t es un vector       
%          columna y u vector columna.   
%       t0: tiempo inicial 
%       tf: tiempo final
%       N : número de pasos 
%       u0: vector columna. Dato inicial
%       b,c,A: coeficientes del tablero de BUTCHER. 
%              A matriz cuadrada de orden s
%              b vector fila de tamaño s
%              c vector columna de tamaño s
%       tol: tolerancia para las iteraciones
%       itmax: Numero máximo de iteraciones
%       
%  Variables de Salida:
%
%       u: matriz de length(u0)·length(t) que contiene la         
%         solución
%       t: vector de tiempos
%       it: Vector que contiene el numero de iteraciones         
%           utilizadas, para cada paso de tiempo en la iteración
%           funcional
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=(T-t0)/N;
t=t0:h:T;
d=length(u0);
u=zeros(d,length(t));
u(:,1)=u0; %En la columna 1 guardamos la solución inicial 
%
%Aplicamos el algoritmo Runge-Kutta
s=length(b); %número de pasos
e=ones(s,1);
I=eye(d);
it=zeros(N,1);
for n=1:N
    tt=t(n)+c*h;% vector columna de tamaño s
    V=kron(e,u(:,n)); %vector columna de tamaño s·d
    error=1;
    while (it(n)<itmax && error>tol)
        it(n)=it(n)+1;
        U=V;
        V=kron(e,u(:,n))+h*kron(A,I)*f(tt,U); %atencion
        error=norm(V-U)/norm(V);        
    end    
    u(:,n+1)=u(:,n)+h*kron(b,I)*f(tt,V);
end
end