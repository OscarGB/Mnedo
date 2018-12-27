function [u,t,no,ha,tol_r]=par_RK(f,tf,t0,u0,bp,bq,p,c,A,rTOL,alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Esta función resuelve el problema de valor inicial
%       u'=f(t,u)
%       u(t0)=u0
% utilizando para la estimación del error local pares encajados p(q) 
% de métodos Runge-Kutta
%
%           [u,t,no,ha,tol_r]=par_RK(f,tf,t0,u0,bp,bq,p,c
%           ,A,rTOL,alpha)
%
%  Variables de Entrada:
%
%       f: vector columna. función que rige el sistema de EDO, 
%          tiene dos argumentos f(t,u) donde t es escalar
%          y u vector columna.   
%       t0: tiempo inicial 
%       tf: tiempo final
%       u0: vector columna. Dato inicial
%       bp,bq,c,A: coeficientes del tablero de BUTCHER.   
%       p: orden del método de avance.
%       rTOL: toleracncia relativa
%       alpha: parámetro de control de paso 
%
%  Variables de Salida:
%
%       u: matriz de length(u0) x length(tiempo) que contiene la solución
%       t: vector de tiempos
%       no: Numero de pasos rechazados
%       ha: tamaño de los pasos aceptados
%       tol_r: valores de la tolerancia utilizados
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u(:,1)=u0;  
r0=10; r1=0.1;
ha=[]; %tamaño de los pasos aceptados

tol_r=[];

s=size(bp,2); %número de pasos de RK
t(1)=t0;
no=0; %pasos fallados


%Selección del paso inicial
tol=rTOL*norm(u0); 
h=min((tf-t0),alpha*(tol/(1+norm(u0)))^(1/(p+1)));

%Empiezo el algoritmo
n=0;
while t(n+1)+h<=tf %para parar cuando llegaos al tiempo final
    n=n+1;
    t(n+1)=t(n)+h;
    % Calculamos las etapas
    K(:,1)=f(t(n),u(:,n));
    for i=2:s
        K(:,i)=f(t(n)+c(i)*h,u(:,n)+h*K(:,1:i-1)*transpose(A(i,1:i-1)));
    end
    
    % Se calcula u_{n+1} 
    u(:,n+1)=u(:,n)+h*K*transpose(bp);
    
    %Calculamos el estimador est
    est= norm(h*K*transpose(bq-bp)); %diferencia entre las aproximaciones
    tol=rTOL*max(norm(u(:,n)),norm(u(:,n+1)));
    %
    %Vemos si aceptamos el valor para u_n o no
    if est >= h*tol %no aceptamos y cambiamos h
      no=no+1;
      %hr(n+1)=h;
      r=max(r1,alpha*(h*tol/est)^(1/p));
      h=h*r;
      n=n-1;
    else %aceptamos pasamos al siguiente paso
      ha(n)=h;%me guardo el h aceptado
      r=min(alpha*(h*tol/est)^(1/(p)),r0);
      h=h*r;%hago el siguiente h un poco más grande
      if (t(n+1)<tf)&(t(n+1)+h>tf)%Ajuste al punto final de integración
            h=tf-t(n+1);
      end     
      tol_r=[tol_r,tol];
  end
  
    
end







end