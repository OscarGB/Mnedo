% Ejemplo 1
%
%Resolvemos mediante Euler aproximaciones de la solución del PVI
%
% 	u'=10*u*(1-u) en [0,2]
%	u(0)=0.01
%
%cada una de ellas con un número de pasos distintos
N=[20,200,2000,20000,200000];
e1=zeros(1,5);
e2=e1;
for n=1:5
    [U,t]=euler(@f,N(n),0,2,0.01);
    s=sol(t);    
    e1(n)=abs(U(N(n)+1)-sol(2));
    e2(n)=max(abs(U-s));
end
%
%En e1 evaluamos los errores finales, en e2 los errores totales 
%que luego represeentamos en gráficas logaritmicas
%
loglog(N,e1)
hold on
loglog(N,e2,'r')