%Ejemplo 2

%Primero representamos las aproximaciones mediante Euler con  %N=100 (azul), N=1000 (rojo) y la solución exacta (amarillo)  %del PVI del ejercicio 1

[U,t]=euler(@f2,100,0,10,[2;3]);
plot(U(1,:),U(2,:))
[U,t]=euler(@f2,1000,0,10,[2;3]);
hold on
plot(U(1,:),U(2,:),'r')
uexac=solexac_aPVI(t);
plot(uexac(1,:),uexac(2,:),'y')

% Ahora identificamos el error según el tamaño de los pasos

N=[100,1000,10000,100000,1000000];
e1=zeros(1,5);
e2=e1;
for n=1:5
    [U,t]=euler(@f2,N(n),0,10,[2;3]);
    uexac=solexac_aPVI(t);    
    e1(n)=max(abs(U(:,N(n)+1)-uexac(:,N(n)+1)));
    e2(n)=max(max(abs(U-uexac)));
end
h=10./N;
loglog(h,e1)
hold on
loglog(h,e2,'r')
