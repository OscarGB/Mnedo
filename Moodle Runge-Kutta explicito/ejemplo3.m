%Runge Kutta Explícitos
%
%ejemplo de euler
%aplicado al ejercicio 1
A=[0];
b=[1];
c=[0];
u0=[2;3];
[U,t]=RKexplicito(@f1,10,0,100,u0,b,c,A);
plot(U(1,:),U(2,:))

%Ejemplo de método de Heun de orden 3
%aplicado al ejercicio 1
A=[0, 0,0;1/3,0,0;0,2/3,0];
b=[1/4,0,3/4];
c=[0;1/3;2/3];
u0=[2;3];
[U,t]=RKexplicito(@f1,10,0,100,u0,b,c,A);
plot(U(1,:),U(2,:))
uexac=solexac_PVI1(t);    
e1=max(abs(U(:,101)-uexac(:,101)))
e2=max(max(abs(U-uexac)))
%
%Aplicado al ejercicio 5
[U5,t]=RKexplicito(@f5,10,0,100,u0,b,c,A);
plot(U5(1,:),U5(2,:))