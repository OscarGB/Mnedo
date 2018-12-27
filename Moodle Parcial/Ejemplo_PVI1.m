

u0=[2;3];
tol=10^(-12);

%Gauss
b=[1];
c=[0.5];
A=[0.5];

[u,t,it]=RKIfuncional(@fPVI1,0,10,1000,u0,b,c,A,100,tol);
plot(u(1,:),u(2,:))
%trapecio

b=[0.5,0.5];
c=[0;1];
A=[0,0;0.5,0.5];

[u,t,it]=RKIfuncional(@fPVI1,0,10,1000,u0,b,c,A,100,tol);
plot(u(1,:),u(2,:))

