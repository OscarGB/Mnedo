%Practica Paso Variable
t0=0;
tf=10;
rTOL=1d-8;
alpha=1;

%Par de Ceschino 2(4)
p_C=2;
A_C= [0,0,0,0; 1/4,0,0,0; 0,1/2,0,0;1,-2,2,0];
c_C=[0;1/4;1/2;1];
bp_C=[1,-2,2,0];
bq_C=[1/6,0,4/6,1/6];

u0=[2; 3];
[u,t,no,ha,tol_r]=par_RK2(@fa,tf,t0,u0,bp_C,bq_C,p_C,c_C,A_C,rTOL,alpha);
no
figure(1); plot(u(1,:),u(2,:),'-');
figure(2);plot(ha);
N=size(ha,2)
PorcentajeRechazados=100*no/(N+no)
Minimoh=min(ha(1:N-1))
Maximoh=max(ha(1:N))
u_exac=solexac_PVI1(t);
figure(3);plot(t,abs(u_exac-u))