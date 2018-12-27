%PVI a)

u0=[2; 3];
t0=0;
tf=10;
N=100;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Adams_Bashforth
beta=1;
[u1,t]=Adams_Bashforth(@fa,tf,t0,N,u0,beta);

 plot(u1(1,:), u1(2,:)); hold all
 
 beta=[-1/2, 3/2];
 [u2,t]=Adams_Bashforth(@fa,tf,t0,N,u0,beta);

 plot(u2(1,:), u2(2,:)); hold all
 
  beta=[5/12, -16/12, 23/12];
 [u3,t]=Adams_Bashforth(@fa,tf,t0,N,u0,beta);

 plot(u3(1,:), u3(2,:)); hold all
 
  beta=[-9/24 37/24 -59/24 55/24];
 [u4,t]=Adams_Bashforth(@fa,tf,t0,N,u0,beta);

 plot(u4(1,:), u4(2,:)); hold all
 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BDF predictor corrector con Adams Bashforth
%
%orden 1
beta=1;
alpha=[-1 1]
[u1,t]=BDF_pc(@fa,tf,t0,N,u0,alpha,beta);
plot(u1(1,:), u1(2,:)); hold all
%
%orden 2
beta=[-1/2, 3/2];
alpha=[1/2 -2 3/2];
[u2,t]=BDF_pc(@fa,tf,t0,N,u0,alpha,beta);
plot(u2(1,:), u2(2,:));
%
%orden 3
beta=[5/12, -16/12 23/12];
alpha=[-1/3 3/2 -3 11/6];
[u3,t]=BDF_pc(@fa,tf,t0,N,u0,alpha,beta);
plot(u3(1,:), u3(2,:));
%
%orden 4
beta=[-9/24 37/24 -59/24 55/24];
alpha=[1/4 -4/3 3 -4 25/12];
[u4,t]=BDF_pc(@fa,tf,t0,N,u0,alpha,beta);
plot(u4(1,:), u4(2,:));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BDF Newton
%
%orden 1
alpha=[-1 1]
[u1,t,it]=BDF_Newton(@fa,@dfa,tf,t0,N,u0,alpha);
plot(u1(1,:), u1(2,:)); hold all
%
%orden 2
alpha=[1/2 -2 3/2];
[u2,t,it]=BDF_Newton(@fa,@dfa,tf,t0,N,u0,alpha);
plot(u2(1,:), u2(2,:));
%
%orden 3
alpha=[-1/3 3/2 -3 11/6];
[u3,t,it]=BDF_Newton(@fa,@dfa,tf,t0,N,u0,alpha);
plot(u3(1,:), u3(2,:));
%
%orden 4
alpha=[1/4 -4/3 3 -4 25/12];
[u4,t,it]=BDF_Newton(@fa,@dfa,tf,t0,N,u0,alpha);
plot(u4(1,:), u4(2,:));