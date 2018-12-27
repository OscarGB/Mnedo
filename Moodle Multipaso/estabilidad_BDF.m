%Regiones de estabilidad multipaso
%
%BDF2
N=1000;
h=2*pi/N;
t=[0:h:2*pi];
x=(3-4*cos(t)+cos(2*t))/2;
y=(4*sin(t)-sin(2*t))/2;
plot(x,y)
%
%BDF3
x=11/6-3*cos(t)+1.5*cos(2*t)-cos(3*t)/3;
y=3*sin(t)-1.5*sin(2*t)+sin(3*t)/3;
hold on
plot(x,y)
plot(0*t,t)