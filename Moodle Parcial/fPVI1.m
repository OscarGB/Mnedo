
function y=f5(t,v)
%t vector columna
%vector columna
  %size(v,1)=s*m, m=2
  %size t = s
for i=1:2:size(v,1)
    y(i:i+1,1)=[-2,1;998,-999]*[v(i);v(i+1)]+[2*sin(t((i+1)/2)); 999*(cos(t((i+1)/2))-sin(t((i+1)/2)))];
end
end
