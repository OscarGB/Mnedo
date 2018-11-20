
function y=f1(t,v)
y = [];
for N=1:1:length(t)
    y = [y; [-2,1;998,-999]*v(2*N-1:2*N)+[2*sin(t(N)); 999*(cos(t(N))-sin(t(N)))]];
end
end
