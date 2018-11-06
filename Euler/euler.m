function y = euler(ld, x, y0)
N = length(x)-1;
y(:,1) = y0;
for n=1:N
    y(:,n+1) = y(:,n)+(x(n+1)-x(n))*feval(ld,x(n),y(:,n));
end

