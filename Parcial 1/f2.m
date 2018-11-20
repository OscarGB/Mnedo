
function y=f2(t,v)
g = 9.81;

% dx/dt = v
% dv/dt = -9.8sen(x)
% v = INTEGRAL(-9.8sen(x)dt)
% v = -9.8tsen(x)
% x = -(9.8/2)t**2sen(x)

% x = v(1) (angulo) 
% v = v(2) (velocidad)

y = [v(2); -g*sin(v(1))];

end
