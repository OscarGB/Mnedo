function y=f2(t,U)
  %size U =2, columna
  %size t = 1
  % y columna de tamaño 2
     y=[-2,1;1,-2]*U+[2*sin(t); 2*(cos(t)-sin(t))];
end