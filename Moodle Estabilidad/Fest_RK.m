function F=Fest_RK(x,y,A,b,c)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calcula la funci�n
%       |F_est(z)|-1=|1+zb^T(I-zA)^{-1}e|-1  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       A,b,c: coeficientes del tablero de BUTCHER. 
%       A: matriz cuadrada de orden s
%       b: vector fila de tama�o s
%       c: vector columna de tama�o s
%
z=x+1i*y;
F=abs(1+z.*b*inv(eye(size(A))-z.*A)*ones(size(c)))-1;

end






