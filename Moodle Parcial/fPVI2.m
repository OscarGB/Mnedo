function y=fpendulo(t,u)
    y(1)=u(2);
    y(2)=-9.8*sin(u(1));
    y=y';
end