function sol=sol(t)
    sol(1,:)=0.01./(0.01+0.99*exp(-10*t));
end