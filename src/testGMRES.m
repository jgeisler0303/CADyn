m1= 1.0;
m2= 2.0;
k1= 3.0;
k2= 10.0;
k12= 6.0;
d1= 0.5;
d2= 0.8;

x0= zeros(4, 1);
x0(1)= 1;

f= @(t,x)[x(3:4); (-d1*sign(x(3))*abs(x(3))^3 - k1*x(1) - k12*(x(1)-x(2))); (-d2*x(4) - k2*x(2) - k12*(x(2)-x(1)))];
M= [m1 0; 0 m2];
[t, x]= ode23s(@(t, x)blkdiag(eye(2), M^-1)*f(t,x), [0 10], x0);

plot(t, x(:, 1))