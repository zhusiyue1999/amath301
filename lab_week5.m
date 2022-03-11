f = @ (p) (p(1) - 2)^2 + (p(2) + 1)^2 + 5*sin(p(1))*sin(p(2)) + 100;
fgrad = @(p)[2*(p(1) - 2) + 5*cos(p(1))*sin(p(2)); 2*(p(2) + 1) + 5*sin(p(1))*cos(p(2))];
fgrad([6; 4]);

p0 = [6;4];
grad = fgrad(p0);
phi = @ (t) p0 - t*grad;
phi(1)

f_of_phi = @(t) f(phi(t));
tmin = fminbnd(f_of_phi, 0, 1);
p1 = phi(tmin)