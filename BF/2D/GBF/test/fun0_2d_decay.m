function res = fun0_2d_decay(x,k)

xk = (x(:,1)*k(:,1)' + x(:,2)*k(:,2)');
sx = (2 + sin(2*pi*x(:,1)).*sin(2*pi*x(:,2)))/0.0001;
cx = (2 + cos(2*pi*x(:,1)).*cos(2*pi*x(:,2)))/0.0001;
kr = sqrt(sx.^2*(k(:,1).^2)' + cx.^2*(k(:,2).^2)');

N1 = -min(k(:,1))-1;
N2 = -min(k(:,2))-1;
decay_k1 = (0.5+0.5*cos(pi()-2*pi()*k(:,1)'/N1));
decay_k2 = (0.5+0.5*cos(pi()-2*pi()*k(:,2)'/N2));
decay_k1(abs(k(:,1)')<N1/2) = 1;
decay_k2(abs(k(:,2)')<N2/2) = 1;

decay_k1 = ones(size(x(:,1)))* decay_k1;
decay_k2 = ones(size(x(:,2)))* decay_k2;

decay = decay_k1.*decay_k2;

tmp = (2*pi)* (xk + kr.*decay);

res = exp(1i*tmp);

end
