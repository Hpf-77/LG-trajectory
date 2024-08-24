%LG光束全空间模拟

function A = LGlight(l,p,w0,lamd)

syms x;
syms y;
syms z;
syms xy;
k = 2*pi/lamd;                                                    %波矢k
zr = (pi * w0^2) / lamd;                                     %瑞利距离
wz = w0 * sqrt(1 + (z/zr)^2);  
R = z + zr^2 / z;

a0 = sqrt(2/pi);
a = sqrt(factorial(p) / ( factorial ( p+abs(l) ) ) ); 
b = 1/wz;             
c = ( sqrt(2*(x^2+y^2)) / wz )^abs(l);  
d = exp( - (x^2 + y^2 ) / wz^2);


fx = laguerreL(p,abs(l),xy);             
xy = (2 * (x^2 + y^2) / wz^2);
f = eval(fx);

gg = exp(1i *(k*(x^2+y^2)/(2*R) + k*z ));
h = exp(-1i *  (2*p + abs(l) +1) *atan2( z,zr) ); 
A = a0*a * b*c*d*gg*f*h;                                   
end

