%LG光束傅里叶变换后光束的 全空间模拟

function A = LGlight_p(l,p,w0,lamd)

syms x;
syms y;
syms z;
syms xy;
k = 2*pi/lamd;                                                %波矢k
zr = (pi * w0^2) / lamd;                                 %瑞利距离
wz = 1/(pi*w0);                                            
R = z + zr^2 / z;

a0 = sqrt(2/pi);
a = sqrt(factorial(p) / ( factorial ( p+abs(l) ) ) );  
b = 1/wz*(-1i)^(abs(l));            
c = ( sqrt(2*(x^2+y^2)) / wz )^abs(l);  
d = exp( - (x^2 + y^2 ) / wz^2*(1+1i*(z/zr)));


fx = laguerreL(p,abs(l),xy);            
xy = (2 * (x^2 + y^2) / wz^2);
f = eval(fx);
 
gg = exp(1i *( k*z )); 
A = a0*a * b*c*d*gg*f;   

end

