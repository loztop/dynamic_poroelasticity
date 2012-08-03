
E=10;
nu=0.3;
       
mu = E / (2.0 * (1.0 + nu));
lambda = E * nu / ((1 + nu) * (1 - 2 * nu));

 C=0.5;

 
 for i=1:10
X1(i)=(i/10)*1.5
val(i) = 2*C*(0.5*lambda + mu)*(1+(1+2*C*X1(i))^(-3))%+2*(2*lambda*C + 4*C*C*lambda*X1(i));
 end
 
 close all;
 plot(X1,val);