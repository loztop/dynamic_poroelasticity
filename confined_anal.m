function[y]=confined_anal(z,t)
%z=0.03;
%t=1;


sigma=0.001*10^3;
H=1*10^3;
h=1;
k=10^(-3);

S=0.0;

for n=1:100
   
    S=S+(((-1)^(n))/((n-0.5)^2))*sin( (n-0.5)*pi*((z/h)-1) ) * exp(-((n-0.5)^2)*t*( ((pi^2)*H*k)/(h^2) ) );
    
        
end


y=(sigma/H)*( (z/h) -1 + (2/pi)*S);



