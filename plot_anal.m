close all;
clear all;

I=200;

for i=1:I

y(i)=confined_anal(0,i*0.01)

end

figure;
plot([1:I],y)