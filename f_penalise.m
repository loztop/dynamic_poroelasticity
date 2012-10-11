J=[0:0.001:1];

f=2.*(J-1-log(J))./((J-1).^2);
close all;
%figure;
%plot(J,f);

k=100000;
pzero=0.5
w=k.*(J-pzero).^(6)
close all;
figure;
plot(J,w);