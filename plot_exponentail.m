close all
time=[0:0.001:0.1]

plot(time,(1-exp(-(time.^(1))./0.01)),'-x')