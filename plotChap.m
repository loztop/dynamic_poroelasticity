clear all;
close all;
 
 
A=load('plot.txt')
 


%p1=A(:,4);
%p2=A(:,5);
%p3=A(:,6);

%figure;
%plot(p1,'-r');
%hold all;
%plot(p2,'--g');
%hold all;
%plot(p3,'xb');

m1=A(:,1);
m2=A(:,2);
m3=A(:,3);

figure;
plot(m1,'-r');
hold all;
plot(m2,'--g');
hold all;
plot(m3,'xb');