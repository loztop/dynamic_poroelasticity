%clear all;
close all;

N=1153;
for i=1:500
str=strcat('cube_',int2str(i),'.tec')
R(i)=importdata(str,' ',4);
x(i,:)=R(i).data(1:N,4);
y(i,:)=R(i).data(1:N,5);
z(i,:)=R(i).data(1:N,6);

max_x(i)=max(x(i,:));
max_y(i)=max(y(i,:));
max_z(i)=max(z(i,:));

end

figure;
plot(max_x,'--');
ylabel('Max height');
xlabel('Time step 1-50 dt=0.1, 51-100 dt=0.002,100-500 dt=0.05');

hold all;
%plot(max_y,'-x');
hold all;
%plot(max_z,'o');
hold all;

M=length(max_x);

time=0;
for j=1:M
            dt=0.1;
   if((j+1)/M>0.1)
      dt=0.02;
   end
   if  ((j+1)/M>0.2)
      dt=0.05;
      
   end 
  
   time=time+dt
          C(j)=time ;      

end

figure;
plot(C,max_x);