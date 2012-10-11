clear all;
close all;

A=load_gmsh('Move_inZero.msh');
a=A.POS;
B=load_gmsh('Move_inXYZ.msh');
b=B.POS;


ars=a(:,1).^2 +a(:,2).^2 +a(:,3).^2 ;
brs=b(:,1).^2 +b(:,2).^2 +b(:,3).^2 ;

[inner_idx,val]=find(((a(:,1).^2 +a(:,2).^2 +a(:,3).^2)<9.01)&((a(:,1).^2 +a(:,2).^2 +a(:,3).^2)>8.99)& (a(:,1).^2<0.001));
[outer_idx,val]=find(((a(:,1).^2 +a(:,2).^2 +a(:,3).^2)>80.99)& (a(:,1).^2<0.001));



close all;
plot3(a(inner_idx,1),a(inner_idx,2),a(inner_idx,3),'x');
hold on;
plot3(b(inner_idx,1),b(inner_idx,2),b(inner_idx,3),'ro');

figure;
plot(a(inner_idx,2),a(inner_idx,3),'x');
hold on;
plot(b(inner_idx,2),b(inner_idx,3),'ro');
hold on;
plot(a(outer_idx,2),a(outer_idx,3),'x');
hold on;
plot(b(outer_idx,2),b(outer_idx,3),'ro');

arin=(mean(ars(inner_idx)))^0.5
arout=(mean(ars(outer_idx)))^0.5

brout=(mean(brs(outer_idx)))^0.5
brin=(mean(brs(inner_idx)))^0.5


vol_true=(4/3)*pi*(arout^3-arin^3)
vol_after=(4/3)*pi*(brout^3-brin^3)


exact_inner_rad_after=((3/(4*pi))*((4/3)*pi*(brout^3)-vol_true))^(1/3)

volume_error=vol_true-vol_after
radius_error=exact_inner_rad_after-brin

