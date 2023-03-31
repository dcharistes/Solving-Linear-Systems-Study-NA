function p = pendatiogonal(~,~,~,~,~)
clc;
n = randi([4,10],1,1);
%n=7;
ee = [round(10*rand(1,n-2))+1];
cc = [round(10*rand(1,n-1))+1];
dd = round(10*rand(1,n))+1;
aa = [round(10*rand(1,n-1))+1];
bb = [round(10*rand(1,n-2))+1];

p = diag(dd,0)+diag(aa,1)+diag(bb,2)+diag(cc,-1)+diag(ee,-2);

end

function [x,psi] = ptransii(~,~,~,~,~,~,~)

e = [0 0 round(10*rand(1,n-2))+1];
c = [0 round(10*rand(1,n-1))+1];
d = round(10*rand(1,n))+1;
a = [round(10*rand(1,n-1))+1 0];
b = [round(10*rand(1,n-2))+1 0 0];

pendatiogonal(e,c,d,a,b);
disp(pendatiogonal);
end


