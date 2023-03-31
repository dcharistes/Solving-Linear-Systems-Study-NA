clc;
n = 10;
e = [0 0 round(10*rand(1,n-2))+1];
c = [0 round(10*rand(1,n-1))+1];
d = round(10*rand(1,n))+1;
a = [round(10*rand(1,n-1))+1 0];
b = [round(10*rand(1,n-2))+1 0 0];

p = pendatiagonal(e,c,d,a,b);
y = round(100*rand(n,1)) + 1;

x_cor = (p\y)';

[x] = PTRANSII(n,e,c,d,a,b,y);


