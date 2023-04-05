
function p = pentadiagonal(ee,cc,dd,bb,aa)
clc;
n=10;
ee = round(10*rand(1,n-2))+1;
cc = round(10*rand(1,n-1))+1;
dd = round(10*rand(1,n))+1;
aa = round(10*rand(1,n-1))+1;
bb = round(10*rand(1,n-2))+1;

p = diag(ee,-2)+diag(cc,-1)+diag(dd,0)+diag(aa,1)+diag(bb,2);


    function [x,psi] = PTRANSII(n,e,c,d,a,b,y)
clc;
n=10;

e = [0 0 ee];
c = [0 cc];
d = dd;
a = [aa 0];
b = [bb 0 0];
y = round(100*rand(n,1)) + 1;


psi(n)=d(n);
s(n)=(c(n)/psi(n));
f(n)=(e(n)/psi(n));
w(n)=(y(n)/psi(n));
r(n-1)=a(n-1);
psi(n-1) = d(n-1) - (s(n)*r(n-1));
s(n-1) = ((c(n-1)- (f(n)*r(n-1)))/psi(n-1)); 
f(n-1) = (e(n-1)/psi(n-1));  
w(n-1) = (y(n-1) - (w(n)*r(n-1)))/psi(n-1);

for i=n-2:-1:3
r(i)=(a(i)-(s(i+2)*b(i)));
psi(i)=(d(i) - (f(i+2)*b(i)) - (s(i+1)*r(i)));
s(i)=((c(i)-(f(i+1)*r(i)))/psi(i));
f(i)=(e(i)/psi(i));
w(i)=((y(i)-(w(i+2)*b(i))-(w(i+1)*r(i)))/psi(i));
end

r(2)=(a(2) - (s(4)*b(2)));
psi(2)=(d(2)-(f(4)*b(2))-(s(3)*r(2)));
s(2)=((c(2)-(f(3)*r(2)))/psi(2));
r(1)=(a(1)-(s(3)*b(1)));
psi(1)=(d(1)-(f(3)*b(1))-(s(2)*r(1)));
w(2)=((y(2)-(w(4)*b(2))-(w(3)*r(2)))/psi(2));
w(1)=((y(1) -(w(3)*b(1))- (w(2)*r(1)))/psi(1));

x(1) = w(1);
x(2) = w(2) -(s(2)*x(1));

%disp("psi:"); 
%disp(psi); 

for i = 3:1:n
x(i) = (w(i) - (s(i)*x(i-1)) - (f(i)*x(i-2)));
end

end

for N=4:35
A= p;
b = randn(N,1);


tic;xM=(A\b)'; tM=toc;
tic;xP=PTRANSII(A,b); tP=toc;
tic;xC=cramer(A,b); tC=toc;
tic;xG= gaussianElimination(A,b); tG=toc;
disp(xP);
eM=norm(xP-xM);
eC=norm(xP-xC);
eG=norm(xP-xG); 

fprintf("TIME WITH MATLAB IS:  %12.10f  AND THE DIFFERENCE BETWEEN PTRANSII AND MATLAB IS:  %20.18f\n",tM,eM);
fprintf("TIME WITH CRAMER IS:  %12.10f  AND THE DIFFERENCE BETWEEN PTRANSII AND CRAMER IS:  %20.18f\n",tC,eC);
fprintf("TIME WITH GAUSSIAN ELIMINATION IS %12.10f AND THE DIFFERENCE BETWEEN PTRANSII AND GAUSSIAN ELIMINATION IS:  %20.18f\n",tG,eG);


end



end