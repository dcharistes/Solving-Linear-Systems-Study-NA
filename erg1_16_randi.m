function erg1_16_randi()
clc;
np=35;
Np=4:np;
[~,cols]=size(Np);
n=10;
sum_tM=zeros(1,cols); avg_tM=zeros(1,cols);
sum_tG=zeros(1,cols); avg_tG=zeros(1,cols);
sum_tC=zeros(1,cols); avg_tC=zeros(1,cols);
sum_tP=zeros(1,cols); avg_tP=zeros(1,cols);

for N=4:np
k=N-3;
for l=1:3


ee = randi([1,n+1],1,N-2);
cc = randi([1,n+1],1,N-1);
dd = randi([1,n+1],1,N);
aa = randi([1,n+1],1,N-1);
bb = randi([1,n+1],1,N-2);

p=pentadiagonal(ee,cc,dd,aa,bb);

y= randi([1,10*n+1],1,N);

b=y';

tic;xM=(p\b)'; tM=toc;
sum_tM(k)= sum_tM(k)+tM;

tic;xP=PTRANSII(N,ee,cc,dd,aa,bb,y); tP=toc;
sum_tP(k)= sum_tP(k)+tP;

tic;xC=cramer(p,b)'; tC=toc;
sum_tC(k)= sum_tC(k)+tC;

tic;xG= gaussianElimination(p,b)'; tG=toc;
sum_tG(k)= sum_tG(k)+tG;

% fprintf("REPETITION %d - THE RESOLUTION TIME OF THE 5-DIAG SYSTEM WITH MATLAB IS:  %12.10f", l, tM);
% fprintf("REPETITION %d - THE RESOLUTION TIME OF THE 5-DIAG SYSTEM WITH PTRANSII IS:  %12.10f", l, tP);
% fprintf("REPETITION %d - THE RESOLUTION TIME OF THE 5-DIAG SYSTEM WITH CRAMER IS:  %12.10f", l, tC);
% fprintf("REPETITION %d - THE RESOLUTION TIME OF THE 5-DIAG SYSTEM WITH GAUSSIAN ELIMINATION IS:  %12.10f", l, tG);

%disp("xP= "); disp(xP);

eM=norm(xP-xM);
eC=norm(xP-xC);
eG=norm(xP-xG); 
end

avg_tM(k)=sum_tM(k)/l;
avg_tP(k)=sum_tP(k)/l;
avg_tC(k)=sum_tC(k)/l;
avg_tG(k)=sum_tG(k)/l;

fprintf("\nDIMENSION %d\n",N);
fprintf("THE AVERAGE RESOLUTION TIME OF THE 5-DIAG SYSTEM WITH CRAMER IS:  %12.10f\n", avg_tC(k));
fprintf("THE AVERAGE RESOLUTION TIME OF THE 5-DIAG SYSTEM WITH GAUSSIAN ELIMINATION IS:  %12.10f\n", avg_tG(k));
fprintf("THE AVERAGE RESOLUTION TIME OF THE 5-DIAG SYSTEM WITH MATLAB IS:  %12.10f\n", avg_tM(k));
fprintf("THE AVERAGE RESOLUTION TIME OF THE 5-DIAG SYSTEM WITH PTRANSII IS:  %12.10f\n", avg_tP(k));



end

subplot(2,1,1);
plot(Np,avg_tC,'m');
hold on;
title('Execution Time Comparison of algorithms for solving linear pentadiagonal systems');
xlabel('Dimension');
ylabel('Execution Time');
axis([4 cols+3 0 0.03]);
grid on;
plot(Np,avg_tG, 'b');
plot(Np,avg_tP, 'g');
plot(Np,avg_tM, 'r');
legend('Cramer', 'Gaussian Elimination', 'PTRANSII', 'Matlab');
legend('Location','northwest');

np=2000;
Np=4:50:np;
[~,cols]=size(Np);
sum2_tM=zeros(1,cols); avg2_tM=zeros(1,cols);
sum2_tP=zeros(1,cols); avg2_tP=zeros(1,cols);

for N=4:50:np

q=((N-4)/50)+1;

for l=1:3


ee = randi([1,n+1],1,N-2);
cc = randi([1,n+1],1,N-1);
dd = randi([1,n+1],1,N);
aa = randi([1,n+1],1,N-1);
bb = randi([1,n+1],1,N-2);

p=pentadiagonal(ee,cc,dd,aa,bb);

y= randi([1,10*n+1],1,N);

b=y';

tic;xM=(p\b)'; tM=toc;
sum2_tM(q)= sum2_tM(q)+tM;

tic;xP=PTRANSII(N,ee,cc,dd,aa,bb,y); tP=toc;
sum2_tP(q)= sum2_tP(q)+tP;

% fprintf("REPETITION %d - THE RESOLUTION TIME OF THE 5-DIAG SYSTEM WITH MATLAB IS:  %12.10f", l, tM);
% fprintf("REPETITION %d - THE RESOLUTION TIME OF THE 5-DIAG SYSTEM WITH PTRANSII IS:  %12.10f", l, tP);
%disp("xP= "); disp(xP);

eM=norm(xP-xM);

end

avg2_tM(q)=sum2_tM(q)/l;
avg2_tP(q)=sum2_tP(q)/l;

fprintf("\nFOR DIMENSION: %d\n",N);
fprintf("THE AVERAGE RESOLUTION TIME OF THE 5-DIAG SYSTEM WITH MATLAB IS:  %12.10f\n", avg2_tM(q));
fprintf("THE AVERAGE RESOLUTION TIME OF THE 5-DIAG SYSTEM WITH PTRANSII IS:  %12.10f\n", avg2_tP(q));

end

subplot(2,1,2);
plot(Np, avg2_tM, 'r');
hold on;
title('Execution Time Comparison of PTRANSII and Matlab for solving linear pentadiagonal systems');
xlabel('Dimension');
ylabel('Execution Time');
axis([4 cols*50 0 0.15]);
grid on;
plot(Np,avg2_tP, 'b');
legend('Matlab', 'PTRANSII');
legend('Location','northwest');
hold off;

end

function p = pentadiagonal(ee,cc,dd,aa,bb)

p = diag(ee,-2)+diag(cc,-1)+diag(dd,0)+diag(aa,1)+diag(bb,2);

end

function [x,psi] = PTRANSII(n,ee,cc,dd,aa,bb,y)

e = [0 0 ee];
c = [0 cc];
d = dd;
a = [aa 0];
b = [bb 0 0];

psi(n)=d(n);
s(n)=c(n)/psi(n);
f(n)=e(n)/psi(n);
w(n)=y(n)/psi(n);
r(n-1)=a(n-1);
psi(n-1) = d(n-1) - (s(n)*r(n-1));
s(n-1) = (c(n-1)- (f(n)*r(n-1)))/psi(n-1); 
f(n-1) = e(n-1)/psi(n-1);  
w(n-1) = (y(n-1) - (w(n)*r(n-1)))/psi(n-1);

for i=n-2:-1:3
r(i)=a(i)-(s(i+2)*b(i));
psi(i)=d(i) - (f(i+2)*b(i)) - (s(i+1)*r(i));
s(i)=(c(i)-(f(i+1)*r(i)))/psi(i);
f(i)=e(i)/psi(i);
w(i)=(y(i)-(w(i+2)*b(i))-(w(i+1)*r(i)))/psi(i);
end

r(2)=a(2) - (s(4)*b(2));
psi(2)=d(2)-(f(4)*b(2))-(s(3)*r(2));
s(2)=(c(2)-(f(3)*r(2)))/psi(2);
r(1)=a(1)-(s(3)*b(1));
psi(1)=d(1)-(f(3)*b(1))-(s(2)*r(1));
w(2)=(y(2)-(w(4)*b(2))-(w(3)*r(2)))/psi(2);
w(1)=(y(1) -(w(3)*b(1))- (w(2)*r(1)))/psi(1);

x(1) = w(1);
x(2) = w(2) -(s(2)*x(1));

%disp("psi:"); 
%disp(psi);

for i = 3:1:n
x(i) = (w(i) - (s(i)*x(i-1)) - (f(i)*x(i-2)));
end

%disp("x:");
%disp(x);

end

function x = cramer(A,b)
	P = det(A); 
	x = zeros(size(b));
	for j = 1:size(b)
        x(j) = det([A(:,1:j-1) b A(:,j+1:end)]) / P;
    end
end

function x = gaussianElimination(A, b)

	[~, n] = size(A);
	Ag = [A b];
	
	for k = 1:n - 1
		[~, m] = max(abs(Ag(k:n, k)));
		C = Ag(k, :);
		Ag(k, :) = Ag(m + k - 1, :);
		Ag(m + k - 1, :) = C;
		if Ag(k, k) == 0
			error('Matrix A is singular');
		end
		for i = k + 1:n
			r = Ag(i, k) / Ag(k, k);
			Ag(i, k:n + 1) = Ag(i, k:n + 1) - r * Ag(k, k: n + 1);
		end
end
	
	x = zeros(n, 1);
	x(n) = Ag(n, n + 1) / Ag(n, n);
	for k = n - 1:-1:1
		x(k) = (Ag(k, n + 1) - Ag(k, k + 1:n) * x(k + 1:n)) / Ag(k, k);
	end
end