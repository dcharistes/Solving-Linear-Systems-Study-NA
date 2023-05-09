function erg1_18
clc;
n1=35;
n2=2000;

fprintf("\nIMPLEMENTATION OF THE FIRST STUDY EXPERIMENT TO SOLVE RANDOM PENTADIAGONAL LINEAR SYSTEMS...\n");
t_start1=tic;

Np=4:n1;
[~,cols]=size(Np);
sum_tM=zeros(1,cols); avg_tM=zeros(1,cols);
sum_tG=zeros(1,cols); avg_tG=zeros(1,cols);
sum_tC=zeros(1,cols); avg_tC=zeros(1,cols);
sum_tP=zeros(1,cols); avg_tP=zeros(1,cols);

for n=4:n1
   k=n-3;
for l=1:10

e = round(10*rand(1,n-2))+1;
c = round(10*rand(1,n-1))+1;
d = round(10*rand(1,n))+1;
a = round(10*rand(1,n-1))+1;
b = round(10*rand(1,n-2))+1;

p=pentadiagonal(e,c,d,a,b);

y = round(100*rand(1,n)) + 1;

z = y';

tic;xM=(p\z)'; tM=toc;
sum_tM(k)= sum_tM(k)+tM;

tic;xP=PTRANSII(n,e,c,d,a,b,y); tP=toc;
sum_tP(k)= sum_tP(k)+tP;

tic;xC=cramer(p,z)'; tC=toc;
sum_tC(k)= sum_tC(k)+tC;

tic;xG= gaussianElimination(p,z)'; tG=toc;
sum_tG(k)= sum_tG(k)+tG;

end

avg_tM(k)=sum_tM(k)/l;
avg_tP(k)=sum_tP(k)/l;
avg_tC(k)=sum_tC(k)/l;
avg_tG(k)=sum_tG(k)/l;

% fprintf("\nDIMENSION %d\n",n); //If uncommented, the average execution time is printed for every dimension 
% fprintf("THE AVERAGE RESOLUTION TIME OF THE 5-DIAG SYSTEM WITH GAUSSIAN ELIMINATION IS:  %12.10f\n", avg_tG(k));
% fprintf("THE AVERAGE RESOLUTION TIME OF THE 5-DIAG SYSTEM WITH CRAMER IS:  %12.10f\n", avg_tC(k));
% fprintf("THE AVERAGE RESOLUTION TIME OF THE 5-DIAG SYSTEM WITH MATLAB IS:  %12.10f\n", avg_tM(k));
% fprintf("THE AVERAGE RESOLUTION TIME OF THE 5-DIAG SYSTEM WITH PTRANSII IS:  %12.10f\n", avg_tP(k));

end
t1=toc(t_start1);
subplot(2,1,1);
plot(Np,avg_tG,'g.-');hold on; grid on;
plot(Np,avg_tC,'c.-');
plot(Np,avg_tM,'r.-');
plot(Np,avg_tP,'b.-');
axis([0,n1+0.2,0,avg_tG(k)+10^-4]);

title('Execution Time Comparison of algorithms for solving linear pentadiagonal systems');
xlabel('Dimension');
ylabel('Execution Time');
legend('GAUSS', 'CRAMER', 'MATLAB', 'PTRANSII','Location','northwest');
hold off;
fprintf("\nTHE IMPLEMENTATION OF THE FIRST STUDY EXPERIMENT IS FINISHED AND THE RESULTS CAN BE SEEN IN THE FIGURE.\n");
fprintf("Execution Time of the first experiment: %12.10f\n",t1);
pause(1);
fprintf("\nIMPLEMENTATION OF THE SECOND STUDY EXPERIMENT TO SOLVE RANDOM PENTADIAGONAL LINEAR SYSTEMS...\n");

Np=4:50:n2;
[~,cols]=size(Np);
sum2_tM=zeros(1,cols); avg2_tM=zeros(1,cols);
sum2_tP=zeros(1,cols); avg2_tP=zeros(1,cols);

t_start2=tic;
for n=4:50:n2
    q=((n-4)/50)+1;
for l=1:10

e = round(10*rand(1,n-2))+1;
c = round(10*rand(1,n-1))+1;
d = round(10*rand(1,n))+1;
a = round(10*rand(1,n-1))+1;
b = round(10*rand(1,n-2))+1;

p=pentadiagonal(e,c,d,a,b);

y = round(100*rand(1,n)) + 1;

z=y';

tic;xM=(p\z)'; tM=toc;
sum2_tM(q)= sum2_tM(q)+tM;

tic;xP=PTRANSII(n,e,c,d,a,b,y); tP=toc;
sum2_tP(q)= sum2_tP(q)+tP;
%xN=norm(xM-xP); disp("xN= "+xN);

end

avg2_tM(q)=sum2_tM(q)/l;
avg2_tP(q)=sum2_tP(q)/l;

% fprintf("\nDIMENSION: %d\n",n); //If uncommented, the average execution time is printed for every dimension 
% fprintf("THE AVERAGE RESOLUTION TIME OF THE 5-DIAG SYSTEM WITH MATLAB IS:  %12.10f\n", avg2_tM(q));
% fprintf("THE AVERAGE RESOLUTION TIME OF THE 5-DIAG SYSTEM WITH PTRANSII IS:  %12.10f\n", avg2_tP(q));


end
t2=toc(t_start2);
subplot(2,1,2);
plot(Np, avg2_tM, 'r.-'); hold on; grid on;
plot(Np,avg2_tP, 'b.-');
axis([0,n2,0,avg2_tM(q)+0.01]);

title('Execution Time Comparison of PTRANSII and Matlab for solving linear pentadiagonal systems');
xlabel('Dimension');
ylabel('Execution Time');
legend('MATLAB', 'PTRANSII','Location','northwest');
hold off;
fprintf("\nTHE IMPLEMENTATION OF THE SECOND STUDY EXPERIMENT IS FINISHED AND THE RESULTS CAN BE SEEN IN THE SECOND FIGURE.\n");
fprintf("Execution Time of the second experiment: %12.10f\n",t2);
end
function p = pentadiagonal(e,c,d,a,b)

p = diag(e,-2)+diag(c,-1)+diag(d)+diag(a,1)+diag(b,2);

end

function [x,psi] = PTRANSII(n,e,c,d,a,b,y)
%clc;
e = [0 0 e];
c = [0 c];
a = [a 0];
b = [b 0 0];

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

%disp("psi:"); //if uncommented, the psi is printed for every dimension(x10).
%disp(psi);

for i = 3:1:n
x(i) = w(i) - (s(i)*x(i-1)) - (f(i)*x(i-2));
end

%disp("x:"); //if uncommented, the x is printed for every dimension(x10). 
%disp(x);

end

function x = cramer(A,b)
	d = det(A); 
	x = zeros(size(b));
    for j = 1:size(b)
		x(j) = det([A(:,1:j-1) b A(:,j+1:end)]) / d;
    end
end

function x = gaussianElimination(A, b)

	[~, n] = size(A);
	Ag = [A b];
	
	for k = 1:n - 1
		[~, j] = max(abs(Ag(k:n, k)));
		C = Ag(k, :);
		Ag(k, :) = Ag(j + k - 1, :);
		Ag(j + k - 1, :) = C;
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
