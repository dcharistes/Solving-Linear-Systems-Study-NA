function untitled
clc;

N=4:35;
sumM = []; 
sumG = [];
sumC = [];
sumP = [];
avgM = [];
avgG = [];
avgC = [];
avgP = [];
for i = 1:size(N, 2) 
    sumM(i) = 0; 
    sumG(i) = 0; 
    sumC(i) = 0; 
    sumP(i) = 0; 
    avgM(i) = 0;
    avgG(i) = 0;
    avgC(i) = 0;
    avgP(i) = 0;
end

for n=4:35
    m=n-3;
for i=1:10


E = round(10*rand(1,n-2))+1;
C = round(10*rand(1,n-1))+1;
D = round(10*rand(1,n))+1;
A = round(10*rand(1,n-1))+1;
B = round(10*rand(1,n-2))+1;

p=pentadiagonal(E,C,D,A,B);

y = round(100*rand(1,n)) + 1;

b = y';

tic;xM=(p\b)'; tM=toc;
sumM(m)= sumM(m)+tM;

tic;xP=PTRANSII(n,E,C,D,A,B,y); tP=toc;
sumP(m)= sumP(m)+tP;

tic;xG= gaussianElimination(p,b)'; tG=toc;
sumG(m)= sumG(m)+tG;

tic;xC=cramer(p,b)'; tC=toc;
sumC(m) = sumC(m)+tC;

end

avgM(m)=sumM(m)/i;
avgP(m)=sumP(m)/i;
avgC(m)=sumC(m)/i;
avgG(m)=sumG(m)/i;

end

subplot(2,1,1);
plot(N,avgG,'g',N,avgC,'C',N,avgM,'r',N,avgP,'b');hold on; grid on;
axis([0,35,0,2*(10^-3)]);
hold off;


N=4:50:2000;
sumM = []; 
sumP = [];
avgM = [];
avgP = [];
for i = 1:size(N, 2) 
    sumM(i) = 0; 
    sumP(i) = 0; 
    avgM(i) = 0;
    avgP(i) = 0;
end

for n=4:50:2000
    m=((n-4)/50)+1;
for i=1:10

E = round(10*rand(1,n-2))+1;
C = round(10*rand(1,n-1))+1;
D = round(10*rand(1,n))+1;
A = round(10*rand(1,n-1))+1;
B = round(10*rand(1,n-2))+1;
p=pentadiagonal(E,C,D,A,B);

y = round(100*rand(1,n)) + 1;

b=y';

tic;xM=(p\b)'; tM=toc;
sumM(m)= sumM(m)+tM;

tic;xP=PTRANSII(n,E,C,D,A,B,y); tP=toc;
sumP(m)= sumP(m)+tP;

end

avgM(m)=sumM(m)/i;
avgP(m)=sumP(m)/i;

end

subplot(2,1,2);
plot(N, avgM, 'r',N,avgP, 'b'); hold on; grid on;
axis([0,2000,0,0.08]);
hold off;

end
function p = pentadiagonal(E,C,D,A,B)
p = diag(E,-2)+diag(C,-1)+diag(D,0)+diag(A,1)+diag(B,2);
end

function [x,psi] = PTRANSII(n,E,C,D,A,B,y)

A(n) = 0;
B(n) = 0;
B(n-1) = 0;
C(1) = 0;
E(1) = 0;
E(2) = 0;

e = [E(1) E(2) E];
c = [C(1) C];
d = D;
a = [A A(n)];
b = [B B(n-1) B(n)];


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

% disp("psi:");
% disp(psi);

for i = 3:1:n
x(i) = w(i) - (s(i)*x(i-1)) - (f(i)*x(i-2));
end

% disp("x:");
% disp(x);

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
