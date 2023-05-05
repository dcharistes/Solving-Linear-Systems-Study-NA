function pithikasNo2
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

for n = 4:35
m = n-3;
for i = 1:10


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

tic;xG= gaussianElimination(p,b)'; tG=toc;
sumG(m)= sumG(m)+tG;

tic;xP=PTRANSII(n,E,C,D,A,B,y); tP=toc;
sumP(m)= sumP(m)+tP;

tic;xC=cramer(p,b)'; tC=toc;
sumC(m) = sumC(m)+tC;

end

avgC(m)=sumC(m)/i;
avgG(m)=sumG(m)/i;
avgM(m)=sumM(m)/i;
avgP(m)=sumP(m)/i;


end

subplot(2,1,1);
plot(N,avgG,'g--',N,avgC,'Co-',N,avgM,'r*-',N,avgP,'b^-');hold on; grid on;
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

n = 4;
m = 1;
while n <= 2000
    m = ((n-4)/50) + 1;
    n = n + 50;
for i=1:10

E = round(10*rand(1,n-2))+1;
C = round(10*rand(1,n-1))+1;
D = round(10*rand(1,n))+1;
A = round(10*rand(1,n-1))+1;
B = round(10*rand(1,n-2))+1;
p=pentadiagonal(E,C,D,A,B);
y = round(100*rand(1,n)) + 1;
b=y';

tic;xP=PTRANSII(n,E,C,D,A,B,y); tP=toc;
sumP(m)= sumP(m)+tP;

tic;xM=(p\b)'; tM=toc;
sumM(m)= sumM(m)+tM;
end

avgM(m)=sumM(m)/i;
avgP(m)=sumP(m)/i;

end

subplot(2,1,2);
plot(N, avgM, 'r*-',N,avgP, 'bo-'); hold on; grid on;
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

psi(n) = d(n);
    sigma(n) = c(n) / psi(n);
    fi(n) = e(n) / psi(n);
    w(n) = y(n) / psi(n);

    ro(n-1) = a(n-1);
    psi(n-1) = d(n-1) - sigma(n) * ro(n-1);
    sigma(n-1) = (c(n-1) - fi(n) * ro(n-1)) / psi(n-1);
    fi(n-1) = e(n-1) / psi(n-1);
    w(n-1) = (y(n-1) - w(n) * ro(n-1)) / psi(n-1);

    for i = n-2 : -1 : 3
        ro(i) = a(i) - sigma(i+2) * b(i);
        psi(i) = d(i) - fi(i+2) * b(i) - sigma(i+1) * ro(i);
        sigma(i) = (c(i) - fi(i+1) * ro(i)) / psi(i);
        fi(i) = e(i) / psi(i);
        w(i) = (y(i) - w(i+2) * b(i) - w(i+1) *ro(i)) / psi(i);
    end

    ro(2) = a(2) - sigma(4) * b(2);
    psi(2) = d(2) - fi(4) * b(2) - sigma(3) * ro(2);
    sigma(2) = (c(2) - fi(3) * ro(2)) / psi(2);
    ro(1) = a(1) - sigma(3) *b(1);
    psi(1) = d(1) - fi(3) * b(1) - sigma(2) * ro(1);
    w(2) = (y(2) - w(4) * b(2) - w(3) * ro(2)) / psi(2);
    w(1) = (y(1) - w(3) * b(1) - w(2) *ro(1)) / psi(1);

    x(1) = w(1);
    x(2) = w(2) - sigma(2) * x(1);

    for i = 3 : 1 : n
        x(i) = w(i) - (sigma(i) * x(i-1)) - fi(i) * x(i-2);
    end

    disp('psi');
    disp(psi);

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

function x = cramer(A,b)
    d = det(A);
    x = zeros(size(b));
    for j = 1:size(b)
        x(j) = det([A(:,1:j-1) b A(:,j+1:end)]) / d;
  end

    end

