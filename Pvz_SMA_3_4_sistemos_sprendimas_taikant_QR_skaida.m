
% QR-skaida 

clc
% A=[1  1  1  1;
%    1 -1 -1  1;
%    2  1 -1  2;
%    3  1  2 -1]
A=[1 2 1 0;
   2 5 0 4;
   14 -8 4 1;
   4 10 0 9]
Ap=A;  % bus naudojama patikrinimui
% b=[2;0;9;7]
b=[-9;-9;-2;-18]
% b=[2 1;0  2;9 3;7 3]
n=size(A,1),  nb=size(b,2)

disp(' tiesioginis zingsnis(atspindziai)')
Q=eye(n);
for i=1:n-1
    z=A(i:n,i);
    zp=zeros(size(z));
    zp(1)=sign(z(1))*norm(z);
    omega=(z-zp); omega=omega/norm(omega);
    Qi=eye(n-i+1)-2*omega*omega'
    A(i:n,:)=Qi*A(i:n,:)
    Q=Q*[eye(i-1),zeros(i-1,n-i+1);zeros(n-i+1,i-1),Qi];
end

Q
Q'*Ap
A

disp('qr skaidos patikrinimas:')
[Qp,Rp]=qr(Ap)


disp('Atvirkstinis zingsnis:')
b1=Q'*b
x=zeros(n,nb);
x(n,:)=b1(n,:)/A(n,n);
for i=n-1:-1:1
    x(i,:)=(b1(i,:)-A(i,i+1:n)*x(i+1:n,:))/A(i,i);
end

x

Ap*x-b
