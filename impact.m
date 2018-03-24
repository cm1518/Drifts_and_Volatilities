function a0 = impact(VAR,N)

%calculates the contemporaneous impact matrix A0
%i.e. specifies the vector of shocks (delta) for GIRFs
%Christiane Baumeister
%January 2008

[P,D]=eig(VAR);
        
W=normrnd(0,1,N,N);
[Q,R]=qr(W);
for i=1:N
    if R(i,i)<0
       Q(:,i)=-Q(:,i);
    end
end
        
a0=P*D.^0.5*Q';          %corresponds to 'rdc' in constant-coefficient case
%a0(:,1)=a0(:,1)/a0(1,1); %normalization on the first variable: QO
%a0(:,2)=a0(:,2)/a0(1,2); %           or on the second: PO
%a0(:,3)=a0(:,3)/a0(1,3);
%a0(:,4)=a0(:,4)/a0(1,4);