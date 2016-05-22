function [w,f,c]=UnconstSingleGD(X,P,eta,eps,nmaxiter,X1)
% Gradient descent with logarithmic distance
% inputs :
%     X: an 2*n matrix containing coordinates of points of the patch in a
%          bipolar format (-k:k)
%     P: an 1*n matrix containing the values at points given in X
%     eta: step parameters used in Gradient Descent
%     eps: threshold used in stop condition: norm(dw)/norm(w)<eps (e.g. 0.01)
%     nmaxiter: maximum # of allowed iterations (e.g. 100)
w=15*randn(3,1);
% w=[0;0;0];
g=@(z)1./(1+exp(-z));
n=size(X,2);
X=[ones(1,n);X];
f=g(w'*X);
cnt=0;
dw=(eta/n)*X*(f-P)';
while (norm(dw)/norm(w)>eps)&&(cnt<nmaxiter)
    %subplot(5,5,m);imshow(reshape(f,11,11)>0.5);m=m+1;
    w=w-dw;
    f=g(w'*X);
    dw=(eta/n)*X*(f-P)';
    cnt=cnt+1;
    %e(cnt)=mean(abs(f-P));
end
X1=[ones(1,size(X1,2));X1];
c=g(w'*X1);
    