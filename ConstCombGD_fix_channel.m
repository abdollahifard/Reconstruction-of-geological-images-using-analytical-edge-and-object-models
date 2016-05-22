function [w,f,c,e]=ConstCombGD_fix_channel(X,P,eta,eps,nmaxiter,X1)
% Gradient descent with logarithmic distance
% inputs :
%     X: an 2*n matrix containing coordinates of points of the patch in a
%          bipolar format (-k:k)
%     P: an 1*n matrix containing the values at points given in X
%     eta: step parameters used in Gradient Descent
%     eps: threshold used in stop condition: norm(dw)/norm(w)<eps (e.g. 0.01)
%     nmaxiter: maximum # of allowed iterations (e.g. 100)

w1=3*randn(3,1);
w2=3*randn(3,1);
w=[w1;w2];
g=@(z)1./(1+exp(-z));
n=size(X,2);
X=[ones(1,n);X];
f1=g(w1'*X);
f2=g(w2'*X);
f=f1.*f2;
cnt=0;
dw1=eta/n*X*((1+exp(-w2'*X))./(1+exp(-w2'*X)+exp((w1'-w2')*X)).*(f-P))';
dw2=eta/n*X*((1+exp(-w1'*X))./(1+exp(-w1'*X)+exp((w2'-w1')*X)).*(f-P))';
dw=[dw1;dw2];
e=0;
while (norm(dw)/norm(w)>eps)&&(cnt<nmaxiter)%
    %subplot(5,5,m);imshow(reshape(f,11,11)>0.5);m=m+1;
    w=w-dw;
    w1=w(1:3);w2=w(4:6);
    f1=g(w(1:3)'*X);
    f2=g(w(4:6)'*X);
    f=f1.*f2;
    dw1=eta/n*X*((1+exp(-w2'*X))./(1+exp(-w2'*X)+exp((w1'-w2')*X)).*(f-P))';
    dw2=eta/n*X*((1+exp(-w1'*X))./(1+exp(-w1'*X)+exp((w2'-w1')*X)).*(f-P))';
    dw=[dw1;dw2];
    cnt=cnt+1;
    e(cnt)=mean(abs(f-P));
end
X1=[ones(1,size(X1,2));X1];
w_t_1=w(1:3);
w_t_2=w(4:6);
c=fix_channel(w_t_1,w_t_2,X1,g);

