% Programe for modified prony developed by praveen

% 10/05/09


function [b,a,mu,p] = smprony1(y,t,p)

%MPRONY Fits a sum of exponential functions to data by a modified Prony

% method.

% This KT method for improvement of prony algorithms

%

% simple prony method


y = y(:); t = t(:); % impose column structure

[n cy]=size(y);

dt=t(2)-t(1);


% form Y


Y=zeros(n-p,p+1);

i=1:(n-p);

for j=1:p+1,

Y(:,j)=y(i+j-1);

end;


% using sigular value to modified Y_mod matrix


% [ x d v]=svd(Y);
% 
% 
% % this part of the code impolyes the method proposed in " Dynamic tracking
% 
% % of low frequncy oscillations with improved prony
% 
% sumation=0;
% 
% for i=1:p+1
% 
% sumation=sumation+ d(i,i);
% 
% k(i)=sumation/sum(diag(d));
% 
% if(k(i) >= 0.9)
% 
% kk=i;
% 
% break;
% 
% end
% 
% end
% 
% kk;
% 
% 
% d_mod=d;
% 
% for i=kk+1:p+1
% 
% d_mod(i,i)=0;
% 
% end
% 
% 
% Y_mod=x*d_mod*v';
% 
% p1=kk;


% starting values



[x1 d1]=eig(Y'*Y);

[l jmin]=min(diag(d1));

c1=x1(:,jmin);

%extract rate constants

b=log(roots( c1(p+1:-1:1) ))/dt;

A=exp(t*b.');

a=A\y;

mu=A*a;


