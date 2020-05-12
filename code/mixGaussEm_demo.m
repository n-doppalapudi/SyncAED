close all; clear;
d = 3;
k = 3;
n = 500;
[X,label] = mixGaussRnd(d,k,n);
plotClass(X,label);

m = floor(n/2);
X1 = X(:,1:m);
y1 = label(:,1:m);
X2 = X(:,(m+1):end);
y2= label(:,(m+1):end);
% train
[z1,model,llh] = mixGaussEm(X1,y1);
figure;
plot(llh);
figure;
plotClass(X1,z1);
% predict
z2 = mixGaussPred(X2,model);
figure;
plotClass(X2,z2);
correct=0;
total=numel(y2);
for i=1:1:total
    if y2(i)==z2(i)
        correct=correct+1;
    end
end
accury=correct/total;

