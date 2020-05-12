function[Untitled,inserted_index]=insertbaddata(Untitled)
e=100;
l=0;
for k=1:79
    Untitled(k+e,:)=0;
inserted_index(k+l)=k+e;
% e=e+1;
% l=l+1;
% Untitled(k+e,:)=0;%Untitled(k+e,:)+0.04;
% inserted_index(k+l)=k+e;
% e=e+1;
% l=l+1;
% Untitled(k+e,:)=0;%Untitled(k+e,:)+0.04;
% inserted_index(k+l)=k+e;
% e=e+1;
% l=l+1;
% Untitled(k+e,:)=0;%Untitled(k+e,:)+0.04;
% inserted_index(k+l)=k+e;
% e=e+1;
% l=l+1;
% Untitled(k+e,:)=0;%Untitled(k+e,:)+0.04;
% inserted_index(k+l)=k+e;
 e=e+100;
end

end
