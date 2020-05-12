function [baddata_index,a,b,c]=zip_estimation_new(Untitled,window_size,handles)

%watchon
format long;
vm_col=1; va_col=2; im_col=3; ia_col=4;

wndo_increase=30;
lent=window_size;
swng=0.0015;%3;
pronythreshold=[-4.27000000000000;-3.75000000000000;-3.43000000000000;-3.24000000000000;-3.07000000000000;-2.95000000000000;-2.83000000000000;-2.75000000000000;-2.67000000000000;-2.62000000000000;-2.55000000000000;-2.50000000000000;-2.45000000000000;-2.40000000000000];

baddata_index=1;
senP=zeros(length(Untitled),1);

%% bad data detection algo

[baddata_index]=Chebyshev_Regression(Untitled,handles);
baddata_index=baddata_index';
if baddata_index(1)==0
    baddata_index(1)=[];
end

msgString = {'ZIP is being computed'}; % insert a white space
msg = strcat(msgString);
set(handles.StatusText,'String',msg);

Untitled(:,5:6)=0;
for i=1:length(baddata_index)
    if baddata_index(1)~=0
        Untitled(baddata_index(i),6)=1;
    end
end
t=(0:length(Untitled(:,1))-1)/30; %60;%
%untitled= Untitled;

bad_data_time=baddata_index/30;

%Untitled(baddata_index(1:length(baddata_index)),:)=[];
thrshold =[16.7000000000000e-05;
16.7000000000000e-05;
16.7000000000000e-05;
16.7000000000000e-05;
16.7000000000000e-05;
16.7000000000000e-05;
16.7000000000000e-05;
16.7000000000000e-05;
16.7000000000000e-05;
16.7000000000000e-05;
16.7000000000000e-05;
16.7000000000000e-05;
16.7000000000000e-05;
16.7000000000000e-05];
wind_diff = 1200;
base_voltage= Untitled(1,vm_col);%7.2;%3.087421562500000e+05; %7.2;%
P=abs(sqrt(3)*1*Untitled(:,vm_col).*Untitled(:,im_col).*cosd(Untitled(:,va_col)-Untitled(:,ia_col)+180));
Q=sqrt(3)*1*Untitled(:,vm_col).*Untitled(:,im_col).*sind(Untitled(:,va_col)-Untitled(:,ia_col));

indx1=1;%20900;%1;

indx2=indx1+lent-1;
iteration=1;
%% calculate window 1 first time
while indx2<length(P)
    sig=[];
    tme=[];
    wndo1=indx1:indx2;
    
    %tme=(0:length(wndo1)-1)/30;
    
    k=1;
    for i=min(wndo1):max(wndo1)
        if Untitled(i,6)==0
            sig(k)=Untitled(i,vm_col);
            tme(k)=(i-1)/30;
            k=k+1;
        end
    end
    
    N=length(sig);
%     count=0;
%     for i=min(wndo1):max(wndo1)  %%%%% what does this do?
%         if Untitled(i,6)==1
%             
%             count=count+1;
%         end
%     end
%     N=N-count;
    M=round(.5*N);
    tme=tme';
   
    scle=max(abs(sig));
    if scle==0
        scle=1;
    end
    sig=sig/scle;
    Y = zeros(M,N-M+1);
    for ctr = 1:M
        Y(ctr,:) = sig(ctr:N - M +ctr);
        
    end
    Rest=(Y*Y')/(N-M+1);
    [L,S,U] = svd(Rest);
    
    
    sing=diag(S);
    P2(1,:)=(20*log10(sing(1:23)'))/(20*log10(sing(1)));
    k=0;
    for ordr = 1:19
        if length(wndo1)>=60 && length(wndo1)<=90
            k=1;
        else if length(wndo1)>=90 && length(wndo1)<=120
                k=2;
            else if length(wndo1)>=120 && length(wndo1)<=150
                    k=3;
                else if length(wndo1)>=150 && length(wndo1)<=180
                        k=4;
                    else if length(wndo1)>=180 && length(wndo1)<=210
                            k=5;
                        else if length(wndo1)>=210 && length(wndo1)<=240
                                k=6;
                            else if length(wndo1)>=240 && length(wndo1)<=270
                                    k=7;
                                else if length(wndo1)>=270 && length(wndo1)<=301
                                        k=8;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        if ((P2(1,ordr)<pronythreshold(k)))
            break
        end
        
    end
    
    ordr=ordr-1;
    [beta,~,~,~]=smprony1(sig,tme,ordr);
    
    omega_list=imag(beta);
    
    modes_fr =omega_list/(2*pi);
    damping =real(beta);
    if ordr==0
        damping=0;
    end
    if isempty (beta)
        ordr=1;
        damping=0;
    end
    
    if ( (ordr>1)||(abs(max(damping))>thrshold(k))||(max(sig)-min(sig)>swng))    %% max-min is also included to check transient
        % (if and elseif) of 'if' statement are false --- then increment the window size by 1
        % move on to next window to check for steady state value
        if iteration>1 % this is done so that once a steady window is obtained and it is incremented by 1 sec again if dynamic is obtained the loop is broken and that window is taken.
            wndo1=indx1:indx2-wndo_increase;
            indx2=indx2-wndo_increase;
            indx1=indx2+1;
            indx2=indx1+lent-1;
            
            break
        end
        
        indx1=indx1+wndo_increase;%indx2+1;
        indx2=indx2+wndo_increase;%lent;
        
        
    else
        
        indx1=indx1;%+wndo_increase;%indx2+1;
        indx2=indx2+wndo_increase;%    lent;
        iteration=iteration+1;
                if length(wndo1)>=300
                     indx1=indx2-wndo_increase+1;%+wndo_increase;%indx2+1;
                     indx2=indx1+lent-1;
                    break
                end
        % break
    end
    
end
iteration=1;


%% calculate window 2 first time
while indx2<length(P)
    
    wndo2=indx1:indx2;
    sig=[];
    tme=[];
    k=1;
    for i=min(wndo2):max(wndo2)
        if Untitled(i,6)==0
            sig(k)=Untitled(i,vm_col);
            tme(k)=(i-1)/30;    %%% check carefully
            k=k+1;
        end
    end
    
    % tme=(0:length(wndo2)-1)/30;  moved in loop %% check carefully
    N=length(sig);
%     count=0;
%     for i=min(wndo2):max(wndo2)
%         if Untitled(i,6)==1
%             count=count+1;
%         end
%     end
%     N=N-count;
    
    M=round(.5*N);
    tme=tme';
    %     [brob,stats] = robustfit(1:(indx2-indx1+1),untitled(wndo1,vm_col));  %% calculate the variance to check the threshold
    %sig=Untitled(wndo2,vm_col); %moved up to remove bad data index and obtain sig witout bad data
    scle=max(abs(sig));
    if scle==0
        scle=1;
    end
    sig=sig/scle;
    Y = zeros(M,N-M+1);
    for ctr = 1:M
        Y(ctr,:) = sig(ctr:N - M +ctr);
        
    end
    Rest=(Y*Y')/(N-M+1);
    [L,S,U] = svd(Rest);
    
    
    sing=diag(S);
    P2(1,:)=(20*log10(sing(1:23)'))/(20*log10(sing(1)));
    k=0;
    for ordr = 1:19
        if length(wndo2)>=60 && length(wndo2)<=90
            k=1;
        else if length(wndo2)>=90 && length(wndo2)<=120
                k=2;
            else if length(wndo2)>=120 && length(wndo2)<=150
                    k=3;
                else if length(wndo2)>=150 && length(wndo2)<=180
                        k=4;
                    else if length(wndo2)>=180 && length(wndo2)<=210
                            k=5;
                        else if length(wndo2)>=210 && length(wndo2)<=240
                                k=6;
                            else if length(wndo2)>=240 && length(wndo2)<=270
                                    k=7;
                                else if length(wndo2)>=270 && length(wndo2)<=301
                                        k=8;
                                    end
                                end
                                
                            end
                        end
                    end
                end
            end
        end
        if ((P2(1,ordr)<pronythreshold(k)))
            break
        end
    end
    
    ordr=ordr-1;
    [beta,~,~,~]=smprony1(sig,tme,ordr);
    
    omega_list=imag(beta);
    
    modes_fr =omega_list/(2*pi);
    damping =real(beta);
    if ordr==0
        damping=0;
    end
    if isempty (beta)
        ordr=1;
        damping=0;
    end
    if ( (ordr>1)||(abs(max(damping))>thrshold(k))||(max(sig)-min(sig)>swng))   %% max-min is also included to check transient
        % (if and elseif) of 'if' statement are false --- then increment the window size by 1
        % move on to next window to check for steady state value
        if iteration>1 % this is done so that once a steady window is obtained and it is incremented by 1 sec again if dynamic is obtained the loop is broken and that window is taken.
            wndo2=indx1:indx2-wndo_increase;
            indx2=indx2-wndo_increase;
            indx1=indx2+1;
            indx2=indx1+lent-1;
            
            break
        end
        
        indx1=indx1+wndo_increase;%indx2+1;
        indx2=indx2+wndo_increase;%lent;
        
        
    else
        
        indx1=indx1;%+wndo_increase;%indx2+1;
        indx2=indx2+wndo_increase;%    lent;
        iteration=iteration+1;
                if length(wndo2)>=300
                    indx1=indx2-wndo_increase+1;%+wndo_increase;%indx2+1;
                     indx2=indx1+lent-1;
                    break
                end
        % break
    end
    
end


%% calculate window 3 now onwards
iteration=1;
compute_zip=0;

a=zeros(1,length(t));
b=zeros(1,length(t));
c=zeros(1,length(t));

while indx2<length(P)
    wndo3=indx1:indx2;
    
    %tme=(0:length(wndo3)-1)/30;
    tme=[];
    sig=[];
    k=1;
    for i=min(wndo3):max(wndo3)
        if Untitled(i,6)==0
            sig(k)=Untitled(i,vm_col);
            tme(k)=(i-1)/30;
            k=k+1;
        end
    end
    
    N=length(sig);

    M=round(.5*N);
    tme=tme';
   
    scle=max(abs(sig));
    if scle==0
        scle=1;
    end
    sig=sig/scle;

    Y = zeros(M,N-M+1);
    for ctr = 1:M
        Y(ctr,:) = sig(ctr:N - M +ctr);
        
    end
    Rest=(Y*Y')/(N-M+1);
    [L,S,U] = svd(Rest);
    
    
    sing=diag(S);
    P2(1,:)=(20*log10(sing(1:23)'))/(20*log10(sing(1)));
    k=0;
    for ordr = 1:19
        if length(wndo3)>=60 && length(wndo3)<=90
            k=1;
        else if length(wndo3)>=90 && length(wndo3)<=120
                k=2;
            else if length(wndo3)>=120 && length(wndo3)<=150
                    k=3;
                else if length(wndo3)>=150 && length(wndo3)<=180
                        k=4;
                    else if length(wndo3)>=180 && length(wndo3)<=210
                            k=5;
                        else if length(wndo3)>=210 && length(wndo3)<=240
                                k=6;
                            else if length(wndo3)>=240 && length(wndo3)<=270
                                    k=7;
                                else if length(wndo3)>=270 && length(wndo3)<=301
                                      k=8;
                                    end
                                 end
                            end
                        end
                    end
                end
            end
        end
        if ((P2(1,ordr)<pronythreshold(k)))
            break
        end
    end
    
    ordr=ordr-1;
    [beta,~,~,~]=smprony1(sig,tme,ordr);

    omega_list=imag(beta);
    
    modes_fr =omega_list/(2*pi);
    damping =real(beta);
    
    if ordr==0 
        damping=0;
    end
    if isempty (beta)
        ordr=1;
        damping=0;
    end
        
        
    if ( (ordr>1)||(abs(max(damping))>thrshold(k))||(max(sig)-min(sig)>swng))     %% max-min is also included to check transient
        % (if and elseif) of 'if' statement are false --- then increment the window size by 1
        % move on to next window to check for steady state value
        if iteration>1 % this is done so that once a steady window is obtained and it is incremented by 1 sec again if dynamic is obtained the loop is broken and that window is taken.
            
            wndo3=indx1:indx2-wndo_increase;
            indx2=indx2-wndo_increase;
            indx1=indx2+1;
            indx2=indx1+lent-1;
            compute_zip=1;
        end
        if compute_zip~=1
        indx1=indx1+wndo_increase;%indx2+1;
        indx2=indx2+wndo_increase;%lent;
        end
        
    else
        
        indx1=indx1;%+wndo_increase;%indx2+1;
        indx2=indx2+wndo_increase;%    lent;
        iteration=iteration+1;
        
                if length(wndo3)>=300
                    indx1=indx2-wndo_increase+1;
                    indx2=indx1+lent-1;
                    compute_zip=1;
                   
                end
       
    end

    if compute_zip==1
        k=1;
        omit_measurement1(1)=0;
        for i=min(wndo1):max(wndo1)
            if Untitled(i,6)==1
                omit_measurement1(k)=i;%Untitled(i,:);
                k=k+1;
            end
        end
      
        
        k=1;
        omit_measurement2(1)=0;
        for i=min(wndo2):max(wndo2)
            if Untitled(i,6)==1
                omit_measurement2(k)=i;%Untitled(i,:);
                k=k+1;
            end
        end
        
        %     if omit_measurement2(1)==0
        %         omit_measurement2=[];
        %     end
        
        
        k=1;
        omit_measurement3(1)=0;
        for i=min(wndo3):max(wndo3)
            if Untitled(i,6)==1
                omit_measurement3(k)=i;%Untitled(i,:);
                k=k+1;
            end
        end

        
        if omit_measurement1(1)==0
            vltg1=(mean(Untitled(wndo1,vm_col)));%-length((omit_measurement1))*mean(Untitled(omit_measurement1,vm_col)))./(length(Untitled(wndo1,vm_col))-length((omit_measurement1)));%omit_measuremenmt1(:,vm_col));
            pwr1=(mean(P(wndo1)));%-length((omit_measurement1))*mean(P(omit_measurement1)))./(length(Untitled(wndo1,vm_col))-length((omit_measurement1)));
            
        else
            vltg1=(length(Untitled(wndo1,vm_col))*mean(Untitled(wndo1,vm_col))-length((omit_measurement1))*mean(Untitled(omit_measurement1,vm_col)))./(length(Untitled(wndo1,vm_col))-length((omit_measurement1)));%omit_measuremenmt1(:,vm_col));
            pwr1=(length(Untitled(wndo1,vm_col))*mean(P(wndo1))-length((omit_measurement1))*mean(P(omit_measurement1)))./(length(Untitled(wndo1,vm_col))-length((omit_measurement1)));
        end
        
        if omit_measurement2(1)==0
            vltg2=(mean(Untitled(wndo2,vm_col)));
            pwr2=(mean(P(wndo2)));
            
        else
            vltg2=(length(Untitled(wndo2,vm_col))*mean(Untitled(wndo2,vm_col))-length((omit_measurement2))*mean(Untitled(omit_measurement2,vm_col)))./(length(Untitled(wndo2,vm_col))-length((omit_measurement2)));%-Untitled(omit_measurement2,vm_col));
            pwr2=(length(Untitled(wndo2,vm_col))*mean(P(wndo2))-length((omit_measurement2))*mean(P(omit_measurement2)))./(length(Untitled(wndo2,vm_col))-length((omit_measurement2)));
        end
        
        if omit_measurement3(1)==0
            vltg3=(mean(Untitled(wndo3,vm_col)));
            pwr3=(mean(P(wndo3)));
            
        else
            
            vltg3=(length(Untitled(wndo3,vm_col))*mean(Untitled(wndo3,vm_col))-length(omit_measurement3)*mean(Untitled(omit_measurement3,vm_col)))./(length(Untitled(wndo3,vm_col))-length((omit_measurement3)));
            pwr3=(length(Untitled(wndo3,vm_col))*mean(P(wndo3))-length(omit_measurement3)*mean(P(omit_measurement3)))./(length(Untitled(wndo3,vm_col))-length((omit_measurement3)));
        end
        senP1=(pwr2-pwr1)/(vltg2-vltg1);
        %                 %% To eliminate spikes in sensitivities
%         if (abs(vltg2-vltg1)/base_voltage) <= 0.0001
%             senP1 = 0;
%         end
%         
%         if (min(wndo2)-max(wndo1)) >=wind_diff
%             senP1=0;
%         end
        
        senP(min(wndo1):max(wndo2))=senP1;
        senP2=(pwr3-pwr2)/(vltg3-vltg2);
%         if (min(wndo3)-max(wndo2)) >=wind_diff
%             senP2=0;
%         end
%         
%         if (abs(vltg3-vltg2)/base_voltage) <= 0.0001
%             senP2 = 0;
%         end
        senP(min(wndo2):max(wndo3))=senP2;
        
        if senP1 > 0 && senP2 >0
            pwr_3=[(P(wndo1));(P(wndo2));(P(wndo3))];
            vltge_3=[(Untitled(wndo1,vm_col));(Untitled(wndo2,vm_col));(Untitled(wndo3,vm_col))];
            vltge_3 = vltge_3./base_voltage; 
            k1=30;k2=6;k3=9;
            k1=1;k2=1;k3=1;
            A_3=[vltge_3.*vltge_3 vltge_3 ones(length(vltge_3),1)];
            A_3(:,1)=A_3(:,1)*k1;
            A_3(:,2)=A_3(:,2)*k2;
            A_3(:,3)=A_3(:,3)*k3;
            A=[];B=[];Aeq=[];beq=[];
            lb=[0;0;0];
            ub=[max(pwr_3);max(pwr_3);max(pwr_3)];
            Blsc = lsqlin(A_3,pwr_3,A,B,Aeq,beq,lb,ub);
            Blsc(1)=Blsc(1)*k1;
            Blsc(2)=Blsc(2)*k2;
            Blsc(3)=Blsc(3)*k3;
            X=Blsc';
            per(1,1) = X(1,1)/(X(1,1)+X(1,2)+X(1,3));
            per(1,2) = X(1,2)/(X(1,1)+X(1,2)+X(1,3));
            per(1,3) = X(1,3)/(X(1,1)+X(1,2)+X(1,3));
           
            a(min(wndo1):max(wndo3))=per(1,1);
            b(min(wndo1):max(wndo3))=per(1,2);
            c(min(wndo1):max(wndo3))=per(1,3);
            
            (abs((mean(P(wndo2))-sum(Blsc)))/(mean(P(wndo2))))*100  %% diff between calc and meas power
%             t(wndo1(1))
%             t(wndo3(end))
        end
        
        
        
        wndo1 = wndo2;  % Window 2 is changed to window 1.
        wndo2 = wndo3;  % Window 3 is changed to window 2.
        
        
        
        
        iteration=1;
        indx1=indx2+1;
        indx2=indx1+lent-1;
        compute_zip=0;
    end
    
end
zz= length(t)-length(a);

for i=1:length(senP)
    if senP(i)<=0
        senP(i)=0;
    end
end


end