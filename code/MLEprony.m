clc;
tic
format long;
vm_col=1; va_col=2; im_col=3; ia_col=4;
%prompt = 'Load the input file ';
% Untitled = input(prompt);

Untitled=xlsread('mlepronytest.xlsx');%PMUparam2;
len=length(Untitled(:,1))
time=0.03333333333333333333333*(1:len)-0.03333333333333333333333;


% 
% HRDATA=[hrdata(:,1:4);hrdata(:,5:8);hrdata(:,9:12)];
 %Untitled=hrdata;
 
 
 untitled= Untitled;
 Untitled=[];

 %% insert bad data
 [Untitled,inserted_baddata]=insertbaddata(untitled(:,2));
 Untitled=[time' Untitled];
%  [Untitled,inserted_baddata]=insertbaddata(Untitled);
%  inserted_baddata=inserted_baddata';
 
%Untitled=SEL2raw;
%Untitled=realpmufabS2;
%Untitled=Untitled1;

lent=60;
swng=3;%210;
%Untitled=[Untitled;Untitled;Untitled;];%Untitled;Untitled;Untitled;];

baddata_index=1;
senP=zeros(length(Untitled),1);

%% bad data detection algo
%[baddata_index,result_DBSCAN]=Chebyshev_Regression_DBSCAN(Untitled);
[baddata_index]=MLE(Untitled,0);
baddata_index=baddata_index';
if baddata_index(1)==0
    baddata_index(1)=[];
end
%% precision and recall
%[precision,recall]=precision_recall(baddata_index,inserted_baddata);

%% flagging bad data index
%result_DBSCAN=result_DBSCAN';
% if baddata_index(1)==0
%     baddata_index(1)=[];
% end
%baddata_index=[];
baddata_index=round((baddata_index*30)+1);
Untitled(:,3:4)=0;
for i=1:length(baddata_index)
    if baddata_index(1)~=0
Untitled(baddata_index(i),3)=1;
    end
end
t=(0:length(Untitled(:,1))-1)/30; %60;%
%untitled= Untitled;

bad_data_time=baddata_index/30;
[precision1,recall1]=precision_recall(baddata_index,inserted_baddata)


%Untitled(baddata_index(1:length(baddata_index)),:)=[];
thrshold = 5e-3;


indx1=1;

indx2=lent;%+indx1;

%% calculate window 1 first time
while indx2<length(Untitled)
    sig=[];
    tme=[];
    wndo1=indx1:indx2;

    %tme=(0:length(wndo1)-1)/30;
    
    k=1;
    for i=min(wndo1):max(wndo1)   
        if Untitled(i,3)==0
    sig(k)=Untitled(i,2);
    tme(k)=(i-1)/30;
    k=k+1;
        end
    end
    
    N=length(wndo1);
    count=0;
    for i=min(wndo1):max(wndo1)  %%%%% what does this do?
        if Untitled(i,3)==1
         
            count=count+1;
        end
    end
    N=N-count;
    M=round(.5*N);
    tme=tme';
%     [brob,stats] = robustfit(1:(indx2-indx1+1),untitled(wndo1,vm_col));  %% calculate the variance to check the threshold
    % sig=Untitled(wndo1,vm_col); moved up
    scle=max(abs(sig));
sig=sig/scle;
Y = zeros(M,N-M+1);
for ctr = 1:M
   Y(ctr,:) = sig(ctr:N - M +ctr);

end
Rest=(Y*Y')/(N-M+1);
[L,S,U] = svd(Rest);


sing=diag(S);
P2(1,:)=(20*log10(sing(1:23)'))/(20*log10(sing(1)));

for ordr = 1:19
    if ((P2(1,ordr)<-4.2700))
         break
    end
end

ordr=ordr-1;
[beta,~,~,~]=smprony1(sig,tme,ordr);

omega_list=imag(beta);

modes_fr =omega_list/(2*pi);
damping =real(beta);

if (((length(ordr)>1) || (abs(max(damping))>0.0001)))%||(max(Untitled(wndo1,vm_col))-min(Untitled(wndo1,vm_col))>swng))    %% max-min is also included to check transient
        % (if and elseif) of 'if' statement are false --- then increment the window size by 1
        % move on to next window to check for steady state value
      
        Untitled(indx1:indx2,4)=1; %%%% it is a transient window.
        indx1=indx2+1;
        indx2=indx2+lent;
else 
        
        
        indx1=indx2+1;
        indx2=indx2+lent;
       
%         break
end

end

for i=1:length(baddata_index)
    if Untitled(baddata_index(i),4)==1
       if(Untitled(baddata_index(i),2))==0
           Untitled(baddata_index(i),3)=1;
       else
           Untitled(baddata_index(i),3)=0;
    
       end
    end
end

baddata_index=[];
k=1;
for i=1:length(Untitled)
    if Untitled(i,3)==1
baddata_index(k)=i;
k=k+1;
    end
end

[precision2,recall2]=precision_recall(baddata_index,inserted_baddata)

% %% calculate window 2 first time
% while indx2<length(P)
%     
%     wndo2=indx1:indx2;
% sig=[];
% tme=[];
%       k=1;
%     for i=min(wndo2):max(wndo2)   
%         if Untitled(i,6)==0
%     sig(k)=Untitled(i,vm_col);
%     tme(k)=(i-1)/30;    %%% check carefully
%     k=k+1;
%         end
%     end
%     
%    % tme=(0:length(wndo2)-1)/30;  moved in loop %% check carefully
%       N=length(wndo2);
%     count=0;
%     for i=min(wndo2):max(wndo2)
%         if Untitled(i,6)==1
%             count=count+1;
%         end
%     end
%     N=N-count;
%   
%     M=round(.5*N);
%     tme=tme';
% %     [brob,stats] = robustfit(1:(indx2-indx1+1),untitled(wndo1,vm_col));  %% calculate the variance to check the threshold
%     %sig=Untitled(wndo2,vm_col); %moved up to remove bad data index and obtain sig witout bad data
%     scle=max(abs(sig));
% sig=sig/scle;
% Y = zeros(M,N-M+1);
% for ctr = 1:M
%    Y(ctr,:) = sig(ctr:N - M +ctr);
% 
% end
% Rest=(Y*Y')/(N-M+1);
% [L,S,U] = svd(Rest);
% 
% 
% sing=diag(S);
% P2(1,:)=(20*log10(sing(1:23)'))/(20*log10(sing(1)));
% 
% for ordr = 1:19
%     if ((P2(1,ordr)<-2.6))
%          break
%     end
% end
% 
% ordr=ordr-1;
% [beta,~,~,~]=smprony1(sig,tme,ordr);
% 
% omega_list=imag(beta);
% 
% modes_fr =omega_list/(2*pi);
% damping =real(beta);
% 
% if (((length(ordr)>1)|| (abs(max(damping))>thrshold))||(max(Untitled(wndo2,vm_col))-min(Untitled(wndo2,vm_col))>swng))    %% max-min is also included to check transient
%         % (if and elseif) of 'if' statement are false --- then increment the window size by 1
%         % move on to next window to check for steady state value
%         indx1=indx2+1;
%         indx2=indx2+lent;   
% else 
%         indx1=indx2+1;
%         indx2=indx2+lent;  
%         break
% end
% 
% end
% 
% 
% %% calculate window 3 now onwards
% a=zeros(1,length(t));
% b=zeros(1,length(t));
% c=zeros(1,length(t));
% 
% while indx2<length(P)
%     wndo3=indx1:indx2;
% 
%     %tme=(0:length(wndo3)-1)/30;
%  tme=[];
%     sig=[];
%       k=1;
%     for i=min(wndo3):max(wndo3)   
%         if Untitled(i,6)==0
%     sig(k)=Untitled(i,vm_col);
%     tme(k)=(i-1)/30;
%     k=k+1;
%         end
%     end
%     
%       N=length(wndo3);
%     count=0;
%     for i=min(wndo3):max(wndo3)
%         if Untitled(i,6)==1
%             count=count+1;
%         end
%     end
%     N=N-count;
%     M=round(.5*N);
%     tme=tme';
% %     [brob,stats] = robustfit(1:(indx2-indx1+1),untitled(wndo1,vm_col));  %% calculate the variance to check the threshold
%    % sig=Untitled(wndo3,vm_col); moved up
%     scle=max(abs(sig));
% sig=sig/scle;
% Y = zeros(M,N-M+1);
% for ctr = 1:M
%    Y(ctr,:) = sig(ctr:N - M +ctr);
% 
% end
% Rest=(Y*Y')/(N-M+1);
% [L,S,U] = svd(Rest);
% 
% 
% sing=diag(S);
% P2(1,:)=(20*log10(sing(1:23)'))/(20*log10(sing(1)));
% 
% for ordr = 1:19
%     if ((P2(1,ordr)<-2.6))
%          break
%     end
% end
% 
% ordr=ordr-1;
% [beta,~,~,~]=smprony1(sig,tme,ordr);
% 
% omega_list=imag(beta);
% 
% modes_fr =omega_list/(2*pi);
% damping =real(beta);
% 
% if (((length(ordr)>1)|| (abs(max(damping))>thrshold))||(max(Untitled(wndo3,vm_col))-min(Untitled(wndo3,vm_col))>swng))
%     
%     indx1=indx2+1;
%                 indx2=indx2+lent;
%                 
% else
%     k=1;
%     omit_measurement1(1)=0;
%     for i=min(wndo1):max(wndo1)
%         if Untitled(i,6)==1
%             omit_measurement1(k)=i;%Untitled(i,:);
%             k=k+1;
%         end
%     end
%     
% %      if omit_measurement1(1)==0
% %         omit_measurement1=[];
% %      end
%      
%      k=1;
%     omit_measurement2(1)=0;
%     for i=min(wndo2):max(wndo2)
%         if Untitled(i,6)==1
%             omit_measurement2(k)=i;%Untitled(i,:);
%             k=k+1;
%         end
%     end
%     
% %     if omit_measurement2(1)==0
% %         omit_measurement2=[];
% %     end
%     
%     
%     k=1;
%      omit_measurement3(1)=0;
%     for i=min(wndo3):max(wndo3)
%         if Untitled(i,6)==1
%             omit_measurement3(k)=i;%Untitled(i,:);
%             k=k+1;
%         end
%     end
%     
% %      if omit_measurement3(1)==0
% %         omit_measurement3=[];
% %      end
% 
% % vltg1=(mean(Untitled(wndo1,vm_col)));
% % pwr1=(mean(P(wndo1)));
% % vltg2=(mean(Untitled(wndo2,vm_col)));
% % pwr2=(mean(P(wndo2)));
% % vltg3=(mean(Untitled(wndo3,vm_col)));
% % pwr3=(mean(P(wndo3)));
%        
%     if omit_measurement1(1)==0
%      vltg1=(mean(Untitled(wndo1,vm_col)));%-length((omit_measurement1))*mean(Untitled(omit_measurement1,vm_col)))./(length(Untitled(wndo1,vm_col))-length((omit_measurement1)));%omit_measuremenmt1(:,vm_col));
%         pwr1=(mean(P(wndo1)));%-length((omit_measurement1))*mean(P(omit_measurement1)))./(length(Untitled(wndo1,vm_col))-length((omit_measurement1)));
%            
%     else
%         vltg1=(length(Untitled(wndo1,vm_col))*mean(Untitled(wndo1,vm_col))-length((omit_measurement1))*mean(Untitled(omit_measurement1,vm_col)))./(length(Untitled(wndo1,vm_col))-length((omit_measurement1)));%omit_measuremenmt1(:,vm_col));
%         pwr1=(length(Untitled(wndo1,vm_col))*mean(P(wndo1))-length((omit_measurement1))*mean(P(omit_measurement1)))./(length(Untitled(wndo1,vm_col))-length((omit_measurement1)));
%     end
%     
%     if omit_measurement2(1)==0
%        vltg2=(mean(Untitled(wndo2,vm_col)));
%        pwr2=(mean(P(wndo2)));  
%         
%     else    
%         vltg2=(length(Untitled(wndo2,vm_col))*mean(Untitled(wndo2,vm_col))-length((omit_measurement2))*mean(Untitled(omit_measurement2,vm_col)))./(length(Untitled(wndo2,vm_col))-length((omit_measurement2)));%-Untitled(omit_measurement2,vm_col));
%         pwr2=(length(Untitled(wndo2,vm_col))*mean(P(wndo2))-length((omit_measurement2))*mean(P(omit_measurement2)))./(length(Untitled(wndo2,vm_col))-length((omit_measurement2)));
%     end    
%     
%      if omit_measurement3(1)==0
%        vltg3=(mean(Untitled(wndo3,vm_col)));
%        pwr3=(mean(P(wndo3)));  
%         
%     else 
%     
%         vltg3=(length(Untitled(wndo3,vm_col))*mean(Untitled(wndo3,vm_col))-length(omit_measurement3)*mean(Untitled(omit_measurement3,vm_col)))./(length(Untitled(wndo3,vm_col))-length((omit_measurement3)));     
%         pwr3=(length(Untitled(wndo3,vm_col))*mean(P(wndo3))-length(omit_measurement3)*mean(P(omit_measurement3)))./(length(Untitled(wndo3,vm_col))-length((omit_measurement3)));
%      end   
%         senP1=(pwr2-pwr1)/(vltg2-vltg1);
% %                 %% To eliminate spikes in sensitivities
%                 if (abs(vltg2-vltg1)/base_voltage) <= 0.0001
%                     senP1 = 0;
%                 end
%                 
%               if (min(wndo2)-max(wndo1)) >=wind_diff
%                   senP1=0;
%               end
%                 
%         senP(min(wndo1):max(wndo2))=senP1;
%         senP2=(pwr3-pwr2)/(vltg3-vltg2);
%         if (min(wndo3)-max(wndo2)) >=wind_diff
%             senP2=0;
%         end
%         
%                 if (abs(vltg3-vltg2)/base_voltage) <= 0.0001
%                     senP2 = 0;
%                 end
%  senP(min(wndo2):max(wndo3))=senP2;
%  
%         if senP1 > 0 && senP2 >0
%             pwr_3=[mean(P(wndo1));mean(P(wndo2));mean(P(wndo3))];
%             vltge_3=[mean(Untitled(wndo1,vm_col));mean(Untitled(wndo2,vm_col));mean(Untitled(wndo3,vm_col))];
%             vltge_3 = vltge_3./base_voltage; % converting voltage into per unit
%             A_3=[vltge_3.*vltge_3 vltge_3 ones(length(vltge_3),1)];
%             
%             A=[];
%             B=[];
%             Aeq=[];
%             beq=[];
%             lb=[0;0;0];
%             ub=[max(pwr_3);max(pwr_3);max(pwr_3)];
%             Blsc = lsqlin(A_3,pwr_3,A,B,Aeq,beq,lb,ub);
% %                         Blsc = lsqlin(A_3,pwr_3);
%             X=Blsc';
%             per(1,1) = X(1,1)/(X(1,1)+X(1,2)+X(1,3));
%             per(1,2) = X(1,2)/(X(1,1)+X(1,2)+X(1,3));
%             per(1,3) = X(1,3)/(X(1,1)+X(1,2)+X(1,3));
%             per
%             a(min(wndo1):max(wndo3))=per(1,1);
%             b(min(wndo1):max(wndo3))=per(1,2);
%             c(min(wndo1):max(wndo3))=per(1,3);
%             
%             (abs((mean(P(wndo2))-sum(Blsc)))/(mean(P(wndo2))))*100  %% diff between calc and meas power
%         t(wndo1(1))
%         t(wndo3(end))
%         end
%         
%         
%         
%         wndo1 = wndo2;  % Window 2 is changed to window 1.
%         wndo2 = wndo3;  % Window 3 is changed to window 2.
%  
%  
%                 
% 
%                 indx1=indx2+1;
%                 indx2=indx2+lent;
%                 
% end
% 
% end
% zz= length(t)-length(a);
% % a(max(wndo3):max(wndo3)+zz)=0;
% % b(max(wndo3):max(wndo3)+zz)=0;
% % c(max(wndo3):max(wndo3)+zz)=0;
% 
% for i=1:length(senP)
% if senP(i)<=0
%     senP(i)=0;
% end
% end
% 
% % k=1;
% % for i=1:length(baddata_index)
% %     if a(baddata_index(i))~=0 && b(baddata_index(i))~=0
% %         error_zip(k)=baddata_index(i);
% %    k=k+1;
% %     end
% % end
% % 
% % for i=1:length(error_zip)
% %     find a(error_zip(i)+1:error_zip(i)+1000)==0
% %         set_a(
% % 
% 
% subplot(4,1,1)
% plot(t,a,t,b,t,c)
% 
% ylabel('ZIP')
% 
% subplot(4,1,2);
% plot(t,Untitled(:,1))
% 
% 
% ylabel('voltage')
% 
% subplot(4,1,3);
% plot(t,P')
% 
% 
% ylabel('Power P')
% 
% subplot(4,1,4);
% plot(t,senP)
% 
% 
% ylabel('sensitivity')
% 
% % %ylabel('Rotor Speed')
% xlabel('Time (s)')
% % grid minor
% Untitled= untitled;
% 
% % plot(t,a,t,b,t,c)
% % xlabel('Time (s)')
% % ylabel ('ZIP parameters')
%  toc