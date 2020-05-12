

DATAFILE='luan.xlsx'; % data set with manually inserted bad data points 
%TESTFILE='Book1.xlsx';
DATA=importdata(DATAFILE);
%TEST_DATA=importdata(TESTFILE);
FIG_FONT_SIZE=12;
% dataset_size=size(DATA,1); %total points in dataset
% test_time=TEST_DATA.data(:,1);
% test_voltage=TEST_DATA.data(:,2);
DATA=DATA.data;
time=DATA(:,1);
voltage=DATA(:,2);
fault_table=DATA(:,3);
fault_table(6541:end)=[];
index_train=find(fault_table<=time(end));
Gtruth=fault_table(index_train);

index_Test_fault=find(fault_table>time(end));
Test_fault=fault_table(index_Test_fault);

Gtruth(~any(Gtruth,2),:)=[];
Gtruth(:,~any(Gtruth,1))=[];
% time=(1:dataset_size)';
% prediction matrix
predict(:,1)=time;
predict(:,2)=1;
predict(:,3)=1;
predict(:,4)=1;

%initialization
window_size=72; % size of the window
moving_size=36; %  moving step of the window

time_size=size(time,1); %size of time scale
ending_point=time_size-window_size; % ending point of the rolling window
result_Chebyshev=zeros(1,2);%store results
result_Regression=zeros(1,2);
result_DBSCAN=zeros(1,2);

result_index_Chebyshev=1;
result_index_Regression=1;
result_index_DBSCAN=1;

Result_union_voltage_mag=zeros(0);
Result_union_voltage_angle=zeros(0);
Result_union_current_mag=zeros(0);
Result_union_current_angle=zeros(0);

%% calculate the distance of outlier/missing data point to thresholds
Distant_Chebyshev=zeros(1,2);
Distant_index_Chebyshev=1;

Distant_Regression=zeros(1,2);
Distant_index_Regression=1;
%%

%voltage=DATA(:,2);
for i=1 : moving_size : ending_point
%%
%% Begin Chebyshev 
%Phase 1 of Chebyshev
% step 1: select the P value 0.1/0.05/0.01  suggestion: 0.01 
p1=0.1;
k1=1/sqrt(p1);
window_time=time(i:window_size+i-1,1);
window_voltage=voltage(i:window_size+i-1,1);
mean1=mean(window_voltage); % mean value within the window
variance1=std(window_voltage);% variance value within the window
%step 4 set the upper and lower bound
ODV_1U=mean1+k1*variance1;
ODV_1L=mean1-k1*variance1;
%step 5 find out the data that are not in the bound
if variance1<1e-4
    continue
end
upper=find(window_voltage>ODV_1U);
lower=find(window_voltage<ODV_1L);%get the index of the values that are beyond the boundary
%Phase 2 of Chebyshev
% step 1         get the trimmed window (eliminate protential outlier in phase 1)
trimmed_window_voltage=window_voltage;
trim=[lower;upper];
trimmed_window_voltage(trim)=[];
% step 2          select the p from 0.01/0.001/0.0001  suggestion : 0.0001
% p2=0.0001;
p2=0.05;
k2=1/sqrt(p2);
mean2=mean(trimmed_window_voltage); % mean value within the window
variance2=std(trimmed_window_voltage);% variance value within the window
ODV_2U=mean2+k2*variance2;
ODV_2L=mean2-k2*variance2;
Cheby_difference=ODV_2U-ODV_2L;
%find the outliers
% outlier_upper=find(window_voltage>ODV_2U);
% outlier_lower=find(window_voltage<ODV_2L);

%%Begin Linear Regression
p = polyfit(window_time,window_voltage,1);   % p returns 2 coefficients fitting r = a_1 * x + a_2
r = p(1) .* window_time + p(2); %regression line
mean_variance=mean(abs(r-window_voltage));
% variance=var(window_voltage);% set the standard deviation 
variance=5*mean_variance;
% variance=0.001;
highThres= p(1) .* window_time + p(2)+variance;
lowThres= p(1) .* window_time + p(2)-variance;
difference=highThres-lowThres;
max_value=max(difference);
%%
for j=1: window_size %detecting bad data
    if Cheby_difference>0.005
    if window_voltage(j)>ODV_2U || window_voltage(j)<ODV_2L
        empty_determine=find(result_Chebyshev==window_time(j));
        index_predict=find(predict(:,1)==window_time(j));
        predict(index_predict,2)=-1;
        if isempty(empty_determine)
%        fprintf('find a ourlier data on time %d \n ', window_time(j) ) 
       result_Chebyshev(result_index_Chebyshev)=window_time(j);
       
 %%
 %get the distance for score use
distance_to_ODV_2U=abs(window_voltage(j)-ODV_2U);
distance_to_ODV_2L=abs(window_voltage(j)-ODV_2L);
Distant_Chebyshev(result_index_Chebyshev)=min(distance_to_ODV_2U,distance_to_ODV_2L) ;
%%
       result_index_Chebyshev=result_index_Chebyshev+1;
       

        end
    end
    end
 %end of detection bad data for Chebyshev
    if  max_value>0.005
if window_voltage(j)<lowThres(j) || window_voltage(j)>highThres(j)
        empty_determine=find(result_Regression==window_time(j));
        index_predict=find(predict(:,1)==window_time(j));
        predict(index_predict,3)=-1;
        if isempty(empty_determine)
%        fprintf('find a outlier data on time %d \n ', window_time(j) ) 
       result_Regression(result_index_Regression)=window_time(j);
              %% get the distance and score
   Q2_high=[window_time(2) highThres(2)];
   Q1_high=[window_time(1) highThres(1)];
   Q2_low=[window_time(2) lowThres(2)];
   Q1_low=[window_time(1) lowThres(1)]; 
   P=[window_time(j) window_voltage(j)];
distance_to_highThres=abs(det([Q2_high-Q1_high;P-Q1_high]))/norm(Q2_high-Q1_high);
distance_to_lowThres=abs(det([Q2_low-Q1_low;P-Q1_low]))/norm(Q2_low-Q1_low);
Distant_Regression(result_index_Regression)=min(distance_to_highThres,distance_to_lowThres) ;       
       %%
       result_index_Regression=result_index_Regression+1;
       

        end
end       
        end
%     
% plot Chebyshev
% if window_time(j)==92.916
% ODV_2LL=zeros(1,72);
% ODV_2UU=zeros(1,72);
% for j=1:72
%  ODV_2LL(j)=ODV_2L; 
%  ODV_2UU(j)=ODV_2U; 
% end
% plot(window_time, ODV_2LL,'-','color', 'b');
% hold on;
% plot(window_time, window_voltage,'-','color', 'g');
% plot(window_time, ODV_2UU,'-','color', 'b');
% % end
% % legend('lowThres_Chebyshev','Voltage','highThres_Chebyshev')
% 
% % % add regression plot in Chebyshev
% plot(window_time, lowThres,'-','color', 'r');
% % plot(window_time, window_voltage,'-','color', 'r');
% plot(window_time, r,'-','color', 'r');
% plot(window_time, highThres,'-','color', 'r');
% legend('lowThres Chebyshev','Voltage','highThres Chebyshev','lowThres Regression','regression line','highThres Regression')

% end 
% % legend('lowThres','Voltage','highThres')
% 

%
%% plot Linear Regression
% if window_time(j)==4.633
% plot(window_time, lowThres,'-','color', 'r');
% hold on;
% plot(window_time, window_voltage,'-','color', 'b');
% plot(window_time, r,'-','color', 'k');
% plot(window_time, highThres,'-','color', 'm');
% legend('lowThres','Voltage', 'regression line','highThres')
% xlabel('Time (s)','Fontsize', FIG_FONT_SIZE,'Fontweight','Bold')
% ylabel('Voltage','Fontsize', FIG_FONT_SIZE,'Fontweight','Bold')
% hold off;
% end
%%
end %end of detection bad data for Linear Regression

%% end of Chebyshev

%%Begin DBSCAN
X(:,1)=window_time;
X(:,2)=window_voltage;
epsilon=0.05;
MinPts=3;
IDX=DBSCAN(X,epsilon,MinPts);



window_index_DBSCAN=find(IDX==0);

window_time_DBSCAN=window_time(window_index_DBSCAN);

window_number_outliers=numel(window_time_DBSCAN);

if window_number_outliers~=0
for k=1:window_number_outliers
 index_predict=find(predict(:,1)==window_time_DBSCAN(k));
        predict(index_predict,4)=-1;
 result_DBSCAN(result_index_DBSCAN)=window_time_DBSCAN(k);
 result_index_DBSCAN=result_index_DBSCAN+1;

% plot
%  if window_time_DBSCAN(k)==4.633
%  PlotClusterinResult(X, IDX);
% title(['DBSCAN Clustering (\epsilon = ' num2str(epsilon) ', MinPts = ' num2str(MinPts) ')']);    
% end
% %
 
end
end

%


%%end DBSCAN
end

result_DBSCAN=unique(result_DBSCAN);
%toc

% if k==1
% Result_union_voltage_mag=union(result_Chebyshev,result_Regression);
% end
% if k==2
% Result_union_voltage_angle=union(result_Chebyshev,result_Regression);
% end
% if k==3
% Result_union_current_mag=union(result_Chebyshev,result_Regression);
% end
% if k==4
% Result_union_current_angle=union(result_Chebyshev,result_Regression);
% end

%%get the results
Total=numel(Gtruth);

inter_Chebyshev=intersect(Gtruth,result_Chebyshev);
inter_Chebyshev=numel(inter_Chebyshev);
num_Chebyshev=numel(result_Chebyshev);

inter_DBSCAN=intersect(Gtruth,result_DBSCAN);
inter_DBSCAN=numel(inter_DBSCAN);
num_DBSCAN=numel(result_DBSCAN);

inter_Regression=intersect(Gtruth,result_Regression);
inter_Regression=numel(inter_Regression);
num_Regression=numel(result_Regression);
%% recall
Recall_Chebyshev=inter_Chebyshev/Total;
Recall_DBSCAN=inter_DBSCAN/Total;
Recall_Regression=inter_Regression/Total;
%%precision
Precision_Chebyshev=inter_Chebyshev/num_Chebyshev;
Precision_DBSCAN=inter_DBSCAN/num_DBSCAN;
Precision_Regression=inter_Regression/num_Regression;
%%false positive
FP_Chebyshev=1-Precision_Chebyshev;
FP_DBSCAN=1-Precision_DBSCAN;
FP_Regression=1-Precision_Regression;

%find the different col
false_positive_Chebyshev=setdiff(result_Chebyshev,Gtruth);
false_positive_Regression=setdiff(result_Regression,Gtruth);
false_positive_DBSCAN=setdiff(result_DBSCAN,Gtruth);
%fail to detection
failed_detection_Chebyshev=setdiff(Gtruth,result_Chebyshev);
failed_detection_Regression=setdiff(Gtruth,result_Regression);
failed_detection_DBSCAN=setdiff(Gtruth,result_DBSCAN);

%vote based
%3
% a=intersect(result_Chebyshev,result_DBSCAN);
% b=intersect(result_Chebyshev,result_Regression);
% c=intersect(result_DBSCAN,result_Regression);
% union1=union(a,b);
% uinion2=union(union1,c);
% 
% inter_union=intersect(Gtruth,uinion2);
% inter_uinion2=numel(inter_union);
% num_union=numel(inter_uinion2);
% Recall_union=inter_uinion2/Total;
% Precision_union=inter_uinion2/num_union;

%SML and MLE
pred1=predict(:,2:4);
[HL,V1,MLE,MAP,VO]=SML_WGS(pred1);

estimat_MLE(:,1)=predict(:,1);
estimat_MLE(:,2)=MLE;

temp_index=find(estimat_MLE(:,2)==-1);
fault_MLE=estimat_MLE(temp_index,1);
inter_MLE=intersect(Gtruth, fault_MLE);
total_MLE=numel(fault_MLE);
result_MLE=numel(inter_MLE);
Recall_MLE=result_MLE/Total;
Precision_MLE=result_MLE/total_MLE;
false_positive_MLE=1-Precision_MLE;

%EM
% test_data=test_voltage;
% test_data(:,2)=test_voltage.^2;
% test_data(:,3)=test_voltage.^3;
% train_data=voltage;
% train_data(:,2)=voltage.^2;
% train_data(:,3)=voltage.^3;
% index_mle=find(MLE==-1);
% training_label=MLE;
% training_label(index_mle)=2;
% [z1,model,llh] = mixGaussEm(train_data',training_label');
% z2 = mixGaussPred(test_data',model);
% z2=z2';
% predict_em(:,1)=test_time;
% predict_em(:,2)=z2;
% tem_index=find(predict_em(:,2)==2);
% fault_em=predict_em(tem_index,1);
% inter_em=intersect(Test_fault,fault_em);
% restult_em=numel(inter_em);
% test_total=numel(Test_fault);
% recall_em_test=restult_em/test_total