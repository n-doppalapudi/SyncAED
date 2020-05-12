function [baddata_index]=Chebyshev_Regression(Untitled)
warning('off','all')
tic
%DATAFILE='6_10_2015_8_16.xlsx';

DATA= Untitled;          %importdata(DATAFILE);
dataset_size=size(DATA,1); %total points in dataset

% time=DATA(:,1);
% voltage=DATA(:,2);
time=(1:dataset_size)';



%initialization
window_size=60; % size of the window
moving_size=30; %  moving step of the window

time_size=size(time,1); %size of time scale
ending_point=time_size-window_size; % ending point of the rolling window
result_Chebyshev=zeros(1,2);%store results
result_Regression=zeros(1,2);
result_index_Chebyshev=1;
result_index_Regression=1;
Result_union_voltage_mag=zeros(0);
Result_union_voltage_angle=zeros(0);
Result_union_current_mag=zeros(0);
Result_union_current_angle=zeros(0);
for k=1:1
    voltage=DATA(:,k);%/7.2;%3.087421562500000e+05;%7.2;
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
        ODV_1U=mean1+2*k1*variance1;
        ODV_1L=mean1-2*k1*variance1;
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
        ODV_2U=mean2+2*k2*variance2;
        ODV_2L=mean2-2*k2*variance2;
        %find the outliers
        % outlier_upper=find(window_voltage>ODV_2U);
        % outlier_lower=find(window_voltage<ODV_2L);
        
        %% Begin Linear Regression
%         p = polyfit(window_time,window_voltage,1);   % p returns 2 coefficients fitting r = a_1 * x + a_2
%         r = p(1) .* window_time + p(2); %regression line
%         mean_variance=mean(abs(r-window_voltage));
%         % variance=var(window_voltage);% set the standard deviation
%         variance=5*mean_variance;
%         % variance=0.001;
%         highThres= p(1) .* window_time + p(2)+variance;
%         lowThres= p(1) .* window_time + p(2)-variance;
%         difference=highThres-lowThres;
%         max_value=max(difference);
        
        
        %%
        for j=1: window_size %detecting bad data
            if  variance1>1e-4
                if window_voltage(j)>ODV_2U || window_voltage(j)<ODV_2L
                    empty_determine=find(result_Chebyshev==window_time(j),1);
                    if isempty(empty_determine)
                        %fprintf('found an outlier data on time %d \n ', window_time(j) ) ;
                        result_Chebyshev(result_index_Chebyshev)=window_time(j);
                        result_index_Chebyshev=result_index_Chebyshev+1;
                    end
                end
            end
            %% end of detection bad data for Chebyshev
%             if  max_value>1e-3
%                 if window_voltage(j)<lowThres(j) || window_voltage(j)>highThres(j)
%                     empty_determine=find(result_Regression==window_time(j),1);
%                     if isempty(empty_determine)
%                         % fprintf('found an outlier data on time %d \n ', window_time(j) ) ;
%                         result_Regression(result_index_Regression)=window_time(j);
%                         result_index_Regression=result_index_Regression+1;
%                     end
%                 end
%             end
%             
        end %end of detection bad data for Linear Regression
        
        %% end of Chebyshev
        
        %%
        %%
    end
    if k==1
        Result_union_voltage_mag=union(result_Chebyshev,result_Regression);
    end
%     if k==2
%         Result_union_voltage_angle=union(result_Chebyshev,result_Regression);
%     end
%     if k==3
%         Result_union_current_mag=union(result_Chebyshev,result_Regression);
%     end
%     if k==4
%         Result_union_current_angle=union(result_Chebyshev,result_Regression);
%     end
%     
% end
% bd1=union(Result_union_voltage_mag,Result_union_voltage_angle);
% bd2=union(Result_union_current_mag,Result_union_current_angle);
baddata_index=  Result_union_voltage_mag;%     union(bd1,bd2);  %   
warning('off','all')
toc
end

