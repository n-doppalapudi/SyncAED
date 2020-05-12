function [baddata_index]=Chebyshev(Untitled)
warning('off','all')
tic
%watchon;
DATA= Untitled;          %importdata(DATAFILE);
dataset_size=size(DATA,1); %total points in dataset

% update status


% time=DATA(:,1);
% voltage=DATA(:,2);
time=(1:dataset_size)';



%initialization
window_size=60; % size of the window
moving_size=60; %  moving step of the window

time_size=size(time,1); %size of time scale
ending_point=time_size-window_size; % ending point of the rolling window
result_Chebyshev=zeros(1,2);%store results
result_Regression=zeros(1,2);
result_index_Chebyshev=1;

Result_union_voltage_mag=zeros(0);

for k=1
    voltage=DATA(:,k);%/Untitled(1,1);%7.2;
    for i=1 : moving_size : ending_point
        %%
        %% Begin Chebyshev
        %Phase 1 of Chebyshev
        % step 1: select the P value 0.1/0.05/0.01  suggestion: 0.01
        p1=0.001;%0.01;
        k1=1/sqrt(p1);
        window_time=time(i:window_size+i-1,1);
        window_voltage=voltage(i:window_size+i-1,1);
        mean1=mean(window_voltage); % mean value within the window
        variance1=std(window_voltage);% variance value within the window
        %step 4 set the upper and lower bound
        ODV_1U=mean1+0.5*k1*variance1;
        ODV_1L=mean1-0.5*k1*variance1;
        %step 5 find out the data that are not in the bound
%         if variance1<1e-4
%             continue
%         end
        upper=find(window_voltage>ODV_1U);
        lower=find(window_voltage<ODV_1L);%get the index of the values that are beyond the boundary
        %Phase 2 of Chebyshev
        % step 1         get the trimmed window (eliminate protential outlier in phase 1)
        trimmed_window_voltage=window_voltage;
        trim=[lower;upper];
        trimmed_window_voltage(trim)=[];
        % step 2          select the p from 0.01/0.001/0.0001  suggestion : 0.0001
        % p2=0.0001;
        p2=0.1;%0.15;%0.5;
        k2=1/sqrt(p2);
        mean2=mean(trimmed_window_voltage); % mean value within the window
        variance2=std(trimmed_window_voltage);% variance value within the window
        ODV_2U=mean2+2*k2*variance2;
        ODV_2L=mean2-2*k2*variance2;
     

        
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

%             
        end 
        
    
      
    end
   % if k==1
        Result_union_voltage_mag=union(result_Chebyshev,result_Regression);
    %end

baddata_index=  Result_union_voltage_mag;%     union(bd1,bd2);  %   
warning('off','all')

k=1;
% for i=1:length(baddata_index)-1
%     if baddata_index(i+1)==(baddata_index(i)+1)
%         fault(k)=round(baddata_index(i)/30);
%     k=k+1;
%     end
% end
% baddata_index=[];
% baddata_index=fault;
toc


end

