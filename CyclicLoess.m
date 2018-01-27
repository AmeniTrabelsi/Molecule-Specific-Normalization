function [data1] = CyclicLoess(data)
    if all(data(:,1)==0)
        data(:,1)=[];
    end
    Area = data;
    TablePerSam_Cur = zeros(size(Area));
    TablePerSam_Ref = zeros(size(Area));
       h = waitbar(1,'normalization now. Please wait...'); 
        for i=1:size(Area,2)-1
            for j=i+1:size(Area,2)
                CurrentSamples=[Area(:,i) Area(:,j)];
                nonan_id = find(~(Area(:,i)==0)&~(Area(:,j)==0));% find(~isnan(Area(:,i))&~isnan(Area(:,j)));
                CurrentSamples_used = CurrentSamples(nonan_id,:);
                M = log2(CurrentSamples_used(:,1)./CurrentSamples_used(:,2));
                A = 1/2*log2(CurrentSamples_used(:,1).*CurrentSamples_used(:,2));
                if length(nonan_id)>4%%The default span of FUNCTION SMOOTH for the moving average is 5
                    MofLoess = smooth(A, M, 0.3, 'loess');
                    MaF = M-MofLoess; %MafterFit
                    TablePerSam_Cur(nonan_id,j) = 2.^(A+MaF/2);
                    TablePerSam_Ref(nonan_id,j) = 2.^(A-MaF/2);
                else
                    TablePerSam_Cur(nonan_id,j) = Area(nonan_id,i);
                    TablePerSam_Ref(nonan_id,j) = Area(nonan_id,j);                    
                end
            end
            Table(i).Sam_cur=TablePerSam_Cur; TablePerSam_Cur = zeros(size(Area));
            Table(i).Sam_Ref=TablePerSam_Ref; TablePerSam_Ref = zeros(size(Area));
            waitbar(i/(size(Area,2)-1));
        end        
        %%%%%%%%%%%%%%%%%%%        
        for i=1:size(Table(1).Sam_cur,2)
            if i>1
                for j=1:i-1
                    Table(i).Sam_cur(:,j)=Table(j).Sam_Ref(:,i);
                end
            end
            Sam_cur = Table(i).Sam_cur;
            if i<size(Table(1).Sam_cur,2)
                Sam_cur(:,i) = [];
            end
            for j = 1:size(Sam_cur,1)
                if length(find(~Sam_cur(j,:)==0))==0
                    TableRearrange(j,i)=0;
                else
                    TableRearrange(j,i)=sum(Sam_cur(j,:))/length(find(~Sam_cur(j,:)==0));
                end
            end 
        end
        close(h);
        data1 = [];
    data1=TableRearrange;