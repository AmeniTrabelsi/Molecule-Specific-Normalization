function [data1] = Quantile(data)
   h = waitbar(1,'normalization now. Please wait...'); 
    if all(data(:,1)==0)
        data(:,1)=[];
    end
    ResultData=quantilenorm(data);
    data1 = [];
    data1=ResultData;
    for i=1:10
    waitbar(i/10);
    end
    close(h);