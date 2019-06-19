function [ temp ] = SAH(temp, vitalslab_hold)

% Matthieu Komorowski - Imperial College London 2017 
% will copy a value in the rows below if the missing values are within the
% hold period for this variable (e.g. 48h for weight, 2h for HR...)
% vitalslab_hold = 2x55 cell (with row1 = strings of names ; row 2 = hold time)


h = waitbar(0,'Initializing waitbar...');

hold=table2array(cell2table(vitalslab_hold(2,:)));
nrow=size(temp,1);
ncol=size(temp,2);

lastcharttime=zeros(1,ncol);
lastvalue=zeros(1,ncol);
oldstayid=temp(1,2);

for i=4:ncol
     waitbar(i/ncol,h,i/ncol*100) 

    for j=1:nrow
        
 
        if oldstayid~=temp(j,2)
            lastcharttime=zeros(1,ncol);
            lastvalue=zeros(1,ncol);
            oldstayid=temp(j,2);
        end
                
        if isnan(temp(j,i))==0
            lastcharttime(i)=temp(j,3);
            lastvalue(i)=temp(j,i);
        end
        
        if j>1
        if isnan(temp(j,i)) &&  temp(j,2)==oldstayid && (temp(j,3)-lastcharttime(i))<=hold(i-3)*3600 %note : hold has 53 cols, temp has 55
            temp(j,i)=lastvalue(i);       
        end
        
        end
    end    
end            


close(h);

end

