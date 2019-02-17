%% AI Clinician Building MIMIC-III dataset 

% (c) Matthieu Komorowski, Imperial College London 2015-2019
% as seen in publication: https://www.nature.com/articles/s41591-018-0213-5

% version 16 Feb 19
% uses the sepsis-3 cohort previously defined
% builds the MIMIC-III dataset

% This code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE

%% ########################################################################
%           INITIAL REFORMAT WITH CHARTEVENTS, LABS AND MECHVENT
% ########################################################################

% gives an array with all unique charttime (1 per row) and all items in columns.
% ################## IMPORTANT !!!!!!!!!!!!!!!!!!
% Here i use -24 -> +48 because that's for the MDP


reformat=NaN(2000000,68);  %final table 
qstime=zeros(100000,4);
winb4=25;   %lower limit for inclusion of data (48h before time flag)
winaft=49;  % upper limit (24h after)
irow=1;  %recording row for summary table
h = waitbar(0,'Initializing waitbar...');

tic
for icustayidrow=1:size(sepsis,1)
    
qst=sepsis.sepsis_time(icustayidrow);%,3); %flag for presumed infection
icustayid=sepsis.icustayid(icustayidrow)-200000;
waitbar(icustayidrow/size(sepsis,1),h,icustayidrow/size(sepsis,1)*100) %moved here to save some time


% CHARTEVENTS
    if icustayid<10000
    temp=ce010(ce010(:,1)==icustayid+200000,:);
    elseif icustayid>=10000 & icustayid<20000
    temp=ce1020(ce1020(:,1)==icustayid+200000,:);
    elseif icustayid>=20000 & icustayid<30000
    temp=ce2030(ce2030(:,1)==icustayid+200000,:);
    elseif icustayid>=30000 && icustayid<40000
    temp=ce3040(ce3040(:,1)==icustayid+200000,:);
    elseif icustayid>=40000 & icustayid<50000
    temp=ce4050(ce4050(:,1)==icustayid+200000,:);
    elseif icustayid>=50000 & icustayid<60000
    temp=ce5060(ce5060(:,1)==icustayid+200000,:);
    elseif icustayid>=60000 & icustayid<70000
    temp=ce6070(ce6070(:,1)==icustayid+200000,:);
    elseif icustayid>=70000 & icustayid<80000
    temp=ce7080(ce7080(:,1)==icustayid+200000,:);
    elseif icustayid>=80000 & icustayid<90000
    temp=ce8090(ce8090(:,1)==icustayid+200000,:);
    elseif icustayid>=90000
    temp=ce90100(ce90100(:,1)==icustayid+200000,:);
    end

ii=temp(:,2)>= qst-(winb4+4)*3600 & temp(:,2)<=qst+(winaft+4)*3600; %time period of interest -4h and +4h
temp=temp(ii,:);   %only time period of interest

%LABEVENTS
ii=labU(:,1)==icustayid+200000;
temp2=labU(ii,:);
ii=temp2(:,2)>= qst-(winb4+4)*3600 & temp2(:,2)<=qst+(winaft+4)*3600; %time period of interest -4h and +4h
temp2=temp2(ii,:);   %only time period of interest

%Mech Vent + ?extubated
ii=MV(:,1)==icustayid+200000;
temp3=MV(ii,:);
ii=temp3(:,2)>= qst-(winb4+4)*3600 & temp3(:,2)<=qst+(winaft+4)*3600; %time period of interest -4h and +4h
temp3=temp3(ii,:);   %only time period of interest

t=unique([temp(:,2);temp2(:,2); temp3(:,2)]);   %list of unique timestamps from all 3 sources / sorted in ascending order

if t
for i=1:numel(t)
    
    %CHARTEVENTS
    ii=temp(:,2)==t(i);
    col=temp(ii,3);
    value=temp(ii,4);  
    reformat(irow,1)=i; %timestep  
    reformat(irow,2)=icustayid;
    reformat(irow,3)=t(i); %charttime
    reformat(irow,3+col)=value; %store available values
      
    %LAB VALUES
    ii=temp2(:,2)==t(i);
    col=temp2(ii,3);
    value=temp2(ii,4);
    reformat(irow,31+col)=value; %store available values
      
    %MV  
    ii=temp3(:,2)==t(i);
    if nansum(ii)>0
    value=temp3(ii,3:4);
      reformat(irow,67:68)=value; %store available values
    else
      reformat(irow,67:68)=NaN;
    end
    
    irow=irow+1;
     
end

qstime(icustayid,1)=qst; %time of sepsis
%HERE I SAVE FIRST and LAST TIMESTAMPS, in QSTIME, for each ICUSTAYID
qstime(icustayid,2)=t(1);  %first timestamp
qstime(icustayid,3)=t(end);  %last timestamp
qstime(icustayid,4)=table2array(demog(demog.icustay_id==icustayid+200000,5)); % discharge time

end

% end
end
toc

close(h);
reformat(irow:end,:)=[];  %delete extra unused rows


%% ########################################################################
%                                   OUTLIERS 
% ########################################################################

%weight
reformat=deloutabove(reformat,5,300);  %delete outlier above a threshold (300 kg), for variable # 5

%HR
reformat=deloutabove(reformat,8,250);

%BP
reformat=deloutabove(reformat,9,300);
reformat=deloutbelow(reformat,10,0);
reformat=deloutabove(reformat,10,200);
reformat=deloutbelow(reformat,11,0);
reformat=deloutabove(reformat,11,200);

%RR
reformat=deloutabove(reformat,12,80);

%SpO2
reformat=deloutabove(reformat,13,150);
ii=reformat(:,13)>100;reformat(ii,13)=100;
reformat=deloutbelow(reformat,13,50);

%temp
ii=reformat(:,14)>90 & isnan(reformat(:,15));reformat(ii,15)=reformat(ii,14);
reformat=deloutabove(reformat,14,90);
reformat=deloutbelow(reformat,14,25);

%interface / is in col 22

% FiO2
reformat=deloutabove(reformat,23,100);
ii=reformat(:,23)<1;reformat(ii,23)=reformat(ii,23)*100;
reformat=deloutbelow(reformat,23,20);
reformat=deloutabove(reformat,24,1.5);

% O2 FLOW
reformat=deloutabove(reformat,25,70);

%PEEP
reformat=deloutbelow(reformat,26,0);
reformat=deloutabove(reformat,26,40);

%TV
reformat=deloutabove(reformat,27,1800);

%MV
reformat=deloutabove(reformat,28,50);

%K+
reformat=deloutbelow(reformat,32,1);
reformat=deloutabove(reformat,32,15);

%Na
reformat=deloutbelow(reformat,33,95);
reformat=deloutabove(reformat,33,178);

%Cl
reformat=deloutbelow(reformat,34,70);
reformat=deloutabove(reformat,34,150);

%Glc
reformat=deloutbelow(reformat,35,1);
reformat=deloutabove(reformat,35,1000);

%Creat
reformat=deloutabove(reformat,37,150);

%Mg
reformat=deloutabove(reformat,38,10);

%Ca
reformat=deloutabove(reformat,39,20);

%ionized Ca
reformat=deloutabove(reformat,40,5);

%CO2
reformat=deloutabove(reformat,41,120);

%SGPT/SGOT
reformat=deloutabove(reformat,42,10000);
reformat=deloutabove(reformat,43,10000);

%Hb/Ht
reformat=deloutabove(reformat,50,20);
reformat=deloutabove(reformat,51,65);

%WBC
reformat=deloutabove(reformat,53,500);

%plt
reformat=deloutabove(reformat,54,2000);

%INR
reformat=deloutabove(reformat,58,20);

%pH
reformat=deloutbelow(reformat,59,6.7);
reformat=deloutabove(reformat,59,8);

%po2
reformat=deloutabove(reformat,60,700);

%pco2
reformat=deloutabove(reformat,61,200);

%BE
reformat=deloutbelow(reformat,62,-50);

%lactate
reformat=deloutabove(reformat,63,30);

% ####################################################################
% some more data manip / imputation from existing values

% estimate GCS from RASS - data from Wesley JAMA 2003
ii=isnan(reformat(:,6))&reformat(:,7)>=0;
reformat(ii,6)=15;
ii=isnan(reformat(:,6))&reformat(:,7)==-1;
reformat(ii,6)=14;
ii=isnan(reformat(:,6))&reformat(:,7)==-2;
reformat(ii,6)=12;
ii=isnan(reformat(:,6))&reformat(:,7)==-3;
reformat(ii,6)=11;
ii=isnan(reformat(:,6))&reformat(:,7)==-4;
reformat(ii,6)=6;
ii=isnan(reformat(:,6))&reformat(:,7)==-5;
reformat(ii,6)=3;


% FiO2
ii=~isnan(reformat(:,23)) & isnan(reformat(:,24));
reformat(ii,24)=reformat(ii,23)./100;
ii=~isnan(reformat(:,24)) & isnan(reformat(:,23));
reformat(ii,23)=reformat(ii,24).*100;


%ESTIMATE FiO2 /// with use of interface / device (cannula, mask, ventilator....)

reformatsah=SAH(reformat,sample_and_hold);  % do SAH first to handle this task

%NO FiO2, YES O2 flow, no interface OR cannula
ii=find(isnan(reformatsah(:,23))&~isnan(reformatsah(:,25))&(reformatsah(:,22)==0|reformatsah(:,22)==2)); 
reformat(ii(reformatsah(ii,25)<=15),23)=70;
reformat(ii(reformatsah(ii,25)<=12),23)=62;
reformat(ii(reformatsah(ii,25)<=10),23)=55;
reformat(ii(reformatsah(ii,25)<=8),23)=50;
reformat(ii(reformatsah(ii,25)<=6),23)=44;
reformat(ii(reformatsah(ii,25)<=5),23)=40;
reformat(ii(reformatsah(ii,25)<=4),23)=36;
reformat(ii(reformatsah(ii,25)<=3),23)=32;
reformat(ii(reformatsah(ii,25)<=2),23)=28;
reformat(ii(reformatsah(ii,25)<=1),23)=24;

%NO FiO2, NO O2 flow, no interface OR cannula
ii=find(isnan(reformatsah(:,23))&isnan(reformatsah(:,25))&(reformatsah(:,22)==0|reformatsah(:,22)==2));  %no fio2 given and o2flow given, no interface OR cannula
reformat(ii,23)=21;

%NO FiO2, YES O2 flow, face mask OR.... OR ventilator (assume it's face mask)
ii=find(isnan(reformatsah(:,23))&~isnan(reformatsah(:,25))&(reformatsah(:,22)==NaN|reformatsah(:,22)==1|reformatsah(:,22)==3|reformatsah(:,22)==4|reformatsah(:,22)==5|reformatsah(:,22)==6|reformatsah(:,22)==9|reformatsah(:,22)==10)); 
reformat(ii(reformatsah(ii,25)<=15),23)=75;
reformat(ii(reformatsah(ii,25)<=12),23)=69;
reformat(ii(reformatsah(ii,25)<=10),23)=66;
reformat(ii(reformatsah(ii,25)<=8),23)=58;
reformat(ii(reformatsah(ii,25)<=6),23)=40;
reformat(ii(reformatsah(ii,25)<=4),23)=36;

%NO FiO2, NO O2 flow, face mask OR ....OR ventilator
ii=find(isnan(reformatsah(:,23))&isnan(reformatsah(:,25))&(reformatsah(:,22)==NaN|reformatsah(:,22)==1|reformatsah(:,22)==3|reformatsah(:,22)==4|reformatsah(:,22)==5|reformatsah(:,22)==6|reformatsah(:,22)==9|reformatsah(:,22)==10));  %no fio2 given and o2flow given, no interface OR cannula
reformat(ii,23)=NaN;

%NO FiO2, YES O2 flow, Non rebreather mask
ii=find(isnan(reformatsah(:,23))&~isnan(reformatsah(:,25))&reformatsah(:,22)==7); 
reformat(ii(reformatsah(ii,25)>=10),23)=90;
reformat(ii(reformatsah(ii,25)>=15),23)=100;
reformat(ii(reformatsah(ii,25)<10),23)=80;
reformat(ii(reformatsah(ii,25)<=8),23)=70;
reformat(ii(reformatsah(ii,25)<=6),23)=60;

%NO FiO2, NO O2 flow, NRM
ii=find(isnan(reformatsah(:,23))&isnan(reformatsah(:,25))&reformatsah(:,22)==7);  %no fio2 given and o2flow given, no interface OR cannula
reformat(ii,23)=NaN;

% update again FiO2 columns
ii=~isnan(reformat(:,23)) & isnan(reformat(:,24));
reformat(ii,24)=reformat(ii,23)./100;
ii=~isnan(reformat(:,24)) & isnan(reformat(:,23));
reformat(ii,23)=reformat(ii,24).*100;

%BP
ii=~isnan(reformat(:,9))&~isnan(reformat(:,10)) & isnan(reformat(:,11));
reformat(ii,11)=(3*reformat(ii,10)-reformat(ii,9))./2;
ii=~isnan(reformat(:,09))&~isnan(reformat(:,11)) & isnan(reformat(:,10));
reformat(ii,10)=(reformat(ii,9)+2*reformat(ii,11))./3;
ii=~isnan(reformat(:,10))&~isnan(reformat(:,11)) & isnan(reformat(:,9));
reformat(ii,9)=3*reformat(ii,10)-2*reformat(ii,11);

%TEMP
%some values recorded in the wrong column
ii=reformat(:,15)>25&reformat(:,15)<45; %tempF close to 37deg??!
reformat(ii,14)=reformat(ii,15);
reformat(ii,15)=NaN;
ii=reformat(:,14)>70;  %tempC > 70?!!! probably degF
reformat(ii,15)=reformat(ii,14);
reformat(ii,14)=NaN;
ii=~isnan(reformat(:,14)) & isnan(reformat(:,15));
reformat(ii,15)=reformat(ii,14)*1.8+32;
ii=~isnan(reformat(:,15)) & isnan(reformat(:,14));
reformat(ii,14)=(reformat(ii,15)-32)./1.8;

% Hb/Ht
ii=~isnan(reformat(:,50)) & isnan(reformat(:,51));
reformat(ii,51)=(reformat(ii,50)*2.862)+1.216;
ii=~isnan(reformat(:,51)) & isnan(reformat(:,50));
reformat(ii,50)=(reformat(ii,51)-1.216)./2.862;

%BILI
ii=~isnan(reformat(:,44)) & isnan(reformat(:,45));
reformat(ii,45)=(reformat(ii,44)*0.6934)-0.1752;
ii=~isnan(reformat(:,45)) & isnan(reformat(:,44));
reformat(ii,44)=(reformat(ii,45)+0.1752)./0.6934;


%% ########################################################################
%                      SAMPLE AND HOLD on RAW DATA
% ########################################################################

reformat=SAH(reformat(:,1:68),sample_and_hold);


%% ########################################################################
%                             DATA COMBINATION
% ########################################################################

tic
     save('D:\BACKUP MIT PC\Data_100219.mat', '-v7.3');
toc


% WARNING: the time window of interest has been defined above (here -48 -> +24)! 

timestep=4;  %resolution of timesteps, in hours
irow=1;
icustayidlist=unique(reformat(:,2));
reformat2=nan(size(reformat,1),84);  %output array
h = waitbar(0,'Initializing waitbar...');
npt=numel(icustayidlist);  %number of patients
% Adding 2 empty cols for future shock index=HR/SBP and P/F
reformat(:,69:70)=NaN(size(reformat,1),2);

tic
for i=1:npt
    
    icustayid=icustayidlist(i);  %1 to 100000, NOT 200 to 300K!
     
        %CHARTEVENTS AND LAB VALUES
        temp=reformat(reformat(:,2)==icustayid,:);   %subtable of interest
        beg=temp(1,3);   %timestamp of first record
    
        % IV FLUID STUFF
        iv=find(inputMV(:,1)==icustayid+200000);   %rows of interest in inputMV
        input=inputMV(iv,:);    %subset of interest
        iv=find(inputCV(:,1)==icustayid+200000);   %rows of interest in inputCV
        input2=inputCV(iv,:);    %subset of interest
        startt=input(:,2); %start of all infusions and boluses
        endt=input(:,3); %end of all infusions and boluses
        rate=input(:,8);  %rate of infusion (is NaN for boluses)
        
        pread=inputpreadm(inputpreadm(:,1)==icustayid+200000,2) ;%preadmission volume
            if ~isempty(pread)             %store the value, if available
                totvol=nansum(pread);
                waitbar(i/npt,h,i/npt*100) %moved here to save some time
            else
                totvol=0;   %if not documented: it's zero
            end
       
        % compute volume of fluid given before start of record!!!
        t0=0;
        t1=beg;
        %input from MV (4 ways to compute)
        infu=  nansum(rate.*(endt-startt).*(endt<=t1&startt>=t0)/3600   +    rate.*(endt-t0).*(startt<=t0&endt<=t1&endt>=t0)/3600 +     rate.*(t1-startt).*(startt>=t0&endt>=t1&startt<=t1)/3600 +      rate.*(t1-t0).*(endt>=t1&startt<=t0)   /3600);
        %all boluses received during this timestep, from inputMV (need to check rate is NaN) and inputCV (simpler):
        bolus=nansum(input(isnan(input(:,6))& input(:,2)>=t0&input(:,2)<=t1,7)) + nansum(input2(input2(:,2)>=t0&input2(:,2)<=t1,5));  
        totvol=nansum([totvol,infu,bolus]); 
            
        %VASOPRESSORS    
        iv=find(vasoMV(:,1)==icustayid+200000);   %rows of interest in vasoMV
        vaso1=vasoMV(iv,:);    %subset of interest
        iv=find(vasoCV(:,1)==icustayid+200000);   %rows of interest in vasoCV
        vaso2=vasoCV(iv,:);    %subset of interest
        startv=vaso1(:,3); %start of VP infusion
        endv=vaso1(:,4); %end of VP infusions
        ratev=vaso1(:,5);  %rate of VP infusion
            

        %DEMOGRAPHICS / gender, age, elixhauser, re-admit, died in hosp?, died within
        %48h of out_time (likely in ICU or soon after), died within 90d after admission?        
        demogi=find(demog.icustay_id==icustayid+200000);        
        dem=[  demog.gender(demogi) ; demog.age(demogi) ;demog.elixhauser(demogi) ; demog.adm_order(demogi)>1 ;  demog.morta_hosp(demogi); abs(demog.dod(demogi)-demog.outtime(demogi))<(24*3600*2); demog.morta_90(demogi) ; (qstime(icustayid,4)-qstime(icustayid,3))/3600];     
        
        
        % URINE OUTPUT
        iu=find(UO(:,1)==icustayid+200000);   %rows of interest in inputMV
        output=UO(iu,:);    %subset of interest
        pread=UOpreadm(UOpreadm(:,1)==icustayid,4) ;%preadmission UO
            if ~isempty(pread)     %store the value, if available
                UOtot=nansum(pread);
            else
                UOtot=0;
            end
        % adding the volume of urine produced before start of recording!    
        UOnow=nansum(output(output(:,2)>=t0&output(:,2)<=t1,4));  %t0 and t1 defined above
        UOtot=nansum([UOtot UOnow]);
    
    
    for j=0:timestep:79 % -28 until +52 = 80 hours in total
        t0=3600*j+ beg;   %left limit of time window
        t1=3600*(j+timestep)+beg;   %right limit of time window
        ii=temp(:,3)>=t0 & temp(:,3)<=t1;  %index of items in this time period
        if sum(ii)>0
            
            
        %ICUSTAY_ID, OUTCOMES, DEMOGRAPHICS
        reformat2(irow,1)=(j/timestep)+1;   %'bloc' = timestep (1,2,3...)
        reformat2(irow,2)=icustayid;        %icustay_ID
        reformat2(irow,3)=3600*j+ beg;      %t0 = lower limit of time window
        reformat2(irow,4:11)=dem;           %demographics and outcomes
            
        
        %CHARTEVENTS and LAB VALUES (+ includes empty cols for shock index and P/F)
        value=temp(ii,:);%records all values in this timestep
                
        if sum(ii)==1   %if only 1 row of values at this timestep
          reformat2(irow,12:78)=value(:,4:end);
        else
          reformat2(irow,12:78)=nanmean(value(:,4:end)); %mean of all available values
          
        end
        
        
        %VASOPRESSORS
            % for CV: dose at timestamps.
            % for MV: 4 possibles cases, each one needing a different way to compute the dose of VP actually administered:
            %----t0---start----end-----t1----
            %----start---t0----end----t1----
            %-----t0---start---t1---end
            %----start---t0----t1---end----

        
        %MV
        v=(endv>=t0&endv<=t1)|(startv>=t0&endv<=t1)|(startv>=t0&startv<=t1)|(startv<=t0&endv>=t1);
        %CV
        v2=vaso2(vaso2(:,3)>=t0&vaso2(:,3)<=t1,4);
        v1=nanmedian([ratev(v); v2]);
        v2=nanmax([ratev(v); v2]);
        if ~isempty(v1)&~isnan(v1)&~isempty(v2)&~isnan(v2)
        reformat2(irow,79)=v1;    %median of dose of VP
        reformat2(irow,80)=v2;    %max dose of VP
        end
        
        %INPUT FLUID
        %input from MV (4 ways to compute)
        infu=  nansum(rate.*(endt-startt).*(endt<=t1&startt>=t0)/3600   +    rate.*(endt-t0).*(startt<=t0&endt<=t1&endt>=t0)/3600 +     rate.*(t1-startt).*(startt>=t0&endt>=t1&startt<=t1)/3600 +      rate.*(t1-t0).*(endt>=t1&startt<=t0)   /3600);
        %all boluses received during this timestep, from inputMV (need to check rate is NaN) and inputCV (simpler):
        bolus=nansum(input(isnan(input(:,6))& input(:,2)>=t0&input(:,2)<=t1,7)) + nansum(input2(input2(:,2)>=t0&input2(:,2)<=t1,5));  
 
        %sum fluid given
        totvol=nansum([totvol,infu,bolus]);
        reformat2(irow,81)=totvol;    %total fluid given
        reformat2(irow,82)=nansum([infu,bolus]);   %fluid given at this step
        
        %UO
        UOnow=nansum(output(output(:,2)>=t0&output(:,2)<=t1,4));  
        UOtot=nansum([UOtot UOnow]);
        reformat2(irow,83)=UOtot;    %total UO
        reformat2(irow,84)=nansum(UOnow);   %UO at this step

        %CUMULATED BALANCE
        reformat2(irow,85)=totvol-UOtot;    %cumulated balance

        irow=irow+1;
        end
    end
end
toc

reformat2(irow:end,:)=[];
close(h);


tic
     save('D:\BACKUP MIT PC\Data_110219.mat', '-v7.3');
toc

%% ########################################################################
%             CONVERT TO TABLE AND KEEP ONLY WANTED VARIABLE
% ########################################################################

dataheaders=[sample_and_hold(1,:) {'Shock_Index' 'PaO2_FiO2'}]; 
dataheaders=regexprep(dataheaders,'['']','');
dataheaders = ['bloc','icustayid','charttime','gender','age','elixhauser','re_admission', 'died_in_hosp', 'died_within_48h_of_out_time','mortality_90d','delay_end_of_record_and_discharge_or_death',...
    dataheaders,  'median_dose_vaso','max_dose_vaso','input_total','input_4hourly','output_total','output_4hourly','cumulated_balance'];

reformat2t=array2table(reformat2);
reformat2t.Properties.VariableNames=dataheaders;

% headers I want to keep
dataheaders5 = {'bloc','icustayid','charttime','gender','age','elixhauser','re_admission', 'died_in_hosp', 'died_within_48h_of_out_time','mortality_90d','delay_end_of_record_and_discharge_or_death','SOFA','SIRS',...
    'Weight_kg','GCS','HR','SysBP','MeanBP','DiaBP','RR','SpO2','Temp_C','FiO2_1','Potassium','Sodium','Chloride','Glucose',...
    'BUN','Creatinine','Magnesium','Calcium','Ionised_Ca','CO2_mEqL','SGOT','SGPT','Total_bili','Albumin','Hb','WBC_count','Platelets_count','PTT','PT','INR',...
    'Arterial_pH','paO2','paCO2','Arterial_BE','HCO3','Arterial_lactate','mechvent','Shock_Index','PaO2_FiO2',...
    'median_dose_vaso','max_dose_vaso','input_total','input_4hourly','output_total','output_4hourly','cumulated_balance'};

ii=find(ismember(reformat2t.Properties.VariableNames,dataheaders5));
reformat3t=reformat2t(:,ii); 


%% SOME DATA MANIP BEFORE IMPUTATION

% CORRECT GENDER
reformat3t.gender=reformat3t.gender-1; 

%CORRECT AGE > 200 yo
ii=reformat3t.age>150*365.25;
reformat3t.age(ii)=91.4*365.25;

% FIX MECHVENT
reformat3t.mechvent(isnan(reformat3t.mechvent))=0;
reformat3t.mechvent(reformat3t.mechvent>0)=1;

% FIX Elixhauser missing values
reformat3t.elixhauser(isnan(reformat3t.elixhauser))=nanmedian(reformat3t.elixhauser);  %use the median value / only a few missing data points 

%vasopressors / no NAN
a=find(ismember(reformat3t.Properties.VariableNames,{'median_dose_vaso'}));
ii=isnan(table2array(reformat3t(:,a)));
reformat3t(ii,a)=array2table(zeros(sum(ii),1));
a=find(ismember(reformat3t.Properties.VariableNames,{'max_dose_vaso'}));
ii=isnan(table2array(reformat3t(:,a)));
reformat3t(ii,a)=array2table(zeros(sum(ii),1));

% check prop of missingness here
miss=array2table(sum(isnan(table2array(reformat3t)))./size(reformat3t,1));
miss.Properties.VariableNames=reformat3t.Properties.VariableNames;
figure; bar(sum(isnan(table2array(reformat3t)))./size(reformat3t,1))

% I fill the values temporarily with zeros, otherwise kNN imp doesnt work
reformat3t.Shock_Index=zeros(size(reformat3t,1),1);
reformat3t.PaO2_FiO2=zeros(size(reformat3t,1),1);


%% ########################################################################
%        HANDLING OF MISSING VALUES & CREATE REFORMAT4T
% ########################################################################

% Do linear interpol where missingness is low (kNN imputation doesnt work if all rows have missing values)
reformat3=table2array(reformat3t);
miss=sum(isnan((reformat3)))./size(reformat3,1);
ii=miss>0&miss<0.05;  %less than 5% missingness
mechventcol=find(ismember(reformat3t.Properties.VariableNames,{'mechvent'}));

for i=11:mechventcol-1 % correct col by col, otherwise it does it wrongly
  if ii(i)==1
    reformat3(:,i)=fixgaps(reformat3(:,i));
  end
end

reformat3t(:,11:mechventcol-1)=array2table(reformat3(:,11:mechventcol-1));

% KNN IMPUTATION -  Done on chunks of 10K records.

mechventcol=find(ismember(reformat3t.Properties.VariableNames,{'mechvent'}));
ref=reformat3(:,11:mechventcol-1);  %columns of interest

tic
for i=1:10000:size(reformat3,1)-9999   %dataset divided in 5K rows chunks (otherwise too large)
    i
    ref(i:i+9999,:)=knnimpute(ref(i:i+9999,:)',1, 'distance','seuclidean')';
end

ref(end-9999:end,:)=knnimpute(ref(end-9999:end,:)',1, 'distance','seuclidean')';  %the last bit is imputed from the last 10K rows

toc

% I paste the data interpolated, but not the demographics and the treatments
reformat3t(:,11:mechventcol-1)=array2table(ref);  

reformat4t=reformat3t;
reformat4=table2array(reformat4t);

%% ########################################################################
%        COMPUTE SOME DERIVED VARIABLES: P/F, Shock Index, SOFA, SIRS...
% ########################################################################


% re-compute P/F with no missing values...
p=find(ismember(reformat4t.Properties.VariableNames,{'paO2'}));
f=find(ismember(reformat4t.Properties.VariableNames,{'FiO2_1'}));
a=find(ismember(reformat4t.Properties.VariableNames,{'PaO2_FiO2'}));
reformat4t(:,a)=array2table(reformat4(:,p)./reformat4(:,f));  

%recompute SHOCK INDEX without NAN and INF
p=find(ismember(reformat4t.Properties.VariableNames,{'HR'}));
f=find(ismember(reformat4t.Properties.VariableNames,{'SysBP'}));
a=find(ismember(reformat4t.Properties.VariableNames,{'Shock_Index'}));
reformat4(:,a)=reformat4(:,p)./reformat4(:,f);  
reformat4(isinf(reformat4(:,a)),a)=NaN;
d=nanmean(reformat4(:,a));
reformat4(isnan(reformat4(:,a)),a)=d;  %replace NaN with average value ~ 0.8
reformat4t(:,a)=array2table(reformat4(:,a));
ii=reformat4t.Shock_Index>=quantile(reformat4t.Shock_Index,0.999); %replace outliers with 99.9th percentile
reformat4t.Shock_Index(ii)=quantile(reformat4t.Shock_Index,0.999);

% SOFA - at each timepoint
% need (in this order):  P/F  MV  PLT  TOT_BILI  MAP  NORAD(max)  GCS  CR  UO
a=zeros(8,1);                              % indices of vars used in SOFA
a(1)=find(ismember(reformat4t.Properties.VariableNames,{'PaO2_FiO2'}));
a(2)=find(ismember(reformat4t.Properties.VariableNames,{'Platelets_count'}));
a(3)=find(ismember(reformat4t.Properties.VariableNames,{'Total_bili'}));
a(4)=find(ismember(reformat4t.Properties.VariableNames,{'MeanBP'}));
a(5)=find(ismember(reformat4t.Properties.VariableNames,{'max_dose_vaso'}));
a(6)=find(ismember(reformat4t.Properties.VariableNames,{'GCS'}));
a(7)=find(ismember(reformat4t.Properties.VariableNames,{'Creatinine'}));
a(8)=find(ismember(reformat4t.Properties.VariableNames,{'output_4hourly'}));
s=table2array(reformat4t(:,a));  

p=[0 1 2 3 4];

s1=[s(:,1)>400 s(:,1)>=300 &s(:,1)<400 s(:,1)>=200 &s(:,1)<300 s(:,1)>=100 &s(:,1)<200 s(:,1)<100 ];   %count of points for all 6 criteria of sofa
s2=[s(:,2)>150 s(:,2)>=100 &s(:,2)<150 s(:,2)>=50 &s(:,2)<100 s(:,2)>=20 &s(:,2)<50 s(:,2)<20 ];
s3=[s(:,3)<1.2 s(:,3)>=1.2 &s(:,3)<2 s(:,3)>=2 &s(:,3)<6 s(:,3)>=6 &s(:,3)<12 s(:,3)>12 ];
s4=[s(:,4)>=70 s(:,4)<70&s(:,4)>=65 s(:,4)<65 s(:,5)>0 &s(:,5)<=0.1 s(:,5)>0.1 ];
s5=[s(:,6)>14 s(:,6)>12 &s(:,6)<=14 s(:,6)>9 &s(:,6)<=12 s(:,6)>5 &s(:,6)<=9 s(:,6)<=5 ];
s6=[s(:,7)<1.2 s(:,7)>=1.2 &s(:,7)<2 s(:,7)>=2 &s(:,7)<3.5 (s(:,7)>=3.5 &s(:,7)<5)|(s(:,8)<84) (s(:,7)>5)|(s(:,8)<34) ];

nrcol=size(reformat4,2);   %nr of variables in data
reformat4(1,nrcol+1:nrcol+7)=0;  
for i=1:size(reformat4,1)  
    t=max(p(s1(i,:)))+max(p(s2(i,:)))+max(p(s3(i,:)))+max(p(s4(i,:)))+max(p(s5(i,:)))+max(p(s6(i,:)));  %SUM OF ALL 6 CRITERIA
    
    if t
    reformat4(i,nrcol+1:nrcol+7)=    [max(p(s1(i,:))) max(p(s2(i,:))) max(p(s3(i,:))) max(p(s4(i,:))) max(p(s5(i,:))) max(p(s6(i,:))) t];
    end
end

% SIRS - at each timepoint |  need: temp HR RR PaCO2 WBC 
a=zeros(5,1); % indices of vars used in SOFA
a(1)=find(ismember(reformat4t.Properties.VariableNames,{'Temp_C'}));
a(2)=find(ismember(reformat4t.Properties.VariableNames,{'HR'}));
a(3)=find(ismember(reformat4t.Properties.VariableNames,{'RR'}));
a(4)=find(ismember(reformat4t.Properties.VariableNames,{'paCO2'}));
a(5)=find(ismember(reformat4t.Properties.VariableNames,{'WBC_count'}));
s=table2array(reformat4t(:,a));  

s1=[s(:,1)>=38| s(:,1)<=36];   %count of points for all criteria of SIRS
s2=[s(:,2)>90 ];
s3=[s(:,3)>=20|s(:,4)<=32];
s4=[s(:,5)>=12| s(:,5)<4];
reformat4(:,nrcol+8)=s1+s2+s3+s4;

% adds 2 cols for SOFA and SIRS, if necessary
if sum(ismember(reformat4t.Properties.VariableNames,{'SIRS'}))== 0
reformat4t(:,end+1:end+2)=array2table(0);
reformat4t.Properties.VariableNames(end-1:end)= {'SOFA','SIRS'};  
end

% more IO corrections
 ii=reformat4t.input_total<0;
 reformat4t.input_total(ii)=0;
 ii=reformat4t.input_4hourly<0;
 reformat4t.input_4hourly(ii)=0;

% records values
reformat4t(:,end-1)=array2table(reformat4(:,end-1));
reformat4t(:,end)=array2table(reformat4(:,end));


%% ########################################################################
%                     CREATE FINAL MIMIC_TABLE
% ########################################################################

MIMICtable = reformat4t;

