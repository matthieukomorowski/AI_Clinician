%% AI Clinician Identifiying MIMIC-III sepsis cohort

% (c) Matthieu Komorowski, Imperial College London 2015-2019
% as seen in publication: https://www.nature.com/articles/s41591-018-0213-5

% version 16 Feb 19
% IDENTIFIES THE COHORT OF PATIENTS WITH SEPSIS in MIMIC-III

% PURPOSE:
% ------------------------------
% This creates a list of icustayIDs of patients who develop sepsis at some point 
% in the ICU. records charttime for onset of sepsis. Uses sepsis3 criteria

% STEPS:
% -------------------------------
% IMPORT DATA FROM CSV FILES
% FLAG PRESUMED INFECTION
% PREPROCESSING
% REFORMAT in 4h time slots
% COMPUTE SOFA at each time step
% FLAG SEPSIS

% note: the process generates the same features as the final MDP dataset, most of which are not used to compute SOFA
% External files required: Reflabs, Refvitals, sample_and_hold (all saved in reference_matrices.mat file)

% This code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE


%% ########################################################################
% IMPORT ALL DATA

tic
abx=table2array(readtable('D:/exportdir/abx.csv'));
culture=table2array(readtable('D:/exportdir/culture.csv'));
microbio=table2array(readtable('D:/exportdir/microbio.csv'));
demog=(readtable('D:/exportdir/demog.csv'));
ce010=table2array(readtable('D:/exportdir/ce010000.csv'));
ce1020=table2array(readtable('D:/exportdir/ce1000020000.csv'));
ce2030=table2array(readtable('D:/exportdir/ce2000030000.csv'));
ce3040=table2array(readtable('D:/exportdir/ce3000040000.csv'));
ce4050=table2array(readtable('D:/exportdir/ce4000050000.csv'));
ce5060=table2array(readtable('D:/exportdir/ce5000060000.csv'));
ce6070=table2array(readtable('D:/exportdir/ce6000070000.csv'));
ce7080=table2array(readtable('D:/exportdir/ce7000080000.csv'));
ce8090=table2array(readtable('D:/exportdir/ce8000090000.csv'));
ce90100=table2array(readtable('D:/exportdir/ce90000100000.csv'));
labU=[ table2array(readtable('D:/exportdir/labs_ce.csv')) ; table2array(readtable('D:/exportdir/labs_le.csv'))  ];
MV=table2array(readtable('D:/exportdir/mechvent.csv'));
inputpreadm=table2array(readtable('D:/exportdir/preadm_fluid.csv'));
inputMV=table2array(readtable('D:/exportdir/fluid_mv.csv'));
inputCV=table2array(readtable('D:/exportdir/fluid_cv.csv'));
vasoMV=table2array(readtable('D:/exportdir/vaso_mv.csv'));
vasoCV=table2array(readtable('D:/exportdir/vaso_cv.csv'));
UOpreadm=table2array(readtable('D:/exportdir/preadm_uo.csv'));
UO=table2array(readtable('D:/exportdir/uo.csv'));
toc

%% ########################################################################
%                       INITIAL DATA MANIPULATIONS
% #########################################################################

ii=isnan(microbio(:,3));  %if charttime is empty but chartdate isn't
microbio(ii,3)=microbio(ii,4);   %copy time
microbio( :,4)=[];    %delete chardate
% Add empty col in microbio (# 3 and #5)
microbio(:,4)=microbio(:,3);
microbio(:,[3 5])=0;
% Combine both tables for micro events
bacterio = [microbio ; culture];
% correct NaNs in DEMOG
demog.morta_90(isnan(demog.morta_90))=0;
demog.morta_hosp(isnan(demog.morta_hosp))=0;
demog.elixhauser(isnan(demog.elixhauser))=0;

% compute normalized rate of infusion
% if we give 100 ml of hypertonic fluid (600 mosm/l) at 100 ml/h (given in 1h) it is 200 ml of NS equivalent
% so the normalized rate of infusion is 200 ml/h (different volume in same duration)
inputMV(:,8)=inputMV(:,7).*inputMV(:,6)./inputMV(:,5);

% fill-in missing ICUSTAY IDs in bacterio
for i=1:size(bacterio,1)
if bacterio(i,3)==0   %if missing icustayid
    o=bacterio(i,4);  %charttime
    subjectid=bacterio(i,1);
    hadmid=bacterio(i,2);
   ii=find(demog.subject_id==subjectid);
   jj=find(demog.subject_id==subjectid & demog.hadm_id==hadmid);
    for j=1:numel(ii)
        if o>=demog.intime(ii(j))-48*3600 & o<=demog.outtime(ii(j))+48*3600
            bacterio(i,3)=demog.icustay_id(ii(j));
        elseif numel(ii)==1   %if we cant confirm from admission and discharge time but there is only 1 admission: it's the one!!
            bacterio(i,3)=demog.icustay_id(ii(j));
        end
    end
end   
end
toc


for i=1:size(bacterio,1)
if bacterio(i,3)==0   %if missing icustayid
    subjectid=bacterio(i,1);
    hadmid=bacterio(i,2);
    jj=find(demog.subject_id==subjectid & demog.hadm_id==hadmid);
    if numel(jj)==1
        bacterio(i,3)=demog.icustay_id(jj);
    end
end
end

% fill-in missing ICUSTAY IDs in ABx

for i=1:size(abx,1)
if isnan(abx(i,2))
    o=abx(i,3);  %time of event
    hadmid=abx(i,1);
    ii=find(demog.hadm_id==hadmid);   %row in table demographics
    for j=1:numel(ii)
        if o>=demog.intime(ii(j))-48*3600 & o<=demog.outtime(ii(j))+48*3600
            abx(i,2)=demog.icustay_id(ii(j));
        elseif numel(ii)==1   %if we cant confirm from admission and discharge time but there is only 1 admission: it's the one!!
            abx(i,2)=demog.icustay_id(ii(j));
        end
    end
end   
end

%% ########################################################################
%    find presumed onset of infection according to sepsis3 guidelines
% ########################################################################

% METHOD:
% I loop through all the ABx given, and as soon as there is a sample present
% within the required time criteria I pick this flag and break the loop.

onset=zeros(100000,3);

for icustayid=1:100000

    ab=abx(abx(:,2)==icustayid+200000,3);   %start time of abx for this icustayid
    bact=bacterio(bacterio(:,3)==icustayid+200000,4);  %time of sample
    subj_bact=bacterio(bacterio(:,3)==icustayid+200000,1);  %subjectid
    
    if ~isempty(ab) & ~isempty(bact)   %if we have data for both: proceed
        
      D = pdist2(ab, bact)/3600;  %pairwise distances btw antibio and cultures, in hours
      
      for i=1:size(D,1)  % looping through all rows of AB given, from early to late
        [M,I] = min(D(i,:));   %minimum distance in this row
        ab1=ab(i);       %timestamp of this value in list of antibio
        bact1=bact(I);      %timestamp in list of cultures
              
        if M<=24 & ab1<=bact1      %if ab was first and delay < 24h
            onset(icustayid,1)=subj_bact(1);   %subject_id
            onset(icustayid,2)=icustayid;       % icustay_id
            onset(icustayid,3)=ab1;     %onset of infection = abx time
              icustayid
            break
        elseif M<=72 & ab1>=bact1    %elseif sample was first and delay < 72h
            onset(icustayid,1)=subj_bact(1);
            onset(icustayid,2)=icustayid;
            onset(icustayid,3)=bact1;       %onset of infection = sample time
            break
        end
      end
    end    
end
toc

%sum of records found
sum(onset(:,3)>0)


%% Replacing item_ids with column numbers from reference tables

% replace itemid in labs with column number
% this will accelerate process later

tic
for i=10001:size(labU,1)
[~,locb]=ismember(Reflabs,labU(i,3));
labU(i,3)=find(max(locb')');
end
toc

% replace itemid in vitals with col number

for i=1:size(ce010,1)
[~,locb]=ismember(Refvitals,ce010(i,3));ce010(i,3)=find(max(locb')');
end
for i=1:size(ce1020,1)
[~,locb]=ismember(Refvitals,ce1020(i,3));ce1020(i,3)=find(max(locb')');
end
for i=1:size(ce2030,1)
[~,locb]=ismember(Refvitals,ce2030(i,3));ce2030(i,3)=find(max(locb')');
end
for i=1:size(ce3040,1)
[~,locb]=ismember(Refvitals,ce3040(i,3));ce3040(i,3)=find(max(locb')');
end
for i=1:size(ce4050,1)
[~,locb]=ismember(Refvitals,ce4050(i,3));ce4050(i,3)=find(max(locb')');
end
for i=1:size(ce5060,1)
[~,locb]=ismember(Refvitals,ce5060(i,3));ce5060(i,3)=find(max(locb')');
end
for i=1:size(ce6070,1)
[~,locb]=ismember(Refvitals,ce6070(i,3));ce6070(i,3)=find(max(locb')');
end
for i=1:size(ce7080,1)
[~,locb]=ismember(Refvitals,ce7080(i,3));ce7080(i,3)=find(max(locb')');
end
for i=1:size(ce8090,1)
[~,locb]=ismember(Refvitals,ce8090(i,3));ce8090(i,3)=find(max(locb')');
end
for i=1:size(ce90100,1)
[~,locb]=ismember(Refvitals,ce90100(i,3));ce90100(i,3)=find(max(locb')');
end



%% ########################################################################
%           INITIAL REFORMAT WITH CHARTEVENTS, LABS AND MECHVENT
% ########################################################################

% gives an array with all unique charttime (1 per row) and all items in columns.
% ################## IMPORTANT !!!!!!!!!!!!!!!!!!
% Here i use -48 -> +24 because that's for sepsis3 cohort defintion!!
% I need different time period for the MDP (-24 -> +48)


reformat=NaN(2000000,68);  %final table 
qstime=zeros(100000,4);
winb4=49;   %lower limit for inclusion of data (48h before time flag)
winaft=25;  % upper limit (24h after)
irow=1;  %recording row for summary table
h = waitbar(0,'Initializing waitbar...');

tic
for icustayid=1:100000
qst=onset(icustayid,3); %flag for presumed infection
if qst>0  % if we have a flag
d1=table2array(demog(demog.icustay_id==icustayid+200000,[11 5])); %age of patient + discharge time

if d1(1)>6574  % if older than 18 years old

    waitbar(icustayid/100000,h,icustayid/1000) %moved here to save some time
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
    reformat(irow,3+col)=value;%(locb(:,1)); %store available values

      
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

qstime(icustayid,1)=qst; %flag for presumed infection / this is time of sepsis if SOFA >=2 for this patient
%HERE I SAVE FIRST and LAST TIMESTAMPS, in QSTIME, for each ICUSTAYID
qstime(icustayid,2)=t(1);  %first timestamp
qstime(icustayid,3)=t(end);  %last timestamp
qstime(icustayid,4)=d1(2); %dischargetime

end
end
end
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

%temp
ii=reformat(:,14)>90 & isnan(reformat(:,15));reformat(ii,15)=reformat(ii,14);
reformat=deloutabove(reformat,14,90);

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
        rate=input(:,8);  %rate of infusion (is NaN for boluses) || corrected for tonicity
        
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
    
    
    for j=0:timestep:79 % -52 until +28 = 80 hours in total
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
        
          % #####################   DISCUSS ADDING STUFF HERE / RANGE, MIN, MAX ETC   ################
        
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


%% ########################################################################
%    CONVERT TO TABLE AND DELETE VARIABLES WITH EXCESSIVE MISSINGNESS
% ########################################################################

dataheaders=[sample_and_hold(1,:) {'Shock_Index' 'PaO2_FiO2'}]; 
dataheaders=regexprep(dataheaders,'['']','');
dataheaders = ['bloc','icustayid','charttime','gender','age','elixhauser','re_admission', 'died_in_hosp', 'died_within_48h_of_out_time','mortality_90d','delay_end_of_record_and_discharge_or_death',...
    dataheaders,  'median_dose_vaso','max_dose_vaso','input_total','input_4hourly','output_total','output_4hourly','cumulated_balance'];

reformat2t=array2table(reformat2);
reformat2t.Properties.VariableNames=dataheaders;
miss=sum(isnan(reformat2))./size(reformat2,1);

% if values have less than 70% missing values (over 30% of values present): I keep them
reformat3t=reformat2t(:,[true(1,11) miss(12:74)<0.70 true(1,11)]) ; 

%% ########################################################################
%             HANDLING OF MISSING VALUES  &  CREATE REFORMAT4T
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

% CORRECT GENDER
reformat4t.gender=reformat4t.gender-1; 

%CORRECT AGE > 200 yo
ii=reformat4t.age>150*365.25;
reformat4t.age(ii)=91.4*365.25;

% FIX MECHVENT
reformat4t.mechvent(isnan(reformat4t.mechvent))=0;
reformat4t.mechvent(reformat4t.mechvent>0)=1;

% FIX Elixhauser missing values
reformat5t.elixhauser(isnan(reformat5t.elixhauser))=nanmedian(reformat5t.elixhauser);  %use the median value / only a few missing data points 

%vasopressors / no NAN
a=find(ismember(reformat4t.Properties.VariableNames,{'median_dose_vaso'}));
ii=isnan(reformat4(:,a));
reformat4t(ii,a)=array2table(zeros(sum(ii),1));
a=find(ismember(reformat4t.Properties.VariableNames,{'max_dose_vaso'}));
ii=isnan(reformat4(:,a));
reformat4t(ii,a)=array2table(zeros(sum(ii),1));

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

% SOFA - at each timepoint
% need (in this order):  P/F  MV  PLT  TOT_BILI  MAP  NORAD(max)  GCS  CR  UO
a=zeros(8,1); % indices of vars used in SOFA
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

% records values
reformat4t(:,end-1)=array2table(reformat4(:,end-1));
reformat4t(:,end)=array2table(reformat4(:,end));


%% ########################################################################
%                            EXCLUSION OF SOME PATIENTS 
% ########################################################################

numel(unique(reformat4t.icustayid))  %count before

% check for patients with extreme UO = outliers = to be deleted (>40 litres of UO per 4h!!)
a=find(reformat4t.output_4hourly>12000);
i=unique(reformat4t.icustayid(a));
i=find(ismember(reformat4t.icustayid,i));
reformat4t(i,:)=[];

% some have bili = 999999
a=find(reformat4t.Total_bili>10000); 
i=unique(reformat4t.icustayid(a));
i=find(ismember(reformat4t.icustayid,i));
reformat4t(i,:)=[];

% check for patients with extreme INTAKE = outliers = to be deleted (>10 litres of intake per 4h!!)
a=find(reformat4t.input_4hourly>10000);
i=unique(reformat4t.icustayid(a));  % 28 ids
i=find(ismember(reformat4t.icustayid,i));
reformat4t(i,:)=[];


% #### exclude early deaths from possible withdrawals ####
% stats per patient
q=reformat4t.bloc==1;
% fence_posts=find(q(:,1)==1);
num_of_trials=numel(unique(reformat4t.icustayid));%size(fence_posts,1);
a=array2table([reformat4t.icustayid reformat4t.mortality_90d reformat4t.max_dose_vaso reformat4t.SOFA]);
a.Properties.VariableNames={'id','mortality_90d','vaso','sofa'};
d=grpstats(a,'id','max');

%finds patients who match our criteria
e=zeros(num_of_trials,1);
for i=1:num_of_trials
    if d.max_mortality_90d(i) ==1
    ii=reformat4t.icustayid==d.id(i) & reformat4t.bloc==d.GroupCount(i);  %last row for this patient
    e(i)=sum((reformat4t.max_dose_vaso(ii)==0 & d.max_vaso(i)>0.3 & reformat4t.SOFA(ii)>=d.max_sofa(i)/2))>0;
    end
end
r=d.id(e==1 & d.GroupCount<20); % ids to be removed
ii=ismember(reformat4t.icustayid,r);
reformat4t(ii,:)=[];

% exclude patients who died in ICU during data collection period
ii=reformat4t.bloc==1&reformat4t.died_within_48h_of_out_time==1& reformat4t.delay_end_of_record_and_discharge_or_death<24;
ii=ismember(icustayidlist,reformat4t.icustayid(ii));
reformat4t(ii,:)=[];

numel(unique(reformat4t.icustayid))   %count after

%% #######################################################################
%                       CREATE SEPSIS COHORT
% ########################################################################

% create array with 1 row per icu admission
% keep only patients with flagged sepsis (max sofa during time period of interest >= 2)
% we assume baseline SOFA of zero (like other publications)

sepsis=zeros(30000,5);
irow=1;

tic
for icustayid=1:100000
    ii=find(ismember(reformat4t.icustayid,icustayid));
    if mod(icustayid,10000)==0;disp([num2str(icustayid/1000), ' %']);end
    if ii
    
         sofa=reformat4t.SOFA(ii);
         sirs=reformat4t.SIRS(ii);
         sepsis(irow,1)=icustayid+200000; 
         sepsis(irow,2)=reformat4t.mortality_90d(ii(1)); % 90-day mortality
         sepsis(irow,3)=max(sofa);
         sepsis(irow,4)=max(sirs);
         sepsis(irow,5)=qstime(icustayid);   %time of onset of sepsis
         irow=irow+1;
    end
end
toc
sepsis(irow:end,:)=[];

sepsis=array2table(sepsis);
sepsis.Properties.VariableNames={'icustayid','morta_90d','max_sofa','max_sirs','sepsis_time'};

% delete all non-sepsis
sepsis(sepsis.max_sofa<2,:)=[];

% final count
size(sepsis,1)  

%save cohort
writetable(sepsis,'sepsis_mimiciii.csv','Delimiter',',');
