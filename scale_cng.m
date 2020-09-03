%use threshold calculated in xval_new.m to get thickness retrievals
%ALGiese

clear
close all

%% CALCULATE THRESHOLD: Added Nov. 12
% load ('cross_val_AC_NAME.mat','rmse','threshold_matrix'); %*** %load rmse, threshold_matrix for each transect
load ('cross_val_B_NAME.mat','rmse','threshold_matrix'); 

threshold_mean=nanmean(threshold_matrix);
z=sort(rmse, 'ascend');

y = round(length(z)/10); 

%avg. thresholds associated with lowest 10% rmse values
for i=1:y 
    rt(i)=mean(threshold_mean(find(rmse==z(i)))); %ranked thersholds 
end

meanRT=mean(rt) %--> 42% for AC, applied to all transects but B (53% for B)

%%
threshold = round(meanRT*100)/100; 

for ag=1 %choose transect to which to apply threshold and calculate thickness retrievals
    
if ag==1, R=readgssi('tsB_hilbert.DZT');    h = .27;  ts = 'S-->B'; load transS_db.mat; end
if ag==2, R=readgssi('tsC_hilbert.DZT');    h = .19;  ts = 'W-->C'; load transW_db.mat; end
if ag==3, R=readgssi('tsA_hilbert.DZT');    h = .19;  ts = 'Q-->A'; load transQ_db.mat; end
if ag==4, R=readgssi('tsE_A_hilbert.DZT');  h = .19; b = 0;    ts = 'ZAa-->Ea'; load transZA_A_db.mat; end 
if ag==5, R=readgssi('tsE_B_hilbert.DZT');  h = .19; b = 799;  ts = 'ZAb-->Eb'; load transZA_B.mat; end
if ag==6, R=readgssi('tsD_A_hilbert.DZT');  h = .19; b = 0;    ts = 'PSa-->Da'; load transPS_A.mat; end
if ag==7, R=readgssi('tsD_B_hilbert.DZT');  h = .19; b = 1000; ts = 'PSb-->Db'; load transPS_B.mat; end

%%gssi & display settings
A = size(R.samp);
vel = .3/sqrt(R.head.epsr); %velocity in debris (m/ns)
vel_air = .3;
oorange = R.head.range/2; %one-way range

t_antenna = h/vel_air; % This is how long it takes (ns) the EM wave to reach the surface from the antenna's elevation.
% h_antenna = t_antenna * vel_air; %Antenna "height" with a wave traveling as slowly as it does when eps = 3 (this is a check only, not used)
bin_antenna = round(t_antenna * 2* R.head.nsamp/R.head.range);

% Preallocation
auc    = zeros(length(R.samp),numel(R.samp(1,:))); %area under curve
auc_cs = zeros(length(R.samp),numel(R.samp(1,:))); %cumulative sum
intfc_bin  = zeros(1,numel(R.samp(1,:))); %interface bin vector
intfc  = zeros(1,numel(R.samp(1,:))); %interface depth vector

% Computation of thicknes retrievals across the entire specified radargram
for j = 1:numel(R.samp(1,:)) %COLUMNS of radargram
    for i = 1:R.head.nsamp %ROWS
        auc(i,j) = R.samp(i,j)-min(R.samp(:,j)); %OR COULD DO ZERO; keep consistent with way threshold was calculated 
    end
    auc_cs(:,j) = cumsum(auc(:,j));

    intfc_bin(j) = find(auc_cs(:,j) >= threshold*auc_cs(end,j),1); %first bin where threshold is exceeded
    intfc(j) = (intfc_bin(j)-bin_antenna) * oorange/R.head.nsamp * vel; 
end
end

intfc(intfc<0)=nan; intfc'
meanINTFC= nanmean(intfc)*100
stdINTFC = nanstd(intfc)*100
medINTFC = nanmedian(intfc)*100
iqrINTFC = prctile(intfc,[25 75])*100

