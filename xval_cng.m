%Cross-Validate: Leave-one-out cross validation from https://www.cs.cmu.edu/~schneide/tut5/node42.html
%Transects A,B,C ... leave out thicker, longitudinal transects (D,E) 
% (b/c all points may be minima; the debris is thicker; and the range is partly 50ns, partly 100ns)
%ALGiese

clear
close all
% ts='AC'; meas = csvread('xval_tsAC.csv');
% %-OR-
ts='B'; meas = csvread('xval_tsB.csv'); 
%% Import %columns: d, z, rg
d =meas(:,1); %dist from start of transects (unit: scans)
z =meas(:,2); %depth (cm)
rg=meas(:,3); %transect
threshold_matrix = nan(length(meas),length(meas)); 
rmse = nan(1,length(meas));

R_A= readgssi('tsA_hilbert.DZT'); R_A.head.epsr = 3; 
R_B= readgssi('tsB_hilbert.DZT'); R_B.head.epsr = 3; 
R_C= readgssi('tsC_hilbert.DZT'); R_C.head.epsr = 3; 

for i=1:length(meas) %*testing point*
    tic
for k=[1:i-1 i+1:length(meas)] %training points (all other points)
    
%load appropriate radargram (1 scan = 1 meter after stacking)
if rg(k)==1, R=R_B; h = .27; m= d(k); end % climber's R to L
if rg(k)==2, R=R_C; h = .19; m= d(k); end % climber's L to R
if rg(k)==3, R=R_A; h = .19; m= d(k); end % climber's R to L
R.samp = R.samp-2^15; 

%%gssi & display settings
A = size(R.samp);
vel = .3/sqrt(R.head.epsr); %velocity in debris (m/ns)
vel_air = .3;
oorange = R.head.range/2; %one-way range

t_antenna = h/vel_air; % This is how long it takes (ns) the EM wave to reach the surface from the antenna's elevation.
% h_antenna = t_antenna * vel_air; %Antenna "height" with a wave traveling as slowly as it does when eps = 3 (this is a check only, not used)
bin_antenna = round(t_antenna * 2* R.head.nsamp/R.head.range);

%linearly interpolate stacked radargrams to 10x columns (dm) (and still 1024 rows) 
S    = size(R.samp);
R.dm = zeros(S(1),S(2)*10);
for j = 1:length(R.samp)
   R.dm(j,:) = interp1(1:10:10*numel(R.samp(1,:)),R.samp(j,:),(1:10*numel(R.samp(1,:))));
end

%% Find threshold associated w/ all z (except that of i)
aucA = zeros(1,R.head.nsamp); %area under curve A (working variable)
bin_meas = round(z(i)/100/vel * R.head.nsamp/R.head.range); %measurement depth --> bin #
if m>length(R.dm), m=length(R.dm); disp(m); end

for n = 1:bin_meas+bin_antenna %ROWS [1 --> bin # that matches gt, with antenna included]
    aucA(n) = R.dm(n,m)-min(R.dm(:,m)); 
end
auc_cs = cumsum(aucA); %cumulative sum to ground-truth measurement

for n = 1:R.head.nsamp %ROWS [1 --> all samples] 100% of area under curve 
    aucA(n) = R.dm(n,m)-min(R.dm(:,m)); 
end

auc_full = cumsum(aucA); %total area under curve

threshold_matrix(k,i) = auc_cs(end)/auc_full(end); % percent of auc (threshold)
end 
    toc
    
%% Compare avg. of that(^) threshold to gt(i)
threshold = nanmean(threshold_matrix(:,i));
if rg(i)==1, R=R_B; h = .27; m= round(d(i)/10); end 
if rg(i)==2, R=R_C; h = .19; m= round(d(i)/10); end
if rg(i)==3, R=R_A; h = .19; m= round(d(i)/10); end
R.samp = R.samp-2^15;  

%preallocate space
aucB    = zeros(1,R.head.nsamp); %area under curve B
auc_cs  = zeros(1,R.head.nsamp); %cumulative sum

s2 = size(R.samp);
if m > s2(2), m = s2(2); disp (s2); end
for n = 1:R.head.nsamp
   aucB(n) = R.samp(n,m)-min(R.samp(:,m)); %should be dist, not index
end
auc_cs = cumsum(aucB(:));

intfc_bin = find(auc_cs >= threshold*auc_cs(end),1); %first bin where threshold is exceeded
intfc = vel*oorange * intfc_bin/R.head.nsamp; %bin --> depth

%avg mismatch for each training point
rmse(i) = sqrt(mean( ((z(i)/100-intfc).^2) ));
rmse(i), i
end
save(['cross_val_',ts,'_NAME.mat'])

