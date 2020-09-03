# GPR_loocv
Code for the LOOCV calculation to determine threshold value for calculating debris thickness retrievals from GPR profiles.

Use the ground-truth thickness measurements in xval_tsB.csv for training the threshold for when the antenna unit is elevated 27 cm and those in xval_tsAC.csv for when the antenna unit is 19 cm.

xval_cng.m is the matlab script that will calculate the thresholds via leave one out cross validation (Arlot et al., 2010) associated with each measurement on the selected transect(s) of Changri Nup Glacier (A&C combined or B, which are seperated by antenna unit elevation).  scale_cng.m calculates the reported threshold from the lowest 10% of RMSEs and applies it to yield a thickness retrieval for the selected profile.

.DZT files collected on Changri Nup Glacier are available at https://glacioclim.osug.fr/, pending publication acceptance. 