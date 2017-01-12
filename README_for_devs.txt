
goGPS notes for developers

----------------------------------------------

To force git to ignore changes to specific files (that are already on the repository:
   git update-index --assume-unchanged <.\pathtofile\filename>
   
It is often important to ignore changes to the stations.crd file, thus it might be useful to run this command in the local repo:
   git update-index --assume-unchanged .\data\stations\stations.crd

----------------------------------------------

