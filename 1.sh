#!/bin/bash
# This little script is written by Quantao for the purpose of preparing input files
#for NAMD based Free Energy Pertubation calculation. 
#Please download the zip file from "Feprepare web server" to the same folder as this script is, also please make sure the 
#parameter_patch2.txt and toppar_modified should be in the folder as per the structure 
#below.

                                  ##### your folder to run the calculation should look like this #####
                                                                #                      
                                                                #
                                                                #
    ###################################################################################################################################
    #                  #              #                  #             #                   #              #             #             #
    #                  #              #                  #             #                   #              #             #             #
#files.zip  toppar_modified.zip  parameter_patch2.txt  1.sh         lambda.txt          README.md       LICENSE       2.1.pbs        2.2.sh
     
unzip files.zip # tar -xvf files.tar.gz
cd files # in case the folder name is random, this cd to the first one.
cp ../toppar_modified.zip ./
unzip toppar_modified.zip
cd complex
tar -xzvf *.tar.gz  # this is the conf files
mkdir complex-input-files
mkdir complex-output-files
mv *.namd complex-input-files/
cd complex-input-files/
############complex NVT script modification ##################################
head -n 13 nvt_equil.namd > output1
cat ../../../parameter_patch2.txt >> output1
tail -n +14 nvt_equil.namd > output2
cat output2 >> output1
# change the numSteps from 500000 to 50000
# change the numMinSteps from 50000 to 5000
mv output1 nvt_equil_modi.namd
# you could manually define the output path to ../complex-output-files
############complex NPT script modification #####################################
head -n 13 npt_equil.namd > output3
cat ../../../parameter_patch2.txt >> output3
tail -n +14 npt_equil.namd > output4
cat output4 >> output3
#mannuall change the numSteps from 500000 to 50000.
mv output3 npt_equil_modi.namd
# you could manually define the output path to ../complex-output-files
###############complex Production script modification #################################
tail -n 3 md_forward_1.namd > tmp2.txt
echo "alchEquilSteps          500
      set numSteps            50000" > tmp2.txt  # overwrite the original Lambda0 Lambda2 format
cat ../../../lambda.txt >> tmp2.txt  # you could adjust this file to have different windows numbers.
cat tmp2.txt >> md_forward_1.namd
#To uncomment lines 137 through 141 to close the one by one simulation.
sed -i '137,141 s/^/#/' md_forward_1.namd    #double check the line numbers.
head -n 123 md_forward_2.namd > outputA
tail -n +124 md_forward_2.namd > outputB
echo "source     ../../fep.tcl" >> outputA
cat outputB >> outputA
mv outputA md_forward_modi.namd
#################### To modify the alchFile section #############################                                                
head -n 13 md_forward_2.namd > output5
cat ../../../parameter_patch2.txt >> output5
tail -n +14 md_forward_2.namd > output6
sed -i '/^alchDecouple/s/off/yes/' output6           # turn on the "alchDecouple"
cat output6 >> output5
mv output5 md_forward_modi.namd


############################################################################
cd ../../ # prepare to cd to solvent folder. ##SWITCH TO SOLVENT
###########################################################################
cd solvent
tar -xzvf *.tar.gz # extract all the conf files
mkdir solvent-input-files
mkidr solvent-output-files
mv *.namd solvent-input-files/
cd solvent-input-files/
############solvent NVT script modification ##################################
head -n 13 nvt_equil.namd > output1
cat ../../../parameter_patch2.txt >> output1
tail -n +14 nvt_equil.namd > output2
cat output2 >> output1
mv output1 nvt_equil_modi.namd
#change margin from 1 to 30.0
# mannually change numSteps from 500000 to 50000
# change the numMinSteps from 50000 to 5000 
# you could manually define the output path to ../solvent-output-files
############solvent NPT script modification #####################################
head -n 13 npt_equil.namd > output3
cat ../../../parameter_patch2.txt >> output3
tail -n +14 npt_equil.namd > output4
cat output4 >> output3
#manually change the numSteps from 500000 to 50000.
mv output3 npt_equil_modi.namd
# you could manually define the output path to ../complex-output-files
###############solvent Production script modification #################################
tail -n 3 md_forward_2.namd > tmp2.txt
echo "alchEquilSteps          500
      set numSteps            50000" > tmp2.txt 
cat ../../../lambda.txt >> tmp2.txt  
cat tmp2.txt >> md_forward_2.namd
#To uncomment lines 139 through 141 to close the one by one simulation.
sed -i '139,141 s/^/#/' md_forward_1.namd        
head -n 123 md_forward_1.namd > outputA
tail -n +124 md_forward_1.namd > outputB
echo "source     ../fep.tcl" >> outputA
cat outputB >> outputA
mv outputA md_forward_1.namd
#################### To modify the alchFile section #########################                                                
head -n 13 md_forward_1.namd > output5
cat parameter_patch.txt >> output5
tail -n +14 md_forward_1.namd > output6       
sed -i '/^alchDecouple/s/off/yes/' output6    
cat output6 >> output5
mv output5 md_forward_modi.namd
grep X ../ionized_solvent.fep > X.txt # validate the fep file is is good.
cat X.txt
echo "If you don't see the X labeled atoms, please go back to check your input"
grep X ../../complex/ionized_complex.fep > X1.txt # validate the fep file is is good.
cat X1.txt
echo "If you don't see the X labeled atoms, please go back to check your input"
