#!/bin/bash

#----------------------------------------------------------------------------------------------
#                           goGPS v0.4.2 beta
#
# Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
#----------------------------------------------------------------------------------------------
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#----------------------------------------------------------------------------------------------

ini_path=${1}

# rover RINEX filename example: UBLX0510.14o
# master RINEX filename example: mila051_15s.14O
# navigation RINEX filename example: brdc0510.14n

MARKER_R=UBLX # rover marker name
MARKER_M=mila # master marker name

SESSION_R=0    # rover session id
SESSION_M=_15s # master session id

EXT_R=o # rover observation id in file extension
EXT_M=O # master observation id in file extension

ID_N=brdc # navigation id in filename
SESSION_N=0  # navigation session id
EXT_N=n   # navigation id in file extension

YEAR=14

DOY_START=29
DOY_END=50

OUT_FOLDER="../data/out/batch"

for ((DOY=DOY_START;DOY<=DOY_END;DOY++))
do
    if (( "$DOY" < 100 )); then
        DOYF=0${DOY}
        DOYF_START=0${DOY_START}
        DOYF_END=0${DOY_END}
    else
        DOYF=${DOY}
        DOYF_START=${DOY_START}
        DOYF_END=${DOY_END}
    fi

    #./CRX2RNX ${MARKER_M}${DOYF}${SESSION_M}.${YEAR}d
    #./CRX2RNX ${MARKER_R}${DOYF}${SESSION_R}.${YEAR}d
    
    SED1="s/${MARKER_M}[0-9][0-9][0-9]${SESSION_M}.${YEAR}${EXT_M}/${MARKER_M}${DOYF}${SESSION_M}.${YEAR}${EXT_M}/g"
    SED2="s/${MARKER_R}[0-9][0-9][0-9]${SESSION_R}.${YEAR}${EXT_R}/${MARKER_R}${DOYF}${SESSION_R}.${YEAR}${EXT_R}/g"
    SED3="s/${ID_N}[0-9][0-9][0-9]${SESSION_N}.${YEAR}${EXT_N}/${ID_N}${DOYF}${SESSION_N}.${YEAR}${EXT_N}/g"
    SED4="s/${MARKER_M}_${MARKER_R}_${YEAR}[0-9][0-9][0-9]/${MARKER_M}_${MARKER_R}_${YEAR}${DOYF}/g"
    SED5="s/mode_user=[0-1]/mode_user=0/g"
     
    if [[ "$OSTYPE" == "linux-gnu" ]]; then
        #Linux
        sed -i ${SED1} $ini_path
        sed -i ${SED2} $ini_path
        sed -i ${SED3} $ini_path
        sed -i ${SED4} global_settings.m
        sed -i ${SED5} goGPS.m
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        #MacOS
        sed -i "" ${SED1} $ini_path
        sed -i "" ${SED2} $ini_path
        sed -i "" ${SED3} $ini_path
        sed -i "" ${SED4} global_settings.m
        sed -i "" ${SED5} goGPS.m
    else
        echo "OS detection failed."
        exit
    fi
    
    sh goGPS_cmd.sh > ${OUT_FOLDER}/${MARKER_M}_${MARKER_R}_${YEAR}${DOYF}_batch_stdout.txt
    
    tail -1 ${OUT_FOLDER}/${MARKER_M}_${MARKER_R}_${YEAR}${DOYF}_position.txt | awk '{print $8,$9,$10}' >> ${OUT_FOLDER}/${MARKER_M}_${MARKER_R}_${YEAR}${DOYF_START}${DOYF_END}_extraction.txt
    
    SED6="s/mode_user=[0-1]/mode_user=1/g"
    
    if [[ "$OSTYPE" == "linux-gnu" ]]; then
        #Linux
        sed -i ${SED6} goGPS.m
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        #MacOS
        sed -i "" ${SED6} goGPS.m
    fi
done
