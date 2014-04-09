#!/bin/sh

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

if [ "$(uname)" == "Darwin" ]; then
    matlab_exec=/Applications/MATLAB_R2010b.app/bin/matlab     #MacOS example
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    matlab_exec=/usr/local/MATLAB/R2013a/bin/matlab    #Linux example
elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
    # Do something under Windows NT platform
    echo "MINGW32_NT platform"
fi

X="goGPS"
echo ${X} > matlab_command.m
cat matlab_command.m
${matlab_exec} -nojvm -nodisplay -nosplash < matlab_command.m
rm matlab_command.m
