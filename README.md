# EEG-Source-Reconstruction-NYC-Head-Model
"EEG Source Reconstruction Using the NYC head model" is a MATLAB software for EEG data analysis that has
been developed at Sharif University of Technology in collaboration with the University of Oslo.

The provided scripts can be used to calculate dipole powers using the eLORETA method, which is available in 
the FieldTrip toolbox [1]. Additionally, it utilizes the NYC head model [2]. 

For more information, please visit https://github.com/MahdiBabaei78/EEG-Source-Reconstruction-NYC-Head-Model. 

"EEG-Source-Reconstruction-NYC-Head-Model" is free software: you can redistribute it and/or modify it 
under the terms of the GNU General Public License as published by the Free Software Foundation, either 
version 3 of the License, or (at your option) any later version.

The functions in the "EEG-Source-Reconstruction-NYC-Head-Model" repository are copyrighted by their 
respective authors: 

- Mahdi Babaei, Sharif University of Technology, Tehran, Iran 			- mahdi78babaei@gmail.com (prefered)
                                                                    - babaei.mahdi@ee.sharif.edu

- Bjørn Erik Juel, University of Oslo, Oslo, Norway 				- b.e.juel@medisin.uio.no



The EEG-Source-Reconstruction-NYC-Head-Model software is a group of functions, which in turn can depend on 
other functions. The release of this repository includes functions from the FieldTrip toolbox that are covered
under their respective licenses. Unauthorized copying and distribution of functions that are not explicitly 
covered by the GPL is not allowed!

#############################

This repository contains new functions and some "modified" functions of an external toolbox, FieldTrip. It is 
used to perform these steps in order:

1. Reading the EEG data
2. Selecting the data of specific events
3. Preparing the New York City head model to be used in the FieldTrip toolbox
4. Calculating the dipoles' powers using the eLORETA function available in the FieldTrip toolbox


NOTES ABOUT EACH STEP:

1 & 2. In these steps, the original FieldTrip functions are used without any modifications.

3. The NYC head model is converted to a compatible version with the FieldTrip functions.

4. The ft_eloreta function in FieldTrip is modified. The modified version of this function is separated into 
two functions to increase the running speed. First, the spatial filters are created once. Then, the powers 
are calculated.



