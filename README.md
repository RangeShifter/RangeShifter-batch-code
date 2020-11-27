# RangeShifter-batch-code <img src="RS_logo.png" align="right" height = 100/>
C++ code for the RangeShifter v2 batch mode application

[RangeShifter](https://rangeshifter.github.io/)
is an eco-evolutionary modelling platform that is becoming 
increasingly used worldwide for both theoretical and applied purposes
[(Bocedi et al. 2014)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12162).

RangeShifter is a spatially-explicit, individual-based simulation platform that 
allows modelling species’ population and range dynamics, such as expansion and shifting, and 
patch connectivity by linking complex local population dynamics and dispersal 
behaviour, while also taking into account inter-individual variability and 
evolutionary processes. RangeShifter is highly flexible in terms of the spatial 
resolution and extent, and regarding the complexity of the considered ecological 
processes. Due to its modular structure, the level of detail in genetics, demographic and 
dispersal processes can be easily adapted to different research questions and 
available data.

## Building RangeShifter
RangeShifter can be built in Windows or Linux. To build the programm you need to set the operating system in the file Version.h.

For Windows: set #define LINUX_CLUSTER 0

For Linux: set #define LINUX_CLUSTER 1

## Usage

Please refer to our [website](https://rangeshifter.github.io/) for more information about the RangeShifter platform and to the RangeShifter User Manual which can be found here: https://github.com/RangeShifter/RangeShifter-software-and-documentation


## References

 - Bocedi G, Palmer SCF, Pe’er G, Heikkinen RK, Matsinos YG, Watts K, Travis JMJ (2014). 
 *RangeShifter: A Platform for Modelling Spatial Eco-Evolutionary Dynamics and 
 Species’ Responses to Environmental Changes.* Methods in Ecology and Evolution 5: 388–96. 
