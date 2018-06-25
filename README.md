# DWD_project
Use double white dwarfs as tracers to model the structure of Milky Way.

SeBa_r105_aa_wdwd.data and SeBa_r105_ag_wdwd.data are input files representing two binary population syntesis models.
DWD_pop_BP_ugriz.C is the code produsing the distribution of DWDs in the Milky Way convolved with the SFH from Boissier & Prantzos 99.

To compile the code: 

rm myutil.o

rm DWD_pop_BP_ugriz

g++ -O2 -W -g -I/usr/include -c myutil.C

g++ -O2 -W -g -I/usr/include -o DWD_pop_BP_ugriz DWD_pop_BP_ugriz.C myutil.o


To run the code:

./DWD_pop_BP_ugriz -f input_filename >& output_filename
