# fortran-xorshift-64-star
Fortran implementation of the xorshift64* random number generator.
Fortran does not support unsigned 64 bit integers and ifort does not suppert 128 bit integers. 
This code therefore implements a workaround by using four signed 32 bit integers to represent one unsigned 64 bit number. 

