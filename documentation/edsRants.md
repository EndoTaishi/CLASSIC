Ed's Rants
============
**To be assimilated into the right spot in the documentation**

# Array Allocation Inside a Function Across the GNU/Intel Compilers
A difference between the GNU and the Intel compilers is that the GNU compiler allows for a function to return an allocatable array. So, if one were to allocate and return an array, the values would work fine after the function call.

On the Intel compiler, however, if the return array is allocated inside the function and returned, after the function call the array would still be unallocated.

If, however, one were to allocate in advance the array to which the function would return the array, then the values are stored properly.

# Precompile Directives
