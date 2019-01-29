This is a prototype project of the joint method for marker-free alignment. The third-part modules will be removed in beta version.

The files is build in fedora 25. To use the exe, please add the bin to the $PATH in the computer.

To run:

jointalign -i Energy.st -a Energy.rawtlt -o fin.xf -n Energy.tlt -t 450 -p 50 

newstack -input Energy.st -xf fin.xf -output ali.mrc #IMOD
or 
mrcstack -i Energy.st -x fin.xf -o ali.mrc #AuTom

mpirun -n 15 volrec_sglm -i ali.mrc -o refine.rec -a Energy.tlt -x xtiltangle.txt  -g -1,2.5,-10,500 -m SART,20,0.2

For the test data, please visit: https://drive.google.com/drive/folders/1ybRGbfgxYs12Gr28BNjQ4my73haWjk9v?usp=sharing
