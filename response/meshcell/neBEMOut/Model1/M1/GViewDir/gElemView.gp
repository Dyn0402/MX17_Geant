set title "neBEM elements in gnuplot VIEWER"
#set pm3d
#set style data pm3d
#set palette model CMY
set hidden3d
set nokey
set xlabel "X"
set ylabel "Y"
set zlabel "Z"
set view 70, 335, 1, 1

splot \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim1.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim2.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim3.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim4.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim5.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim6.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim7.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim8.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim9.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim10.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim11.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim12.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim13.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim14.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim15.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim16.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim17.out' w p, \
 'neBEMOut/Model1/M1/GViewDir/gpElemOnPrim18.out' w p

pause-1