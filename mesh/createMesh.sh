# nel caso 3D, il parametro di numero di elementi, si intende 1/4 del numero di elemnti che vuoi lungo un quarto di circonferenza che taglia la sphera passando per il centro.
# vuoi 12 elementi? metti ns = 3 

if [ $1 == 2 ];
then
  echo “Generating circular mesh”

  ./parroth-opt -i generateQuarterCircle.i radius=$2 ns=$3
  rm -rf quarterCircle_0000.xda
  rm -rf quarterCircle_0001.xda
  rm -rf quarterCircle_0001_mesh.xda
  mv quarterCircle_0000_mesh.xda quarterCircle_${2}_${3}.xda

fi

if [ $1 == 3 ];
then
  echo “Generating spherical mesh”

  ./parroth-opt -i generateOctSphere.i radius=$2 ns=$3
  rm -rf octSphere_0000.xda
  rm -rf octSphere_0001.xda
  rm -rf octSphere_0001_mesh.xda
  mv octSphere_0000_mesh.xda octSphere_${2}_${3}.xda

fi
