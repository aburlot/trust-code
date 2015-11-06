#!/bin/bash
echo $0 $*
dirmodel=$1
nommodel=`basename $2 .prep`
if [  -d  $dirmodel/../prepro/trans_geom ]
then
   cd $dirmodel/../prepro/trans_geom
   dist=$dirmodel/..
   # Pb si etude Xprepro differente de Etude TRUST
   # Donc on fixe a ETUDE_TRUST
   [ "$ETUDE_TRUST" != "" ] && dist=$ETUDE_TRUST
   echo $dist
   echo "*****************************"
   echo "Copie des fichiers sous $dist"
   echo "*****************************"
   cp -f ${dirmodel}/$2 $dist
   is2D=0 && is2DAxi=0
   for file in `ls *_2D*.geom 2>/dev/null`
   do
      is2D=1 && [ ${file%Axi.geom} != $file ] && is2DAxi=1
      echo on recupere $file
      cp -f $file $dist/${nommodel}_$file
      geoms=$geoms" "${nommodel}_$file
   done
   if [ $is2D = 0 ]
   then
      echo on recupere $file
      for file in `ls *.geom 2>/dev/null`
      do
         echo on recupere $file
         cp -f $file $dist/${nommodel}_$file
         geoms=$geoms" "${nommodel}_$file
      done
   fi
   # On recupere aussi les .geos
   for file in `ls *.geos 2>/dev/null`
   do
      echo on recupere $file
      cp -f $file $dist/${nommodel}_$file
   done   
   # on regarde si Def_Suz existe si oui il faut copier Suz_Def
   for fileD in `ls Def_Suz_def* 2>/dev/null`
   do
      file=${fileD#Def_}
      echo on recupere $file
      cp -f $file $dist/${nommodel}_$file
      Suzs=$Suzs" "${nommodel}_$file 2>/dev/null
   done
   for file in `ls *.med 2>/dev/null`
   do
      echo on recupere $file
      cp -f $file $dist/$file
   done

   cd ../..
   NOMCAS=$nommodel
   file=$dist/$NOMCAS.mesh
   echo $ECHO_OPTS "# ############################################ #
# File to save and include into your data file #
# either by copy-paste, or by the instruction: #
# lire_fichier name_of_this_file ;             #
# It defines dimension, domains and read mesh  #
# generated by Xprepro                         #
# ############################################ #
" > $file
   # 2D ou 3D ?
   if [ $is2D = 1 ]
   then
      echo $ECHO_OPTS "dimension 2" >> $file
   else
      echo $ECHO_OPTS "dimension 3" >> $file
   fi
   # Axi 2D ?
   [ $is2DAxi = 1 ] && echo $ECHO_OPTS "Bidim_axi" >> $file
   
   [ ${#geoms} = 0 ] && $TRUST_ROOT/bin/IHM/Erreur.ihm "Files containing meshes not generated !"
   n=0
   for geom in $geoms
   do
      let n=$n+1
      dom=dom_pb$n && [ $is2D = 1 ] && dom=dom_$n
      echo $ECHO_OPTS "export domaine $dom" >> $file
   done
   echo $ECHO_OPTS "# DEBUT MAILLAGE #" >> $file
   n=0
   for geom in $geoms
   do
      let n=$n+1
      dom=dom_pb$n && [ $is2D = 1 ] && dom=dom_$n
      echo $ECHO_OPTS "lire_fichier $dom $geom" >> $file 
   done
   n=0
   # Ajout d'un test pour prevenir que les Sous Zones
   # restent 3D et ne sont pas coupees par 3Dto2D (a developper?)
   if [ $is2D = 1 ]
   then
      echo $ECHO_OPTS "# FIN MAILLAGE #" >> $file
      echo "# Remarque:Les Sous Zones restent 3D et ne sont pas encore coupees par 3Dto2D. #\n" >> $file
   else
      # Les sous zones
      cat prepro/trans_geom/Def_Suz_def* >> $file
      echo $ECHO_OPTS "# FIN MAILLAGE #\n" >> $file
   fi
   echo $ECHO_OPTS "# DEBUT DECOUPAGE" >> $file
   n=0
   for geom in $geoms
   do
      let n=$n+1
      dom=dom_pb$n && [ $is2D = 1 ] && dom=dom_$n
      echo $ECHO_OPTS "Decouper $dom
{
   Partitionneur metis { Nb_parts 2 }
   Larg_joint 1
   Nom_Zones DOM$n
}" >> $file
   done
   echo $ECHO_OPTS "Fin
FIN DECOUPAGE #" >> $file

   echo $ECHO_OPTS "\n# DEBUT LECTURE" >> $file
   n=0
   for geom in $geoms
   do
      let n=$n+1
      dom=dom_pb$n && [ $is2D = 1 ] && dom=dom_$n
      echo $ECHO_OPTS "Scatter DOM"$n".Zones "$dom >> $file
   done
   n=0
   # Ajout les Sous Zones uniquement en 3D
   if [ $is2D = 0 ]
   then
      # La ligne suivante n'est pas suffisante pour le //
      cat prepro/trans_geom/Def_Suz_def* >> $file
   fi
   echo $ECHO_OPTS "1,$ s? Suz_def_pb? "$nommodel"_Suz_def_pb?g\nw" | ed $file 1>/dev/null 2>&1
   echo $ECHO_OPTS "FIN LECTURE #" >> $file
   echo $ECHO_OPTS "\n#" >> $file
   $TRUST_Awk '/Lecture du fichier/ {print "---------------------------\nBOUNDARY NAMES OF "$NF":\n---------------------------"}' prepro/trans_geom/TRUST.log >> $file 
   if (( `grep -i "bord conserve" prepro/trans_geom/TRUST.log | wc -l` == 0 ))
   then
      $TRUST_Awk '/commence a la face/ && !/coupe_2D/ && !/coupe_Axi2D/ {print $3}' prepro/trans_geom/TRUST.log | sort -u >> $file 
   else
      $TRUST_Awk '/bord conserve/ && !/coupe_2D/ && !/coupe_Axi2D/ {print $3}' prepro/trans_geom/TRUST.log | sort -u >> $file 
   fi
   echo $ECHO_OPTS "#" >> $file   
   $EDITEUR $file
   
   cd $dirmodel
   `dirname $0`/nettoie ${dirmodel}
else
   echo prepro must be run before!!!
fi
