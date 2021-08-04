#!/bin/bash

for ref in {17,19,38}
do
  # put /opt before eagle
  map_file=/opt/Eagle_v2.4.1/tables/genetic_map_hg${ref}_withX.txt.gz
  for chrom in {1..23}
  do
    echo -e 'chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)' > /opt/genetic_maps_eagle/hg${ref}/genetic_map_hg${ref}_chr${chrom}_withX.txt
    cat $map_file | gunzip | awk '{if ( ($1==CHROM) ) print $0}' CHROM=${chrom} >> /opt/genetic_maps_eagle/hg${ref}/genetic_map_hg${ref}_chr${chrom}_withX.txt
    gzip /opt/genetic_maps_eagle/hg${ref}/genetic_map_hg${ref}_chr${chrom}_withX.txt

    # echo -e 'pos\tchr\tcM' > /opt/genetic_maps_shapeit/hg${ref}/genetic_map_hg${ref}_chr${chrom}_withX.txt
    cat /opt/genetic_maps_eagle/hg${ref}/genetic_map_hg${ref}_chr${chrom}_withX.txt.gz | gunzip | awk '{print $2,"\t",$1,"\t",$4}' > /opt/genetic_maps_shapeit/hg${ref}/genetic_map_hg${ref}_chr${chrom}_withX.txt
    sed 's/position/~~/g; s/Genetic_Map(cM)/cM/g; s/~~/pos/g' /opt/genetic_maps_shapeit/hg${ref}/genetic_map_hg${ref}_chr${chrom}_withX.txt | gzip > /opt/genetic_maps_shapeit/hg${ref}/genetic_map_hg${ref}_chr${chrom}_withX.txt.gz
    rm /opt/genetic_maps_shapeit/hg${ref}/genetic_map_hg${ref}_chr${chrom}_withX.txt
  done
done

