BASEDIR=~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results

zcat $BASEDIR/ld/all_snps.beagle.ld.gz | sed -n '1,2983657p' > $BASEDIR/ld/output_chromosomes/OZ060699.beagle.ld.gz
zcat $BASEDIR/ld/all_snps.beagle.ld.gz | sed -n '2983658,5475162p' > $BASEDIR/ld/output_chromosomes/OZ060700.beagle.ld.gz
zcat $BASEDIR/ld/all_snps.beagle.ld.gz | sed -n '5475163,8225777p' > $BASEDIR/ld/output_chromosomes/OZ060701.beagle.ld.gz
zcat $BASEDIR/ld/all_snps.beagle.ld.gz | sed -n '8225778,10425822p' > $BASEDIR/ld/output_chromosomes/OZ060702.beagle.ld.gz
zcat $BASEDIR/ld/all_snps.beagle.ld.gz | sed -n '10425823,12076294p' > $BASEDIR/ld/output_chromosomes/OZ060703.beagle.ld.gz
zcat $BASEDIR/ld/all_snps.beagle.ld.gz | sed -n '12076295,14046564p' > $BASEDIR/ld/output_chromosomes/OZ060704.beagle.ld.gz
zcat $BASEDIR/ld/all_snps.beagle.ld.gz | sed -n '14046565,16137657p' > $BASEDIR/ld/output_chromosomes/OZ060705.beagle.ld.gz
zcat $BASEDIR/ld/all_snps.beagle.ld.gz | sed -n '16137658,18145831p' > $BASEDIR/ld/output_chromosomes/OZ060706.beagle.ld.gz
zcat $BASEDIR/ld/all_snps.beagle.ld.gz | sed -n '18145832,19585886p' > $BASEDIR/ld/output_chromosomes/OZ060707.beagle.ld.gz
zcat $BASEDIR/ld/all_snps.beagle.ld.gz | sed -n '19585887,21460703p' > $BASEDIR/ld/output_chromosomes/OZ060708.beagle.ld.gz
zcat $BASEDIR/ld/all_snps.beagle.ld.gz | sed -n '21460704,23201800p' > $BASEDIR/ld/output_chromosomes/OZ060709.beagle.ld.gz
zcat $BASEDIR/ld/all_snps.beagle.ld.gz | sed -n '23201801,24510522p' > $BASEDIR/ld/output_chromosomes/OZ060710.beagle.ld.gz
zcat $BASEDIR/ld/all_snps.beagle.ld.gz | sed -n '24510523,26190954p' > $BASEDIR/ld/output_chromosomes/OZ060711.beagle.ld.gz
zcat $BASEDIR/ld/all_snps.beagle.ld.gz | sed -n '26190955,27581866p' > $BASEDIR/ld/output_chromosomes/OZ060712.beagle.ld.gz
zcat $BASEDIR/ld/all_snps.beagle.ld.gz | sed -n '27581867,27582048p' > $BASEDIR/ld/output_chromosomes/OZ060713.beagle.ld.gz

