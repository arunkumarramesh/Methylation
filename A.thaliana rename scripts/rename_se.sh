#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=100000
#SBATCH -J rename_se
#SBATCH -o rename_se.out
#SBATCH -e rename_se.err

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data/

cat SRR797874_trim.fq.gz SRR3300169_trim.fq.gz > 7125_trim.fq.gz
cat SRR771664_trim.fq.gz SRR3300083_trim.fq.gz > 6897_trim.fq.gz
cat SRR771666_trim.fq.gz SRR3300121_trim.fq.gz > 6987_trim.fq.gz
cat SRR771668_trim.fq.gz SRR771669_trim.fq.gz > 6990_trim.fq.gz
cat SRR771671_trim.fq.gz SRR3300072_trim.fq.gz SRR3300073_trim.fq.gz SRR3300074_trim.fq.gz > 6680_trim.fq.gz
cat SRR771672_trim.fq.gz SRR3300131_trim.fq.gz SRR3300132_trim.fq.gz SRR3300133_trim.fq.gz > 6994_trim.fq.gz
cat SRR771673_trim.fq.gz SRR771674_trim.fq.gz > 9759_trim.fq.gz
cat SRR771679_trim.fq.gz SRR3301099_trim.fq.gz SRR3301100_trim.fq.gz SRR3301101_trim.fq.gz > 9761_trim.fq.gz
cat SRR771682_trim.fq.gz SRR3300137_trim.fq.gz SRR3300138_trim.fq.gz > 7026_trim.fq.gz
cat SRR771683_trim.fq.gz SRR771684_trim.fq.gz > 6903_trim.fq.gz
cat SRR771688_trim.fq.gz SRR771689_trim.fq.gz > 7031_trim.fq.gz
cat SRR771690_trim.fq.gz SRR771691_trim.fq.gz > 7036_trim.fq.gz
cat SRR771696_trim.fq.gz SRR3299928_trim.fq.gz SRR3299929_trim.fq.gz > 5726_trim.fq.gz
cat SRR771700_trim.fq.gz SRR771701_trim.fq.gz > 6911_trim.fq.gz
cat SRR771704_trim.fq.gz SRR3300155_trim.fq.gz > 7096_trim.fq.gz
cat SRR771706_trim.fq.gz SRR3300092_trim.fq.gz > 6915_trim.fq.gz
cat SRR771713_trim.fq.gz SRR3300170_trim.fq.gz > 7130_trim.fq.gz
cat SRR771714_trim.fq.gz SRR3300171_trim.fq.gz > 7138_trim.fq.gz
cat SRR771721_trim.fq.gz SRR3300095_trim.fq.gz > 6921_trim.fq.gz
cat SRR771722_trim.fq.gz SRR3299880_trim.fq.gz > 430_trim.fq.gz
cat SRR771729_trim.fq.gz SRR771730_trim.fq.gz > 7162_trim.fq.gz
cat SRR771734_trim.fq.gz SRR3300187_trim.fq.gz > 7181_trim.fq.gz
cat SRR771752_trim.fq.gz SRR3300197_trim.fq.gz > 7236_trim.fq.gz
cat SRR771756_trim.fq.gz SRR3300204_trim.fq.gz > 7252_trim.fq.gz
cat SRR771768_trim.fq.gz SRR771769_trim.fq.gz SRR3300379_trim.fq.gz > 8354_trim.fq.gz
cat SRR771770_trim.fq.gz SRR3300218_trim.fq.gz > 7298_trim.fq.gz
cat SRR771775_trim.fq.gz SRR3300219_trim.fq.gz SRR3300220_trim.fq.gz > 7305_trim.fq.gz
cat SRR771776_trim.fq.gz SRR771777_trim.fq.gz > 6951_trim.fq.gz
cat SRR771781_trim.fq.gz SRR771782_trim.fq.gz > 7314_trim.fq.gz
cat SRR771789_trim.fq.gz SRR3300306_trim.fq.gz SRR3300307_trim.fq.gz > 7514_trim.fq.gz
cat SRR771792_trim.fq.gz SRR771793_trim.fq.gz > 7332_trim.fq.gz
cat SRR771796_trim.fq.gz SRR3300247_trim.fq.gz SRR3300248_trim.fq.gz > 7337_trim.fq.gz
cat SRR771801_trim.fq.gz SRR3300249_trim.fq.gz SRR3300250_trim.fq.gz > 7342_trim.fq.gz
cat SRR771809_trim.fq.gz SRR771810_trim.fq.gz > 7378_trim.fq.gz
cat SRR4295280_trim.fq.gz SRR4295281_trim.fq.gz > 108_trim.fq.gz
cat SRR4295282_trim.fq.gz SRR4295283_trim.fq.gz > 428_trim.fq.gz
cat SRR4295284_trim.fq.gz SRR4295285_trim.fq.gz > 801_trim.fq.gz
cat SRR4295288_trim.fq.gz SRR4295289_trim.fq.gz > 1872_trim.fq.gz
cat SRR4295296_trim.fq.gz SRR4295297_trim.fq.gz > 5165_trim.fq.gz
cat SRR4295303_trim.fq.gz SRR4295304_trim.fq.gz > 5741_trim.fq.gz
cat SRR4295305_trim.fq.gz SRR4295306_trim.fq.gz > 5779_trim.fq.gz
cat SRR4295313_trim.fq.gz SRR4295314_trim.fq.gz SRR4295315_trim.fq.gz SRR4295316_trim.fq.gz SRR4295317_trim.fq.gz SRR4295318_trim.fq.gz > 6217_trim.fq.gz
cat SRR4295319_trim.fq.gz SRR4295320_trim.fq.gz > 6296_trim.fq.gz
cat SRR4295339_trim.fq.gz SRR4295340_trim.fq.gz > 7081_trim.fq.gz
cat SRR4295342_trim.fq.gz SRR4295343_trim.fq.gz SRR4295344_trim.fq.gz SRR4295345_trim.fq.gz SRR4295346_trim.fq.gz SRR4295347_trim.fq.gz SRR4295348_trim.fq.gz > 7213_trim.fq.gz
cat SRR4295352_trim.fq.gz SRR4295353_trim.fq.gz SRR4295354_trim.fq.gz SRR4295355_trim.fq.gz SRR4295356_trim.fq.gz SRR4295357_trim.fq.gz SRR4295358_trim.fq.gz SRR4295359_trim.fq.gz SRR4295360_trim.fq.gz SRR4295361_trim.fq.gz SRR4295362_trim.fq.gz SRR4295363_trim.fq.gz SRR4295364_trim.fq.gz SRR4295365_trim.fq.gz SRR4295366_trim.fq.gz SRR4295367_trim.fq.gz SRR4295368_trim.fq.gz SRR4295369_trim.fq.gz SRR4295370_trim.fq.gz SRR4295371_trim.fq.gz > 7417_trim.fq.gz
cat SRR4295374_trim.fq.gz SRR4295375_trim.fq.gz > 7717_trim.fq.gz
cat SRR4295384_trim.fq.gz SRR4295385_trim.fq.gz > 9067_trim.fq.gz
cat SRR4295386_trim.fq.gz SRR4295387_trim.fq.gz SRR4295388_trim.fq.gz SRR4295389_trim.fq.gz SRR4295390_trim.fq.gz SRR4295391_trim.fq.gz SRR4295392_trim.fq.gz SRR4295393_trim.fq.gz > 9078_trim.fq.gz
cat SRR4295397_trim.fq.gz SRR4295398_trim.fq.gz > 9102_trim.fq.gz
cat SRR4295399_trim.fq.gz SRR4295400_trim.fq.gz > 9121_trim.fq.gz
cat SRR4295401_trim.fq.gz SRR4295402_trim.fq.gz > 9128_trim.fq.gz
cat SRR4295404_trim.fq.gz SRR4295405_trim.fq.gz > 9133_trim.fq.gz
cat SRR4295406_trim.fq.gz SRR4295407_trim.fq.gz SRR4295408_trim.fq.gz SRR4295409_trim.fq.gz SRR4295410_trim.fq.gz SRR4295411_trim.fq.gz SRR4295412_trim.fq.gz SRR4295413_trim.fq.gz SRR4295414_trim.fq.gz SRR4295415_trim.fq.gz SRR4295416_trim.fq.gz SRR4295417_trim.fq.gz SRR4295418_trim.fq.gz SRR4295419_trim.fq.gz > 9622_trim.fq.gz
cat SRR4295422_trim.fq.gz SRR4295423_trim.fq.gz SRR4295424_trim.fq.gz SRR4295425_trim.fq.gz SRR4295426_trim.fq.gz SRR4295427_trim.fq.gz SRR4295428_trim.fq.gz SRR4295429_trim.fq.gz SRR4295430_trim.fq.gz SRR4295431_trim.fq.gz > 9642_trim.fq.gz
cat SRR4295433_trim.fq.gz SRR4295434_trim.fq.gz > 9665_trim.fq.gz
cat SRR4295435_trim.fq.gz SRR4295436_trim.fq.gz > 9689_trim.fq.gz
cat SRR4295438_trim.fq.gz SRR4295439_trim.fq.gz SRR4295440_trim.fq.gz SRR4295441_trim.fq.gz SRR4295442_trim.fq.gz SRR4295443_trim.fq.gz SRR4295444_trim.fq.gz SRR4295445_trim.fq.gz SRR4295446_trim.fq.gz SRR4295447_trim.fq.gz > 9725_trim.fq.gz
cat SRR4295450_trim.fq.gz SRR4295451_trim.fq.gz SRR4295452_trim.fq.gz SRR4295453_trim.fq.gz SRR4295454_trim.fq.gz SRR4295455_trim.fq.gz > 9742_trim.fq.gz
cat SRR4295471_trim.fq.gz SRR4295472_trim.fq.gz > 9988_trim.fq.gz
cat SRR4295473_trim.fq.gz SRR4295474_trim.fq.gz > 9993_trim.fq.gz
cat SRR4295475_trim.fq.gz SRR4295476_trim.fq.gz > 9995_trim.fq.gz
cat SRR4295477_trim.fq.gz SRR4295478_trim.fq.gz > 10000_trim.fq.gz
cat SRR4295481_trim.fq.gz SRR4295482_trim.fq.gz SRR4295483_trim.fq.gz SRR4295484_trim.fq.gz SRR4295485_trim.fq.gz SRR4295486_trim.fq.gz SRR4295487_trim.fq.gz SRR4295488_trim.fq.gz SRR4295489_trim.fq.gz SRR4295490_trim.fq.gz SRR4295491_trim.fq.gz SRR4295492_trim.fq.gz > 10005_trim.fq.gz
cat SRR4295494_trim.fq.gz SRR4295495_trim.fq.gz SRR4295496_trim.fq.gz SRR4295497_trim.fq.gz SRR4295498_trim.fq.gz SRR4295499_trim.fq.gz > 10019_trim.fq.gz
cat SRR3299775_trim.fq.gz SRR3299776_trim.fq.gz SRR3299777_trim.fq.gz > 10001_trim.fq.gz
cat SRR3299778_trim.fq.gz SRR3299779_trim.fq.gz SRR3299780_trim.fq.gz SRR3299781_trim.fq.gz > 10004_trim.fq.gz
cat SRR3299782_trim.fq.gz SRR3299783_trim.fq.gz SRR3299784_trim.fq.gz > 10006_trim.fq.gz
cat SRR3299786_trim.fq.gz SRR3299787_trim.fq.gz SRR3299788_trim.fq.gz > 10010_trim.fq.gz
cat SRR3299791_trim.fq.gz SRR3299792_trim.fq.gz > 10013_trim.fq.gz
cat SRR3299793_trim.fq.gz SRR3299794_trim.fq.gz > 10014_trim.fq.gz
cat SRR3299796_trim.fq.gz SRR3299797_trim.fq.gz > 10015_trim.fq.gz
cat SRR3299798_trim.fq.gz SRR3299799_trim.fq.gz SRR3299800_trim.fq.gz SRR3299801_trim.fq.gz > 10016_trim.fq.gz
cat SRR3299802_trim.fq.gz SRR3299803_trim.fq.gz SRR3299804_trim.fq.gz > 10017_trim.fq.gz
cat SRR3299805_trim.fq.gz SRR3299806_trim.fq.gz SRR3299807_trim.fq.gz > 10018_trim.fq.gz
cat SRR3299808_trim.fq.gz SRR3299809_trim.fq.gz SRR3299810_trim.fq.gz > 10020_trim.fq.gz
cat SRR3299811_trim.fq.gz SRR3299812_trim.fq.gz SRR3299813_trim.fq.gz > 10021_trim.fq.gz
cat SRR3299814_trim.fq.gz SRR3299815_trim.fq.gz SRR3299816_trim.fq.gz > 10022_trim.fq.gz
cat SRR3299817_trim.fq.gz SRR3299818_trim.fq.gz SRR3299819_trim.fq.gz > 10023_trim.fq.gz
cat SRR3299820_trim.fq.gz SRR3299821_trim.fq.gz SRR3299822_trim.fq.gz > 10024_trim.fq.gz
cat SRR3299823_trim.fq.gz SRR3299824_trim.fq.gz SRR3299825_trim.fq.gz > 10026_trim.fq.gz
cat SRR3299826_trim.fq.gz SRR3299827_trim.fq.gz SRR3299828_trim.fq.gz > 1002_trim.fq.gz
cat SRR3299829_trim.fq.gz SRR3299830_trim.fq.gz SRR3299831_trim.fq.gz > 1066_trim.fq.gz
cat SRR3299832_trim.fq.gz SRR3299833_trim.fq.gz SRR3299834_trim.fq.gz > 1158_trim.fq.gz
cat SRR3299835_trim.fq.gz SRR3299836_trim.fq.gz > 1435_trim.fq.gz
cat SRR3299840_trim.fq.gz SRR3299841_trim.fq.gz > 1829_trim.fq.gz
cat SRR3299843_trim.fq.gz SRR3299844_trim.fq.gz SRR3299845_trim.fq.gz > 1890_trim.fq.gz
cat SRR3299847_trim.fq.gz SRR3299848_trim.fq.gz > 1925_trim.fq.gz
cat SRR3299849_trim.fq.gz SRR3299850_trim.fq.gz > 1954_trim.fq.gz
cat SRR3299852_trim.fq.gz SRR3299853_trim.fq.gz SRR3299854_trim.fq.gz > 2016_trim.fq.gz
cat SRR3299855_trim.fq.gz SRR3299856_trim.fq.gz > 2202_trim.fq.gz
cat SRR3299857_trim.fq.gz SRR3299858_trim.fq.gz SRR3299859_trim.fq.gz > 2278_trim.fq.gz
cat SRR3299860_trim.fq.gz SRR3299861_trim.fq.gz SRR3299862_trim.fq.gz > 2317_trim.fq.gz
cat SRR3299863_trim.fq.gz SRR3299864_trim.fq.gz > 265_trim.fq.gz
cat SRR3299868_trim.fq.gz SRR3299869_trim.fq.gz > 351_trim.fq.gz
cat SRR3299870_trim.fq.gz SRR3299871_trim.fq.gz > 403_trim.fq.gz
cat SRR3299873_trim.fq.gz SRR3299874_trim.fq.gz SRR3299875_trim.fq.gz > 410_trim.fq.gz
cat SRR3299877_trim.fq.gz SRR3299878_trim.fq.gz SRR3299879_trim.fq.gz > 424_trim.fq.gz
cat SRR3299881_trim.fq.gz SRR3299882_trim.fq.gz > 4807_trim.fq.gz
cat SRR3299883_trim.fq.gz SRR3299884_trim.fq.gz > 4884_trim.fq.gz
cat SRR3299885_trim.fq.gz SRR3299886_trim.fq.gz > 4931_trim.fq.gz
cat SRR3299890_trim.fq.gz SRR3299891_trim.fq.gz > 5023_trim.fq.gz
cat SRR3299892_trim.fq.gz SRR3299893_trim.fq.gz > 506_trim.fq.gz
cat SRR3299894_trim.fq.gz SRR3299895_trim.fq.gz > 5104_trim.fq.gz
cat SRR3299897_trim.fq.gz SRR3299898_trim.fq.gz SRR3299899_trim.fq.gz > 5151_trim.fq.gz
cat SRR3299900_trim.fq.gz SRR3299901_trim.fq.gz SRR3299902_trim.fq.gz SRR3299903_trim.fq.gz SRR3299904_trim.fq.gz > 5210_trim.fq.gz
cat SRR3299905_trim.fq.gz SRR3299906_trim.fq.gz SRR3299907_trim.fq.gz SRR3299908_trim.fq.gz SRR3299909_trim.fq.gz > 5236_trim.fq.gz
cat SRR3299910_trim.fq.gz SRR3299911_trim.fq.gz > 5249_trim.fq.gz
cat SRR3299912_trim.fq.gz SRR3299913_trim.fq.gz SRR3299914_trim.fq.gz > 5276_trim.fq.gz
cat SRR3299915_trim.fq.gz SRR3299916_trim.fq.gz > 5349_trim.fq.gz
cat SRR3299917_trim.fq.gz SRR3299918_trim.fq.gz > 5353_trim.fq.gz
cat SRR3299919_trim.fq.gz SRR3299920_trim.fq.gz SRR3299921_trim.fq.gz > 5577_trim.fq.gz
cat SRR3299922_trim.fq.gz SRR3299923_trim.fq.gz SRR3299924_trim.fq.gz SRR3299925_trim.fq.gz > 5644_trim.fq.gz
cat SRR3299926_trim.fq.gz SRR3299927_trim.fq.gz > 5717_trim.fq.gz
cat SRR3299932_trim.fq.gz SRR3299933_trim.fq.gz > 5772_trim.fq.gz
cat SRR3299938_trim.fq.gz SRR3299939_trim.fq.gz SRR3299940_trim.fq.gz SRR3299941_trim.fq.gz SRR3299942_trim.fq.gz SRR3299943_trim.fq.gz > 5837_trim.fq.gz
cat SRR3299944_trim.fq.gz SRR3299945_trim.fq.gz SRR3299946_trim.fq.gz > 5860_trim.fq.gz
cat SRR3299948_trim.fq.gz SRR3299949_trim.fq.gz SRR3299950_trim.fq.gz > 5890_trim.fq.gz
cat SRR3299952_trim.fq.gz SRR3299953_trim.fq.gz SRR3299954_trim.fq.gz > 5893_trim.fq.gz
cat SRR3299955_trim.fq.gz SRR3299956_trim.fq.gz > 5907_trim.fq.gz
cat SRR3299960_trim.fq.gz SRR3299961_trim.fq.gz SRR3299962_trim.fq.gz > 5950_trim.fq.gz
cat SRR3299966_trim.fq.gz SRR3299967_trim.fq.gz > 5993_trim.fq.gz
cat SRR3299968_trim.fq.gz SRR3299969_trim.fq.gz > 6009_trim.fq.gz
cat SRR3299970_trim.fq.gz SRR3299971_trim.fq.gz > 6017_trim.fq.gz
cat SRR3299972_trim.fq.gz SRR3299973_trim.fq.gz > 6019_trim.fq.gz
cat SRR3299974_trim.fq.gz SRR3299975_trim.fq.gz > 6020_trim.fq.gz
cat SRR3299976_trim.fq.gz SRR3299977_trim.fq.gz > 6024_trim.fq.gz
cat SRR3299978_trim.fq.gz SRR3299979_trim.fq.gz > 6025_trim.fq.gz
cat SRR3299980_trim.fq.gz SRR3299981_trim.fq.gz > 6035_trim.fq.gz
cat SRR3299982_trim.fq.gz SRR3299983_trim.fq.gz SRR3299984_trim.fq.gz > 6038_trim.fq.gz
cat SRR3299985_trim.fq.gz SRR3299986_trim.fq.gz SRR3299987_trim.fq.gz > 6041_trim.fq.gz
cat SRR3299988_trim.fq.gz SRR3299989_trim.fq.gz SRR3299990_trim.fq.gz > 6064_trim.fq.gz
cat SRR3299991_trim.fq.gz SRR3299992_trim.fq.gz > 6073_trim.fq.gz
cat SRR3299993_trim.fq.gz SRR3299994_trim.fq.gz > 6086_trim.fq.gz
cat SRR3299995_trim.fq.gz SRR3299996_trim.fq.gz > 6099_trim.fq.gz
cat SRR3299997_trim.fq.gz SRR3299998_trim.fq.gz > 6106_trim.fq.gz
cat SRR3299999_trim.fq.gz SRR3300000_trim.fq.gz SRR3300001_trim.fq.gz > 6107_trim.fq.gz
cat SRR3300002_trim.fq.gz SRR3300003_trim.fq.gz > 6108_trim.fq.gz
cat SRR3300004_trim.fq.gz SRR3300005_trim.fq.gz > 6111_trim.fq.gz
cat SRR3300006_trim.fq.gz SRR3300007_trim.fq.gz > 6133_trim.fq.gz
cat SRR3300008_trim.fq.gz SRR3300009_trim.fq.gz > 6137_trim.fq.gz
cat SRR3300010_trim.fq.gz SRR3300011_trim.fq.gz > 6138_trim.fq.gz
cat SRR3300012_trim.fq.gz SRR3300013_trim.fq.gz SRR3300014_trim.fq.gz SRR3300015_trim.fq.gz > 6140_trim.fq.gz
cat SRR3300016_trim.fq.gz SRR3300017_trim.fq.gz SRR3300018_trim.fq.gz SRR3300019_trim.fq.gz > 6145_trim.fq.gz
cat SRR3300020_trim.fq.gz SRR3300021_trim.fq.gz > 6148_trim.fq.gz
cat SRR3300022_trim.fq.gz SRR3300023_trim.fq.gz SRR3300024_trim.fq.gz > 6154_trim.fq.gz
cat SRR3300025_trim.fq.gz SRR3300026_trim.fq.gz SRR3300027_trim.fq.gz SRR3300028_trim.fq.gz > 6163_trim.fq.gz
cat SRR3300029_trim.fq.gz SRR3300030_trim.fq.gz > 6166_trim.fq.gz
cat SRR3300031_trim.fq.gz SRR3300032_trim.fq.gz SRR3300033_trim.fq.gz > 6169_trim.fq.gz
cat SRR3300034_trim.fq.gz SRR3300035_trim.fq.gz > 6171_trim.fq.gz
cat SRR3300036_trim.fq.gz SRR3300037_trim.fq.gz > 6174_trim.fq.gz
cat SRR3300038_trim.fq.gz SRR3300039_trim.fq.gz > 6177_trim.fq.gz
cat SRR3300040_trim.fq.gz SRR3300041_trim.fq.gz > 6195_trim.fq.gz
cat SRR3300042_trim.fq.gz SRR3300043_trim.fq.gz > 6198_trim.fq.gz
cat SRR3300044_trim.fq.gz SRR3300045_trim.fq.gz > 6201_trim.fq.gz
cat SRR3300046_trim.fq.gz SRR3300047_trim.fq.gz > 6203_trim.fq.gz
cat SRR3300048_trim.fq.gz SRR3300049_trim.fq.gz > 6218_trim.fq.gz
cat SRR3300050_trim.fq.gz SRR3300051_trim.fq.gz > 6220_trim.fq.gz
cat SRR3300052_trim.fq.gz SRR3300053_trim.fq.gz > 6241_trim.fq.gz
cat SRR3300054_trim.fq.gz SRR3300055_trim.fq.gz > 6255_trim.fq.gz
cat SRR3300056_trim.fq.gz SRR3300057_trim.fq.gz > 6268_trim.fq.gz
cat SRR3300058_trim.fq.gz SRR3300059_trim.fq.gz > 6284_trim.fq.gz
cat SRR3300063_trim.fq.gz SRR3300064_trim.fq.gz SRR3300065_trim.fq.gz > 6390_trim.fq.gz
cat SRR3300066_trim.fq.gz SRR3300067_trim.fq.gz > 6434_trim.fq.gz
cat SRR3300068_trim.fq.gz SRR3300069_trim.fq.gz SRR3300070_trim.fq.gz SRR3300071_trim.fq.gz > 6445_trim.fq.gz
cat SRR3300078_trim.fq.gz SRR3300079_trim.fq.gz > 6744_trim.fq.gz
cat SRR3300080_trim.fq.gz SRR3300081_trim.fq.gz SRR3300082_trim.fq.gz > 6830_trim.fq.gz
cat SRR3300084_trim.fq.gz SRR3300085_trim.fq.gz SRR3300086_trim.fq.gz > 6900_trim.fq.gz
cat SRR3300087_trim.fq.gz SRR3300088_trim.fq.gz > 6906_trim.fq.gz
cat SRR3300089_trim.fq.gz SRR3300090_trim.fq.gz SRR3300091_trim.fq.gz > 6907_trim.fq.gz
cat SRR3300093_trim.fq.gz SRR3300094_trim.fq.gz > 6917_trim.fq.gz
cat SRR3300098_trim.fq.gz SRR3300099_trim.fq.gz SRR3300100_trim.fq.gz > 6937_trim.fq.gz
cat SRR3300103_trim.fq.gz SRR3300104_trim.fq.gz > 6944_trim.fq.gz
cat SRR3300105_trim.fq.gz SRR3300106_trim.fq.gz > 6956_trim.fq.gz
cat SRR3300108_trim.fq.gz SRR3300109_trim.fq.gz > 6959_trim.fq.gz
cat SRR3300110_trim.fq.gz SRR3300111_trim.fq.gz SRR3300112_trim.fq.gz > 6970_trim.fq.gz
cat SRR3300114_trim.fq.gz SRR3300115_trim.fq.gz SRR3300116_trim.fq.gz > 6971_trim.fq.gz
cat SRR3300117_trim.fq.gz SRR3300118_trim.fq.gz > 6974_trim.fq.gz
cat SRR3300119_trim.fq.gz SRR3300120_trim.fq.gz > 6979_trim.fq.gz
cat SRR3300122_trim.fq.gz SRR3300123_trim.fq.gz SRR3300124_trim.fq.gz SRR3300125_trim.fq.gz > 6989_trim.fq.gz
cat SRR3300126_trim.fq.gz SRR3300127_trim.fq.gz SRR3300128_trim.fq.gz SRR3300129_trim.fq.gz SRR3300130_trim.fq.gz > 6992_trim.fq.gz
cat SRR3300134_trim.fq.gz SRR3300135_trim.fq.gz SRR3300136_trim.fq.gz > 7013_trim.fq.gz
cat SRR3300139_trim.fq.gz SRR3300140_trim.fq.gz SRR3300141_trim.fq.gz > 7028_trim.fq.gz
cat SRR3300143_trim.fq.gz SRR3300144_trim.fq.gz SRR3300145_trim.fq.gz > 7058_trim.fq.gz
cat SRR3300148_trim.fq.gz SRR3300149_trim.fq.gz > 7067_trim.fq.gz
cat SRR3300150_trim.fq.gz SRR3300151_trim.fq.gz > 7068_trim.fq.gz
cat SRR3300153_trim.fq.gz SRR3300154_trim.fq.gz > 7094_trim.fq.gz
cat SRR3300156_trim.fq.gz SRR3300157_trim.fq.gz SRR3300158_trim.fq.gz > 7102_trim.fq.gz
cat SRR3300159_trim.fq.gz SRR3300160_trim.fq.gz SRR3300161_trim.fq.gz > 7103_trim.fq.gz
cat SRR3300162_trim.fq.gz SRR3300163_trim.fq.gz SRR3300164_trim.fq.gz > 7106_trim.fq.gz
cat SRR3300165_trim.fq.gz SRR3300166_trim.fq.gz > 7107_trim.fq.gz
cat SRR3300167_trim.fq.gz SRR3300168_trim.fq.gz > 7119_trim.fq.gz
cat SRR3300172_trim.fq.gz SRR3300173_trim.fq.gz > 7158_trim.fq.gz
cat SRR3300174_trim.fq.gz SRR3300175_trim.fq.gz > 7161_trim.fq.gz
cat SRR3300176_trim.fq.gz SRR3300177_trim.fq.gz SRR3300178_trim.fq.gz > 7164_trim.fq.gz
cat SRR3300179_trim.fq.gz SRR3300180_trim.fq.gz > 7165_trim.fq.gz
cat SRR3300181_trim.fq.gz SRR3300182_trim.fq.gz SRR3300183_trim.fq.gz > 7167_trim.fq.gz
cat SRR3300184_trim.fq.gz SRR3300185_trim.fq.gz SRR3300186_trim.fq.gz > 7169_trim.fq.gz
cat SRR3300188_trim.fq.gz SRR3300189_trim.fq.gz SRR3300190_trim.fq.gz > 7186_trim.fq.gz
cat SRR3300191_trim.fq.gz SRR3300192_trim.fq.gz SRR3300193_trim.fq.gz > 7202_trim.fq.gz
cat SRR3300194_trim.fq.gz SRR3300195_trim.fq.gz SRR3300196_trim.fq.gz > 7218_trim.fq.gz
cat SRR3300198_trim.fq.gz SRR3300199_trim.fq.gz > 7244_trim.fq.gz
cat SRR3300200_trim.fq.gz SRR3300201_trim.fq.gz SRR3300202_trim.fq.gz SRR3300203_trim.fq.gz > 7250_trim.fq.gz
cat SRR3300205_trim.fq.gz SRR3300206_trim.fq.gz SRR3300207_trim.fq.gz > 7263_trim.fq.gz
cat SRR3300208_trim.fq.gz SRR3300209_trim.fq.gz > 7273_trim.fq.gz
cat SRR3300210_trim.fq.gz SRR3300211_trim.fq.gz > 7282_trim.fq.gz
cat SRR3300213_trim.fq.gz SRR3300214_trim.fq.gz SRR3300215_trim.fq.gz > 7288_trim.fq.gz
cat SRR3300216_trim.fq.gz SRR3300217_trim.fq.gz > 728_trim.fq.gz
cat SRR3300222_trim.fq.gz SRR3300223_trim.fq.gz SRR3300224_trim.fq.gz > 7307_trim.fq.gz
cat SRR3300225_trim.fq.gz SRR3300226_trim.fq.gz SRR3300227_trim.fq.gz SRR3300228_trim.fq.gz > 7308_trim.fq.gz
cat SRR3300229_trim.fq.gz SRR3300230_trim.fq.gz > 7316_trim.fq.gz
cat SRR3300231_trim.fq.gz SRR3300232_trim.fq.gz SRR3300233_trim.fq.gz > 7320_trim.fq.gz
cat SRR3300234_trim.fq.gz SRR3300235_trim.fq.gz > 7322_trim.fq.gz
cat SRR3300237_trim.fq.gz SRR3300238_trim.fq.gz SRR3300239_trim.fq.gz > 7327_trim.fq.gz
cat SRR3300240_trim.fq.gz SRR3300241_trim.fq.gz SRR3300242_trim.fq.gz SRR3300243_trim.fq.gz > 7328_trim.fq.gz
cat SRR3300244_trim.fq.gz SRR3300245_trim.fq.gz SRR3300246_trim.fq.gz > 7330_trim.fq.gz
cat SRR3300251_trim.fq.gz SRR3300252_trim.fq.gz > 7346_trim.fq.gz
cat SRR3300253_trim.fq.gz SRR3300254_trim.fq.gz > 7350_trim.fq.gz
cat SRR3300255_trim.fq.gz SRR3300256_trim.fq.gz > 7353_trim.fq.gz
cat SRR3300257_trim.fq.gz SRR3300258_trim.fq.gz SRR3300259_trim.fq.gz > 7354_trim.fq.gz
cat SRR3300260_trim.fq.gz SRR3300261_trim.fq.gz SRR3300262_trim.fq.gz > 7355_trim.fq.gz
cat SRR3300263_trim.fq.gz SRR3300264_trim.fq.gz SRR3300265_trim.fq.gz SRR3300266_trim.fq.gz SRR3300267_trim.fq.gz > 7373_trim.fq.gz
cat SRR3300268_trim.fq.gz SRR3300269_trim.fq.gz SRR3300270_trim.fq.gz > 7377_trim.fq.gz
cat SRR3300271_trim.fq.gz SRR3300272_trim.fq.gz SRR3300273_trim.fq.gz > 7384_trim.fq.gz
cat SRR3300274_trim.fq.gz SRR3300275_trim.fq.gz SRR3300276_trim.fq.gz SRR3300277_trim.fq.gz > 7396_trim.fq.gz
cat SRR3300278_trim.fq.gz SRR3300279_trim.fq.gz SRR3300280_trim.fq.gz SRR3300281_trim.fq.gz SRR3300282_trim.fq.gz > 7404_trim.fq.gz
cat SRR3300283_trim.fq.gz SRR3300284_trim.fq.gz SRR3300285_trim.fq.gz SRR3300286_trim.fq.gz > 7413_trim.fq.gz
cat SRR3300287_trim.fq.gz SRR3300288_trim.fq.gz > 7415_trim.fq.gz
cat SRR3300289_trim.fq.gz SRR3300290_trim.fq.gz > 7418_trim.fq.gz
cat SRR3300291_trim.fq.gz SRR3300292_trim.fq.gz > 7424_trim.fq.gz
cat SRR3300293_trim.fq.gz SRR3300294_trim.fq.gz SRR3300295_trim.fq.gz SRR3300296_trim.fq.gz SRR3300297_trim.fq.gz > 7427_trim.fq.gz
cat SRR3300298_trim.fq.gz SRR3300299_trim.fq.gz SRR3300300_trim.fq.gz SRR3300301_trim.fq.gz > 7458_trim.fq.gz
cat SRR3300304_trim.fq.gz SRR3300305_trim.fq.gz > 7477_trim.fq.gz
cat SRR3300308_trim.fq.gz SRR3300309_trim.fq.gz SRR3300310_trim.fq.gz SRR3300311_trim.fq.gz > 7515_trim.fq.gz
cat SRR3300315_trim.fq.gz SRR3300316_trim.fq.gz SRR3300317_trim.fq.gz > 763_trim.fq.gz
cat SRR3300318_trim.fq.gz SRR3300319_trim.fq.gz SRR3300320_trim.fq.gz SRR3300321_trim.fq.gz > 768_trim.fq.gz
cat SRR3300322_trim.fq.gz SRR3300323_trim.fq.gz > 7917_trim.fq.gz
cat SRR3300325_trim.fq.gz SRR3300326_trim.fq.gz SRR3300327_trim.fq.gz > 7947_trim.fq.gz
cat SRR3300328_trim.fq.gz SRR3300329_trim.fq.gz SRR3300330_trim.fq.gz SRR3300331_trim.fq.gz > 8077_trim.fq.gz
cat SRR3300332_trim.fq.gz SRR3300333_trim.fq.gz SRR3300334_trim.fq.gz > 8132_trim.fq.gz
cat SRR3300335_trim.fq.gz SRR3300336_trim.fq.gz > 8231_trim.fq.gz
cat SRR3300337_trim.fq.gz SRR3300338_trim.fq.gz SRR3300339_trim.fq.gz > 8233_trim.fq.gz
cat SRR3300344_trim.fq.gz SRR3300345_trim.fq.gz SRR3300346_trim.fq.gz > 8239_trim.fq.gz
cat SRR3300348_trim.fq.gz SRR3300349_trim.fq.gz > 8244_trim.fq.gz
cat SRR3300352_trim.fq.gz SRR3300353_trim.fq.gz > 8247_trim.fq.gz
cat SRR3300354_trim.fq.gz SRR3300355_trim.fq.gz SRR3300356_trim.fq.gz SRR3300357_trim.fq.gz SRR3300358_trim.fq.gz > 8259_trim.fq.gz
cat SRR3300360_trim.fq.gz SRR3300361_trim.fq.gz > 8266_trim.fq.gz
cat SRR3300362_trim.fq.gz SRR3300363_trim.fq.gz > 8275_trim.fq.gz
cat SRR3300367_trim.fq.gz SRR3300368_trim.fq.gz > 8290_trim.fq.gz
cat SRR3300369_trim.fq.gz SRR3300370_trim.fq.gz > 8326_trim.fq.gz
cat SRR3300371_trim.fq.gz SRR3300372_trim.fq.gz SRR3300373_trim.fq.gz > 8337_trim.fq.gz
cat SRR3300374_trim.fq.gz SRR3300375_trim.fq.gz SRR3300376_trim.fq.gz > 8343_trim.fq.gz
cat SRR3300377_trim.fq.gz SRR3300378_trim.fq.gz > 8353_trim.fq.gz
cat SRR3300386_trim.fq.gz SRR3300387_trim.fq.gz > 8386_trim.fq.gz
cat SRR3300392_trim.fq.gz SRR3300393_trim.fq.gz > 8428_trim.fq.gz
cat SRR3300395_trim.fq.gz SRR3300396_trim.fq.gz > 8472_trim.fq.gz
cat SRR3300397_trim.fq.gz SRR3300398_trim.fq.gz > 8584_trim.fq.gz
cat SRR3300402_trim.fq.gz SRR3300403_trim.fq.gz SRR3300404_trim.fq.gz > 9045_trim.fq.gz
cat SRR3300406_trim.fq.gz SRR3300407_trim.fq.gz > 9069_trim.fq.gz
cat SRR3300411_trim.fq.gz SRR3300412_trim.fq.gz SRR3300413_trim.fq.gz > 9075_trim.fq.gz
cat SRR3300415_trim.fq.gz SRR3300416_trim.fq.gz SRR3300417_trim.fq.gz > 9081_trim.fq.gz
cat SRR3300419_trim.fq.gz SRR3300420_trim.fq.gz > 9084_trim.fq.gz
cat SRR3300421_trim.fq.gz SRR3300422_trim.fq.gz > 9089_trim.fq.gz
cat SRR3300423_trim.fq.gz SRR3300424_trim.fq.gz SRR3300425_trim.fq.gz > 9091_trim.fq.gz
cat SRR3300426_trim.fq.gz SRR3300427_trim.fq.gz > 9095_trim.fq.gz
cat SRR3300429_trim.fq.gz SRR3300430_trim.fq.gz > 9099_trim.fq.gz
cat SRR3300432_trim.fq.gz SRR3300433_trim.fq.gz SRR3300434_trim.fq.gz > 9103_trim.fq.gz
cat SRR3300435_trim.fq.gz SRR3300436_trim.fq.gz > 9104_trim.fq.gz
cat SRR3300437_trim.fq.gz SRR3300438_trim.fq.gz SRR3300439_trim.fq.gz SRR3300440_trim.fq.gz > 9105_trim.fq.gz
cat SRR3300442_trim.fq.gz SRR3300443_trim.fq.gz SRR3300444_trim.fq.gz > 9106_trim.fq.gz
cat SRR3300446_trim.fq.gz SRR3300447_trim.fq.gz > 9111_trim.fq.gz
cat SRR3300451_trim.fq.gz SRR3300452_trim.fq.gz SRR3300453_trim.fq.gz > 9114_trim.fq.gz
cat SRR3300455_trim.fq.gz SRR3300456_trim.fq.gz SRR3300457_trim.fq.gz > 9115_trim.fq.gz
cat SRR3300458_trim.fq.gz SRR3300459_trim.fq.gz SRR3300460_trim.fq.gz SRR3300461_trim.fq.gz > 9125_trim.fq.gz
cat SRR3300462_trim.fq.gz SRR3300463_trim.fq.gz > 9131_trim.fq.gz
cat SRR3300464_trim.fq.gz SRR3300465_trim.fq.gz SRR3300466_trim.fq.gz > 9134_trim.fq.gz
cat SRR3300468_trim.fq.gz SRR3300469_trim.fq.gz SRR3300470_trim.fq.gz > 9298_trim.fq.gz
cat SRR3300472_trim.fq.gz SRR3300473_trim.fq.gz > 9312_trim.fq.gz
cat SRR3300474_trim.fq.gz SRR3300475_trim.fq.gz > 9314_trim.fq.gz
cat SRR3300476_trim.fq.gz SRR3300477_trim.fq.gz SRR3300478_trim.fq.gz > 9321_trim.fq.gz
cat SRR3300479_trim.fq.gz SRR3300480_trim.fq.gz > 932_trim.fq.gz
cat SRR3300481_trim.fq.gz SRR3300482_trim.fq.gz > 9386_trim.fq.gz
cat SRR3300483_trim.fq.gz SRR3300484_trim.fq.gz > 9388_trim.fq.gz
cat SRR3300485_trim.fq.gz SRR3300486_trim.fq.gz > 9404_trim.fq.gz
cat SRR3300487_trim.fq.gz SRR3300488_trim.fq.gz > 9405_trim.fq.gz
cat SRR3300489_trim.fq.gz SRR3300490_trim.fq.gz > 9407_trim.fq.gz
cat SRR3300491_trim.fq.gz SRR3300492_trim.fq.gz > 9408_trim.fq.gz
cat SRR3300493_trim.fq.gz SRR3300494_trim.fq.gz > 9409_trim.fq.gz
cat SRR3300495_trim.fq.gz SRR3300496_trim.fq.gz SRR3300497_trim.fq.gz > 9421_trim.fq.gz
cat SRR3300498_trim.fq.gz SRR3300499_trim.fq.gz SRR3300500_trim.fq.gz > 9433_trim.fq.gz
cat SRR3300501_trim.fq.gz SRR3300502_trim.fq.gz > 9450_trim.fq.gz
cat SRR3300503_trim.fq.gz SRR3300504_trim.fq.gz > 9451_trim.fq.gz
cat SRR3300505_trim.fq.gz SRR3300506_trim.fq.gz SRR3300507_trim.fq.gz SRR3300508_trim.fq.gz > 9503_trim.fq.gz
cat SRR3300509_trim.fq.gz SRR3300510_trim.fq.gz SRR3300511_trim.fq.gz > 9506_trim.fq.gz
cat SRR3300512_trim.fq.gz SRR3300513_trim.fq.gz SRR3300514_trim.fq.gz > 9507_trim.fq.gz
cat SRR3300515_trim.fq.gz SRR3300516_trim.fq.gz SRR3300517_trim.fq.gz > 9508_trim.fq.gz
cat SRR3300518_trim.fq.gz SRR3300519_trim.fq.gz SRR3300520_trim.fq.gz > 9509_trim.fq.gz
cat SRR3300521_trim.fq.gz SRR3300522_trim.fq.gz SRR3300523_trim.fq.gz > 9510_trim.fq.gz
cat SRR3300524_trim.fq.gz SRR3300525_trim.fq.gz SRR3300526_trim.fq.gz > 9511_trim.fq.gz
cat SRR3300527_trim.fq.gz SRR3300528_trim.fq.gz SRR3300529_trim.fq.gz > 9512_trim.fq.gz
cat SRR3300530_trim.fq.gz SRR3300531_trim.fq.gz SRR3300532_trim.fq.gz > 9513_trim.fq.gz
cat SRR3300533_trim.fq.gz SRR3300534_trim.fq.gz SRR3300535_trim.fq.gz > 9514_trim.fq.gz
cat SRR3300536_trim.fq.gz SRR3300537_trim.fq.gz SRR3300538_trim.fq.gz > 9515_trim.fq.gz
cat SRR3300539_trim.fq.gz SRR3300540_trim.fq.gz SRR3300541_trim.fq.gz > 9516_trim.fq.gz
cat SRR3300542_trim.fq.gz SRR3300543_trim.fq.gz > 9517_trim.fq.gz
cat SRR3300544_trim.fq.gz SRR3300545_trim.fq.gz SRR3300546_trim.fq.gz > 9518_trim.fq.gz
cat SRR3300547_trim.fq.gz SRR3300548_trim.fq.gz SRR3300549_trim.fq.gz > 9519_trim.fq.gz
cat SRR3300550_trim.fq.gz SRR3300551_trim.fq.gz SRR3300552_trim.fq.gz > 9520_trim.fq.gz
cat SRR3300553_trim.fq.gz SRR3300554_trim.fq.gz SRR3300555_trim.fq.gz > 9521_trim.fq.gz
cat SRR3300556_trim.fq.gz SRR3300557_trim.fq.gz SRR3300558_trim.fq.gz > 9522_trim.fq.gz
cat SRR3300559_trim.fq.gz SRR3300560_trim.fq.gz SRR3300561_trim.fq.gz > 9523_trim.fq.gz
cat SRR3300562_trim.fq.gz SRR3300563_trim.fq.gz SRR3300564_trim.fq.gz > 9524_trim.fq.gz
cat SRR3300565_trim.fq.gz SRR3300566_trim.fq.gz SRR3300567_trim.fq.gz > 9525_trim.fq.gz
cat SRR3300568_trim.fq.gz SRR3300569_trim.fq.gz SRR3300570_trim.fq.gz > 9526_trim.fq.gz
cat SRR3300571_trim.fq.gz SRR3300572_trim.fq.gz SRR3300573_trim.fq.gz > 9527_trim.fq.gz
cat SRR3300574_trim.fq.gz SRR3300575_trim.fq.gz SRR3300576_trim.fq.gz > 9529_trim.fq.gz
cat SRR3300577_trim.fq.gz SRR3300578_trim.fq.gz SRR3300579_trim.fq.gz > 9530_trim.fq.gz
cat SRR3300580_trim.fq.gz SRR3300581_trim.fq.gz SRR3300582_trim.fq.gz > 9531_trim.fq.gz
cat SRR3300583_trim.fq.gz SRR3300584_trim.fq.gz SRR3300585_trim.fq.gz > 9532_trim.fq.gz
cat SRR3300586_trim.fq.gz SRR3300587_trim.fq.gz SRR3300588_trim.fq.gz > 9533_trim.fq.gz
cat SRR3300589_trim.fq.gz SRR3300590_trim.fq.gz SRR3300591_trim.fq.gz > 9534_trim.fq.gz
cat SRR3300592_trim.fq.gz SRR3300593_trim.fq.gz SRR3300594_trim.fq.gz > 9535_trim.fq.gz
cat SRR3300595_trim.fq.gz SRR3300596_trim.fq.gz SRR3300597_trim.fq.gz > 9537_trim.fq.gz
cat SRR3300598_trim.fq.gz SRR3300599_trim.fq.gz SRR3300600_trim.fq.gz > 9538_trim.fq.gz
cat SRR3300601_trim.fq.gz SRR3300602_trim.fq.gz > 9540_trim.fq.gz
cat SRR3300603_trim.fq.gz SRR3300604_trim.fq.gz SRR3300605_trim.fq.gz > 9541_trim.fq.gz
cat SRR3300606_trim.fq.gz SRR3300607_trim.fq.gz > 9542_trim.fq.gz
cat SRR3300608_trim.fq.gz SRR3300609_trim.fq.gz SRR3300610_trim.fq.gz > 9543_trim.fq.gz
cat SRR3300611_trim.fq.gz SRR3300612_trim.fq.gz SRR3300613_trim.fq.gz > 9544_trim.fq.gz
cat SRR3300614_trim.fq.gz SRR3300615_trim.fq.gz SRR3300616_trim.fq.gz > 9545_trim.fq.gz
cat SRR3300617_trim.fq.gz SRR3300618_trim.fq.gz SRR3300619_trim.fq.gz > 9546_trim.fq.gz
cat SRR3300620_trim.fq.gz SRR3300621_trim.fq.gz SRR3300622_trim.fq.gz > 9547_trim.fq.gz
cat SRR3300623_trim.fq.gz SRR3300624_trim.fq.gz SRR3300625_trim.fq.gz > 9548_trim.fq.gz
cat SRR3300626_trim.fq.gz SRR3300627_trim.fq.gz SRR3300628_trim.fq.gz > 9549_trim.fq.gz
cat SRR3300629_trim.fq.gz SRR3300630_trim.fq.gz SRR3300631_trim.fq.gz > 9550_trim.fq.gz
cat SRR3300632_trim.fq.gz SRR3300633_trim.fq.gz SRR3300634_trim.fq.gz > 9551_trim.fq.gz
cat SRR3300635_trim.fq.gz SRR3300636_trim.fq.gz SRR3300637_trim.fq.gz > 9552_trim.fq.gz
cat SRR3300638_trim.fq.gz SRR3300639_trim.fq.gz SRR3300640_trim.fq.gz > 9553_trim.fq.gz
cat SRR3300641_trim.fq.gz SRR3300642_trim.fq.gz SRR3300643_trim.fq.gz > 9554_trim.fq.gz
cat SRR3300644_trim.fq.gz SRR3300645_trim.fq.gz SRR3300646_trim.fq.gz > 9556_trim.fq.gz
cat SRR3300647_trim.fq.gz SRR3300648_trim.fq.gz SRR3300649_trim.fq.gz > 9557_trim.fq.gz
cat SRR3300650_trim.fq.gz SRR3300651_trim.fq.gz SRR3300652_trim.fq.gz > 9558_trim.fq.gz
cat SRR3300653_trim.fq.gz SRR3300654_trim.fq.gz SRR3300655_trim.fq.gz > 9559_trim.fq.gz
cat SRR3300656_trim.fq.gz SRR3300657_trim.fq.gz SRR3300658_trim.fq.gz > 9560_trim.fq.gz
cat SRR3300659_trim.fq.gz SRR3300660_trim.fq.gz SRR3300661_trim.fq.gz > 9561_trim.fq.gz
cat SRR3300662_trim.fq.gz SRR3300663_trim.fq.gz SRR3300664_trim.fq.gz > 9562_trim.fq.gz
cat SRR3300665_trim.fq.gz SRR3300666_trim.fq.gz SRR3300667_trim.fq.gz > 9563_trim.fq.gz
cat SRR3300668_trim.fq.gz SRR3300669_trim.fq.gz SRR3300670_trim.fq.gz > 9564_trim.fq.gz
cat SRR3300671_trim.fq.gz SRR3300672_trim.fq.gz SRR3300673_trim.fq.gz > 9565_trim.fq.gz
cat SRR3300674_trim.fq.gz SRR3300675_trim.fq.gz SRR3300676_trim.fq.gz > 9566_trim.fq.gz
cat SRR3300677_trim.fq.gz SRR3300678_trim.fq.gz SRR3300679_trim.fq.gz > 9567_trim.fq.gz
cat SRR3300680_trim.fq.gz SRR3300681_trim.fq.gz SRR3300682_trim.fq.gz > 9568_trim.fq.gz
cat SRR3300683_trim.fq.gz SRR3300684_trim.fq.gz > 9569_trim.fq.gz
cat SRR3300685_trim.fq.gz SRR3300686_trim.fq.gz SRR3300687_trim.fq.gz > 9570_trim.fq.gz
cat SRR3300688_trim.fq.gz SRR3300689_trim.fq.gz SRR3300690_trim.fq.gz > 9572_trim.fq.gz
cat SRR3300691_trim.fq.gz SRR3300692_trim.fq.gz SRR3300693_trim.fq.gz > 9573_trim.fq.gz
cat SRR3300694_trim.fq.gz SRR3300695_trim.fq.gz SRR3300696_trim.fq.gz > 9574_trim.fq.gz
cat SRR3300697_trim.fq.gz SRR3300698_trim.fq.gz SRR3300699_trim.fq.gz > 9575_trim.fq.gz
cat SRR3300700_trim.fq.gz SRR3300701_trim.fq.gz SRR3300702_trim.fq.gz > 9576_trim.fq.gz
cat SRR3300703_trim.fq.gz SRR3300704_trim.fq.gz SRR3300705_trim.fq.gz > 9577_trim.fq.gz
cat SRR3300706_trim.fq.gz SRR3300707_trim.fq.gz SRR3300708_trim.fq.gz > 9578_trim.fq.gz
cat SRR3300709_trim.fq.gz SRR3300710_trim.fq.gz > 9579_trim.fq.gz
cat SRR3300711_trim.fq.gz SRR3300712_trim.fq.gz SRR3300713_trim.fq.gz > 9580_trim.fq.gz
cat SRR3300714_trim.fq.gz SRR3300715_trim.fq.gz SRR3300716_trim.fq.gz > 9581_trim.fq.gz
cat SRR3300717_trim.fq.gz SRR3300718_trim.fq.gz SRR3300719_trim.fq.gz > 9582_trim.fq.gz
cat SRR3300720_trim.fq.gz SRR3300721_trim.fq.gz > 9583_trim.fq.gz
cat SRR3300722_trim.fq.gz SRR3300723_trim.fq.gz SRR3300724_trim.fq.gz > 9584_trim.fq.gz
cat SRR3300725_trim.fq.gz SRR3300726_trim.fq.gz SRR3300727_trim.fq.gz > 9585_trim.fq.gz
cat SRR3300728_trim.fq.gz SRR3300729_trim.fq.gz SRR3300730_trim.fq.gz > 9586_trim.fq.gz
cat SRR3300731_trim.fq.gz SRR3300732_trim.fq.gz SRR3300733_trim.fq.gz > 9587_trim.fq.gz
cat SRR3300734_trim.fq.gz SRR3300735_trim.fq.gz SRR3300736_trim.fq.gz > 9588_trim.fq.gz
cat SRR3300737_trim.fq.gz SRR3300738_trim.fq.gz SRR3300739_trim.fq.gz > 9589_trim.fq.gz
cat SRR3300740_trim.fq.gz SRR3300741_trim.fq.gz SRR3300742_trim.fq.gz > 9590_trim.fq.gz
cat SRR3300743_trim.fq.gz SRR3300744_trim.fq.gz SRR3300745_trim.fq.gz > 9591_trim.fq.gz
cat SRR3300746_trim.fq.gz SRR3300747_trim.fq.gz SRR3300748_trim.fq.gz > 9592_trim.fq.gz
cat SRR3300749_trim.fq.gz SRR3300750_trim.fq.gz SRR3300751_trim.fq.gz > 9593_trim.fq.gz
cat SRR3300752_trim.fq.gz SRR3300753_trim.fq.gz SRR3300754_trim.fq.gz > 9594_trim.fq.gz
cat SRR3300755_trim.fq.gz SRR3300756_trim.fq.gz SRR3300757_trim.fq.gz > 9595_trim.fq.gz
cat SRR3300758_trim.fq.gz SRR3300759_trim.fq.gz > 9596_trim.fq.gz
cat SRR3300760_trim.fq.gz SRR3300761_trim.fq.gz SRR3300762_trim.fq.gz > 9597_trim.fq.gz
cat SRR3300763_trim.fq.gz SRR3300764_trim.fq.gz SRR3300765_trim.fq.gz > 9598_trim.fq.gz
cat SRR3300766_trim.fq.gz SRR3300767_trim.fq.gz SRR3300768_trim.fq.gz > 9599_trim.fq.gz
cat SRR3300769_trim.fq.gz SRR3300770_trim.fq.gz SRR3300771_trim.fq.gz > 9600_trim.fq.gz
cat SRR3300772_trim.fq.gz SRR3300773_trim.fq.gz SRR3300774_trim.fq.gz > 9601_trim.fq.gz
cat SRR3300775_trim.fq.gz SRR3300776_trim.fq.gz SRR3300777_trim.fq.gz SRR3300778_trim.fq.gz SRR3300779_trim.fq.gz > 9606_trim.fq.gz
cat SRR3300780_trim.fq.gz SRR3300781_trim.fq.gz SRR3300782_trim.fq.gz SRR3300783_trim.fq.gz SRR3300784_trim.fq.gz SRR3300785_trim.fq.gz SRR3300786_trim.fq.gz SRR3300787_trim.fq.gz > 9607_trim.fq.gz
cat SRR3300788_trim.fq.gz SRR3300789_trim.fq.gz > 9608_trim.fq.gz
cat SRR3300790_trim.fq.gz SRR3300791_trim.fq.gz > 9609_trim.fq.gz
cat SRR3300792_trim.fq.gz SRR3300793_trim.fq.gz SRR3300794_trim.fq.gz SRR3300795_trim.fq.gz > 9610_trim.fq.gz
cat SRR3300796_trim.fq.gz SRR3300797_trim.fq.gz > 9611_trim.fq.gz
cat SRR3300798_trim.fq.gz SRR3300799_trim.fq.gz > 9612_trim.fq.gz
cat SRR3300800_trim.fq.gz SRR3300801_trim.fq.gz SRR3300802_trim.fq.gz > 9613_trim.fq.gz
cat SRR3300803_trim.fq.gz SRR3300804_trim.fq.gz SRR3300805_trim.fq.gz > 9615_trim.fq.gz
cat SRR3300806_trim.fq.gz SRR3300807_trim.fq.gz > 9616_trim.fq.gz
cat SRR3300808_trim.fq.gz SRR3300809_trim.fq.gz SRR3300810_trim.fq.gz > 9617_trim.fq.gz
cat SRR3300811_trim.fq.gz SRR3300812_trim.fq.gz SRR3300813_trim.fq.gz SRR3300814_trim.fq.gz > 9619_trim.fq.gz
cat SRR3300815_trim.fq.gz SRR3300816_trim.fq.gz SRR3300817_trim.fq.gz SRR3300818_trim.fq.gz > 9620_trim.fq.gz
cat SRR3300819_trim.fq.gz SRR3300820_trim.fq.gz > 9621_trim.fq.gz
cat SRR3300821_trim.fq.gz SRR3300822_trim.fq.gz > 9623_trim.fq.gz
cat SRR3300823_trim.fq.gz SRR3300824_trim.fq.gz > 9625_trim.fq.gz
cat SRR3300825_trim.fq.gz SRR3300826_trim.fq.gz SRR3300827_trim.fq.gz SRR3300828_trim.fq.gz > 9626_trim.fq.gz
cat SRR3300829_trim.fq.gz SRR3300830_trim.fq.gz SRR3300831_trim.fq.gz > 9627_trim.fq.gz
cat SRR3300832_trim.fq.gz SRR3300833_trim.fq.gz SRR3300834_trim.fq.gz SRR3300835_trim.fq.gz > 9628_trim.fq.gz
cat SRR3300836_trim.fq.gz SRR3300837_trim.fq.gz SRR3300838_trim.fq.gz SRR3300839_trim.fq.gz > 9629_trim.fq.gz
cat SRR3300840_trim.fq.gz SRR3300841_trim.fq.gz SRR3300842_trim.fq.gz > 9631_trim.fq.gz
cat SRR3300843_trim.fq.gz SRR3300844_trim.fq.gz SRR3300845_trim.fq.gz SRR3300846_trim.fq.gz SRR3300847_trim.fq.gz > 9632_trim.fq.gz
cat SRR3300848_trim.fq.gz SRR3300849_trim.fq.gz > 9633_trim.fq.gz
cat SRR3300850_trim.fq.gz SRR3300851_trim.fq.gz SRR3300852_trim.fq.gz > 9634_trim.fq.gz
cat SRR3300853_trim.fq.gz SRR3300854_trim.fq.gz > 9635_trim.fq.gz
cat SRR3300855_trim.fq.gz SRR3300856_trim.fq.gz SRR3300857_trim.fq.gz SRR3300858_trim.fq.gz > 9636_trim.fq.gz
cat SRR3300859_trim.fq.gz SRR3300860_trim.fq.gz SRR3300861_trim.fq.gz SRR3300862_trim.fq.gz > 9637_trim.fq.gz
cat SRR3300863_trim.fq.gz SRR3300864_trim.fq.gz SRR3300865_trim.fq.gz SRR3300866_trim.fq.gz SRR3300867_trim.fq.gz > 9638_trim.fq.gz
cat SRR3300868_trim.fq.gz SRR3300869_trim.fq.gz SRR3300870_trim.fq.gz SRR3300871_trim.fq.gz SRR3300872_trim.fq.gz SRR3300873_trim.fq.gz SRR3300874_trim.fq.gz > 9639_trim.fq.gz
cat SRR3300875_trim.fq.gz SRR3300876_trim.fq.gz > 9640_trim.fq.gz
cat SRR3300877_trim.fq.gz SRR3300878_trim.fq.gz > 9641_trim.fq.gz
cat SRR3300879_trim.fq.gz SRR3300880_trim.fq.gz SRR3300881_trim.fq.gz SRR3300882_trim.fq.gz > 9643_trim.fq.gz
cat SRR3300883_trim.fq.gz SRR3300884_trim.fq.gz > 9644_trim.fq.gz
cat SRR3300885_trim.fq.gz SRR3300886_trim.fq.gz SRR3300887_trim.fq.gz > 9645_trim.fq.gz
cat SRR3300888_trim.fq.gz SRR3300889_trim.fq.gz > 9646_trim.fq.gz
cat SRR3300890_trim.fq.gz SRR3300891_trim.fq.gz > 9647_trim.fq.gz
cat SRR3300892_trim.fq.gz SRR3300893_trim.fq.gz SRR3300894_trim.fq.gz > 9648_trim.fq.gz
cat SRR3300895_trim.fq.gz SRR3300896_trim.fq.gz > 9649_trim.fq.gz
cat SRR3300897_trim.fq.gz SRR3300898_trim.fq.gz > 9650_trim.fq.gz
cat SRR3300899_trim.fq.gz SRR3300900_trim.fq.gz SRR3300901_trim.fq.gz SRR3300902_trim.fq.gz > 9651_trim.fq.gz
cat SRR3300903_trim.fq.gz SRR3300904_trim.fq.gz SRR3300905_trim.fq.gz SRR3300906_trim.fq.gz > 9652_trim.fq.gz
cat SRR3300907_trim.fq.gz SRR3300908_trim.fq.gz > 9653_trim.fq.gz
cat SRR3300909_trim.fq.gz SRR3300910_trim.fq.gz > 9654_trim.fq.gz
cat SRR3300911_trim.fq.gz SRR3300912_trim.fq.gz SRR3300913_trim.fq.gz > 9655_trim.fq.gz
cat SRR3300914_trim.fq.gz SRR3300915_trim.fq.gz SRR3300916_trim.fq.gz > 9656_trim.fq.gz
cat SRR3300917_trim.fq.gz SRR3300918_trim.fq.gz > 9657_trim.fq.gz
cat SRR3300919_trim.fq.gz SRR3300920_trim.fq.gz SRR3300921_trim.fq.gz SRR3300922_trim.fq.gz > 9658_trim.fq.gz
cat SRR3300923_trim.fq.gz SRR3300924_trim.fq.gz SRR3300925_trim.fq.gz > 9659_trim.fq.gz
cat SRR3300926_trim.fq.gz SRR3300927_trim.fq.gz SRR3300928_trim.fq.gz > 9660_trim.fq.gz
cat SRR3300929_trim.fq.gz SRR3300930_trim.fq.gz SRR3300931_trim.fq.gz > 9661_trim.fq.gz
cat SRR3300932_trim.fq.gz SRR3300933_trim.fq.gz > 9662_trim.fq.gz
cat SRR3300934_trim.fq.gz SRR3300935_trim.fq.gz SRR3300936_trim.fq.gz > 9668_trim.fq.gz
cat SRR3300937_trim.fq.gz SRR3300938_trim.fq.gz SRR3300939_trim.fq.gz > 9671_trim.fq.gz
cat SRR3300940_trim.fq.gz SRR3300941_trim.fq.gz SRR3300942_trim.fq.gz > 9676_trim.fq.gz
cat SRR3300943_trim.fq.gz SRR3300944_trim.fq.gz > 9678_trim.fq.gz
cat SRR3300945_trim.fq.gz SRR3300946_trim.fq.gz SRR3300947_trim.fq.gz > 9679_trim.fq.gz
cat SRR3300948_trim.fq.gz SRR3300949_trim.fq.gz > 9684_trim.fq.gz
cat SRR3300950_trim.fq.gz SRR3300951_trim.fq.gz SRR3300952_trim.fq.gz > 9690_trim.fq.gz
cat SRR3300953_trim.fq.gz SRR3300954_trim.fq.gz SRR3300955_trim.fq.gz > 9693_trim.fq.gz
cat SRR3300958_trim.fq.gz SRR3300959_trim.fq.gz SRR3300960_trim.fq.gz > 9696_trim.fq.gz
cat SRR3300961_trim.fq.gz SRR3300962_trim.fq.gz SRR3300963_trim.fq.gz SRR3300964_trim.fq.gz > 9697_trim.fq.gz
cat SRR3300965_trim.fq.gz SRR3300966_trim.fq.gz SRR3300967_trim.fq.gz SRR3300968_trim.fq.gz > 9698_trim.fq.gz
cat SRR3300969_trim.fq.gz SRR3300970_trim.fq.gz SRR3300971_trim.fq.gz SRR3300972_trim.fq.gz > 9699_trim.fq.gz
cat SRR3300973_trim.fq.gz SRR3300974_trim.fq.gz SRR3300975_trim.fq.gz SRR3300976_trim.fq.gz SRR3300977_trim.fq.gz > 9700_trim.fq.gz
cat SRR3300978_trim.fq.gz SRR3300979_trim.fq.gz > 9701_trim.fq.gz
cat SRR3300980_trim.fq.gz SRR3300981_trim.fq.gz > 9702_trim.fq.gz
cat SRR3300982_trim.fq.gz SRR3300983_trim.fq.gz > 9703_trim.fq.gz
cat SRR3300984_trim.fq.gz SRR3300985_trim.fq.gz SRR3300986_trim.fq.gz > 9704_trim.fq.gz
cat SRR3300987_trim.fq.gz SRR3300988_trim.fq.gz SRR3300989_trim.fq.gz SRR3300990_trim.fq.gz SRR3300991_trim.fq.gz > 9705_trim.fq.gz
cat SRR3300992_trim.fq.gz SRR3300993_trim.fq.gz SRR3300994_trim.fq.gz SRR3300995_trim.fq.gz SRR3300996_trim.fq.gz > 9706_trim.fq.gz
cat SRR3300997_trim.fq.gz SRR3300998_trim.fq.gz SRR3300999_trim.fq.gz > 9707_trim.fq.gz
cat SRR3301000_trim.fq.gz SRR3301001_trim.fq.gz > 9708_trim.fq.gz
cat SRR3301002_trim.fq.gz SRR3301003_trim.fq.gz > 9709_trim.fq.gz
cat SRR3301004_trim.fq.gz SRR3301005_trim.fq.gz SRR3301006_trim.fq.gz > 9710_trim.fq.gz
cat SRR3301007_trim.fq.gz SRR3301008_trim.fq.gz > 9713_trim.fq.gz
cat SRR3301009_trim.fq.gz SRR3301010_trim.fq.gz SRR3301011_trim.fq.gz SRR3301012_trim.fq.gz > 9714_trim.fq.gz
cat SRR3301013_trim.fq.gz SRR3301014_trim.fq.gz SRR3301015_trim.fq.gz > 9716_trim.fq.gz
cat SRR3301016_trim.fq.gz SRR3301017_trim.fq.gz SRR3301018_trim.fq.gz SRR3301019_trim.fq.gz > 9717_trim.fq.gz
cat SRR3301020_trim.fq.gz SRR3301021_trim.fq.gz SRR3301022_trim.fq.gz SRR3301023_trim.fq.gz > 9718_trim.fq.gz
cat SRR3301024_trim.fq.gz SRR3301025_trim.fq.gz > 9719_trim.fq.gz
cat SRR3301026_trim.fq.gz SRR3301027_trim.fq.gz SRR3301028_trim.fq.gz SRR3301029_trim.fq.gz SRR3301030_trim.fq.gz SRR3301031_trim.fq.gz > 9721_trim.fq.gz
cat SRR3301032_trim.fq.gz SRR3301033_trim.fq.gz > 9722_trim.fq.gz
cat SRR3301034_trim.fq.gz SRR3301035_trim.fq.gz SRR3301036_trim.fq.gz SRR3301037_trim.fq.gz > 9723_trim.fq.gz
cat SRR3301038_trim.fq.gz SRR3301039_trim.fq.gz SRR3301040_trim.fq.gz > 9726_trim.fq.gz
cat SRR3301041_trim.fq.gz SRR3301042_trim.fq.gz > 9727_trim.fq.gz
cat SRR3301043_trim.fq.gz SRR3301044_trim.fq.gz > 9728_trim.fq.gz
cat SRR3301045_trim.fq.gz SRR3301046_trim.fq.gz SRR3301047_trim.fq.gz SRR3301048_trim.fq.gz > 9729_trim.fq.gz
cat SRR3301049_trim.fq.gz SRR3301050_trim.fq.gz > 9730_trim.fq.gz
cat SRR3301051_trim.fq.gz SRR3301052_trim.fq.gz SRR3301053_trim.fq.gz SRR3301054_trim.fq.gz > 9732_trim.fq.gz
cat SRR3301055_trim.fq.gz SRR3301056_trim.fq.gz > 9733_trim.fq.gz
cat SRR3301057_trim.fq.gz SRR3301058_trim.fq.gz SRR3301059_trim.fq.gz SRR3301060_trim.fq.gz > 9734_trim.fq.gz
cat SRR3301061_trim.fq.gz SRR3301062_trim.fq.gz SRR3301063_trim.fq.gz > 9736_trim.fq.gz
cat SRR3301064_trim.fq.gz SRR3301065_trim.fq.gz SRR3301066_trim.fq.gz SRR3301067_trim.fq.gz > 9737_trim.fq.gz
cat SRR3301068_trim.fq.gz SRR3301069_trim.fq.gz SRR3301070_trim.fq.gz > 9738_trim.fq.gz
cat SRR3301071_trim.fq.gz SRR3301072_trim.fq.gz > 9739_trim.fq.gz
cat SRR3301073_trim.fq.gz SRR3301074_trim.fq.gz > 9741_trim.fq.gz
cat SRR3301075_trim.fq.gz SRR3301076_trim.fq.gz > 9743_trim.fq.gz
cat SRR3301077_trim.fq.gz SRR3301078_trim.fq.gz > 9744_trim.fq.gz
cat SRR3301079_trim.fq.gz SRR3301080_trim.fq.gz > 9745_trim.fq.gz
cat SRR3301081_trim.fq.gz SRR3301082_trim.fq.gz > 9746_trim.fq.gz
cat SRR3301083_trim.fq.gz SRR3301084_trim.fq.gz SRR3301085_trim.fq.gz SRR3301086_trim.fq.gz > 9747_trim.fq.gz
cat SRR3301087_trim.fq.gz SRR3301088_trim.fq.gz > 9749_trim.fq.gz
cat SRR3301089_trim.fq.gz SRR3301090_trim.fq.gz > 9754_trim.fq.gz
cat SRR3301091_trim.fq.gz SRR3301092_trim.fq.gz SRR3301093_trim.fq.gz SRR3301094_trim.fq.gz > 9755_trim.fq.gz
cat SRR3301095_trim.fq.gz SRR3301096_trim.fq.gz SRR3301097_trim.fq.gz SRR3301098_trim.fq.gz > 9756_trim.fq.gz
cat SRR3301102_trim.fq.gz SRR3301103_trim.fq.gz > 9768_trim.fq.gz
cat SRR3301104_trim.fq.gz SRR3301105_trim.fq.gz > 9769_trim.fq.gz
cat SRR3301106_trim.fq.gz SRR3301107_trim.fq.gz SRR3301108_trim.fq.gz > 9770_trim.fq.gz
cat SRR3301109_trim.fq.gz SRR3301110_trim.fq.gz > 9771_trim.fq.gz
cat SRR3301111_trim.fq.gz SRR3301112_trim.fq.gz > 9772_trim.fq.gz
cat SRR3301113_trim.fq.gz SRR3301114_trim.fq.gz SRR3301115_trim.fq.gz SRR3301116_trim.fq.gz > 9774_trim.fq.gz
cat SRR3301117_trim.fq.gz SRR3301118_trim.fq.gz > 9775_trim.fq.gz
cat SRR3301119_trim.fq.gz SRR3301120_trim.fq.gz SRR3301121_trim.fq.gz > 9776_trim.fq.gz
cat SRR3301122_trim.fq.gz SRR3301123_trim.fq.gz > 9777_trim.fq.gz
cat SRR3301124_trim.fq.gz SRR3301125_trim.fq.gz SRR3301126_trim.fq.gz > 9778_trim.fq.gz
cat SRR3301127_trim.fq.gz SRR3301128_trim.fq.gz > 9779_trim.fq.gz
cat SRR3301129_trim.fq.gz SRR3301130_trim.fq.gz > 9780_trim.fq.gz
cat SRR3301131_trim.fq.gz SRR3301132_trim.fq.gz > 9781_trim.fq.gz
cat SRR3301133_trim.fq.gz SRR3301134_trim.fq.gz > 9782_trim.fq.gz
cat SRR3301137_trim.fq.gz SRR3301138_trim.fq.gz SRR3301139_trim.fq.gz > 9784_trim.fq.gz
cat SRR3301140_trim.fq.gz SRR3301141_trim.fq.gz SRR3301142_trim.fq.gz > 9785_trim.fq.gz
cat SRR3301143_trim.fq.gz SRR3301144_trim.fq.gz SRR3301145_trim.fq.gz > 9786_trim.fq.gz
cat SRR3301146_trim.fq.gz SRR3301147_trim.fq.gz SRR3301148_trim.fq.gz > 9787_trim.fq.gz
cat SRR3301149_trim.fq.gz SRR3301150_trim.fq.gz > 9788_trim.fq.gz
cat SRR3301151_trim.fq.gz SRR3301152_trim.fq.gz > 9789_trim.fq.gz
cat SRR3301153_trim.fq.gz SRR3301154_trim.fq.gz SRR3301155_trim.fq.gz > 9790_trim.fq.gz
cat SRR3301156_trim.fq.gz SRR3301157_trim.fq.gz SRR3301158_trim.fq.gz > 9791_trim.fq.gz
cat SRR3301159_trim.fq.gz SRR3301160_trim.fq.gz SRR3301161_trim.fq.gz > 9792_trim.fq.gz
cat SRR3301162_trim.fq.gz SRR3301163_trim.fq.gz SRR3301164_trim.fq.gz > 9793_trim.fq.gz
cat SRR3301167_trim.fq.gz SRR3301168_trim.fq.gz > 9795_trim.fq.gz
cat SRR3301169_trim.fq.gz SRR3301170_trim.fq.gz SRR3301171_trim.fq.gz > 9796_trim.fq.gz
cat SRR3301172_trim.fq.gz SRR3301173_trim.fq.gz > 9797_trim.fq.gz
cat SRR3301174_trim.fq.gz SRR3301175_trim.fq.gz SRR3301176_trim.fq.gz > 9798_trim.fq.gz
cat SRR3301177_trim.fq.gz SRR3301178_trim.fq.gz > 9799_trim.fq.gz
cat SRR3301179_trim.fq.gz SRR3301180_trim.fq.gz > 9800_trim.fq.gz
cat SRR3301181_trim.fq.gz SRR3301182_trim.fq.gz SRR3301183_trim.fq.gz > 9801_trim.fq.gz
cat SRR3301184_trim.fq.gz SRR3301185_trim.fq.gz SRR3301186_trim.fq.gz > 9802_trim.fq.gz
cat SRR3301187_trim.fq.gz SRR3301188_trim.fq.gz SRR3301189_trim.fq.gz > 9803_trim.fq.gz
cat SRR3301190_trim.fq.gz SRR3301191_trim.fq.gz SRR3301192_trim.fq.gz > 9804_trim.fq.gz
cat SRR3301193_trim.fq.gz SRR3301194_trim.fq.gz > 9805_trim.fq.gz
cat SRR3301195_trim.fq.gz SRR3301196_trim.fq.gz > 9806_trim.fq.gz
cat SRR3301197_trim.fq.gz SRR3301198_trim.fq.gz > 9807_trim.fq.gz
cat SRR3301222_trim.fq.gz SRR3301223_trim.fq.gz > 9815_trim.fq.gz
cat SRR3301228_trim.fq.gz SRR3301229_trim.fq.gz > 9817_trim.fq.gz
cat SRR3301230_trim.fq.gz SRR3301231_trim.fq.gz > 9819_trim.fq.gz
cat SRR3301232_trim.fq.gz SRR3301233_trim.fq.gz > 9820_trim.fq.gz
cat SRR3301234_trim.fq.gz SRR3301235_trim.fq.gz > 9821_trim.fq.gz
cat SRR3301236_trim.fq.gz SRR3301237_trim.fq.gz SRR3301238_trim.fq.gz > 9822_trim.fq.gz
cat SRR3301239_trim.fq.gz SRR3301240_trim.fq.gz > 9823_trim.fq.gz
cat SRR3301241_trim.fq.gz SRR3301242_trim.fq.gz > 9824_trim.fq.gz
cat SRR3301243_trim.fq.gz SRR3301244_trim.fq.gz > 9825_trim.fq.gz
cat SRR3301245_trim.fq.gz SRR3301246_trim.fq.gz SRR3301247_trim.fq.gz > 9826_trim.fq.gz
cat SRR3301248_trim.fq.gz SRR3301249_trim.fq.gz > 9827_trim.fq.gz
cat SRR3301250_trim.fq.gz SRR3301251_trim.fq.gz > 9829_trim.fq.gz
cat SRR3301252_trim.fq.gz SRR3301253_trim.fq.gz > 9830_trim.fq.gz
cat SRR3301254_trim.fq.gz SRR3301255_trim.fq.gz > 9831_trim.fq.gz
cat SRR3301256_trim.fq.gz SRR3301257_trim.fq.gz SRR3301258_trim.fq.gz > 9832_trim.fq.gz
cat SRR3301259_trim.fq.gz SRR3301260_trim.fq.gz > 9833_trim.fq.gz
cat SRR3301261_trim.fq.gz SRR3301262_trim.fq.gz SRR3301263_trim.fq.gz > 9834_trim.fq.gz
cat SRR3301264_trim.fq.gz SRR3301265_trim.fq.gz SRR3301266_trim.fq.gz > 9835_trim.fq.gz
cat SRR3301267_trim.fq.gz SRR3301268_trim.fq.gz > 9836_trim.fq.gz
cat SRR3301269_trim.fq.gz SRR3301270_trim.fq.gz SRR3301271_trim.fq.gz > 9837_trim.fq.gz
cat SRR3301272_trim.fq.gz SRR3301273_trim.fq.gz > 9838_trim.fq.gz
cat SRR3301274_trim.fq.gz SRR3301275_trim.fq.gz > 9839_trim.fq.gz
cat SRR3301276_trim.fq.gz SRR3301277_trim.fq.gz SRR3301278_trim.fq.gz > 9840_trim.fq.gz
cat SRR3301279_trim.fq.gz SRR3301280_trim.fq.gz > 9841_trim.fq.gz
cat SRR3301281_trim.fq.gz SRR3301282_trim.fq.gz SRR3301283_trim.fq.gz SRR3301284_trim.fq.gz SRR3301285_trim.fq.gz SRR3301286_trim.fq.gz > 9843_trim.fq.gz
cat SRR3301287_trim.fq.gz SRR3301288_trim.fq.gz SRR3301289_trim.fq.gz SRR3301290_trim.fq.gz SRR3301291_trim.fq.gz SRR3301292_trim.fq.gz > 9844_trim.fq.gz
cat SRR3301293_trim.fq.gz SRR3301294_trim.fq.gz > 9845_trim.fq.gz
cat SRR3301295_trim.fq.gz SRR3301296_trim.fq.gz SRR3301297_trim.fq.gz SRR3301298_trim.fq.gz > 9846_trim.fq.gz
cat SRR3301299_trim.fq.gz SRR3301300_trim.fq.gz SRR3301301_trim.fq.gz SRR3301302_trim.fq.gz SRR3301303_trim.fq.gz > 9847_trim.fq.gz
cat SRR3301304_trim.fq.gz SRR3301305_trim.fq.gz > 9848_trim.fq.gz
cat SRR3301306_trim.fq.gz SRR3301307_trim.fq.gz > 9849_trim.fq.gz
cat SRR3301308_trim.fq.gz SRR3301309_trim.fq.gz SRR3301310_trim.fq.gz > 9850_trim.fq.gz
cat SRR3301311_trim.fq.gz SRR3301312_trim.fq.gz > 9851_trim.fq.gz
cat SRR3301313_trim.fq.gz SRR3301314_trim.fq.gz > 9852_trim.fq.gz
cat SRR3301315_trim.fq.gz SRR3301316_trim.fq.gz SRR3301317_trim.fq.gz SRR3301318_trim.fq.gz > 9853_trim.fq.gz
cat SRR3301319_trim.fq.gz SRR3301320_trim.fq.gz > 9854_trim.fq.gz
cat SRR3301321_trim.fq.gz SRR3301322_trim.fq.gz > 9855_trim.fq.gz
cat SRR3301323_trim.fq.gz SRR3301324_trim.fq.gz SRR3301325_trim.fq.gz > 9856_trim.fq.gz
cat SRR3301326_trim.fq.gz SRR3301327_trim.fq.gz SRR3301328_trim.fq.gz > 9857_trim.fq.gz
cat SRR3301329_trim.fq.gz SRR3301330_trim.fq.gz > 9858_trim.fq.gz
cat SRR3301331_trim.fq.gz SRR3301332_trim.fq.gz SRR3301333_trim.fq.gz > 9859_trim.fq.gz
cat SRR3301334_trim.fq.gz SRR3301335_trim.fq.gz > 9860_trim.fq.gz
cat SRR3301336_trim.fq.gz SRR3301337_trim.fq.gz > 9861_trim.fq.gz
cat SRR3301338_trim.fq.gz SRR3301339_trim.fq.gz > 9862_trim.fq.gz
cat SRR3301340_trim.fq.gz SRR3301341_trim.fq.gz > 9864_trim.fq.gz
cat SRR3301342_trim.fq.gz SRR3301343_trim.fq.gz SRR3301344_trim.fq.gz > 9866_trim.fq.gz
cat SRR3301345_trim.fq.gz SRR3301346_trim.fq.gz > 9867_trim.fq.gz
cat SRR3301347_trim.fq.gz SRR3301348_trim.fq.gz > 9868_trim.fq.gz
cat SRR3301349_trim.fq.gz SRR3301350_trim.fq.gz > 9869_trim.fq.gz
cat SRR3301351_trim.fq.gz SRR3301352_trim.fq.gz SRR3301353_trim.fq.gz SRR3301354_trim.fq.gz > 9870_trim.fq.gz
cat SRR3301355_trim.fq.gz SRR3301356_trim.fq.gz > 9871_trim.fq.gz
cat SRR3301357_trim.fq.gz SRR3301358_trim.fq.gz > 9872_trim.fq.gz
cat SRR3301359_trim.fq.gz SRR3301360_trim.fq.gz SRR3301361_trim.fq.gz > 9873_trim.fq.gz
cat SRR3301362_trim.fq.gz SRR3301363_trim.fq.gz > 9874_trim.fq.gz
cat SRR3301364_trim.fq.gz SRR3301365_trim.fq.gz > 9875_trim.fq.gz
cat SRR3301366_trim.fq.gz SRR3301367_trim.fq.gz SRR3301368_trim.fq.gz > 9876_trim.fq.gz
cat SRR3301369_trim.fq.gz SRR3301370_trim.fq.gz > 9877_trim.fq.gz
cat SRR3301371_trim.fq.gz SRR3301372_trim.fq.gz SRR3301373_trim.fq.gz > 9878_trim.fq.gz
cat SRR3301374_trim.fq.gz SRR3301375_trim.fq.gz SRR3301376_trim.fq.gz > 9879_trim.fq.gz
cat SRR3301377_trim.fq.gz SRR3301378_trim.fq.gz > 9880_trim.fq.gz
cat SRR3301379_trim.fq.gz SRR3301380_trim.fq.gz > 9881_trim.fq.gz
cat SRR3301381_trim.fq.gz SRR3301382_trim.fq.gz > 9882_trim.fq.gz
cat SRR3301383_trim.fq.gz SRR3301384_trim.fq.gz SRR3301385_trim.fq.gz > 9883_trim.fq.gz
cat SRR3301386_trim.fq.gz SRR3301387_trim.fq.gz SRR3301388_trim.fq.gz > 9885_trim.fq.gz
cat SRR3301389_trim.fq.gz SRR3301390_trim.fq.gz > 9886_trim.fq.gz
cat SRR3301391_trim.fq.gz SRR3301392_trim.fq.gz SRR3301393_trim.fq.gz SRR3301394_trim.fq.gz > 9888_trim.fq.gz
cat SRR3301395_trim.fq.gz SRR3301396_trim.fq.gz > 9890_trim.fq.gz
cat SRR3301397_trim.fq.gz SRR3301398_trim.fq.gz > 9891_trim.fq.gz
cat SRR3301399_trim.fq.gz SRR3301400_trim.fq.gz > 9892_trim.fq.gz
cat SRR3301401_trim.fq.gz SRR3301402_trim.fq.gz SRR3301403_trim.fq.gz > 9893_trim.fq.gz
cat SRR3301404_trim.fq.gz SRR3301405_trim.fq.gz SRR3301406_trim.fq.gz > 9894_trim.fq.gz
cat SRR3301407_trim.fq.gz SRR3301408_trim.fq.gz > 9895_trim.fq.gz
cat SRR3301409_trim.fq.gz SRR3301410_trim.fq.gz > 9897_trim.fq.gz
cat SRR3301411_trim.fq.gz SRR3301412_trim.fq.gz SRR3301413_trim.fq.gz > 9898_trim.fq.gz
cat SRR3301414_trim.fq.gz SRR3301415_trim.fq.gz SRR3301416_trim.fq.gz SRR3301417_trim.fq.gz SRR3301418_trim.fq.gz > 9899_trim.fq.gz
cat SRR3301419_trim.fq.gz SRR3301420_trim.fq.gz SRR3301421_trim.fq.gz > 9900_trim.fq.gz
cat SRR3301422_trim.fq.gz SRR3301423_trim.fq.gz SRR3301424_trim.fq.gz > 9901_trim.fq.gz
cat SRR3301425_trim.fq.gz SRR3301426_trim.fq.gz SRR3301427_trim.fq.gz SRR3301428_trim.fq.gz > 9902_trim.fq.gz
cat SRR3301429_trim.fq.gz SRR3301430_trim.fq.gz SRR3301431_trim.fq.gz > 9903_trim.fq.gz
cat SRR3301432_trim.fq.gz SRR3301433_trim.fq.gz > 9904_trim.fq.gz
cat SRR3301434_trim.fq.gz SRR3301435_trim.fq.gz SRR3301436_trim.fq.gz SRR3301437_trim.fq.gz > 9905_trim.fq.gz
cat SRR3301438_trim.fq.gz SRR3301439_trim.fq.gz > 9906_trim.fq.gz
cat SRR3301442_trim.fq.gz SRR3301443_trim.fq.gz SRR3301444_trim.fq.gz SRR3301445_trim.fq.gz > 9910_trim.fq.gz
cat SRR3301446_trim.fq.gz SRR3301447_trim.fq.gz > 9913_trim.fq.gz
cat SRR3301451_trim.fq.gz SRR3301452_trim.fq.gz > 9915_trim.fq.gz
cat SRR3301453_trim.fq.gz SRR3301454_trim.fq.gz SRR3301455_trim.fq.gz > 9916_trim.fq.gz
cat SRR3301456_trim.fq.gz SRR3301457_trim.fq.gz > 9917_trim.fq.gz
cat SRR3301458_trim.fq.gz SRR3301459_trim.fq.gz SRR3301460_trim.fq.gz > 9919_trim.fq.gz
cat SRR3301462_trim.fq.gz SRR3301463_trim.fq.gz > 9920_trim.fq.gz
cat SRR3301464_trim.fq.gz SRR3301465_trim.fq.gz > 9921_trim.fq.gz
cat SRR3301468_trim.fq.gz SRR3301469_trim.fq.gz > 9923_trim.fq.gz
cat SRR3301470_trim.fq.gz SRR3301471_trim.fq.gz > 9924_trim.fq.gz
cat SRR3301472_trim.fq.gz SRR3301473_trim.fq.gz > 9926_trim.fq.gz
cat SRR3301474_trim.fq.gz SRR3301475_trim.fq.gz SRR3301476_trim.fq.gz SRR3301477_trim.fq.gz SRR3301478_trim.fq.gz > 9928_trim.fq.gz
cat SRR3301479_trim.fq.gz SRR3301480_trim.fq.gz > 9929_trim.fq.gz
cat SRR3301481_trim.fq.gz SRR3301482_trim.fq.gz > 9931_trim.fq.gz
cat SRR3301483_trim.fq.gz SRR3301484_trim.fq.gz > 9932_trim.fq.gz
cat SRR3301485_trim.fq.gz SRR3301486_trim.fq.gz > 9933_trim.fq.gz
cat SRR3301487_trim.fq.gz SRR3301488_trim.fq.gz > 9935_trim.fq.gz
cat SRR3301489_trim.fq.gz SRR3301490_trim.fq.gz > 9936_trim.fq.gz
cat SRR3301492_trim.fq.gz SRR3301493_trim.fq.gz SRR3301494_trim.fq.gz > 9937_trim.fq.gz
cat SRR3301495_trim.fq.gz SRR3301496_trim.fq.gz SRR3301497_trim.fq.gz SRR3301498_trim.fq.gz > 9938_trim.fq.gz
cat SRR3301499_trim.fq.gz SRR3301500_trim.fq.gz SRR3301501_trim.fq.gz SRR3301502_trim.fq.gz > 9939_trim.fq.gz
cat SRR3301503_trim.fq.gz SRR3301504_trim.fq.gz > 9940_trim.fq.gz
cat SRR3301507_trim.fq.gz SRR3301508_trim.fq.gz SRR3301509_trim.fq.gz > 9942_trim.fq.gz
cat SRR3301510_trim.fq.gz SRR3301511_trim.fq.gz SRR3301512_trim.fq.gz > 9944_trim.fq.gz
cat SRR3301513_trim.fq.gz SRR3301514_trim.fq.gz SRR3301515_trim.fq.gz > 9945_trim.fq.gz
cat SRR3301516_trim.fq.gz SRR3301517_trim.fq.gz SRR3301518_trim.fq.gz > 9946_trim.fq.gz
cat SRR3301519_trim.fq.gz SRR3301520_trim.fq.gz SRR3301521_trim.fq.gz > 9947_trim.fq.gz
cat SRR3301522_trim.fq.gz SRR3301523_trim.fq.gz SRR3301524_trim.fq.gz > 9948_trim.fq.gz
cat SRR3301527_trim.fq.gz SRR3301528_trim.fq.gz > 9950_trim.fq.gz
cat SRR3301529_trim.fq.gz SRR3301530_trim.fq.gz SRR3301531_trim.fq.gz SRR3301532_trim.fq.gz > 9951_trim.fq.gz
cat SRR3301533_trim.fq.gz SRR3301534_trim.fq.gz > 9953_trim.fq.gz
cat SRR3301537_trim.fq.gz SRR3301538_trim.fq.gz SRR3301539_trim.fq.gz SRR3301540_trim.fq.gz > 9956_trim.fq.gz
cat SRR3301541_trim.fq.gz SRR3301542_trim.fq.gz > 9957_trim.fq.gz
cat SRR3301543_trim.fq.gz SRR3301544_trim.fq.gz SRR3301545_trim.fq.gz SRR3301546_trim.fq.gz > 9958_trim.fq.gz
cat SRR3301549_trim.fq.gz SRR3301550_trim.fq.gz > 9962_trim.fq.gz
cat SRR3301551_trim.fq.gz SRR3301552_trim.fq.gz > 9963_trim.fq.gz
cat SRR3301553_trim.fq.gz SRR3301554_trim.fq.gz SRR3301555_trim.fq.gz > 9964_trim.fq.gz
cat SRR3301558_trim.fq.gz SRR3301559_trim.fq.gz > 9966_trim.fq.gz
cat SRR3301560_trim.fq.gz SRR3301561_trim.fq.gz SRR3301562_trim.fq.gz > 9969_trim.fq.gz
cat SRR3301564_trim.fq.gz SRR3301565_trim.fq.gz SRR3301566_trim.fq.gz > 9973_trim.fq.gz
cat SRR3301567_trim.fq.gz SRR3301568_trim.fq.gz SRR3301569_trim.fq.gz SRR3301570_trim.fq.gz > 9976_trim.fq.gz
cat SRR3301571_trim.fq.gz SRR3301572_trim.fq.gz SRR3301573_trim.fq.gz > 997_trim.fq.gz
cat SRR3301574_trim.fq.gz SRR3301575_trim.fq.gz SRR3301576_trim.fq.gz > 9980_trim.fq.gz
cat SRR3301577_trim.fq.gz SRR3301578_trim.fq.gz > 9981_trim.fq.gz
cat SRR3301579_trim.fq.gz SRR3301580_trim.fq.gz > 9983_trim.fq.gz
cat SRR3301581_trim.fq.gz SRR3301582_trim.fq.gz > 9985_trim.fq.gz
cat SRR3301583_trim.fq.gz SRR3301584_trim.fq.gz > 9986_trim.fq.gz
cat SRR3301585_trim.fq.gz SRR3301586_trim.fq.gz SRR3301587_trim.fq.gz > 9987_trim.fq.gz
cat SRR3301588_trim.fq.gz SRR3301589_trim.fq.gz > 9989_trim.fq.gz
cat SRR3301590_trim.fq.gz SRR3301591_trim.fq.gz SRR3301592_trim.fq.gz > 9990_trim.fq.gz
cat SRR3301593_trim.fq.gz SRR3301594_trim.fq.gz > 9991_trim.fq.gz
cat SRR3301597_trim.fq.gz SRR3301598_trim.fq.gz SRR3301599_trim.fq.gz > 9996_trim.fq.gz
cat SRR3301600_trim.fq.gz SRR3301601_trim.fq.gz SRR3301602_trim.fq.gz > 9997_trim.fq.gz
