chr_num=$1

bcftools view input/gnomAD/SNP/rawdata/gnomad.genomes.v4.0.sites.chr${chr_num}.vcf.bgz | awk -F '\t' 'BEGIN { OFS = FS }
$0!~/#/{   af_afr = "-"
    ac_afr = "-"
    an_afr = "-"
    af_eas = "-"
    ac_eas = "-"
    an_eas = "-"
    af_nfe = "-"
    ac_nfe = "-"
    an_nfe = "-"
    split($8, info_arr, ";")
    for (i in info_arr) {
        if (info_arr[i] ~ /^AC_afr=/) {
            split(info_arr[i], tmp, "=")
            ac_afr = tmp[2]
        }
        else if (info_arr[i] ~ /^AF_afr=/) {
            split(info_arr[i], tmp, "=")
            af_afr = tmp[2]
        }
        else if (info_arr[i] ~ /^AN_afr=/) {
            split(info_arr[i], tmp, "=")
            an_afr = tmp[2]
        }
        else if (info_arr[i] ~ /^AC_eas=/) {
            split(info_arr[i], tmp, "=")
            ac_eas = tmp[2]
        }
        else if (info_arr[i] ~ /^AF_eas=/) {
            split(info_arr[i], tmp, "=")
            af_eas = tmp[2]
        }
        else if (info_arr[i] ~ /^AN_eas=/) {
            split(info_arr[i], tmp, "=")
            an_eas = tmp[2]
        }
        else if (info_arr[i] ~ /^AC_nfe=/) {
            split(info_arr[i], tmp, "=")
            ac_nfe = tmp[2]
        }
        else if (info_arr[i] ~ /^AF_nfe=/) {
            split(info_arr[i], tmp, "=")
            af_nfe = tmp[2]
        }
        else if (info_arr[i] ~ /^AN_nfe=/) {
            split(info_arr[i], tmp, "=")
            an_nfe = tmp[2]
        }
    }
    if(af_afr > 0.01 || af_eas > 0.01 || af_nfe > 0.01){
        print $1":"$2":"$4":"$5"\t"$3"\t"ac_afr"\t"an_afr"\t"af_afr"\t"ac_eas"\t"an_eas"\t"af_eas"\t"ac_nfe"\t"an_nfe"\t"af_nfe
    }
}' - > input/gnomAD/SNP/tmp/extracted_snps_chr${chr_num}.txt
