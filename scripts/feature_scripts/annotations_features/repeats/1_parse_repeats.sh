#!/bin/bash
module load tabix

#gunzip -d SL5.0.fa.mod.out.gz
#bgzip -d SL5.0.repeat.gff3.gz # I compressed this

# A script to document the parsing I used
# Change the Helitron to DNA/Helitron
awk '{if ($11 == "Helitron") $11 = "DNA/Helitron"; print}' SL5.0.fa.mod.out > SL5.0.fa.mod.out2

# Separate into family and class
awk '{
  if ($11 ~ "/") {
    split($11, arr, "/");
    $11 = arr[1];
    $12 = arr[2];
  } else {
    $12 = "other";
  }

  print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15$16;
}' SL5.0.fa.mod.out2 > SL5.0.fa.mod.out3

# We want to keep DNA and LTR families. The thing is LTR has an 'unknown' family
# We want to differentiate that with the other families which will be removed. So I will name it ltrunknown
awk '{if ($12 == "unknown") $12 = "ltrunknown"; print}' SL5.0.fa.mod.out3 > SL5.0.fa.mod.out4

# Now, we want to make any family that is not DNA or LTR be other.
# These are EnSpm_CACTA, hAT, MuDR_Mutator, PIF_Harbinger
# The condition has to include all lines where $11 is not DNA or LTR
awk '{
  if (!($11 ~ /DNA/ || $11 ~ /LTR/)) {
    $12 = "other";
  }
}
{
  print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15$16;
}' SL5.0.fa.mod.out4 > SL5.0.fa.mod.out5

# Now I wish to change the value of the 11th column to 'Unknown' if the value does not contain one of the following:DNA, LTR, Simple_repeat, MITE
awk '{
  if (!($11 ~ /DNA/ || $11 ~ /LTR/ || $11 ~ /Simple_repeat/ || $11 ~ /MITE/)) {
    $11 = "Unknown";
  }
}
{
  print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15$16;
}' SL5.0.fa.mod.out5 > SL5.0.fa.mod.out6

# Header is a bit fucked, whatever

# Now to extract usefull information, I need column numbers 1, 5, 6, 7, 11, 12
# Adding header
touch parsed_repeats.mod
echo -e "#SW_Score\tchr\tstart\tend\tclass\tfamily" > parsed_repeats.mod
awk '{print $5"\t"$6"\t"$7"\t"$11"\t"$12}' SL5.0.fa.mod.out6 >> parsed_repeats.mod

# Split to chromosomes
awk 'NR > 4 {print > "parsed_repeats_chr_"$1".mod"}' parsed_repeats.mod

rm SL5.0.fa.mod.out2 SL5.0.fa.mod.out3 SL5.0.fa.mod.out4 SL5.0.fa.mod.out5 SL5.0.fa.mod.out6 parsed_repeats.mod