
#The following script demonstrates running the rate4site program to calculate conservation scores for a given protein sequence.
#INPUT: multiple sequence alignment. This MSA was generated for the seed sequence: input_seed.fasta
#OUTPUT: 

#Rate4site program can be downloaded from the following website:
# http://www.tau.ac.il/~itaymay/cp/rate4site.html

#Path of the binary:
CMD_MAYROSE=/home/yohan/workspace/positional_bias/experiments/20110808_conservation_profile/src/rate4site/rate4site.3.2.source_temp/sourceMar09/rate4site_fast

#Id of the reference sequence against which residue position-specific conservation scores will be calculated.
ID_REF_SEQUENCE="ref|NP_570292.1|_1_342"


FNAME_MSA=msa.aln
FNAME_OUTPUT=conservation_scores.txt


$CMD_MAYROSE -a $ID_REF_SEQUENCE -s $FNAME_MSA -o $FNAME_OUTPUT
