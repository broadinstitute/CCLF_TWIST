import argparse
import pysam

############################################
# Annotate vcf file with mutect version used
# Inputs:
#  See args
# Outputs:
#	${vcf_file}
# mcadosch 9/17
############################################

parser = argparse.ArgumentParser(description='Annoatate MuTect Version')
parser.add_argument('--vcf_file', type=str, required=True,
                    help='vcf file')
parser.add_argument('--mutect_version', type=str, required=True,
                    help='MuTect version used')
args = parser.parse_args()

### Read the input file
in_vcf=pysam.VariantFile(args.vcf_file,"r")

### Add the MuTect_version field to header
in_vcf.header.formats.add("MuTect_version",".","String","MuTect version used for variant call")

### Create output vcf file, with MuTect version field
out_vcf = pysam.VariantFile(args.vcf_file, 'w', header=in_vcf.header)

with open(args.vcf_file, "a") as out:
    for variant in in_vcf:
        for sample in variant.samples:
            variant.samples[sample]['MuTect_version'] = args.mutect_version
            out.write(variant)
