import argparse
import numpy as np
import pandas as pd

############################################
# Determine if sample is mouse-contaminated or not
# Contaminated if more than 8 snps with mouse allelic fraction >= 75%
# mcadosch 01/18
############################################


#######################################################################################################
# 0. Parse command line inputs
#######################################################################################################
parser = argparse.ArgumentParser(description='Determine mousecontamination')
parser.add_argument('--mouse_qc_file', type=str, required=True,
                    help='result from mouse qc')

args = parser.parse_args()
mouse_qc_file                 = args.mouse_qc_file

#######################################################################################################
# 1. Analysis
#######################################################################################################
def determine_contamination_status(file):
    """Determine contamination status. If more than 8 snps with mouse AF >= 75%, then sample is mouse-contaminated
    """
    data = pd.read_table(file)
    out_file = open("mouse_qc_status.txt", "w")
    if (data[data['prelim_call']=="mouse_contamination"].shape[0] >= 8):
        out_file.write("contaminated")
    else:
        out_file.write("not_contaminated")
    out_file.close()
    return


####################################################################################################
# MAIN
####################################################################################################
if __name__ == '__main__':
    determine_contamination_status(mouse_qc_file)
