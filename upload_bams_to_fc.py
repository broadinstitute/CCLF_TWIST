
##################################################################
### Upload BAM to Google Cloud
### M. Cadosch
### 12/17
##################################################################
import os

##################################################################
# EDIT: Params
##################################################################
sn_id = "SN0176239"
tsca_id = "TSCA50"
date = "20190701"
##################################################################

target_dir_path = "/xchip/clf/seq_data/process_for_fc/%s__%s_%s"%(tsca_id.lower(), date, sn_id)

print("Uploading BAM files to Google Cloud...")
cmd = "gsutil -m cp -r %s/* gs://fc-35446f22-ea37-483a-bd6c-5e9fc56851ff/seq_data/%s/ "%(target_dir_path, tsca_id)
print(cmd)
res = os.system(cmd)
print(res)
