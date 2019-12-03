import os
import matplotlib
matplotlib.use('Agg')
import argparse
import pandas as pd
from jinja2 import Environment, FileSystemLoader

############################################
# Create batch report
# Inputs:
#   See args
# Outputs:
# 
# 
# mcadosch 1/18
############################################


#######################################################################################################
# 0. Parse command line inputs
#######################################################################################################
parser = argparse.ArgumentParser(description='Plot Results from Somatic CNV Calls')
parser.add_argument('--tsca_id', type=str, required=True,
                    help='Batch ID')
parser.add_argument('--sample_ids', type=str, required=True, nargs='+', 
                    help='List of sample ids')
parser.add_argument('--external_ids', type=str, required=True, nargs='+', 
                    help='List of external validation ids')
parser.add_argument('--participant_ids', type=str, required=True, nargs='+', 
                    help='List of participant ids')
parser.add_argument('--mean_interval_cvg', type=float, required=True, nargs='+',
                    help='List of mean_interval_cvg files')
parser.add_argument('--mouse_qc_results', type=str, required=True, nargs='+',
                    help='List of mouse_qc_results .txt files')
parser.add_argument('--fng_should_match_file', type=str, required=True,
                    help='Fingerprinting file of samples that should match but do not')
parser.add_argument('--fng_should_not_match_file', type=str, required=True,
                    help='Fingerprinting file of samples that should not match but do')
parser.add_argument('--report_template', type=str, required=True,
                    help='Report html template')



args = parser.parse_args()

def assert_inputs(sample_ids, external_ids, participant_ids, mean_interval_cvg, mouse_qc_results):
    """Ensure that all inputs contain information for the same number of samples
    Returns:
        - True / False if inputs are correct
    """
    print( "Number of external_ids: %s \n Number of sample_ids: %s \n Number of mean interval coverage: %s \
                     \n Number of mouseq QCs results: %s \n Number of participant_ids: %s" \
                        %(len(external_ids), len(sample_ids), len(mean_interval_cvg), len(mouse_qc_results), len(participant_ids)))
    
    if len(external_ids) == len(sample_ids) == len(mean_interval_cvg) == len(mouse_qc_results) == len(participant_ids):
        print("Same number of sample_ids, external_ids, mean_interval_cvg, and mouse_qc_results results passed. Inputs are correct...")
        return True
    else:
        print("Inputs have different lengths. One of the inputs is wrong...")
        return False

def prepare_fingerprinting_results(fng_should_match_file, fng_should_not_match_file, external_ids, participant_ids):
    """Prepare tables to display fingerprinting results
    """
    fng_should_match = pd.read_table(fng_should_match_file)
    fng_should_not_match = pd.read_table(fng_should_not_match_file)

    # Only keep snps with 50 reads or more
    fng_should_match = fng_should_match[fng_should_match["n_powered_snps"] >= 50]
    fng_should_not_match = fng_should_not_match[fng_should_not_match["n_powered_snps"] >= 50]

    ### Further filter samples that should not match. Some samples refer to the same participant, but previous code fails to filter those out
    samples = pd.DataFrame({"participant": participant_ids, "external_id_validation": external_ids})
    # Use code to cross-reference rows
    fng_should_not_match["code"] = fng_should_not_match.index
    
    # Queries and database
    fng_should_not_match_queries = fng_should_not_match[["sample.query", "code"]]
    fng_should_not_match_database = fng_should_not_match[["sample.database", "code"]]
    
    # Add participant_id to queries and database
    fng_should_not_match_queries_with_participant_id = \
        pd.merge(fng_should_not_match_queries, samples, left_on="sample.query", right_on="external_id_validation", how="left")
    fng_should_not_match_database_with_participant_id = \
        pd.merge(fng_should_not_match_database, samples, left_on="sample.database", right_on="external_id_validation", how="left")

    # Merge on code, and check if participants are the same
    fng_should_not_match_query_database = \
        pd.merge(fng_should_not_match_queries_with_participant_id, fng_should_not_match_database_with_participant_id, \
                on="code", suffixes=["_query", "_database"])
    fng_should_not_match_query_database['same_participant'] = \
        fng_should_not_match_query_database['participant_query'] == fng_should_not_match_query_database['participant_database']

    participant_ids_dont_match = \
        fng_should_not_match_query_database.loc[fng_should_not_match_query_database["same_participant"] == False, "code"].values

    fng_should_not_match_filtered = fng_should_not_match[fng_should_not_match["code"].isin(participant_ids_dont_match)]
    ### End of filtering

    # Dictionaries to display as tables in HTML
    samples_should_match = fng_should_match.to_dict(orient="records")
    samples_should_not_match = fng_should_not_match_filtered.to_dict(orient="records")

    return samples_should_match, samples_should_not_match

def create_report(report_template, sample_ids, external_ids, mean_interval_cvg, mouse_qc_results, tsca_id, samples_should_match, samples_should_not_match):
    """Create batch report
    """
    # Gather data into pandas df
    df = pd.DataFrame({ "sample_id": sample_ids, "external_id": external_ids, \
                        "mean_interval_cvg": mean_interval_cvg, "mouse_qc_results": mouse_qc_results})

    # Threshold for min average coverage for QC
    THRESHOLD = 50
    df['depth_of_cov_qc_results'] = df['mean_interval_cvg'] > THRESHOLD

    # Turn into dictionary for displaying as a table in HTML
    sample_results = df.to_dict(orient="records")

    # Render template
    THIS_DIR = os.path.dirname(os.path.abspath(__file__))
    env = Environment(loader=FileSystemLoader(THIS_DIR),trim_blocks=True)
    # template = env.get_template('/TSCA/report_template.html')
    # template = env.get_template('report_template.html')
    template  = env.get_template(report_template)
    report_html = template.render(title='CCLF batch report for %s'%tsca_id, \
                                    samples=sample_results, tsca_id=tsca_id, \
                                    samples_should_match=samples_should_match, samples_should_not_match=samples_should_not_match)
    
    # to save the results
    with open("report.html", "w") as f:
        f.write(report_html)
    return

####################################################################################################
# MAIN
####################################################################################################
if __name__ == '__main__':
    print("Running create_batch_report.py")
    inputs_correct = assert_inputs(args.sample_ids, args.external_ids, args.participant_ids, args.mean_interval_cvg, args.mouse_qc_results)
    if not inputs_correct:
        import sys; sys.exit()

    samples_should_match, samples_should_not_match = \
        prepare_fingerprinting_results(args.fng_should_match_file, args.fng_should_not_match_file, args.external_ids, args.participant_ids)

    create_report(args.report_template, args.sample_ids, args.external_ids, args.mean_interval_cvg, args.mouse_qc_results, args.tsca_id, \
                samples_should_match, samples_should_not_match)


