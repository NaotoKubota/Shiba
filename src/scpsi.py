import warnings
warnings.simplefilter('ignore')
import argparse

def get_args():
    ## Get arguments from command line

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "PSI calculation for alternative splicing events in scRNA-seq data"
    )

    parser.add_argument("junctions", type = str, help = "A bed file of Junction read counts generated by bam2junc.sh")
    parser.add_argument("event", type = str, help = "Directory that contains text files of alternative splicing events generated by gtf2event.py")
    parser.add_argument("output", type = str, help = "Directory for output files")
    parser.add_argument("-p", "--num-process", type = int, help = "Number of processors to use", default = 1)
    parser.add_argument("-f", "--fdr", type = float, help = "FDR for detecting differential events", default = 0.05)
    parser.add_argument("-d", "--psi", type = float, help = "Threshold of delta PSI for detecting differential events", default = 0.1)
    parser.add_argument("-r", "--reference", type = str, help = "Reference group for detecting differential events")
    parser.add_argument("-a", "--alternative", type = str, help = "Alternative group for detecting differential events")
    parser.add_argument("-m", "--minimum-reads", type = int, help = "Minumum value of total reads for each junction for detecting differential events", default = 10)
    parser.add_argument("--onlypsi", help = "Just calculate PSI for each sample, not perform statistical tests", action = 'store_true')
    parser.add_argument("--excel", help = "Make result files in excel format", action = 'store_true')

    args = parser.parse_args()

    return(args)


def main():
    ## Main

    args = get_args()
    junction_path = args.junctions
    event_path = args.event
    output_path = args.output
    num_process = args.num_process
    FDR = args.fdr
    dPSI = args.psi
    reference = args.reference
    alternative = args.alternative
    minimum_reads = args.minimum_reads
    onlypsi = args.onlypsi
    excel = args.excel

    # Load modules
    import sys
    import os
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
    sys.path.append(parent_dir)
    import time
    import pandas as pd
    from styleframe import StyleFrame, Styler, utils
    from lib import shibalib

    # Read events
    print("Load event files....", file = sys.stdout)
    event_df_dict = shibalib.read_events_sc(event_path)

    # Read junction count
    print("Load junction files....", file = sys.stdout)
    junc_df = shibalib.read_junctions(junction_path)

    # Dictionary for junction read counts from each sample
    junc_dict_all = shibalib.junc_dict(junc_df)

    # Get sample list
    sample_list = shibalib.make_sample_list(junc_df)
    # print("Samples: " + str(sample_list), file = sys.stdout)

    # Get group list
    if onlypsi == False:

        try:

            print(reference + " vs. " + alternative, file = sys.stdout)
            group_list = [reference, alternative]
            junc_group_df = junc_df[["chr", "start", "end", "ID"] + group_list]
            junc_dict_group = shibalib.junc_dict(junc_group_df)

        except KeyError:

            print("Error: " + reference + " or " + alternative + " is not in the sample list", file = sys.stderr)
            sys.exit()

    # Set of junctions
    junc_set = shibalib.make_junc_set(junc_df)

    ################################################################################
    # Skipped exon
    print("PSI for each skipped exon....", file = sys.stdout)

    event_for_analysis_df = shibalib.event_for_analysis_se(event_df_dict["SE"], junc_set)

    if onlypsi == False:

        psi_table_df = shibalib.make_psi_table_group(
            group_list, event_for_analysis_df, junc_dict_group, shibalib.se,
            shibalib.col_se, num_process, minimum_reads
        )

    else:

        psi_table_df = shibalib.make_psi_table_sample(
            sample_list, event_for_analysis_df, junc_dict_all, shibalib.se,
            shibalib.col_se, num_process, minimum_reads
        )

    if onlypsi == False:

        SE_df = shibalib.diff_event(
            event_for_analysis_df, psi_table_df, junc_dict_all,
            False, group_list, sample_list,
            shibalib.diff_se, shibalib.se_ind, num_process,
            FDR, dPSI, False, False
        )

    else:

        SE_df, output_mtx_SE_df = shibalib.make_psi_mtx(psi_table_df)

    ##############################################################################
    # Alternative five ss
    print("PSI for each alternative five prime ss....", file = sys.stdout)

    event_for_analysis_df = shibalib.event_for_analysis_five_three_afe_ale(event_df_dict["FIVE"], junc_set)

    if onlypsi == False:

        psi_table_df = shibalib.make_psi_table_group(
            group_list, event_for_analysis_df, junc_dict_group, shibalib.five_three_afe_ale,
            shibalib.col_five_three_afe_ale, num_process, minimum_reads
        )

    else:

        psi_table_df = shibalib.make_psi_table_sample(
            sample_list, event_for_analysis_df, junc_dict_all, shibalib.five_three_afe_ale,
            shibalib.col_five_three_afe_ale, num_process, minimum_reads
        )

    if onlypsi == False:

        FIVE_df = shibalib.diff_event(
            event_for_analysis_df, psi_table_df, junc_dict_all,
            False, group_list, sample_list,
            shibalib.diff_five_three_afe_ale, shibalib.five_three_afe_ale_ind, num_process,
            FDR, dPSI, False, False
        )

    else:

        FIVE_df, output_mtx_FIVE_df = shibalib.make_psi_mtx(psi_table_df)

    ##############################################################################
    # Alternative three ss
    print("PSI for each alternative three prime ss....", file = sys.stdout)

    event_for_analysis_df = shibalib.event_for_analysis_five_three_afe_ale(event_df_dict["THREE"], junc_set)

    if onlypsi == False:

        psi_table_df = shibalib.make_psi_table_group(
            group_list, event_for_analysis_df, junc_dict_group, shibalib.five_three_afe_ale,
            shibalib.col_five_three_afe_ale, num_process, minimum_reads
        )

    else:

        psi_table_df = shibalib.make_psi_table_sample(
            sample_list, event_for_analysis_df, junc_dict_all, shibalib.five_three_afe_ale,
            shibalib.col_five_three_afe_ale, num_process, minimum_reads
        )

    if onlypsi == False:

        THREE_df = shibalib.diff_event(
            event_for_analysis_df, psi_table_df, junc_dict_all,
            False, group_list, sample_list,
            shibalib.diff_five_three_afe_ale, shibalib.five_three_afe_ale_ind, num_process,
            FDR, dPSI, False, False
        )

    else:

        THREE_df, output_mtx_THREE_df = shibalib.make_psi_mtx(psi_table_df)

    ##############################################################################
    # Mutually exclusive exons
    print("PSI for each mutually exclusive exons....", file = sys.stdout)

    event_for_analysis_df = shibalib.event_for_analysis_mxe(event_df_dict["MXE"], junc_set)

    if onlypsi == False:

        psi_table_df = shibalib.make_psi_table_group(
            group_list, event_for_analysis_df, junc_dict_group, shibalib.mxe,
            shibalib.col_mxe, num_process, minimum_reads
        )

    else:

        psi_table_df = shibalib.make_psi_table_sample(
            sample_list, event_for_analysis_df, junc_dict_all, shibalib.mxe,
            shibalib.col_mxe, num_process, minimum_reads
        )

    if onlypsi == False:

        MXE_df = shibalib.diff_event(
            event_for_analysis_df, psi_table_df, junc_dict_all,
            False, group_list, sample_list,
            shibalib.diff_mxe, shibalib.mxe_ind, num_process,
            FDR, dPSI, False, False
        )

    else:

        MXE_df, output_mtx_MXE_df = shibalib.make_psi_mtx(psi_table_df)

    ################################################################################
    # Multiple skipped exons
    print("PSI for each multiple skipped exons....", file = sys.stdout)

    event_for_analysis_df = shibalib.event_for_analysis_mse(event_df_dict["MSE"], junc_set)

    if onlypsi == False:

        psi_table_df = shibalib.make_psi_table_group(
            group_list, event_for_analysis_df, junc_dict_group, shibalib.mse,
            shibalib.col_mse, num_process, minimum_reads
        )

    else:

        psi_table_df = shibalib.make_psi_table_sample(
            sample_list, event_for_analysis_df, junc_dict_all, shibalib.mse,
            shibalib.col_mse, num_process, minimum_reads
        )

    if onlypsi == False:

        MSE_df = shibalib.diff_event(
            event_for_analysis_df, psi_table_df, junc_dict_all,
            False, group_list, sample_list,
            shibalib.diff_mse, shibalib.mse_ind, num_process,
            FDR, dPSI, False, False
        )

    else:

        MSE_df, output_mtx_MSE_df = shibalib.make_psi_mtx(psi_table_df)

    ################################################################################
    # Alternative first exon
    print("PSI for each alternative first exon....", file = sys.stdout)

    event_for_analysis_df = shibalib.event_for_analysis_five_three_afe_ale(event_df_dict["AFE"], junc_set)

    if onlypsi == False:

        psi_table_df = shibalib.make_psi_table_group(
            group_list, event_for_analysis_df, junc_dict_group, shibalib.five_three_afe_ale,
            shibalib.col_five_three_afe_ale, num_process, minimum_reads
        )

    else:

        psi_table_df = shibalib.make_psi_table_sample(
            sample_list, event_for_analysis_df, junc_dict_all, shibalib.five_three_afe_ale,
            shibalib.col_five_three_afe_ale, num_process, minimum_reads
        )

    if onlypsi == False:

        AFE_df = shibalib.diff_event(
            event_for_analysis_df, psi_table_df, junc_dict_all,
            False, group_list, sample_list,
            shibalib.diff_five_three_afe_ale, shibalib.five_three_afe_ale_ind, num_process,
            FDR, dPSI, False, False
        )

    else:

        AFE_df, output_mtx_AFE_df = shibalib.make_psi_mtx(psi_table_df)

    ################################################################################
    # Alternative last exon
    print("PSI for each alternative last exon....", file = sys.stdout)

    event_for_analysis_df = shibalib.event_for_analysis_five_three_afe_ale(event_df_dict["ALE"], junc_set)

    if onlypsi == False:

        psi_table_df = shibalib.make_psi_table_group(
            group_list, event_for_analysis_df, junc_dict_group, shibalib.five_three_afe_ale,
            shibalib.col_five_three_afe_ale, num_process, minimum_reads
        )

    else:

        psi_table_df = shibalib.make_psi_table_sample(
            sample_list, event_for_analysis_df, junc_dict_all, shibalib.five_three_afe_ale,
            shibalib.col_five_three_afe_ale, num_process, minimum_reads
        )

    if onlypsi == False:

        ALE_df = shibalib.diff_event(
            event_for_analysis_df, psi_table_df, junc_dict_all,
            False, group_list, sample_list,
            shibalib.diff_five_three_afe_ale, shibalib.five_three_afe_ale_ind, num_process,
            FDR, dPSI, False, False
        )

    else:

        ALE_df, output_mtx_ALE_df = shibalib.make_psi_mtx(psi_table_df)

    ################################################################################
    # Make output directory
    os.makedirs(output_path, exist_ok = True)

    if onlypsi: # Simple PSI matrix

        simple_psi_df = pd.concat(

            [output_mtx_SE_df, output_mtx_FIVE_df, output_mtx_THREE_df, output_mtx_MXE_df, output_mtx_MSE_df, output_mtx_AFE_df, output_mtx_ALE_df],

        )

        simple_psi_df.to_csv(

            output_path + "/PSI_matrix.txt",
            sep = "\t",
            index = False,

        )

        del simple_psi_df

    else: # Summary file (up-regulated and down-regulated exons)

        event_counter_SE = shibalib.EventCounter(SE_df, dPSI)
        event_counts_SE = event_counter_SE.count_all_events()
        event_counter_FIVE = shibalib.EventCounter(FIVE_df, dPSI)
        event_counts_FIVE = event_counter_FIVE.count_all_events()
        event_counter_THREE = shibalib.EventCounter(THREE_df, dPSI)
        event_counts_THREE = event_counter_THREE.count_all_events()
        event_counter_MXE = shibalib.EventCounter(MXE_df, dPSI)
        event_counts_MXE = event_counter_MXE.count_all_events()
        event_counter_MSE = shibalib.EventCounter(MSE_df, dPSI)
        event_counts_MSE = event_counter_MSE.count_all_events()
        event_counter_AFE = shibalib.EventCounter(AFE_df, dPSI)
        event_counts_AFE = event_counter_AFE.count_all_events()
        event_counter_ALE = shibalib.EventCounter(ALE_df, dPSI)
        event_counts_ALE = event_counter_ALE.count_all_events()

        summary_l = [
            ["SE", "Up", "annotated", event_counts_SE["up_annotated_num"]],
            ["SE", "Down", "annotated", event_counts_SE["down_annotated_num"]],
            ["SE", "Up", "unannotated", event_counts_SE["up_unannotated_num"]],
            ["SE", "Down", "unannotated", event_counts_SE["down_unannotated_num"]],
            ["FIVE", "Up", "annotated", event_counts_FIVE["up_annotated_num"]],
            ["FIVE", "Down", "annotated", event_counts_FIVE["down_annotated_num"]],
            ["FIVE", "Up", "unannotated", event_counts_FIVE["up_unannotated_num"]],
            ["FIVE", "Down", "unannotated", event_counts_FIVE["down_unannotated_num"]],
            ["THREE", "Up", "annotated", event_counts_THREE["up_annotated_num"]],
            ["THREE", "Down", "annotated", event_counts_THREE["down_annotated_num"]],
            ["THREE", "Up", "unannotated", event_counts_THREE["up_unannotated_num"]],
            ["THREE", "Down", "unannotated", event_counts_THREE["down_unannotated_num"]],
            ["MXE", "Up", "annotated", event_counts_MXE["up_annotated_num"]],
            ["MXE", "Down", "annotated", event_counts_MXE["down_annotated_num"]],
            ["MXE", "Up", "unannotated", event_counts_MXE["up_unannotated_num"]],
            ["MXE", "Down", "unannotated", event_counts_MXE["down_unannotated_num"]],
            ["MSE", "Up", "annotated", event_counts_MSE["up_annotated_num"]],
            ["MSE", "Down", "annotated", event_counts_MSE["down_annotated_num"]],
            ["MSE", "Up", "unannotated", event_counts_MSE["up_unannotated_num"]],
            ["MSE", "Down", "unannotated", event_counts_MSE["down_unannotated_num"]],
            ["AFE", "Up", "annotated", event_counts_AFE["up_annotated_num"]],
            ["AFE", "Down", "annotated", event_counts_AFE["down_annotated_num"]],
            ["AFE", "Up", "unannotated", event_counts_AFE["up_unannotated_num"]],
            ["AFE", "Down", "unannotated", event_counts_AFE["down_unannotated_num"]],
            ["ALE", "Up", "annotated", event_counts_ALE["up_annotated_num"]],
            ["ALE", "Down", "annotated", event_counts_ALE["down_annotated_num"]],
            ["ALE", "Up", "unannotated", event_counts_ALE["up_unannotated_num"]],
            ["ALE", "Down", "unannotated", event_counts_ALE["down_unannotated_num"]],
        ]

        summary_df = pd.DataFrame(

            summary_l,
            columns = ["AS", "Direction", "Label", "Number"]

        )

        summary_df.to_csv(

            output_path + "/summary.txt",
            sep = "\t",
            index = False,

        )

    # Write to a file

    SE_df.to_csv(

        output_path + "/PSI_SE.txt",
        sep = "\t",
        index = False

    )

    FIVE_df.to_csv(

        output_path + "/PSI_FIVE.txt",
        sep = "\t",
        index = False

    )

    THREE_df.to_csv(

        output_path + "/PSI_THREE.txt",
        sep = "\t",
        index = False

    )

    MXE_df.to_csv(

        output_path + "/PSI_MXE.txt",
        sep = "\t",
        index = False

    )

    MSE_df.to_csv(

        output_path + "/PSI_MSE.txt",
        sep = "\t",
        index = False

    )

    AFE_df.to_csv(

        output_path + "/PSI_AFE.txt",
        sep = "\t",
        index = False

    )

    ALE_df.to_csv(

        output_path + "/PSI_ALE.txt",
        sep = "\t",
        index = False

    )

    ################################################################################
    # Excel file
    if excel:

        print("Export to an excel file....", file = sys.stdout)
        shibalib.save_excel_sc(output_path, SE_df, FIVE_df, THREE_df, MXE_df, MSE_df, AFE_df, ALE_df)

    ################################################################################

    print("Output file: " + output_path, file = sys.stdout)


if __name__ == '__main__':

    main()
