import os


def check_progression(res_path, input_path):
    all_done = [ result.split('_aggregation.txt')[0] for result in os.listdir(res_path) ]
    subset_done = [ result for result in all_done if result in os.listdir(input_path) ]
    prcnt = round(len(subset_done)/len(os.listdir(input_path))*100)
    print(f'{prcnt}% completed for {input_path.split("/")[-1]}')

if __name__ == '__main__':
    res_path = '/media/savvy/DATA3/savvy/project_2018/computational_mutagenesis/FINAL_RESULTS'
    # MM_input_path = '/media/savvy/DATA3/savvy/project_2018/computational_mutagenesis/MM_CHAP_CLT'
    # HG_input_path = '/media/savvy/DATA3/savvy/project_2018/computational_mutagenesis/HG_CHAP_CLT'
    # MM_input_path = '/media/savvy/DATA3/savvy/project_2018/computational_mutagenesis/MM_OTHERS'
    # HG_input_path = '/media/savvy/DATA3/savvy/project_2018/computational_mutagenesis/HG_OTHERS'
    RERUN_path = '/media/savvy/DATA3/savvy/project_2018/computational_mutagenesis/RERUN'

    # check_progression(res_path, MM_input_path)
    # check_progression(res_path, HG_input_path)
    check_progression(res_path, RERUN_path)
