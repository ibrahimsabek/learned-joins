import os, sys
import collections
from parsed_benchmark_output import ParsedLearnedSortMergeBenchmarkOutput

class LearnedSortMergeJoinAnalyzer:

    def analyze_datasets_sort_results(self, datasets_path):
        files = os.listdir(datasets_path)

        results_dict = {}

        for file in files:
            #print(datasets_path+file)
            parsed_obj = ParsedLearnedSortMergeBenchmarkOutput(datasets_path+file)
            parsed_obj.parse_benchmark_parameters()
            parsed_obj.parse_benchmark_output()

            curr_result_configs = {}
            curr_result_configs["datasets_r"] = parsed_obj.datasets_r
            curr_result_configs["datasets_s"] = parsed_obj.datasets_s
            curr_result_configs["threads"] = parsed_obj.threads
            curr_result_configs["use_avxsort_for_sorting_minor_bckts"] = parsed_obj.use_avxsort_for_sorting_minor_bckts
            curr_result_configs["ls_default_threshold"] = parsed_obj.ls_default_threshold
            curr_result_configs["ls_default_arch"] = parsed_obj.ls_default_arch
            curr_result_configs["ls_imv_avx"] = parsed_obj.ls_imv_avx
            curr_result_configs["ls_prefetch_minor_bckt_sizes_off"] = parsed_obj.ls_prefetch_minor_bckt_sizes_off
            curr_result_configs["ls_prefetch_slopes_intercepts_minor_bckts"] = parsed_obj.ls_prefetch_slopes_intercepts_minor_bckts
            curr_result_configs["ls_simdstate"] = parsed_obj.ls_simdstate
            curr_result_configs["ls_pdis"] = parsed_obj.ls_pdis

            results_dict[parsed_obj.total_algorithm_time_in_ms] = curr_result_configs

        sorted_results_dict = collections.OrderedDict(sorted(results_dict.items()))
        for key, value in sorted_results_dict.items():
            print(f'total_algorithm_time_in_ms: {key}, r_ds: {value["datasets_r"]}, s_ds: {value["datasets_s"]}, th: {value["threads"]}, avxmb: {value["use_avxsort_for_sorting_minor_bckts"]}, lsdt: {value["ls_default_threshold"]}, lsda: {value["ls_default_arch"]}, lsimv: {value["ls_imv_avx"]}, lsps: {value["ls_prefetch_minor_bckt_sizes_off"]}, lspsi: {value["ls_prefetch_slopes_intercepts_minor_bckts"]}, lsss: {value["ls_simdstate"]}, lspdis: {value["ls_pdis"]}')

    def analyze_datasets_get_min(self, datasets_path):
        files = os.listdir(datasets_path)

        i = 0
        min_total_algorithm_time_in_ms = 0
        min_r_dataset = ''
        min_s_dataset = ''
        min_threads = 0
        min_use_avxsort_for_sorting_minor_bckts = 0 
        min_ls_default_threshold = 0
        min_ls_default_arch = 0
        min_ls_imv_avx = 0
        min_ls_prefetch_minor_bckt_sizes_off = 0
        min_ls_prefetch_slopes_intercepts_minor_bckts = 0
        min_ls_simdstate = 0
        min_ls_pdis = 0

        for file in files:
            #print(datasets_path+file)
            parsed_obj = ParsedLearnedSortMergeBenchmarkOutput(datasets_path+file)
            parsed_obj.parse_benchmark_parameters()
            parsed_obj.parse_benchmark_output()

            if i == 0:
                min_total_algorithm_time_in_ms = parsed_obj.total_algorithm_time_in_ms
                min_r_dataset = parsed_obj.datasets_r
                min_s_dataset = parsed_obj.datasets_s
                min_threads = parsed_obj.threads
                min_use_avxsort_for_sorting_minor_bckts = parsed_obj.use_avxsort_for_sorting_minor_bckts 
                min_ls_default_threshold = parsed_obj.ls_default_threshold
                min_ls_default_arch = parsed_obj.ls_default_arch
                min_ls_imv_avx = parsed_obj.ls_imv_avx
                min_ls_prefetch_minor_bckt_sizes_off = parsed_obj.ls_prefetch_minor_bckt_sizes_off
                min_ls_prefetch_slopes_intercepts_minor_bckts = parsed_obj.ls_prefetch_slopes_intercepts_minor_bckts
                min_ls_simdstate = parsed_obj.ls_simdstate
                min_ls_pdis = parsed_obj.ls_pdis                
            else:
                if(parsed_obj.total_algorithm_time_in_ms < min_total_algorithm_time_in_ms):
                    min_total_algorithm_time_in_ms = parsed_obj.total_algorithm_time_in_ms
                    min_r_dataset = parsed_obj.datasets_r
                    min_s_dataset = parsed_obj.datasets_s
                    min_threads = parsed_obj.threads
                    min_use_avxsort_for_sorting_minor_bckts = parsed_obj.use_avxsort_for_sorting_minor_bckts 
                    min_ls_default_threshold = parsed_obj.ls_default_threshold
                    min_ls_default_arch = parsed_obj.ls_default_arch
                    min_ls_imv_avx = parsed_obj.ls_imv_avx
                    min_ls_prefetch_minor_bckt_sizes_off = parsed_obj.ls_prefetch_minor_bckt_sizes_off
                    min_ls_prefetch_slopes_intercepts_minor_bckts = parsed_obj.ls_prefetch_slopes_intercepts_minor_bckts
                    min_ls_simdstate = parsed_obj.ls_simdstate
                    min_ls_pdis = parsed_obj.ls_pdis

            i = i + 1

        print(f'min_total_algorithm_time_in_ms: {min_total_algorithm_time_in_ms}')
        print(f'min_r_dataset: {min_r_dataset}')
        print(f'min_s_dataset: {min_s_dataset}')
        print(f'min_threads: {min_threads}')
        print(f'min_use_avxsort_for_sorting_minor_bckts: {min_use_avxsort_for_sorting_minor_bckts}')
        print(f'min_ls_default_threshold: {min_ls_default_threshold}')
        print(f'min_ls_default_arch: {min_ls_default_arch}')
        print(f'min_ls_imv_avx: {min_ls_imv_avx}')
        print(f'min_ls_prefetch_minor_bckt_sizes_off: {min_ls_prefetch_minor_bckt_sizes_off}')
        print(f'min_ls_prefetch_slopes_intercepts_minor_bckts: {min_ls_prefetch_slopes_intercepts_minor_bckts}')
        print(f'min_ls_simdstate: {min_ls_simdstate}')
        print(f'min_ls_pdis: {min_ls_pdis}')


analyzer = LearnedSortMergeJoinAnalyzer()

#print(f'learned_sort_merge_join_tuning_unique/v0/')
#analyzer.analyze_datasets_get_min("/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_unique/v0/")
#print(f'learned_sort_merge_join_tuning_unique/v1/')
#analyzer.analyze_datasets_get_min("/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_unique/v1/")
#print(f'learned_sort_merge_join_tuning_unique/v2/')
#analyzer.analyze_datasets_get_min("/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_unique/v2/")
#print(f'learned_sort_merge_join_tuning_unique/v3/')
#analyzer.analyze_datasets_get_min("/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_unique/v3/")
#print(f'learned_sort_merge_join_tuning_unique/v4/')
#analyzer.analyze_datasets_get_min("/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_unique/v4/")
#print(f'learned_sort_merge_join_tuning_unique/v8/')
#analyzer.analyze_datasets_get_min("/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_unique/v8/")

#print(f'learned_sort_merge_join_tuning_lognormal/v0/')
#analyzer.analyze_datasets_get_min("/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_lognormal/v0/")
#print(f'learned_sort_merge_join_tuning_lognormal/v1/')
#analyzer.analyze_datasets_get_min("/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_lognormal/v1/")
#print(f'learned_sort_merge_join_tuning_lognormal/v2/')
#analyzer.analyze_datasets_get_min("/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_lognormal/v2/")
#print(f'learned_sort_merge_join_tuning_lognormal/v3/')
#analyzer.analyze_datasets_get_min("/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_lognormal/v3/")
#print(f'learned_sort_merge_join_tuning_lognormal/v4/')
#analyzer.analyze_datasets_get_min("/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_lognormal/v4/")
#print(f'learned_sort_merge_join_tuning_lognormal/v8/')
#analyzer.analyze_datasets_get_min("/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_lognormal/v8/")

print(f'sj_with_learned_uniform/32E6_128E6/non_imv')
analyzer.analyze_datasets_sort_results("/spinning/sabek/learned_join_results/sj_with_learned_uniform/32E6_128E6/non_imv/")
print(f'sj_with_learned_uniform/32E6_128E6/imv')
analyzer.analyze_datasets_sort_results("/spinning/sabek/learned_join_results/sj_with_learned_uniform/32E6_128E6/imv/")
print(f'sj_with_learned_uniform/32E6_640E6/non_imv')
analyzer.analyze_datasets_sort_results("/spinning/sabek/learned_join_results/sj_with_learned_uniform/32E6_640E6/non_imv/")
print(f'sj_with_learned_uniform/32E6_640E6/imv')
analyzer.analyze_datasets_sort_results("/spinning/sabek/learned_join_results/sj_with_learned_uniform/32E6_640E6/imv/")
print(f'sj_with_learned_uniform/128E6_640E6/non_imv')
analyzer.analyze_datasets_sort_results("/spinning/sabek/learned_join_results/sj_with_learned_uniform/128E6_640E6/non_imv/")
print(f'sj_with_learned_uniform/128E6_640E6/imv')
analyzer.analyze_datasets_sort_results("/spinning/sabek/learned_join_results/sj_with_learned_uniform/128E6_640E6/imv/")


#print(f'sj_with_learned_sosd_books_800M_uint64/')
#analyzer.analyze_datasets_sort_results("/spinning/sabek/learned_join_results/sj_with_learned_sosd_books_800M_uint64/")

#print(f'sj_with_learned_sosd_fb_200M_uint64/')
#analyzer.analyze_datasets_sort_results("/spinning/sabek/learned_join_results/sj_with_learned_sosd_fb_200M_uint64/")

#print(f'sj_with_learned_sosd_osm_cellids_800M_uint64/')
#analyzer.analyze_datasets_sort_results("/spinning/sabek/learned_join_results/sj_with_learned_sosd_osm_cellids_800M_uint64/")

#print(f'sj_with_learned_sosd_wiki_ts_200M_uint64/')
#analyzer.analyze_datasets_sort_results("/spinning/sabek/learned_join_results/sj_with_learned_sosd_wiki_ts_200M_uint64/")
