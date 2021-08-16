import os, sys

from parsed_benchmark_output import ParsedNonPartitionJoinBenchmarkOutput
from parsed_benchmark_output import ParsedRadixJoinBenchmarkOutput
from parsed_benchmark_output import ParsedETHSortMergeBenchmarkOutput
from parsed_benchmark_output import ParsedLearnedSortMergeBenchmarkOutput

class AlgosVsRuntimeFor1DatasetPreparer:
    
    def prepare(self, output_path, non_partition_join_1_path,
                      radix_join_1_path,
                      eth_sort_merge_join_1_path,
                      eth_sort_merge_with_learned_join_1_path,
                      learned_sort_merge_join_1_path):

        non_partition_join_1_parsed_obj = ParsedNonPartitionJoinBenchmarkOutput(non_partition_join_1_path)
        non_partition_join_1_parsed_obj.parse_benchmark_parameters()
        non_partition_join_1_parsed_obj.parse_benchmark_output()

        radix_join_1_parsed_obj = ParsedRadixJoinBenchmarkOutput(radix_join_1_path)
        radix_join_1_parsed_obj.parse_benchmark_parameters()
        radix_join_1_parsed_obj.parse_benchmark_output()

        eth_sort_merge_join_1_parsed_obj = ParsedETHSortMergeBenchmarkOutput(eth_sort_merge_join_1_path)
        eth_sort_merge_join_1_parsed_obj.parse_benchmark_parameters()
        eth_sort_merge_join_1_parsed_obj.parse_benchmark_output()

        eth_sort_merge_with_learned_join_1_parsed_obj = ParsedETHSortMergeBenchmarkOutput(eth_sort_merge_with_learned_join_1_path)
        eth_sort_merge_with_learned_join_1_parsed_obj.parse_benchmark_parameters()
        eth_sort_merge_with_learned_join_1_parsed_obj.parse_benchmark_output()

        learned_sort_merge_join_1_parsed_obj = ParsedLearnedSortMergeBenchmarkOutput(learned_sort_merge_join_1_path)
        learned_sort_merge_join_1_parsed_obj.parse_benchmark_parameters()
        learned_sort_merge_join_1_parsed_obj.parse_benchmark_output()

        lines = []
        lines.append("#Algos vs Runtime\n")

        lines.append(f'NPJ {non_partition_join_1_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'RJ {radix_join_1_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'SMJ {eth_sort_merge_join_1_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'LSMJ-1 {eth_sort_merge_with_learned_join_1_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'LSMJ-2 {learned_sort_merge_join_1_parsed_obj.total_algorithm_time_usec}\n')

        with open(output_path, "w") as output_file:
            output_file.writelines(lines); 


preparer = AlgosVsRuntimeFor1DatasetPreparer()

print(f'AlgosVsRuntimeFor1DatasetPreparer')
preparer.prepare("/spinning/sabek/learned_join_plots/algos_vs_runtime_for_1_dataset/algos_vs_runtime_for_1_dataset.dat",
                 "/spinning/sabek/learned_join_results/non_partition_join_tuning_unique/v8/non_partition_join_tuning_v8_th_64_bs_10_pd_32.csv",
                 "/spinning/sabek/learned_join_results/radix_join_tuning_unique/v8/radix_join_tuning_v8_th_64_rb_12_np_1.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v8/typical/sort_merge_join_tuning_v8_th_64_ls_0_ft_128_lsavx_0.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v8/learned/sort_merge_join_tuning_v8_th_64_ls_1_ft_128_lsavx_0.csv",
                 "/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_unique/v8/learned_sort_merge_join_tuning_v8_th_64_avxmb_1_lsdt_30000_lsda_200000_lsimv_0_lsps_0_lspsi_0_lsss_5_lspdis_320.csv")
