# A script for parsing profiled results by different join algorithms 

import os, sys

from parsed_benchmark_output import ParsedNonPartitionJoinBenchmarkOutput
from parsed_benchmark_output import ParsedRadixJoinBenchmarkOutput
from parsed_benchmark_output import ParsedETHSortMergeBenchmarkOutput
from parsed_benchmark_output import ParsedLearnedSortMergeBenchmarkOutput

class AlgosVsRuntimeFor5DatasetsPreparer:
    
    def prepare(self, output_path, non_partition_join_1_path, non_partition_join_2_path, non_partition_join_3_path, non_partition_join_4_path, non_partition_join_5_path, non_partition_join_6_path,
                      radix_join_1_path, radix_join_2_path, radix_join_3_path, radix_join_4_path, radix_join_5_path, radix_join_6_path,
                      eth_sort_merge_join_1_path, eth_sort_merge_join_2_path, eth_sort_merge_join_3_path, eth_sort_merge_join_4_path, eth_sort_merge_join_5_path, eth_sort_merge_join_6_path,
                      eth_sort_merge_with_learned_join_1_path, eth_sort_merge_with_learned_join_2_path, eth_sort_merge_with_learned_join_3_path, eth_sort_merge_with_learned_join_4_path, eth_sort_merge_with_learned_join_5_path, eth_sort_merge_with_learned_join_6_path,
                      learned_sort_merge_join_1_path, learned_sort_merge_join_2_path, learned_sort_merge_join_3_path, learned_sort_merge_join_4_path, learned_sort_merge_join_5_path, learned_sort_merge_join_6_path):

        non_partition_join_1_parsed_obj = ParsedNonPartitionJoinBenchmarkOutput(non_partition_join_1_path)
        non_partition_join_1_parsed_obj.parse_benchmark_parameters()
        non_partition_join_1_parsed_obj.parse_benchmark_output()
        non_partition_join_2_parsed_obj = ParsedNonPartitionJoinBenchmarkOutput(non_partition_join_2_path)
        non_partition_join_2_parsed_obj.parse_benchmark_parameters()
        non_partition_join_2_parsed_obj.parse_benchmark_output()
        non_partition_join_3_parsed_obj = ParsedNonPartitionJoinBenchmarkOutput(non_partition_join_3_path)
        non_partition_join_3_parsed_obj.parse_benchmark_parameters()
        non_partition_join_3_parsed_obj.parse_benchmark_output()
        non_partition_join_4_parsed_obj = ParsedNonPartitionJoinBenchmarkOutput(non_partition_join_4_path)
        non_partition_join_4_parsed_obj.parse_benchmark_parameters()
        non_partition_join_4_parsed_obj.parse_benchmark_output()
        non_partition_join_5_parsed_obj = ParsedNonPartitionJoinBenchmarkOutput(non_partition_join_5_path)
        non_partition_join_5_parsed_obj.parse_benchmark_parameters()
        non_partition_join_5_parsed_obj.parse_benchmark_output()
        non_partition_join_6_parsed_obj = ParsedNonPartitionJoinBenchmarkOutput(non_partition_join_6_path)
        non_partition_join_6_parsed_obj.parse_benchmark_parameters()
        non_partition_join_6_parsed_obj.parse_benchmark_output()

        radix_join_1_parsed_obj = ParsedRadixJoinBenchmarkOutput(radix_join_1_path)
        radix_join_1_parsed_obj.parse_benchmark_parameters()
        radix_join_1_parsed_obj.parse_benchmark_output()
        radix_join_2_parsed_obj = ParsedRadixJoinBenchmarkOutput(radix_join_2_path)
        radix_join_2_parsed_obj.parse_benchmark_parameters()
        radix_join_2_parsed_obj.parse_benchmark_output()
        radix_join_3_parsed_obj = ParsedRadixJoinBenchmarkOutput(radix_join_3_path)
        radix_join_3_parsed_obj.parse_benchmark_parameters()
        radix_join_3_parsed_obj.parse_benchmark_output()
        radix_join_4_parsed_obj = ParsedRadixJoinBenchmarkOutput(radix_join_4_path)
        radix_join_4_parsed_obj.parse_benchmark_parameters()
        radix_join_4_parsed_obj.parse_benchmark_output()
        radix_join_5_parsed_obj = ParsedRadixJoinBenchmarkOutput(radix_join_5_path)
        radix_join_5_parsed_obj.parse_benchmark_parameters()
        radix_join_5_parsed_obj.parse_benchmark_output()
        radix_join_6_parsed_obj = ParsedRadixJoinBenchmarkOutput(radix_join_6_path)
        radix_join_6_parsed_obj.parse_benchmark_parameters()
        radix_join_6_parsed_obj.parse_benchmark_output()

        eth_sort_merge_join_1_parsed_obj = ParsedETHSortMergeBenchmarkOutput(eth_sort_merge_join_1_path)
        eth_sort_merge_join_1_parsed_obj.parse_benchmark_parameters()
        eth_sort_merge_join_1_parsed_obj.parse_benchmark_output()
        eth_sort_merge_join_2_parsed_obj = ParsedETHSortMergeBenchmarkOutput(eth_sort_merge_join_2_path)
        eth_sort_merge_join_2_parsed_obj.parse_benchmark_parameters()
        eth_sort_merge_join_2_parsed_obj.parse_benchmark_output()
        eth_sort_merge_join_3_parsed_obj = ParsedETHSortMergeBenchmarkOutput(eth_sort_merge_join_3_path)
        eth_sort_merge_join_3_parsed_obj.parse_benchmark_parameters()
        eth_sort_merge_join_3_parsed_obj.parse_benchmark_output()
        eth_sort_merge_join_4_parsed_obj = ParsedETHSortMergeBenchmarkOutput(eth_sort_merge_join_4_path)
        eth_sort_merge_join_4_parsed_obj.parse_benchmark_parameters()
        eth_sort_merge_join_4_parsed_obj.parse_benchmark_output()
        eth_sort_merge_join_5_parsed_obj = ParsedETHSortMergeBenchmarkOutput(eth_sort_merge_join_5_path)
        eth_sort_merge_join_5_parsed_obj.parse_benchmark_parameters()
        eth_sort_merge_join_5_parsed_obj.parse_benchmark_output()
        eth_sort_merge_join_6_parsed_obj = ParsedETHSortMergeBenchmarkOutput(eth_sort_merge_join_6_path)
        eth_sort_merge_join_6_parsed_obj.parse_benchmark_parameters()
        eth_sort_merge_join_6_parsed_obj.parse_benchmark_output()

        eth_sort_merge_with_learned_join_1_parsed_obj = ParsedETHSortMergeBenchmarkOutput(eth_sort_merge_with_learned_join_1_path)
        eth_sort_merge_with_learned_join_1_parsed_obj.parse_benchmark_parameters()
        eth_sort_merge_with_learned_join_1_parsed_obj.parse_benchmark_output()
        eth_sort_merge_with_learned_join_2_parsed_obj = ParsedETHSortMergeBenchmarkOutput(eth_sort_merge_with_learned_join_2_path)
        eth_sort_merge_with_learned_join_2_parsed_obj.parse_benchmark_parameters()
        eth_sort_merge_with_learned_join_2_parsed_obj.parse_benchmark_output()
        eth_sort_merge_with_learned_join_3_parsed_obj = ParsedETHSortMergeBenchmarkOutput(eth_sort_merge_with_learned_join_3_path)
        eth_sort_merge_with_learned_join_3_parsed_obj.parse_benchmark_parameters()
        eth_sort_merge_with_learned_join_3_parsed_obj.parse_benchmark_output()
        eth_sort_merge_with_learned_join_4_parsed_obj = ParsedETHSortMergeBenchmarkOutput(eth_sort_merge_with_learned_join_4_path)
        eth_sort_merge_with_learned_join_4_parsed_obj.parse_benchmark_parameters()
        eth_sort_merge_with_learned_join_4_parsed_obj.parse_benchmark_output()
        eth_sort_merge_with_learned_join_5_parsed_obj = ParsedETHSortMergeBenchmarkOutput(eth_sort_merge_with_learned_join_5_path)
        eth_sort_merge_with_learned_join_5_parsed_obj.parse_benchmark_parameters()
        eth_sort_merge_with_learned_join_5_parsed_obj.parse_benchmark_output()
        eth_sort_merge_with_learned_join_6_parsed_obj = ParsedETHSortMergeBenchmarkOutput(eth_sort_merge_with_learned_join_6_path)
        eth_sort_merge_with_learned_join_6_parsed_obj.parse_benchmark_parameters()
        eth_sort_merge_with_learned_join_6_parsed_obj.parse_benchmark_output()

        learned_sort_merge_join_1_parsed_obj = ParsedLearnedSortMergeBenchmarkOutput(learned_sort_merge_join_1_path)
        learned_sort_merge_join_1_parsed_obj.parse_benchmark_parameters()
        learned_sort_merge_join_1_parsed_obj.parse_benchmark_output()
        learned_sort_merge_join_2_parsed_obj = ParsedLearnedSortMergeBenchmarkOutput(learned_sort_merge_join_2_path)
        learned_sort_merge_join_2_parsed_obj.parse_benchmark_parameters()
        learned_sort_merge_join_2_parsed_obj.parse_benchmark_output()
        learned_sort_merge_join_3_parsed_obj = ParsedLearnedSortMergeBenchmarkOutput(learned_sort_merge_join_3_path)
        learned_sort_merge_join_3_parsed_obj.parse_benchmark_parameters()
        learned_sort_merge_join_3_parsed_obj.parse_benchmark_output()
        learned_sort_merge_join_4_parsed_obj = ParsedLearnedSortMergeBenchmarkOutput(learned_sort_merge_join_4_path)
        learned_sort_merge_join_4_parsed_obj.parse_benchmark_parameters()
        learned_sort_merge_join_4_parsed_obj.parse_benchmark_output()
        learned_sort_merge_join_5_parsed_obj = ParsedLearnedSortMergeBenchmarkOutput(learned_sort_merge_join_5_path)
        learned_sort_merge_join_5_parsed_obj.parse_benchmark_parameters()
        learned_sort_merge_join_5_parsed_obj.parse_benchmark_output()
        learned_sort_merge_join_6_parsed_obj = ParsedLearnedSortMergeBenchmarkOutput(learned_sort_merge_join_6_path)
        learned_sort_merge_join_6_parsed_obj.parse_benchmark_parameters()
        learned_sort_merge_join_6_parsed_obj.parse_benchmark_output()

        lines = []
        lines.append("#Algos vs Runtime\n")

        lines.append("\n#NPJ\n")
        lines.append(f'16M {non_partition_join_1_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'32M {non_partition_join_2_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'128M {non_partition_join_3_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'640M {non_partition_join_4_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'1664M {non_partition_join_5_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'1920M {non_partition_join_6_parsed_obj.total_algorithm_time_usec}\n')

        lines.append("\n#RJ\n")       
        lines.append(f'16M {radix_join_1_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'32M {radix_join_2_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'128M {radix_join_3_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'640M {radix_join_4_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'1664M {radix_join_5_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'1920M {radix_join_6_parsed_obj.total_algorithm_time_usec}\n')

        lines.append("\n#SMJ\n")       
        lines.append(f'16M {eth_sort_merge_join_1_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'32M {eth_sort_merge_join_2_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'128M {eth_sort_merge_join_3_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'640M {eth_sort_merge_join_4_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'1664M {eth_sort_merge_join_5_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'1920M {eth_sort_merge_join_6_parsed_obj.total_algorithm_time_usec}\n')

        lines.append("\n#LSMJ-1\n")       
        lines.append(f'16M {eth_sort_merge_with_learned_join_1_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'32M {eth_sort_merge_with_learned_join_2_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'128M {eth_sort_merge_with_learned_join_3_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'640M {eth_sort_merge_with_learned_join_4_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'1664M {eth_sort_merge_with_learned_join_5_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'1920M {eth_sort_merge_with_learned_join_6_parsed_obj.total_algorithm_time_usec}\n')

        lines.append("\n#LSMJ-2\n")
        lines.append(f'16M {learned_sort_merge_join_1_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'32M {learned_sort_merge_join_2_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'128M {learned_sort_merge_join_3_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'640M {learned_sort_merge_join_4_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'1664M {learned_sort_merge_join_5_parsed_obj.total_algorithm_time_usec}\n')
        lines.append(f'1920M {learned_sort_merge_join_6_parsed_obj.total_algorithm_time_usec}\n')

        with open(output_path, "w") as output_file:
            output_file.writelines(lines); 


preparer = AlgosVsRuntimeFor5DatasetsPreparer()

print(f'AlgosVsRuntimeFor5DatasetsPreparer')
preparer.prepare("/spinning/sabek/learned_join_plots/algos_vs_runtime_for_5_datasets_unique/algos_vs_runtime_for_5_datasets_unique.dat",
                 "/spinning/sabek/learned_join_results/non_partition_join_tuning_unique/v0/non_partition_join_tuning_v0_th_64_bs_2_pd_128.csv",
                 "/spinning/sabek/learned_join_results/non_partition_join_tuning_unique/v1/non_partition_join_tuning_v1_th_64_bs_2_pd_128.csv",
                 "/spinning/sabek/learned_join_results/non_partition_join_tuning_unique/v2/non_partition_join_tuning_v2_th_64_bs_2_pd_128.csv",
                 "/spinning/sabek/learned_join_results/non_partition_join_tuning_unique/v4/non_partition_join_tuning_v4_th_64_bs_2_pd_128.csv",
                 "/spinning/sabek/learned_join_results/non_partition_join_tuning_unique/v7/non_partition_join_tuning_v7_th_64_bs_2_pd_128.csv",
                 "/spinning/sabek/learned_join_results/non_partition_join_tuning_unique/v8/non_partition_join_tuning_v8_th_64_bs_10_pd_32.csv",
                 "/spinning/sabek/learned_join_results/radix_join_tuning_unique/v0/radix_join_tuning_v0_th_64_rb_12_np_1.csv",
                 "/spinning/sabek/learned_join_results/radix_join_tuning_unique/v1/radix_join_tuning_v1_th_64_rb_12_np_1.csv",
                 "/spinning/sabek/learned_join_results/radix_join_tuning_unique/v2/radix_join_tuning_v2_th_64_rb_12_np_1.csv",
                 "/spinning/sabek/learned_join_results/radix_join_tuning_unique/v4/radix_join_tuning_v4_th_64_rb_12_np_1.csv",
                 "/spinning/sabek/learned_join_results/radix_join_tuning_unique/v7/radix_join_tuning_v7_th_64_rb_12_np_1.csv",
                 "/spinning/sabek/learned_join_results/radix_join_tuning_unique/v8/radix_join_tuning_v8_th_64_rb_12_np_1.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v0/typical/sort_merge_join_tuning_v0_th_64_ls_0_ft_128_lsavx_0.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v1/typical/sort_merge_join_tuning_v1_th_64_ls_0_ft_128_lsavx_0.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v2/typical/sort_merge_join_tuning_v2_th_64_ls_0_ft_128_lsavx_0.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v4/typical/sort_merge_join_tuning_v4_th_64_ls_0_ft_128_lsavx_0.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v7/typical/sort_merge_join_tuning_v7_th_64_ls_0_ft_1024_lsavx_0.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v8/typical/sort_merge_join_tuning_v8_th_64_ls_0_ft_1024_lsavx_0.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v0/learned/sort_merge_join_tuning_v0_th_64_ls_1_ft_128_lsavx_0.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v1/learned/sort_merge_join_tuning_v1_th_64_ls_1_ft_128_lsavx_0.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v2/learned/sort_merge_join_tuning_v2_th_64_ls_1_ft_128_lsavx_0.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v4/learned/sort_merge_join_tuning_v4_th_64_ls_1_ft_128_lsavx_0.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v7/learned/sort_merge_join_tuning_v7_th_64_ls_1_ft_1024_lsavx_0.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v8/learned/sort_merge_join_tuning_v8_th_64_ls_1_ft_1024_lsavx_0.csv", 
                 "/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_unique/v0/learned_sort_merge_join_tuning_v0_th_64_avxmb_1_lsdt_30000_lsda_100000_lsimv_0_lsps_0_lspsi_0_lsss_5_lspdis_320.csv", 
                 "/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_unique/v1/learned_sort_merge_join_tuning_v1_th_64_avxmb_1_lsdt_30000_lsda_100000_lsimv_0_lsps_0_lspsi_0_lsss_5_lspdis_320.csv",
                 "/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_unique/v2/learned_sort_merge_join_tuning_v2_th_64_avxmb_1_lsdt_30000_lsda_100000_lsimv_0_lsps_0_lspsi_0_lsss_5_lspdis_320.csv",
                 "/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_unique/v4/learned_sort_merge_join_tuning_v4_th_64_avxmb_1_lsdt_30000_lsda_200000_lsimv_0_lsps_0_lspsi_0_lsss_5_lspdis_320.csv",
                 "/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_unique/v7/learned_sort_merge_join_tuning_v7_th_64_avxmb_1_lsdt_30000_lsda_200000_lsimv_0_lsps_0_lspsi_0_lsss_5_lspdis_320.csv",
                 "/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_unique/v8/learned_sort_merge_join_tuning_v8_th_64_avxmb_1_lsdt_30000_lsda_200000_lsimv_0_lsps_0_lspsi_0_lsss_5_lspdis_320.csv")

preparer.prepare("/spinning/sabek/learned_join_plots/algos_vs_runtime_for_5_datasets_lognormal/algos_vs_runtime_for_5_datasets_lognormal.dat",
                 "/spinning/sabek/learned_join_results/non_partition_join_tuning_lognormal/v0/non_partition_join_tuning_v0_th_64_bs_10_pd_128.csv",
                 "/spinning/sabek/learned_join_results/non_partition_join_tuning_lognormal/v1/non_partition_join_tuning_v1_th_64_bs_10_pd_128.csv",
                 "/spinning/sabek/learned_join_results/non_partition_join_tuning_lognormal/v2/non_partition_join_tuning_v2_th_64_bs_10_pd_128.csv",
                 "/spinning/sabek/learned_join_results/non_partition_join_tuning_lognormal/v4/non_partition_join_tuning_v4_th_64_bs_10_pd_128_tmp.csv",
                 "/spinning/sabek/learned_join_results/non_partition_join_tuning_lognormal/v7/non_partition_join_tuning_v7_th_64_bs_10_pd_128_tmp.csv",
                 "/spinning/sabek/learned_join_results/non_partition_join_tuning_lognormal/v8/non_partition_join_tuning_v8_th_64_bs_10_pd_128_tmp.csv",
                 "/spinning/sabek/learned_join_results/radix_join_tuning_lognormal/v0/radix_join_tuning_v0_th_64_rb_10_np_1.csv",
                 "/spinning/sabek/learned_join_results/radix_join_tuning_lognormal/v1/radix_join_tuning_v1_th_64_rb_12_np_1.csv",
                 "/spinning/sabek/learned_join_results/radix_join_tuning_lognormal/v2/radix_join_tuning_v2_th_64_rb_12_np_1.csv",
                 "/spinning/sabek/learned_join_results/radix_join_tuning_lognormal/v4/radix_join_tuning_v4_th_64_rb_12_np_1.csv",
                 "/spinning/sabek/learned_join_results/radix_join_tuning_lognormal/v7/radix_join_tuning_v7_th_64_rb_12_np_1.csv",
                 "/spinning/sabek/learned_join_results/radix_join_tuning_lognormal/v8/radix_join_tuning_v8_th_64_rb_12_np_1.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_lognormal/v0/typical/sort_merge_join_tuning_v0_th_64_ls_0_ft_128_lsavx_0.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_lognormal/v1/typical/sort_merge_join_tuning_v1_th_64_ls_0_ft_128_lsavx_0.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_lognormal/v2/typical/sort_merge_join_tuning_v2_th_64_ls_0_ft_128_lsavx_0.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_lognormal/v4/typical/sort_merge_join_tuning_v4_th_64_ls_0_ft_128_lsavx_0_tmp.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_lognormal/v7/typical/sort_merge_join_tuning_v7_th_64_ls_0_ft_1024_lsavx_0_tmp.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_lognormal/v8/typical/sort_merge_join_tuning_v8_th_64_ls_0_ft_1024_lsavx_0_tmp.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_lognormal/v0/learned/sort_merge_join_tuning_v0_th_64_ls_1_ft_128_lsavx_0.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_lognormal/v1/learned/sort_merge_join_tuning_v1_th_64_ls_1_ft_128_lsavx_0.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_lognormal/v2/learned/sort_merge_join_tuning_v2_th_64_ls_1_ft_128_lsavx_0.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_lognormal/v4/learned/sort_merge_join_tuning_v4_th_64_ls_1_ft_128_lsavx_0_tmp.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_lognormal/v7/learned/sort_merge_join_tuning_v7_th_64_ls_1_ft_1024_lsavx_0_tmp.csv",
                 "/spinning/sabek/learned_join_results/sort_merge_join_tuning_lognormal/v8/learned/sort_merge_join_tuning_v8_th_64_ls_1_ft_1024_lsavx_0_tmp.csv", 
                 "/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_lognormal/v0/learned_sort_merge_join_tuning_v0_th_64_avxmb_1_lsdt_60000_lsda_50000_lsimv_0_lsps_0_lspsi_0_lsss_5_lspdis_320.csv", 
                 "/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_lognormal/v1/learned_sort_merge_join_tuning_v1_th_64_avxmb_1_lsdt_60000_lsda_50000_lsimv_0_lsps_0_lspsi_0_lsss_5_lspdis_320.csv",
                 "/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_lognormal/v2/learned_sort_merge_join_tuning_v2_th_64_avxmb_1_lsdt_60000_lsda_50000_lsimv_0_lsps_0_lspsi_0_lsss_5_lspdis_320.csv",
                 "/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_lognormal/v4/learned_sort_merge_join_tuning_v4_th_64_avxmb_1_lsdt_60000_lsda_50000_lsimv_0_lsps_0_lspsi_0_lsss_5_lspdis_320.csv",
                 "/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_lognormal/v7/learned_sort_merge_join_tuning_v7_th_64_avxmb_1_lsdt_60000_lsda_50000_lsimv_0_lsps_0_lspsi_0_lsss_5_lspdis_320.csv",
                 "/spinning/sabek/learned_join_results/learned_sort_merge_join_tuning_lognormal/v8/learned_sort_merge_join_tuning_v8_th_64_avxmb_1_lsdt_60000_lsda_50000_lsimv_0_lsps_0_lspsi_0_lsss_5_lspdis_320.csv")                 