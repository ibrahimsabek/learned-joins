# A script for parsing parameters used in experiments and profiled results by partitioned hash join 

import os, sys
from parsed_benchmark_output import ParsedRadixJoinBenchmarkOutput

class RadixJoinAnalyzer:
    
    def analyze_datasets(self, datasets_path):
        files = os.listdir(datasets_path)

        i = 0
        min_total_algorithm_time_usec = 0
        min_threads = 0
        min_num_radix_bits = 0 
        min_num_passes = 0


        for file in files:
            #print(datasets_path+file)
            parsed_obj = ParsedRadixJoinBenchmarkOutput(datasets_path+file)
            parsed_obj.parse_benchmark_parameters()
            parsed_obj.parse_benchmark_output()

            if i == 0:
                min_total_algorithm_time_usec = parsed_obj.total_algorithm_time_usec
                min_threads = parsed_obj.threads
                min_num_radix_bits = parsed_obj.num_radix_bits
                min_num_passes = parsed_obj.num_passes              
            else:
                if(parsed_obj.total_algorithm_time_usec < min_total_algorithm_time_usec):
                    min_total_algorithm_time_usec = parsed_obj.total_algorithm_time_usec
                    min_threads = parsed_obj.threads
                    min_num_radix_bits = parsed_obj.num_radix_bits
                    min_num_passes = parsed_obj.num_passes

            i = i + 1

        print(f'min_total_algorithm_time_usec: {min_total_algorithm_time_usec}')
        print(f'min_threads: {min_threads}')
        print(f'min_num_radix_bits: {min_num_radix_bits}')
        print(f'min_num_passes: {min_num_passes}')

analyzer = RadixJoinAnalyzer()

print(f'radix_join_tuning_unique/v0/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/radix_join_tuning_unique/v0/")
print(f'radix_join_tuning_unique/v1/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/radix_join_tuning_unique/v1/")
print(f'radix_join_tuning_unique/v2/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/radix_join_tuning_unique/v2/")
print(f'radix_join_tuning_unique/v3/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/radix_join_tuning_unique/v3/")
print(f'radix_join_tuning_unique/v4/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/radix_join_tuning_unique/v4/")
print(f'radix_join_tuning_unique/v8/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/radix_join_tuning_unique/v8/")

print(f'radix_join_tuning_lognormal/v0/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/radix_join_tuning_lognormal/v0/")
print(f'radix_join_tuning_lognormal/v1/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/radix_join_tuning_lognormal/v1/")
print(f'radix_join_tuning_lognormal/v2/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/radix_join_tuning_lognormal/v2/")
print(f'radix_join_tuning_lognormal/v3/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/radix_join_tuning_lognormal/v3/")
print(f'radix_join_tuning_lognormal/v4/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/radix_join_tuning_lognormal/v4/")
print(f'radix_join_tuning_lognormal/v8/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/radix_join_tuning_lognormal/v8/")
