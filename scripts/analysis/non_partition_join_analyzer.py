# A script for parsing parameters used in experiments and profiled results by non-partition hash join 

import os, sys
from parsed_benchmark_output import ParsedNonPartitionJoinBenchmarkOutput

class NonPartitionJoinAnalyzer:
    
    def analyze_datasets(self, datasets_path):
        files = os.listdir(datasets_path)

        i = 0
        min_total_algorithm_time_usec = 0
        min_threads = 0
        min_bucket_size = 0 
        min_prefetch_distance = 0

        for file in files:
            #print(datasets_path+file)
            parsed_obj = ParsedNonPartitionJoinBenchmarkOutput(datasets_path+file)
            parsed_obj.parse_benchmark_parameters()
            parsed_obj.parse_benchmark_output()

            if i == 0:
                min_total_algorithm_time_usec = parsed_obj.total_algorithm_time_usec
                min_threads = parsed_obj.threads
                min_bucket_size = parsed_obj.bucket_size
                min_prefetch_distance = parsed_obj.prefetch_distance              
            else:
                if(parsed_obj.total_algorithm_time_usec < min_total_algorithm_time_usec):
                    min_total_algorithm_time_usec = parsed_obj.total_algorithm_time_usec
                    min_threads = parsed_obj.threads
                    min_bucket_size = parsed_obj.bucket_size
                    min_prefetch_distance = parsed_obj.prefetch_distance
            i = i + 1

        print(f'min_total_algorithm_time_usec: {min_total_algorithm_time_usec}')
        print(f'min_threads: {min_threads}')
        print(f'min_bucket_size: {min_bucket_size}')
        print(f'min_prefetch_distance: {min_prefetch_distance}')

analyzer = NonPartitionJoinAnalyzer()

print(f'non_partition_join_tuning_unique/v0/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/non_partition_join_tuning_unique/v0/")
print(f'non_partition_join_tuning_unique/v1/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/non_partition_join_tuning_unique/v1/")
print(f'non_partition_join_tuning_unique/v2/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/non_partition_join_tuning_unique/v2/")
print(f'non_partition_join_tuning_unique/v3/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/non_partition_join_tuning_unique/v3/")
print(f'non_partition_join_tuning_unique/v4/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/non_partition_join_tuning_unique/v4/")
print(f'non_partition_join_tuning_unique/v8/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/non_partition_join_tuning_unique/v8/")

print(f'non_partition_join_tuning_lognormal/v0/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/non_partition_join_tuning_lognormal/v0/")
print(f'non_partition_join_tuning_lognormal/v1/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/non_partition_join_tuning_lognormal/v1/")
print(f'non_partition_join_tuning_lognormal/v2/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/non_partition_join_tuning_lognormal/v2/")
