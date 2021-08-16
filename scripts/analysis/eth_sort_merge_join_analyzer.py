import os, sys
from parsed_benchmark_output import ParsedETHSortMergeBenchmarkOutput

class ETHSortMergeJoinAnalyzer:
    
    def analyze_datasets(self, datasets_path):
        files = os.listdir(datasets_path)

        i = 0
        min_total_algorithm_time_usec = 0
        min_threads = 0
        min_use_learned_sort = 0 
        min_fanout_per_thread = 0
        min_use_avxsort_as_std_sort = 0        

        for file in files:
            #print(datasets_path+file)
            parsed_obj = ParsedETHSortMergeBenchmarkOutput(datasets_path+file)
            parsed_obj.parse_benchmark_parameters()
            parsed_obj.parse_benchmark_output()

            if i == 0:
                min_total_algorithm_time_usec = parsed_obj.total_algorithm_time_usec
                min_threads = parsed_obj.threads
                min_use_learned_sort = parsed_obj.use_learned_sort
                min_fanout_per_thread = parsed_obj.fanout_per_thread
                min_use_avxsort_as_std_sort = parsed_obj.use_avxsort_as_std_sort
              
            else:
                if(parsed_obj.total_algorithm_time_usec < min_total_algorithm_time_usec):
                    min_total_algorithm_time_usec = parsed_obj.total_algorithm_time_usec
                    min_threads = parsed_obj.threads
                    min_use_learned_sort = parsed_obj.use_learned_sort
                    min_fanout_per_thread = parsed_obj.fanout_per_thread
                    min_use_avxsort_as_std_sort = parsed_obj.use_avxsort_as_std_sort

            i = i + 1

        print(f'min_total_algorithm_time_usec: {min_total_algorithm_time_usec}')
        print(f'min_threads: {min_threads}')
        print(f'min_use_learned_sort: {min_use_learned_sort}')
        print(f'min_fanout_per_thread: {min_fanout_per_thread}')
        print(f'min_use_avxsort_as_std_sort: {min_use_avxsort_as_std_sort}')

analyzer = ETHSortMergeJoinAnalyzer()

print(f'sort_merge_join_tuning_unique/v0/typical/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v0/typical/")
print(f'sort_merge_join_tuning_unique/v1/typical/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v1/typical/")
print(f'sort_merge_join_tuning_unique/v2/typical/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v2/typical/")
print(f'sort_merge_join_tuning_unique/v3/typical/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v3/typical/")
print(f'sort_merge_join_tuning_unique/v4/typical/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v4/typical/")
print(f'sort_merge_join_tuning_unique/v8/typical/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v8/typical/")

print(f'sort_merge_join_tuning_lognormal/v0/typical/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/sort_merge_join_tuning_lognormal/v0/typical/")
print(f'sort_merge_join_tuning_lognormal/v1/typical/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/sort_merge_join_tuning_lognormal/v1/typical/")
print(f'sort_merge_join_tuning_lognormal/v2/typical/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/sort_merge_join_tuning_lognormal/v2/typical/")


print(f'sort_merge_join_tuning_unique/v0/learned/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v0/learned/")
print(f'sort_merge_join_tuning_unique/v1/learned/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v1/learned/")
print(f'sort_merge_join_tuning_unique/v2/learned/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v2/learned/")
print(f'sort_merge_join_tuning_unique/v3/learned/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v3/learned/")
print(f'sort_merge_join_tuning_unique/v4/learned/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v4/learned/")
print(f'sort_merge_join_tuning_unique/v8/learned/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/sort_merge_join_tuning_unique/v8/learned/")


print(f'sort_merge_join_tuning_lognormal/v0/learned/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/sort_merge_join_tuning_lognormal/v0/learned/")
print(f'sort_merge_join_tuning_lognormal/v1/learned/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/sort_merge_join_tuning_lognormal/v1/learned/")
print(f'sort_merge_join_tuning_lognormal/v2/learned/')
analyzer.analyze_datasets("/spinning/sabek/learned_join_results/sort_merge_join_tuning_lognormal/v2/learned/")
