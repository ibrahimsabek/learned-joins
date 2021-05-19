import csv
import os, sys

class ParsedLearnedSortMergeBenchmarkOutput:
    def __init__(self, benchmark_output_path_in, num_parsed_items_in=1):
        self.benchmark_output_path = benchmark_output_path_in
        self.num_parsed_items = num_parsed_items_in

    def parse_benchmark_parameters(self):
        parts = os.path.split(self.benchmark_output_path)
        file_name = parts[1]
        configs = file_name[24:-4]
        configs_arr = configs.split("_")
        self.datasets_r = configs_arr[0]
        self.datasets_s = configs_arr[1]
        self.threads = configs_arr[3]
        self.use_avxsort_for_sorting_minor_bckts = configs_arr[5] 
        self.ls_default_threshold = configs_arr[7]
        self.ls_default_arch = configs_arr[9]
        self.ls_imv_avx = configs_arr[11]
        self.ls_prefetch_minor_bckt_sizes_off = configs_arr[13]
        self.ls_prefetch_slopes_intercepts_minor_bckts = configs_arr[15]
        self.ls_simdstate = configs_arr[17]
        self.ls_pdis = configs_arr[19]
        print(self.datasets_r, self.datasets_s,  self.ls_simdstate, self.ls_pdis)

    def parse_benchmark_output(self):
        keys = []
        values = [] 
        with open(self.benchmark_output_path, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            i = 0
            for row in reader:
                if i == 0:
                    keys.extend(row)
                else:
                    values.append(row)
                i = i + 1

        if len(values) > 1 :
            printf('Not supported multiple items for now')
        else:
            for j in range(len(values[0])):
                if keys[j] == 'learned_model_in_ms':
                    self.learned_model_in_ms = values[0][j]          
                elif keys[j] == 'partition_in_ms':
                    self.partition_in_ms = values[0][j]
                elif keys[j] == 'sorting_in_ms':
                    self.sorting_in_ms = values[0][j]
                elif keys[j] == 'join_in_ms':
                    self.join_in_ms = values[0][j]
                elif keys[j] == 'learned_model_Throughput_in_mtuples_per_sec':
                    self.learned_model_Throughput_in_mtuples_per_sec = values[0][j]
                elif keys[j] == 'learned_model_Cycles':
                    self.learned_model_Cycles = values[0][j]
                elif keys[j] == 'learned_model_LLC_misses':
                    self.learned_model_LLC_misses = values[0][j]
                elif keys[j] == 'learned_model_L1_misses':
                    self.learned_model_L1_misses = values[0][j]
                elif keys[j] == 'learned_model_Instructions':
                    self.learned_model_Instructions = values[0][j]
                elif keys[j] == 'learned_model_Branch_misses':
                    self.learned_model_Branch_misses = values[0][j]
                elif keys[j] == 'learned_model_Task_clock':
                    self.learned_model_Task_clock = float(values[0][j])
                elif keys[j] == 'partition_Throughput_in_mtuples_per_sec':
                    self.partition_Throughput_in_mtuples_per_sec = float(values[0][j])
                elif keys[j] == 'partition_Cycles':
                    self.partition_Cycles = float(values[0][j])
                elif keys[j] == 'partition_LLC_misses':
                    self.partition_LLC_misses = float(values[0][j])
                elif keys[j] == 'partition_L1_misses':
                    self.partition_L1_misses = float(values[0][j])
                elif keys[j] == 'partition_Instructions':
                    self.partition_Instructions = float(values[0][j])
                elif keys[j] == 'partition_Branch_misses':
                    self.partition_Branch_misses = float(values[0][j])                    
                elif keys[j] == 'partition_Task_clock':
                    self.partition_Task_clock = float(values[0][j])
                elif keys[j] == 'sorting_Throughput_in_mtuples_per_sec':
                    self.sorting_Throughput_in_mtuples_per_sec = float(values[0][j])
                elif keys[j] == 'sorting_Cycles':
                    self.sorting_Cycles = float(values[0][j])                    
                elif keys[j] == 'sorting_LLC_misses':
                    self.sorting_LLC_misses = float(values[0][j])                
                elif keys[j] == 'sorting_L1_misses':
                    self.sorting_L1_misses = float(values[0][j])                
                elif keys[j] == 'sorting_Instructions':
                    self.sorting_Instructions = float(values[0][j])                
                elif keys[j] == 'sorting_Branch_misses':
                    self.sorting_Branch_misses = float(values[0][j])                
                elif keys[j] == 'sorting_Task_clock':
                    self.sorting_Task_clock = float(values[0][j])                
                elif keys[j] == 'join_Throughput_in_mtuples_per_sec':
                    self.join_Throughput_in_mtuples_per_sec = float(values[0][j])                
                elif keys[j] == 'join_Cycles':
                    self.join_Cycles = float(values[0][j])                
                elif keys[j] == 'join_LLC_misses':
                    self.join_LLC_misses = float(values[0][j])       
                elif keys[j] == 'join_L1_misses':
                    self.join_L1_misses = float(values[0][j])                
                elif keys[j] == 'join_Instructions':
                    self.join_Instructions = float(values[0][j])                
                elif keys[j] == 'join_Branch_misses':
                    self.join_Branch_misses = float(values[0][j])                
                elif keys[j] == 'join_Task_clock':
                    self.join_Task_clock = float(values[0][j]) 

class ParsedETHSortMergeBenchmarkOutput:
    def __init__(self, benchmark_output_path_in, num_parsed_items_in=1):
        self.benchmark_output_path = benchmark_output_path_in
        self.num_parsed_items = num_parsed_items_in

    def parse_benchmark_parameters(self):
        parts = os.path.split(self.benchmark_output_path)
        file_name = parts[1]
        configs = file_name[23:-4]
        configs_arr = configs.split("_")
        self.datasets = configs_arr[0]
        self.threads = configs_arr[2]
        self.use_learned_sort = configs_arr[4] 
        self.fanout_per_thread = configs_arr[6]
        self.use_avxsort_as_std_sort = configs_arr[8]


    def parse_benchmark_output(self):
        keys = []
        values = [] 
        with open(self.benchmark_output_path, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            i = 0
            for row in reader:
                if i == 0:
                    keys.extend(row)
                else:
                    values.append(row)
                i = i + 1

        if len(values) > 1 :
            printf('Not supported multiple items for now')
        else:
            for j in range(len(values[0])):
                if keys[j] == 'name':
                    self.name = values[0][j]          
                elif keys[j] == 'iterations':
                    self.iterations = values[0][j]
                elif keys[j] == 'real_time':
                    self.real_time = values[0][j]
                elif keys[j] == 'cpu_time':
                    self.cpu_time = values[0][j]
                elif keys[j] == 'time_unit':
                    self.time_unit = values[0][j]
                elif keys[j] == 'bytes_per_second':
                    self.bytes_per_second = values[0][j]
                elif keys[j] == 'items_per_second':
                    self.items_per_second = values[0][j]
                elif keys[j] == 'label':
                    self.label = values[0][j]
                elif keys[j] == 'error_occurred':
                    self.error_occurred = values[0][j]
                elif keys[j] == 'error_message':
                    self.error_message = values[0][j]
                elif keys[j] == 'num_join_results':
                    self.num_join_results = float(values[0][j])
                elif keys[j] == 'total_algorithm_time_usec':
                    self.total_algorithm_time_usec = float(values[0][j])
                elif keys[j] == 'total_joining_branch_misses':
                    self.total_joining_branch_misses = float(values[0][j])
                elif keys[j] == 'total_joining_cpus':
                    self.total_joining_cpus = float(values[0][j])
                elif keys[j] == 'total_joining_cycles':
                    self.total_joining_cycles = float(values[0][j])
                elif keys[j] == 'total_joining_ghz':
                    self.total_joining_ghz = float(values[0][j])
                elif keys[j] == 'total_joining_instructions':
                    self.total_joining_instructions = float(values[0][j])                    
                elif keys[j] == 'total_joining_instructions_per_cycle':
                    self.total_joining_instructions_per_cycle = float(values[0][j])
                elif keys[j] == 'total_joining_l1_misses':
                    self.total_joining_l1_misses = float(values[0][j])
                elif keys[j] == 'total_joining_llc_misses':
                    self.total_joining_llc_misses = float(values[0][j])                    
                elif keys[j] == 'total_joining_task_clock':
                    self.total_joining_task_clock = float(values[0][j])                
                elif keys[j] == 'total_joining_time_usec':
                    self.total_joining_time_usec = float(values[0][j])                
                elif keys[j] == 'total_partitioning_branch_misses':
                    self.total_partitioning_branch_misses = float(values[0][j])                
                elif keys[j] == 'total_partitioning_cpus':
                    self.total_partitioning_cpus = float(values[0][j])                
                elif keys[j] == 'total_partitioning_cycles':
                    self.total_partitioning_cycles = float(values[0][j])                
                elif keys[j] == 'total_partitioning_ghz':
                    self.total_partitioning_ghz = float(values[0][j])                
                elif keys[j] == 'total_partitioning_instructions':
                    self.total_partitioning_instructions = float(values[0][j])                
                elif keys[j] == 'total_partitioning_instructions_per_cycle':
                    self.total_partitioning_instructions_per_cycle = float(values[0][j])                
                elif keys[j] == 'total_partitioning_l1_misses':
                    self.total_partitioning_l1_misses = float(values[0][j])                
                elif keys[j] == 'total_partitioning_llc_misses':
                    self.total_partitioning_llc_misses = float(values[0][j])                
                elif keys[j] == 'total_partitioning_task_clock':
                    self.total_partitioning_task_clock = float(values[0][j])                
                elif keys[j] == 'total_partitioning_time_usec':
                    self.total_partitioning_time_usec = float(values[0][j])                
                elif keys[j] == 'total_merging_branch_misses':
                    self.total_merging_branch_misses = float(values[0][j])                
                elif keys[j] == 'total_merging_cpus':
                    self.total_merging_cpus = float(values[0][j])                
                elif keys[j] == 'total_merging_cycles':
                    self.total_merging_cycles = float(values[0][j])                
                elif keys[j] == 'total_merging_ghz':
                    self.total_merging_ghz = float(values[0][j])                
                elif keys[j] == 'total_merging_instructions':
                    self.total_merging_instructions = float(values[0][j])                
                elif keys[j] == 'total_merging_instructions_per_cycle':
                    self.total_merging_instructions_per_cycle = float(values[0][j])                
                elif keys[j] == 'total_merging_l1_misses':
                    self.total_merging_l1_misses = float(values[0][j])                
                elif keys[j] == 'total_merging_llc_misses':
                    self.total_merging_llc_misses = float(values[0][j])                
                elif keys[j] == 'total_merging_time_usec':
                    self.total_merging_time_usec = float(values[0][j])                
                elif keys[j] == 'total_sorting_branch_misses':
                    self.total_sorting_branch_misses = float(values[0][j])    
                elif keys[j] == 'total_sorting_cpus':
                    self.total_sorting_cpus = float(values[0][j])    
                elif keys[j] == 'total_sorting_cycles':
                    self.total_sorting_cycles = float(values[0][j])    
                elif keys[j] == 'total_sorting_ghz':
                    self.total_sorting_ghz = float(values[0][j])    
                elif keys[j] == 'total_sorting_instructions':
                    self.total_sorting_instructions = float(values[0][j])
                elif keys[j] == 'total_sorting_instructions_per_cycle':
                    self.total_sorting_instructions_per_cycle = float(values[0][j])
                elif keys[j] == 'total_sorting_l1_misses':
                    self.total_sorting_l1_misses = float(values[0][j])
                elif keys[j] == 'total_sorting_llc_misses':
                    self.total_sorting_llc_misses = float(values[0][j])
                elif keys[j] == 'total_sorting_task_clock':
                    self.total_sorting_task_clock = float(values[0][j])
                elif keys[j] == 'total_sorting_time_usec':
                    self.total_sorting_time_usec = float(values[0][j])

class ParsedNonPartitionJoinBenchmarkOutput:
    def __init__(self, benchmark_output_path_in, num_parsed_items_in=1):
        self.benchmark_output_path = benchmark_output_path_in
        self.num_parsed_items = num_parsed_items_in

    def parse_benchmark_parameters(self):
        parts = os.path.split(self.benchmark_output_path)
        file_name = parts[1]
        configs = file_name[26:-4]
        configs_arr = configs.split("_")
        self.datasets = configs_arr[0]
        self.threads = configs_arr[2]
        self.bucket_size = configs_arr[4] 
        self.prefetch_distance = configs_arr[6]


    def parse_benchmark_output(self):
        keys = []
        values = [] 
        with open(self.benchmark_output_path, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            i = 0
            for row in reader:
                if i == 0:
                    keys.extend(row)
                else:
                    values.append(row)
                i = i + 1

        if len(values) > 1 :
            printf('Not supported multiple items for now')
        else:
            for j in range(len(values[0])):
                if keys[j] == 'name':
                    self.name = values[0][j]          
                elif keys[j] == 'iterations':
                    self.iterations = values[0][j]
                elif keys[j] == 'real_time':
                    self.real_time = values[0][j]
                elif keys[j] == 'cpu_time':
                    self.cpu_time = values[0][j]
                elif keys[j] == 'time_unit':
                    self.time_unit = values[0][j]
                elif keys[j] == 'bytes_per_second':
                    self.bytes_per_second = values[0][j]
                elif keys[j] == 'items_per_second':
                    self.items_per_second = values[0][j]
                elif keys[j] == 'label':
                    self.label = values[0][j]
                elif keys[j] == 'error_occurred':
                    self.error_occurred = values[0][j]
                elif keys[j] == 'error_message':
                    self.error_message = values[0][j]
                elif keys[j] == 'num_join_results':
                    self.num_join_results = float(values[0][j])
                elif keys[j] == 'total_algorithm_time_usec':
                    self.total_algorithm_time_usec = float(values[0][j])
                elif keys[j] == 'total_probing_branch_misses':
                    self.total_probing_branch_misses = float(values[0][j])
                elif keys[j] == 'total_probing_cpus':
                    self.total_probing_cpus = float(values[0][j])
                elif keys[j] == 'total_probing_cycles':
                    self.total_probing_cycles = float(values[0][j])
                elif keys[j] == 'total_probing_ghz':
                    self.total_probing_ghz = float(values[0][j])
                elif keys[j] == 'total_probing_instructions':
                    self.total_probing_instructions = float(values[0][j])                    
                elif keys[j] == 'total_probing_instructions_per_cycle':
                    self.total_probing_instructions_per_cycle = float(values[0][j])
                elif keys[j] == 'total_probing_l1_misses':
                    self.total_probing_l1_misses = float(values[0][j])
                elif keys[j] == 'total_probing_llc_misses':
                    self.total_probing_llc_misses = float(values[0][j])                    
                elif keys[j] == 'total_probing_task_clock':
                    self.total_probing_task_clock = float(values[0][j])                
                elif keys[j] == 'total_probing_time_usec':
                    self.total_probing_time_usec = float(values[0][j])                
                elif keys[j] == 'total_building_branch_misses':
                    self.total_building_branch_misses = float(values[0][j])                
                elif keys[j] == 'total_building_cpus':
                    self.total_building_cpus = float(values[0][j])                
                elif keys[j] == 'total_building_cycles':
                    self.total_building_cycles = float(values[0][j])                
                elif keys[j] == 'total_building_ghz':
                    self.total_building_ghz = float(values[0][j])                
                elif keys[j] == 'total_building_instructions':
                    self.total_building_instructions = float(values[0][j])                
                elif keys[j] == 'total_building_instructions_per_cycle':
                    self.total_building_instructions_per_cycle = float(values[0][j])                
                elif keys[j] == 'total_building_l1_misses':
                    self.total_building_l1_misses = float(values[0][j])                
                elif keys[j] == 'total_building_llc_misses':
                    self.total_building_llc_misses = float(values[0][j])                
                elif keys[j] == 'total_building_task_clock':
                    self.total_building_task_clock = float(values[0][j])                
                elif keys[j] == 'total_building_time_usec':
                    self.total_building_time_usec = float(values[0][j])                


class ParsedRadixJoinBenchmarkOutput:
    def __init__(self, benchmark_output_path_in, num_parsed_items_in=1):
        self.benchmark_output_path = benchmark_output_path_in
        self.num_parsed_items = num_parsed_items_in
        
    def parse_benchmark_parameters(self):
        parts = os.path.split(self.benchmark_output_path)
        file_name = parts[1]
        configs = file_name[18:-4]
        configs_arr = configs.split("_")
        self.datasets = configs_arr[0]
        self.threads = configs_arr[2]
        self.num_radix_bits = configs_arr[4] 
        self.num_passes = configs_arr[6]

    def parse_benchmark_output(self):
        keys = []
        values = [] 
        with open(self.benchmark_output_path, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            i = 0
            for row in reader:
                if i == 0:
                    keys.extend(row)
                else:
                    values.append(row)
                i = i + 1

        if len(values) > 1 :
            printf('Not supported multiple items for now')
        else:
            for j in range(len(values[0])):
                if keys[j] == 'name':
                    self.name = values[0][j]          
                elif keys[j] == 'iterations':
                    self.iterations = values[0][j]
                elif keys[j] == 'real_time':
                    self.real_time = values[0][j]
                elif keys[j] == 'cpu_time':
                    self.cpu_time = values[0][j]
                elif keys[j] == 'time_unit':
                    self.time_unit = values[0][j]
                elif keys[j] == 'bytes_per_second':
                    self.bytes_per_second = values[0][j]
                elif keys[j] == 'items_per_second':
                    self.items_per_second = values[0][j]
                elif keys[j] == 'label':
                    self.label = values[0][j]
                elif keys[j] == 'error_occurred':
                    self.error_occurred = values[0][j]
                elif keys[j] == 'error_message':
                    self.error_message = values[0][j]
                elif keys[j] == 'num_join_results':
                    self.num_join_results = float(values[0][j])
                elif keys[j] == 'total_algorithm_time_usec':
                    self.total_algorithm_time_usec = float(values[0][j])
                elif keys[j] == 'total_joining_branch_misses':
                    self.total_joining_branch_misses = float(values[0][j])
                elif keys[j] == 'total_joining_cpus':
                    self.total_joining_cpus = float(values[0][j])
                elif keys[j] == 'total_joining_cycles':
                    self.total_joining_cycles = float(values[0][j])
                elif keys[j] == 'total_joining_ghz':
                    self.total_joining_ghz = float(values[0][j])
                elif keys[j] == 'total_joining_instructions':
                    self.total_joining_instructions = float(values[0][j])                    
                elif keys[j] == 'total_joining_instructions_per_cycle':
                    self.total_joining_instructions_per_cycle = float(values[0][j])
                elif keys[j] == 'total_joining_l1_misses':
                    self.total_joining_l1_misses = float(values[0][j])
                elif keys[j] == 'total_joining_llc_misses':
                    self.total_joining_llc_misses = float(values[0][j])                    
                elif keys[j] == 'total_joining_task_clock':
                    self.total_joining_task_clock = float(values[0][j])                
                elif keys[j] == 'total_joining_time_usec':
                    self.total_joining_time_usec = float(values[0][j])                
                elif keys[j] == 'total_partitioning_branch_misses':
                    self.total_partitioning_branch_misses = float(values[0][j])                
                elif keys[j] == 'total_partitioning_cpus':
                    self.total_partitioning_cpus = float(values[0][j])                
                elif keys[j] == 'total_partitioning_cycles':
                    self.total_partitioning_cycles = float(values[0][j])                
                elif keys[j] == 'total_partitioning_ghz':
                    self.total_partitioning_ghz = float(values[0][j])                
                elif keys[j] == 'total_partitioning_instructions':
                    self.total_partitioning_instructions = float(values[0][j])                
                elif keys[j] == 'total_partitioning_instructions_per_cycle':
                    self.total_partitioning_instructions_per_cycle = float(values[0][j])                
                elif keys[j] == 'total_partitioning_l1_misses':
                    self.total_partitioning_l1_misses = float(values[0][j])                
                elif keys[j] == 'total_partitioning_llc_misses':
                    self.total_partitioning_llc_misses = float(values[0][j])                
                elif keys[j] == 'total_partitioning_task_clock':
                    self.total_partitioning_task_clock = float(values[0][j])                
                elif keys[j] == 'total_partitioning_time_usec':
                    self.total_partitioning_time_usec = float(values[0][j])                
