import sys
import numpy as np
import struct
import os

def prepare_binary_format(input_file_path, key_type_str, payload_type_str, should_be_unique, should_downsample, downsample_ratio, output_key_file_path):
    # read and parse the input data
    if os.path.exists(input_file_path):
        print(input_file_path, " exists and binarizing it should start now ...")
    else:
        print("File does not exist!")
        sys.exit(2)

    key_type = np.uint64
    
    if key_type_str == 'uint64':
        key_type = np.uint64
    elif key_type_str == 'uint32':
        key_type = np.uint32

    keys_list = []
    first_line = 0
    file_input = open(input_file_path, 'r') 
    for line in file_input: 
        if first_line == 1:
            if key_type_str == 'uint64':
                keys_list.append((long)(line.split()[0]))
            elif key_type_str == 'uint32':
                keys_list.append((int)(line.split()[0]))
        first_line = 1
      
    file_input.close()

    if should_be_unique == 0:
        #sort the key data
        keys_list.sort()

    keys_np = np.array(keys_list).astype(key_type)

    if should_be_unique == 1:
        keys_np_unique = np.unique(keys_np, return_index=False, return_inverse=False)
    else: 
        keys_np_unique = keys_np

    # sample the sorted key data
    if should_downsample:
        downsample_ratio_arr = [0.50, 0.25, 0.20, 0.10, 0.05, 0.01]
        downsample_ratio_slice = [2, 4, 5, 10, 20, 100]

        selected_slice = 100
        for i in range(len(downsample_ratio_arr)):
            if downsample_ratio_arr[i] == downsample_ratio:
                selected_slice = downsample_ratio_slice[i]
                break

        keys_np_unique_list = keys_np_unique.tolist()
        sampled_keys_list = keys_np_unique_list[::selected_slice]
        #sampled_keys_list = keys_list[::selected_slice]  
        sampled_keys_list.append(keys_np_unique_list[-1])

        sampled_keys_arr = np.array(sampled_keys_list).astype(key_type)
        print("Count of sampled unique keys: ", len(sampled_keys_arr))
        with open(output_key_file_path, "wb") as f:
            f.write(struct.pack("Q", len(sampled_keys_arr)))
            sampled_keys_arr.tofile(f)
    else:
        #keys_arr = np.array(keys_list).astype(key_type)
        keys_arr = keys_np_unique
        print("Count of unique keys: ", len(keys_arr))
        with open(output_key_file_path, "wb") as f:
            f.write(struct.pack("Q", len(keys_arr)))
            keys_arr.tofile(f)

        #with open('/spinning/sabek/learned_join_datasets/tmp.txt', 'w') as f:
        #    f.write("#KEY, VAL\n")
        #    for item in keys_arr:
        #        f.write(f'{item} {item}\n')

    print(" Finished the binarization process")

    
# You should endup your file key filename, which will be input to the RMI library, with key type, e.g. _uint32, _uint64
#prepare_binary_format("/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v1_int_int_1000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v1_int_int_1000000_key_uint32")
#prepare_binary_format("/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v2_int_int_100000000.txt", 'uint32', 'uint32', 0, 1, 0.50, "/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v2_int_int_100000000_key_uint32")
#prepare_binary_format("/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v3_int_int_100000000.txt", 'uint32', 'uint32', 0, 0, 0.50, "/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v3_int_int_100000000_key_uint32")

#prepare_binary_format("/spinning/sabek/learned_join_datasets/r_UNIQUE_v2_uint32_uint32_32000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/r_UNIQUE_v2_uint32_uint32_32000000_key_uint32")
#prepare_binary_format("/spinning/sabek/learned_join_datasets/s_UNIQUE_v2_uint32_uint32_32000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/s_UNIQUE_v2_uint32_uint32_32000000_key_uint32")

#prepare_binary_format("/spinning/sabek/learned_join_datasets/r_UNIQUE_v3_uint32_uint32_128000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/r_UNIQUE_v3_uint32_uint32_128000000_key_uint32")
#prepare_binary_format("/spinning/sabek/learned_join_datasets/s_UNIQUE_v3_uint32_uint32_128000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/s_UNIQUE_v3_uint32_uint32_128000000_key_uint32")

#prepare_binary_format("/spinning/sabek/learned_join_datasets/r_UNIQUE_v5_uint32_uint32_640000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/r_UNIQUE_v5_uint32_uint32_640000000_key_uint32")
#prepare_binary_format("/spinning/sabek/learned_join_datasets/s_UNIQUE_v5_uint32_uint32_640000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/s_UNIQUE_v5_uint32_uint32_640000000_key_uint32")

#prepare_binary_format("/spinning/sabek/learned_join_datasets/r_UNIQUE_v8_uint32_uint32_1664000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/r_UNIQUE_v8_uint32_uint32_1664000000_key_uint32")
#prepare_binary_format("/spinning/sabek/learned_join_datasets/s_UNIQUE_v8_uint32_uint32_1664000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/s_UNIQUE_v8_uint32_uint32_1664000000_key_uint32")

#prepare_binary_format("/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000_key_uint32")
#prepare_binary_format("/spinning/sabek/learned_join_datasets/s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000_key_uint32")

#prepare_binary_format("/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000_key_uint32")
#prepare_binary_format("/spinning/sabek/learned_join_datasets/s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000_key_uint32")

#prepare_binary_format("/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000_key_uint32")
#prepare_binary_format("/spinning/sabek/learned_join_datasets/s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000_key_uint32")

#prepare_binary_format("/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000_key_uint32")
#prepare_binary_format("/spinning/sabek/learned_join_datasets/s_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/s_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000_key_uint32")


prepare_binary_format("/spinning/sabek/learned_join_datasets/r_SEQ_HOLE_v1_uint32_uint32_16000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/r_SEQ_HOLE_v1_uint32_uint32_16000000_key_uint32")
prepare_binary_format("/spinning/sabek/learned_join_datasets/s_SEQ_HOLE_v1_uint32_uint32_16000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/s_SEQ_HOLE_v1_uint32_uint32_16000000_key_uint32")

prepare_binary_format("/spinning/sabek/learned_join_datasets/r_SEQ_HOLE_v2_uint32_uint32_32000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/r_SEQ_HOLE_v2_uint32_uint32_32000000_key_uint32")
prepare_binary_format("/spinning/sabek/learned_join_datasets/s_SEQ_HOLE_v2_uint32_uint32_32000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/s_SEQ_HOLE_v2_uint32_uint32_32000000_key_uint32")

prepare_binary_format("/spinning/sabek/learned_join_datasets/r_SEQ_HOLE_v3_uint32_uint32_128000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/r_SEQ_HOLE_v3_uint32_uint32_128000000_key_uint32")
prepare_binary_format("/spinning/sabek/learned_join_datasets/s_SEQ_HOLE_v3_uint32_uint32_128000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/s_SEQ_HOLE_v3_uint32_uint32_128000000_key_uint32")

prepare_binary_format("/spinning/sabek/learned_join_datasets/r_SEQ_HOLE_v5_uint32_uint32_640000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/r_SEQ_HOLE_v5_uint32_uint32_640000000_key_uint32")
prepare_binary_format("/spinning/sabek/learned_join_datasets/s_SEQ_HOLE_v5_uint32_uint32_640000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/s_SEQ_HOLE_v5_uint32_uint32_640000000_key_uint32")

prepare_binary_format("/spinning/sabek/learned_join_datasets/r_SEQ_HOLE_v8_uint32_uint32_1664000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/r_SEQ_HOLE_v8_uint32_uint32_1664000000_key_uint32")
prepare_binary_format("/spinning/sabek/learned_join_datasets/s_SEQ_HOLE_v8_uint32_uint32_1664000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/s_SEQ_HOLE_v8_uint32_uint32_1664000000_key_uint32")


prepare_binary_format("/spinning/sabek/learned_join_datasets/r_UNIFORM_v1_uint32_uint32_16000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/r_UNIFORM_v1_uint32_uint32_16000000_key_uint32")
prepare_binary_format("/spinning/sabek/learned_join_datasets/s_UNIFORM_v1_uint32_uint32_16000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/s_UNIFORM_v1_uint32_uint32_16000000_key_uint32")

prepare_binary_format("/spinning/sabek/learned_join_datasets/r_UNIFORM_v2_uint32_uint32_32000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/r_UNIFORM_v2_uint32_uint32_32000000_key_uint32")
prepare_binary_format("/spinning/sabek/learned_join_datasets/s_UNIFORM_v2_uint32_uint32_32000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/s_UNIFORM_v2_uint32_uint32_32000000_key_uint32")

prepare_binary_format("/spinning/sabek/learned_join_datasets/r_UNIFORM_v3_uint32_uint32_128000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/r_UNIFORM_v3_uint32_uint32_128000000_key_uint32")
prepare_binary_format("/spinning/sabek/learned_join_datasets/s_UNIFORM_v3_uint32_uint32_128000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/s_UNIFORM_v3_uint32_uint32_128000000_key_uint32")

prepare_binary_format("/spinning/sabek/learned_join_datasets/r_UNIFORM_v5_uint32_uint32_640000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/r_UNIFORM_v5_uint32_uint32_640000000_key_uint32")
prepare_binary_format("/spinning/sabek/learned_join_datasets/s_UNIFORM_v5_uint32_uint32_640000000.txt", 'uint32', 'uint32', 1, 0, 0.01, "/spinning/sabek/learned_join_datasets/s_UNIFORM_v5_uint32_uint32_640000000_key_uint32")
