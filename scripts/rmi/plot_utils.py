import matplotlib.pyplot as plt
import numpy as np
import struct
import os
import sys

def plt_x_y_values(input_file_path, output_folder_path, output_fig_file_name):
    
    # read and parse the input data
    if os.path.exists(input_file_path):
        print(input_file_path, " exists!")
    else:
        print("File does not exist!")
        sys.exit(2)

    keys_list = []
    visits_list = []    
    with open(input_file_path, "r") as f:
        zeros_count = 0
        for line in f: 
            curr_key = (np.uint64)(line.split()[0])
            curr_visit = (np.uint64)(line.split()[1])
#            if curr_visit > 10:
#                print("This key is frequent: ", curr_key, " its visits are: ", curr_visit)
#            elif curr_visit == 0:
#                zeros_count = zeros_count + 1
            keys_list.append(curr_key) 
            visits_list.append(curr_visit)
#        print("Zeros count is ", zeros_count)
    keys_arr = np.array(keys_list).astype(np.uint64)
    visits_arr = np.array(visits_list).astype(np.uint64)

    print(" Starting the plotting process ...")

    plt.plot(keys_arr, visits_arr, 'ro')
    #plt.bar(keys_arr, visits_arr)
    plt.savefig(output_folder_path + output_fig_file_name)
    plt.close()

    print(" Finished the plotting process")       
    
def plt_histogram(is_cdf_plot, input_file_path, key_type_str, is_rmi_keys_only, output_folder_path, output_fig_file_name):

    # read and parse the input data
    if os.path.exists(input_file_path):
        print(input_file_path, " exists!")
    else:
        print("File does not exist!")
        sys.exit(2)

    key_type = np.uint64

    if key_type_str == 'uint64':
        key_type = np.uint64
    elif key_type_str == 'uint32':
        key_type = np.uint32

    if is_rmi_keys_only:
        print("Loading the binary file ....")
        with open(input_file_path, "rb") as f:
            keys_count = struct.unpack("Q", f.read(8))
            print("Count of expected elements to be loaded: ", keys_count[0])

            f.seek(8, os.SEEK_SET)
            keys_arr = np.fromfile(f, dtype=key_type)

    else:
        print("Loading the .txt input data file ....")
        keys_list = []
        first_line = 0
        with open(input_file_path, "r") as f:
            for line in f: 
                if first_line == 1:
                    if key_type_str == 'uint64':
                        keys_list.append((long)(line.split()[0]))
                    elif key_type_str == 'uint32':
                        keys_list.append((int)(line.split()[0]))
                first_line = 1

            keys_arr = np.array(keys_list).astype(key_type)

    print("Count of loaded elements: ", len(keys_arr))

    if is_cdf_plot:
        print("Plotting the CDF now ....")
        n, bins, patches = plt.hist(x=keys_arr, bins='auto', density=True, cumulative=True, label='CDF',
                histtype='step', alpha=0.7, rwidth=0.85, color='k')
        plt.grid(axis='y', alpha=0.75)
        plt.xlabel('Key')
        plt.ylabel('Probability')
        #plt.title(output_fig_file_name)
    else:
        print("Plotting the histogram now ....")
        n, bins, patches = plt.hist(x=keys_arr, bins='auto', color='#0504aa',
                                    alpha=0.7, rwidth=0.85)
        plt.grid(axis='y', alpha=0.75)
        plt.xlabel('Key')
        plt.ylabel('Frequency')
        #plt.title(output_fig_file_name)
        maxfreq = n.max()
        #plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)

    plt.savefig(output_folder_path + output_fig_file_name)

    print(" Finished the plotting process")

#LOGNORMAL
############
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v1_uint32_uint32_segma_1_16000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_LOGNORMAL_v1_uint32_uint32_segma_1_16000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v9_uint32_uint32_segma_1_1920000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_LOGNORMAL_v9_uint32_uint32_segma_1_1920000000_hist.png")    
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_LOGNORMAL_v1_uint32_uint32_segma_1_16000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_LOGNORMAL_v1_uint32_uint32_segma_1_16000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_LOGNORMAL_v9_uint32_uint32_segma_1_1920000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_LOGNORMAL_v9_uint32_uint32_segma_1_1920000000_hist.png")    

#UNIQUE
############
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_UNIQUE_v1_uint32_uint32_16000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_UNIQUE_v1_uint32_uint32_16000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_UNIQUE_v2_uint32_uint32_32000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_UNIQUE_v2_uint32_uint32_32000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_UNIQUE_v3_uint32_uint32_128000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_UNIQUE_v3_uint32_uint32_128000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_UNIQUE_v5_uint32_uint32_640000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_UNIQUE_v5_uint32_uint32_640000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_UNIQUE_v8_uint32_uint32_1664000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_UNIQUE_v8_uint32_uint32_1664000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_UNIQUE_v9_uint32_uint32_1920000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_UNIQUE_v9_uint32_uint32_1920000000_hist.png")    
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_UNIQUE_v1_uint32_uint32_16000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_UNIQUE_v1_uint32_uint32_16000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_UNIQUE_v2_uint32_uint32_32000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_UNIQUE_v2_uint32_uint32_32000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_UNIQUE_v3_uint32_uint32_128000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_UNIQUE_v3_uint32_uint32_128000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_UNIQUE_v5_uint32_uint32_640000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_UNIQUE_v5_uint32_uint32_640000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_UNIQUE_v8_uint32_uint32_1664000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_UNIQUE_v8_uint32_uint32_1664000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_UNIQUE_v9_uint32_uint32_1920000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_UNIQUE_v9_uint32_uint32_1920000000_hist.png")    

#SEQ_HOLE
############
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_SEQ_HOLE_v1_uint32_uint32_16000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_SEQ_HOLE_v1_uint32_uint32_16000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_SEQ_HOLE_v2_uint32_uint32_32000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_SEQ_HOLE_v2_uint32_uint32_32000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_SEQ_HOLE_v3_uint32_uint32_128000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_SEQ_HOLE_v3_uint32_uint32_128000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_SEQ_HOLE_v5_uint32_uint32_640000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_SEQ_HOLE_v5_uint32_uint32_640000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_SEQ_HOLE_v8_uint32_uint32_1664000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_SEQ_HOLE_v8_uint32_uint32_1664000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_SEQ_HOLE_v9_uint32_uint32_1920000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_SEQ_HOLE_v9_uint32_uint32_1920000000_hist.png")    
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_SEQ_HOLE_v1_uint32_uint32_16000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_SEQ_HOLE_v1_uint32_uint32_16000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_SEQ_HOLE_v2_uint32_uint32_32000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_SEQ_HOLE_v2_uint32_uint32_32000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_SEQ_HOLE_v3_uint32_uint32_128000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_SEQ_HOLE_v3_uint32_uint32_128000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_SEQ_HOLE_v5_uint32_uint32_640000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_SEQ_HOLE_v5_uint32_uint32_640000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_SEQ_HOLE_v8_uint32_uint32_1664000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_SEQ_HOLE_v8_uint32_uint32_1664000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_SEQ_HOLE_v9_uint32_uint32_1920000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_SEQ_HOLE_v9_uint32_uint32_1920000000_hist.png")    

#UNIFORM
############
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_UNIFORM_v1_uint32_uint32_16000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_UNIFORM_v1_uint32_uint32_16000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_UNIFORM_v2_uint32_uint32_32000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_UNIFORM_v2_uint32_uint32_32000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_UNIFORM_v3_uint32_uint32_128000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_UNIFORM_v3_uint32_uint32_128000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_UNIFORM_v5_uint32_uint32_640000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_UNIFORM_v5_uint32_uint32_640000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_UNIFORM_v8_uint32_uint32_1664000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_UNIFORM_v8_uint32_uint32_1664000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/r_UNIFORM_v9_uint32_uint32_1920000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_UNIFORM_v9_uint32_uint32_1920000000_hist.png")    
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_UNIFORM_v1_uint32_uint32_16000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_UNIFORM_v1_uint32_uint32_16000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_UNIFORM_v2_uint32_uint32_32000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_UNIFORM_v2_uint32_uint32_32000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_UNIFORM_v3_uint32_uint32_128000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_UNIFORM_v3_uint32_uint32_128000000_hist.png")
plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_UNIFORM_v5_uint32_uint32_640000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_UNIFORM_v5_uint32_uint32_640000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_UNIFORM_v8_uint32_uint32_1664000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_UNIFORM_v8_uint32_uint32_1664000000_hist.png")
#plt_histogram(0, "/spinning/sabek/learned_join_datasets/s_UNIFORM_v9_uint32_uint32_1920000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "s_UNIFORM_v9_uint32_uint32_1920000000_hist.png")    

#SOSD
############
plt_histogram(0, "/spinning/sabek/learned_join_datasets_sosd/fb_200M_uint64.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "fb_200M_uint64_hist.png")
plt_histogram(1, "/spinning/sabek/learned_join_datasets_sosd/fb_200M_uint64.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "fb_200M_uint64_cdf.png")




#plt_histogram(0, "../learned_join_data/r_LOGNORMAL_v1_int_int_1000000_key_uint32", 'uint32', 1, "../learned_join_plots/", "r_LOGNORMAL_v1_int_int_1000000_key_uint32_hist.png")
#plt_histogram(0, "../learned_join_data/r_LOGNORMAL_v1_int_int_1000000.txt", 'uint32', 0, "../learned_join_plots/", "r_LOGNORMAL_v1_int_int_1000000_hist.png")
#plt_histogram(0, "../learned_join_data/s_LOGNORMAL_v1_int_int_1000000.txt", 'uint32', 0, "../learned_join_plots/", "s_LOGNORMAL_v1_int_int_1000000_hist.png")

#plt_histogram(0, "../learned_join_data/r_LOGNORMAL_v2_int_int_100000000_key_uint32", 'uint32', 1, "../learned_join_plots/", "r_LOGNORMAL_v2_int_int_100000000_key_uint32_hist.png")
#plt_histogram(0, "../learned_join_data/r_LOGNORMAL_v2_int_int_100000000.txt", 'uint32', 0, "../learned_join_plots/", "r_LOGNORMAL_v2_int_int_100000000_hist.png")
#plt_histogram(0, "../learned_join_data/s_LOGNORMAL_v2_int_int_100000000.txt", 'uint32', 0, "../learned_join_plots/", "s_LOGNORMAL_v2_int_int_100000000_hist.png")
#plt_histogram(1, "/spinning/sabek/learned_join_datasets/r_LOGNORMAL_v1_uint32_uint32_segma_1_16000000.txt", 'uint32', 0, "/spinning/sabek/learned_join_plots/", "r_LOGNORMAL_v1_uint32_uint32_segma_1_16000000_cdf.png")
#plt_histogram(1, "../learned_join_data/r_LOGNORMAL_v2_int_int_100000000_key_uint32", 'uint32', 1, "../learned_join_plots/", "r_LOGNORMAL_v2_int_int_100000000_key_uint32_cdf.png")
#plt_histogram(1, "../learned_join_data/r_LOGNORMAL_v2_int_int_100000000.txt", 'uint32', 0, "../learned_join_plots/", "r_LOGNORMAL_v2_int_int_100000000_cdf.png")
#plt_histogram(1, "../learned_join_data/s_LOGNORMAL_v2_int_int_100000000.txt", 'uint32', 0, "../learned_join_plots/", "s_LOGNORMAL_v2_int_int_100000000_cdf.png")


#plt_x_y_values("/Users/ibrahimsabek/Downloads/r_LOGNORMAL_v2_int_int_100000000_eth_build_visits.txt", "/Users/ibrahimsabek/Downloads/",  "r_LOGNORMAL_v2_int_int_100000000_eth_build_visits_1.png")
#plt_x_y_values("/Users/ibrahimsabek/Downloads/r_LOGNORMAL_v2_int_int_100000000_eth_probe_visits.txt", "/Users/ibrahimsabek/Downloads/",  "r_LOGNORMAL_v2_int_int_100000000_eth_probe_visits_1.png")
#plt_x_y_values("/Users/ibrahimsabek/Downloads/r_LOGNORMAL_v2_int_int_100000000_rmi_build_visits_2.txt", "/Users/ibrahimsabek/Downloads/",  "r_LOGNORMAL_v2_int_int_100000000_rmi_build_visits_1.png")
#plt_x_y_values("/Users/ibrahimsabek/Downloads/r_LOGNORMAL_v2_int_int_100000000_rmi_probe_visits_2.txt", "/Users/ibrahimsabek/Downloads/",  "r_LOGNORMAL_v2_int_int_100000000_rmi_probe_visits_1.png")
#plt_x_y_values("../learned_join_debug/r_LOGNORMAL_v1_int_int_1000000_eth_build_visits.txt", "../learned_join_plots/",  "r_LOGNORMAL_v1_int_int_1000000_eth_build_visits.png")
#plt_x_y_values("../learned_join_debug/r_LOGNORMAL_v1_int_int_1000000_eth_probe_visits.txt", "../learned_join_plots/",  "r_LOGNORMAL_v1_int_int_1000000_eth_probe_visits.png")

#plt_x_y_values("../learned_join_debug/r_LOGNORMAL_v2_int_int_100000000_eth_build_visits.txt", "../learned_join_plots/",  "r_LOGNORMAL_v2_int_int_100000000_eth_build_visits.png")
#plt_x_y_values("../learned_join_debug/r_LOGNORMAL_v2_int_int_100000000_eth_probe_visits.txt", "../learned_join_plots/",  "r_LOGNORMAL_v2_int_int_100000000_eth_probe_visits.png")
#plt_x_y_values("../learned_join_debug/r_LOGNORMAL_v2_int_int_100000000_rmi_build_visits.txt", "../learned_join_plots/",  "r_LOGNORMAL_v2_int_int_100000000_rmi_build_visits.png")
#plt_x_y_values("../learned_join_debug/r_LOGNORMAL_v2_int_int_100000000_rmi_probe_visits.txt", "../learned_join_plots/",  "r_LOGNORMAL_v2_int_int_100000000_rmi_probe_visits.png")

#plt_x_y_values("../learned_join_debug/r_LOGNORMAL_v2_int_int_100000000_eth_build_keys_hashes.txt", "../learned_join_plots/",  "r_LOGNORMAL_v2_int_int_100000000_eth_build_keys_hashes.png")
#plt_x_y_values("../learned_join_debug/r_LOGNORMAL_v2_int_int_100000000_eth_probe_keys_hashes.txt", "../learned_join_plots/",  "r_LOGNORMAL_v2_int_int_100000000_eth_probe_keys_hashes.png")
