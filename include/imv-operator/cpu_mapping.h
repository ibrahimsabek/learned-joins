/**
 * @file    cpu_mapping.h
 * @author  Cagri Balkesen <cagri.balkesen@inf.ethz.ch>
 * @date    Tue May 22 16:35:12 2012
 * @version $Id: cpu_mapping.h 4548 2013-12-07 16:05:16Z bcagri $
 * 
 * @brief  Provides cpu mapping utility function.
 * 
 * 
 */
#ifndef CPU_MAPPING_H
#define CPU_MAPPING_H

/** 
 * @defgroup cpumapping CPU mapping tool
 * @{
 */

/** 
 * if the custom cpu mapping file exists, logical to physical mappings are
 * initialized from that file, otherwise it will be round-robin 
 */
#ifndef CUSTOM_CPU_MAPPING
#define CUSTOM_CPU_MAPPING "cpu-mapping.txt"
#endif

/**
 * Returns SMT aware logical to physical CPU mapping for a given thread id.
 */
int get_cpu_id(int thread_id);

/** 
 * Returns the NUMA id of the given thread id returned from get_cpu_id(int)
 * 
 * @param mytid 
 * 
 * @return 
 */
int
get_numa_id(int mytid);

/** 
 * Returns number of NUMA regions.
 * 
 * @return 
 */
int
get_num_numa_regions(void);

/**
 * Returns the NUMA-node id of a given memory address
 */
//int
//get_numa_node_of_address(void * ptr);

/** @} */

#endif /* CPU_MAPPING_H */
