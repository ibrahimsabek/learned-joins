/* @version $Id: cpu_mapping.c 4548 2013-12-07 16:05:16Z bcagri $ */

#include <stdio.h>  /* FILE, fopen */
#include <stdlib.h> /* exit, perror */
#include <unistd.h> /* sysconf */
#include <numaif.h> /* get_mempolicy() */

#include "imv-operator/cpu_mapping.h"
#include "imv-operator/prefetch.hpp"
/** \internal
 * @{
 */

#define MAX_NODES 512

static int inited = 0;
static int max_cpus;
static int node_mapping[MAX_NODES];

/**
 * Initializes the cpu mapping from the file defined by CUSTOM_CPU_MAPPING.
 * The mapping used for our machine Intel L5520 is = "8 0 1 2 3 8 9 10 11".
 */
static int init_mappings_from_file() {
  FILE* cfg;
  int i, j, begin, end;

  cfg = fopen("cpu-mapping.txt", "r");
  if (cfg != NULL) {
    if (fscanf(cfg, "%d", &max_cpus) <= 0) {
      perror("Could not parse input!\n");
    }

    if (max_cpus >= 0) {  // 4 0,1,2,3,
      for (i = 0; i < max_cpus; i++) {
        if (fscanf(cfg, "%d,", &node_mapping[i]) <= 0) {
          perror("Could not parse input!\n");
        }
      }
    } else {  // -4 0-3,
      max_cpus = -max_cpus;
      for (i = 0; i < max_cpus;) {
        if (fscanf(cfg, "%d-%d,", &begin, &end) <= 0) {
          perror("Could not parse input!\n");
        } else {
          for (j = begin; j <= end; ++j) {
            node_mapping[i] = j;
            i++;
          }
        }
      }
    }

    printf("all cpus,  %d, there are, ", max_cpus);
    for (i = 0; i < max_cpus; ++i) {
      printf("%d\t", node_mapping[i]);
    }
    puts("end");

    fclose(cfg);
    return 1;
  }
  perror("Could not open input!\n");
  fclose(cfg);
  /* perror("Custom cpu mapping file not found!\n"); */
  return 0;
}

/**
 * Try custom cpu mapping file first, if does not exist then round-robin
 * initialization among available CPUs reported by the system.
 */
static void init_mappings() {
  if (init_mappings_from_file() == 0) {
    int i;

    max_cpus = sysconf(_SC_NPROCESSORS_ONLN);
    for (i = 0; i < max_cpus; i++) {
      node_mapping[i] = i;
    }
  }
}

/** @} */

/**
 * Returns SMT aware logical to physical CPU mapping for a given thread id.
 */
int get_cpu_id(int thread_id) {
  if (!inited) {
    init_mappings();
    inited = 1;
  }

  return node_mapping[thread_id % max_cpus];
}

/* TODO: These are just place-holder implementations. */
/**
 * Topology of Intel E5-4640
 node 0 cpus: 0 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60
 node 1 cpus: 1 5 9 13 17 21 25 29 33 37 41 45 49 53 57 61
 node 2 cpus: 2 6 10 14 18 22 26 30 34 38 42 46 50 54 58 62
 node 3 cpus: 3 7 11 15 19 23 27 31 35 39 43 47 51 55 59 63
*/
#define INTEL_E5 1

static int numa[][16] = {
    {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60},
    {1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61},
    {2, 6, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62},
    {3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63}};

int get_numa_id(int mytid) {
#if INTEL_E5
  int ret = 0;

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 16; j++)
      if (numa[i][j] == mytid) {
        ret = i;
        break;
      }

  return ret;
#else
  return 0;
#endif
}

int get_num_numa_regions(void) {
/* TODO: FIXME automate it from the system config. */
#if INTEL_E5
  return 4;
#else
  return 1;
#endif
}

//int get_numa_node_of_address(void* ptr) {
//  int numa_node = -1;
//#if !KNL
//  get_mempolicy(&numa_node, NULL, 0, ptr, MPOL_F_NODE | MPOL_F_ADDR);
//#endif
//  return numa_node;
//}
