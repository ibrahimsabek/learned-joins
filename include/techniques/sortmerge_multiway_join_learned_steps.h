#pragma once

#include "techniques/sortmerge_multiway_join_eth_steps.h"
#include "utils/learned_sort_for_sort_merge.h"

template<typename KeyType, typename PayloadType, typename JoinThreadType>
class LearnedSortMergeMultiwayJoinSteps : public ETHSortMergeMultiwayJoinSteps<KeyType, PayloadType, JoinThreadType>{
 public:
  
  void partition_phase(Relation<KeyType, PayloadType> *** relRparts, Relation<KeyType, PayloadType> *** relSparts, JoinThreadType * in_args) 
  {
    //NOTE: Not implemented
 /*     
    ETHSortMergeMultiwayJoinThread<KeyType, PayloadType> * args = (ETHSortMergeMultiwayJoinThread<KeyType, PayloadType> *) in_args;
    const int PARTFANOUT = args->joincfg->PARTFANOUT;
    const int NRADIXBITS = log2(PARTFANOUT);

    Relation<KeyType, PayloadType> ** partsR = (Relation<KeyType, PayloadType> **)
            alloc_aligned(PARTFANOUT * sizeof(Relation<KeyType, PayloadType> *));
    Relation<KeyType, PayloadType> ** partsS = (Relation<KeyType, PayloadType> **)
            alloc_aligned(PARTFANOUT * sizeof(Relation<KeyType, PayloadType> *));

    // note: only free prels[0] when releasing memory
    Relation<KeyType, PayloadType> * prels = (Relation<KeyType, PayloadType> *) malloc(2 * PARTFANOUT * sizeof(Relation<KeyType, PayloadType>));
    for(int i = 0; i < PARTFANOUT; i++) {
        partsR[i] = prels + i;
        partsS[i] = prels + PARTFANOUT + i;
    }

    Relation<KeyType, PayloadType> relR, relS;
    Relation<KeyType, PayloadType> tmpR, tmpS;
    relR.tuples     = args->relR;
    relR.num_tuples = args->numR;
    relS.tuples     = args->relS;
    relS.num_tuples = args->numS;
    tmpR.tuples     = args->tmp_partR;
    tmpR.num_tuples = args->numR;
    tmpS.tuples     = args->tmp_partS;
    tmpS.num_tuples = args->numS;

    // after partitioning tmpR, tmpS holds the partitioned data
    int bitshift = ceil(log2(relR.num_tuples * args->nthreads)) - 1;
    if(args->nthreads == 1)
        bitshift = bitshift - NRADIXBITS + 1;
    else {
#if SKEW_HANDLING
        // NOTE: Special to skew handling code, we must set the radix bits in a
           way that the MSB-Radix partitioning results in range partitioning
           assuming keys are dense. 
        bitshift = bitshift - NRADIXBITS + 1;
#else
        bitshift = bitshift - NRADIXBITS - 1;
#endif
    }

    //DEBUGMSG(1, "[INFO ] bitshift = %d\n", bitshift);

    SortMergeJoinSteps<KeyType, PayloadType, JoinThreadType>::partition_relation_optimized_V2(partsR, &relR, &tmpR, NRADIXBITS, bitshift);
    SortMergeJoinSteps<KeyType, PayloadType, JoinThreadType>::partition_relation_optimized_V2(partsS, &relS, &tmpS, NRADIXBITS, bitshift);

    // return parts
    *relRparts = partsR;
    *relSparts = partsS;
    */
  }
  
  void sorting_phase(Relation<KeyType, PayloadType> ** relRparts, Relation<KeyType, PayloadType> ** relSparts, JoinThreadType * in_args) 
  {
    LearnedSortMergeMultiwayJoinThread<KeyType, PayloadType> * args = (LearnedSortMergeMultiwayJoinThread<KeyType, PayloadType> *) in_args;  
    const int PARTFANOUT = args->joincfg->PARTFANOUT;
    const int learnedsortflag = args->joincfg->LEARNEDSORT;

    int32_t my_tid = args->my_tid;

    args->threadrelchunks[my_tid] = (RelationPair<KeyType, PayloadType> *)
                                  alloc_aligned(PARTFANOUT * sizeof(RelationPair<KeyType, PayloadType>));

    // Sorting R relation
    Tuple<KeyType, PayloadType> * inptrR  = args->tmp_partR;
    Tuple<KeyType, PayloadType> * outptrR = args->tmp_sortR;

    if(learnedsortflag)
    {
#ifdef USE_LEARNED_SORT_FOR_SORT_MERGE
        printf("here for sorting \n");
        learned_sort_for_sort_merge::sort_avx(outptrR, args->rmi_r, args->NUM_MINOR_BCKT_PER_MAJOR_BCKT_r, args->MINOR_BCKTS_OFFSET_r, args->TOT_NUM_MINOR_BCKTS_r, 
                                    args->INPUT_SZ_r, inptrR, args->numR - args->tmp_total_repeatedKeysCountsR, 
                                    args->tmp_minor_bckts_r, args->tmp_minor_bckt_sizes_r,
                                    args->tmp_spill_bucket_r, args->sorted_spill_bucket_r, 
                                    args->tmp_total_repeatedKeysCountsR, args->tmp_repeatedKeysPredictedRanksR, 
                                    args->tmp_repeatedKeysPredictedRanksCountsR, args->tmp_repeatedKeysCountsR, my_tid, 0);
#else
        int64_t * inputR  = (int64_t*)(inptrR);
        int64_t * outputR = (int64_t*)(outptrR);

#ifdef USE_LEARNED_SORT_AVX
        learned_sort::sort_avx(inputR, inputR + args->numR);
#else
        learned_sort::sort(inputR, inputR + args->numR, my_tid, 0);
#endif
        inptrR = (Tuple<KeyType, PayloadType> *)(outputR);
        outptrR = (Tuple<KeyType, PayloadType> *)(inputR);
#endif
    }
    else
        avxsort_tuples<KeyType, PayloadType>(&inptrR, &outptrR, args->numR);

    args->threadrelchunks[my_tid][0].R.tuples     = outptrR;
    args->threadrelchunks[my_tid][0].R.num_tuples = args->numR;


    // Sorting S relation
   Tuple<KeyType, PayloadType> * inptrS  = args->tmp_partS;
   Tuple<KeyType, PayloadType> * outptrS = args->tmp_sortS;

    if(learnedsortflag)
    {
#ifdef USE_LEARNED_SORT_FOR_SORT_MERGE
        learned_sort_for_sort_merge::sort_avx(outptrS, args->rmi_s, args->NUM_MINOR_BCKT_PER_MAJOR_BCKT_s, args->MINOR_BCKTS_OFFSET_s, args->TOT_NUM_MINOR_BCKTS_s, 
                                    args->INPUT_SZ_s, inptrS, args->numS - args->tmp_total_repeatedKeysCountsS, 
                                    args->tmp_minor_bckts_s, args->tmp_minor_bckt_sizes_s,
                                    args->tmp_spill_bucket_s, args->sorted_spill_bucket_s, 
                                    args->tmp_total_repeatedKeysCountsS, args->tmp_repeatedKeysPredictedRanksS, 
                                    args->tmp_repeatedKeysPredictedRanksCountsS, args->tmp_repeatedKeysCountsS, my_tid, 0);
#else
        int64_t * inputS  = (int64_t*)(inptrS);
        int64_t * outputS = (int64_t*)(outptrS);

#ifdef USE_LEARNED_SORT_AVX
        learned_sort::sort_avx(inputS, inputS + args->numS);
#else
        learned_sort::sort(inputS, inputS + args->numS, my_tid, 0);
#endif
        inptrS = (Tuple<KeyType, PayloadType> *)(outputS);
        outptrS = (Tuple<KeyType, PayloadType> *)(inputS);
#endif
    }
    else
        avxsort_tuples<KeyType, PayloadType>(&inptrS, &outptrS, args->numS);

    args->threadrelchunks[my_tid][0].S.tuples = outptrS;
    args->threadrelchunks[my_tid][0].S.num_tuples = args->numS;
  }

  void multiwaymerge_phase(int numaregionid, Relation<KeyType, PayloadType> ** relRparts, Relation<KeyType, PayloadType> ** relSparts, JoinThreadType * in_args,
                                             Relation<KeyType, PayloadType> * mergedRelR, Relation<KeyType, PayloadType> * mergedRelS) 
  {
  //NOTE: Not implemented    
/*      
    ETHSortMergeMultiwayJoinThread<KeyType, PayloadType> * args = (ETHSortMergeMultiwayJoinThread<KeyType, PayloadType> *) in_args;  
    const int PARTFANOUT = args->joincfg->PARTFANOUT;
    const int scalarmergeflag = args->joincfg->SCALARMERGE;

    int32_t my_tid = args->my_tid;
    uint64_t mergeRtotal = 0, mergeStotal = 0;
    Tuple<KeyType, PayloadType> * tmpoutR = NULL;
    Tuple<KeyType, PayloadType> * tmpoutS = NULL;

    if(args->nthreads == 1) {
        // single threaded execution; no multi-way merge.
        for(int i = 0; i < PARTFANOUT; i ++) {
            RelationPair<KeyType, PayloadType> * rels = & args->threadrelchunks[my_tid][i];
            mergeRtotal += rels->R.num_tuples;
            mergeStotal += rels->S.num_tuples;

            // evaluate join between each sorted part 
            relRparts[i]->tuples = rels->R.tuples;
            relRparts[i]->num_tuples = rels->R.num_tuples;
            relSparts[i]->tuples = rels->S.tuples;
            relSparts[i]->num_tuples = rels->S.num_tuples;
        }
    }
    else 
    {
        uint32_t       j;
        const uint32_t perthread   = PARTFANOUT / args->nthreads;

        // multi-threaded execution
        // merge remote relations and bring to local memory
        const uint32_t start = my_tid * perthread;
        const uint32_t end = start + perthread;

        Relation<KeyType, PayloadType> * Rparts[PARTFANOUT];
        Relation<KeyType, PayloadType> * Sparts[PARTFANOUT];

        // compute the size of merged relations to be stored locally
        uint32_t f = 0;
        for(j = start; j < end; j ++) {
            for(int i = 0; i < args->nthreads; i ++) {
                //uint32_t tid = (my_tid + i) % args->nthreads;
                uint32_t tid = get_numa_shuffle_strategy(my_tid, i, args->nthreads);
                //printf("SHUF %d %d --> %d\n", i, my_tid, tid);
                RelationPair<KeyType, PayloadType> * rels = & args->threadrelchunks[tid][j];
                //fprintf(stdout, "TID=%d Part-%d-size = %d\n", my_tid, f, rels->S.num_tuples);
                Rparts[f] = & rels->R;
                Sparts[f] = & rels->S;
                f++;

                mergeRtotal += rels->R.num_tuples;
                mergeStotal += rels->S.num_tuples;
            }
        }

        // allocate memory at local node for temporary merge results
        tmpoutR = (Tuple<KeyType, PayloadType> *) alloc_aligned_threadlocal(mergeRtotal*sizeof(Tuple<KeyType, PayloadType>));
        tmpoutS = (Tuple<KeyType, PayloadType> *) alloc_aligned_threadlocal(mergeStotal*sizeof(Tuple<KeyType, PayloadType>));

        // returned merged relations, only if nthreads > 1 
        mergedRelR->tuples = tmpoutR;
        mergedRelR->num_tuples = mergeRtotal;
        mergedRelS->tuples = tmpoutS;
        mergedRelS->num_tuples = mergeStotal;

        // determine the L3 cache-size per thread
        // int nnuma = get_num_numa_regions();

        // active number of threads in the current NUMA-region: 
        int active_nthreads_in_numa = get_num_active_threads_in_numa(numaregionid);

        // index of the current thread in its NUMA-region: 
        int numatidx = get_thread_index_in_numa(my_tid);

        // get the exclusive part of the merge buffer for the current thread 
        int bufsz_thr = (args->joincfg->MWAYMERGEBUFFERSIZE/active_nthreads_in_numa)
                / sizeof(Tuple<KeyType, PayloadType>);
        Tuple<KeyType, PayloadType> * mergebuf = args->sharedmergebuffer[numaregionid]
                                                     + (numatidx * bufsz_thr);


#ifdef SKEW_HANDLING
        // We decompose the merge tasks into number of merge tasks and
        //   add them to a task queue if the total merge size > expected
        //   average merge size + some threashold(currently 10%). Task
        //   queue is local to NUMA-region. TODO: we should also consider task
        //   stealing among NUMA-regions as an extension. We also
        //   consider skew only in S for the time being.
        //
        taskqueue_t * mergetaskqueue = args->numa_taskqueues[numaregionid];
        uint64_t DECOMPOSE_THRESHOLD = args->numS * SKEW_DECOMPOSE_MARGIN;
        int rv;
        if(mergeStotal > DECOMPOSE_THRESHOLD) 
        {
            // 1) Break-down the merge task and add to a task queue;
            //   first determine how many merge tasks we will create. 
            int nmergetasks = (mergeStotal + DECOMPOSE_THRESHOLD - 1) 
                              / DECOMPOSE_THRESHOLD;
            KeyType * merge_splitters = (KeyType *)
                                         malloc(sizeof(KeyType) * nmergetasks);
            Relation<KeyType, PayloadType> * heavyhits = 0;
            // 2) Find splitters for range part. to different merge tasks.
            balanced_sorted_partition(Sparts,PARTFANOUT,SKEW_DECOMPOSE_SAMPLES,
                                      merge_splitters, nmergetasks, my_tid, 
                                      &heavyhits);

            // 3) Create and add merge tasks to the merge task queue. Also,
               return the largest merge task without adding to the queue. 
            MergeTask<KeyType, PayloadType> * largest_task = 0;
            create_and_add_merge_tasks(mergetaskqueue,
                                       &largest_task, 
                                       merge_splitters,
                                       Sparts,
                                       tmpoutS,
                                       nmergetasks,
                                       PARTFANOUT);
            free(merge_splitters);

            // 5) Check if heavy hitters exist in the largest task, if yes, then
            //   find the boundaries of the heavy hitters and create specialized
            //   'direct copy' tasks for those.
            detect_heavy_hitters_and_add_tasks(heavyhits,
                                               mergetaskqueue,
                                               largest_task,
                                               my_tid,
                                               PARTFANOUT);
            free(heavyhits->tuples);
            free(heavyhits);

        }
        else 
#endif
        {
            // Just add the merge tasks for rel-S to the merge task queue. 
            MergeTask<KeyType, PayloadType> * mwtask = create_mergetask(PARTFANOUT);
            for(int p = 0; p < PARTFANOUT; p++){
                mwtask->runstomerge[p] = Sparts[p];
            }
            mwtask->totaltuples = mergeStotal;
            mwtask->output = tmpoutS;
            //DEBUGMSG(1, "[---- TID=%d adding merge task with %d tuples.\n",
            //         my_tid, mwtask->totaltuples);
            taskqueue_addfront(mergetaskqueue, mwtask);
        }
        // Wait until all threads complete skew handling logic. 
        BARRIER_ARRIVE(args->barrier, rv);

        // Now do the multi-way merging for the relation R. 
        if(scalarmergeflag){
            scalar_multiway_merge<KeyType, PayloadType>(tmpoutR, Rparts, PARTFANOUT, 
                                  mergebuf, bufsz_thr);
        }
        else {
            avx_multiway_merge<KeyType, PayloadType>(tmpoutR, Rparts, PARTFANOUT, 
                               mergebuf, bufsz_thr);
        }

#ifdef SKEW_HANDLING
        // 6) Grab merge tasks from the task queue and merge as usual or
        //    directly copy if it is a heavy-hitter. 
        MergeTask<KeyType, PayloadType> * mwtask = 0;
        while((mwtask = (MergeTask<KeyType, PayloadType>*)taskqueue_getfront(mergetaskqueue)))
        {
            //DEBUGMSG(1, " TID=%d GOT TASK = %d tuples; nruns=%d, is-HH=%d \n",
            //         my_tid, mwtask->totaltuples, 
            //         mwtask->numruns, mwtask->isheavyhitter);

            if(mwtask->totaltuples > 0)
            {
                if(mwtask->isheavyhitter){
                    // if heavy-hitter then not merged, directly copied 
                    do_fast_memcpy(mwtask);
                }
                else {
                    // Now do the multi-way merging
                    if(scalarmergeflag)
                    {
                        // do normal scalar multi-way merge
                        scalar_multiway_merge<KeyType, PayloadType>(mwtask->output, 
                                              mwtask->runstomerge, 
                                              mwtask->numruns, 
                                              mergebuf, bufsz_thr);
                    }
                    else {
                        avx_multiway_merge<KeyType, PayloadType>(mwtask->output, 
                                           mwtask->runstomerge, 
                                           mwtask->numruns, 
                                           mergebuf, bufsz_thr);
                    }
                }
            }
            // clean-up. 
            free_mergetask(mwtask);
        }
#endif

    }
*/      
  }

  void mergejoin_phase(Relation<KeyType, PayloadType> ** relRparts, Relation<KeyType, PayloadType> ** relSparts,
                       Relation<KeyType, PayloadType> * mergedRelR, Relation<KeyType, PayloadType> * mergedRelS, JoinThreadType * in_args) 
  {
    LearnedSortMergeMultiwayJoinThread<KeyType, PayloadType> * args = (LearnedSortMergeMultiwayJoinThread<KeyType, PayloadType> *) in_args;  
   
    void * chainedbuf = NULL;
    uint64_t nresults = 0;

    Tuple<KeyType, PayloadType> * rtuples = args->threadrelchunks[args->my_tid][0].R.tuples;
    Tuple<KeyType, PayloadType> * stuples = args->threadrelchunks[args->my_tid][0].S.tuples;


    printf("my_tid %d R num_elems %ld min %ld %ld max %ld %ld\n", args->my_tid, args->threadrelchunks[args->my_tid][0].R.num_tuples, 
        rtuples[0].key, rtuples[0].payload, rtuples[args->threadrelchunks[args->my_tid][0].R.num_tuples - 1].key, rtuples[args->threadrelchunks[args->my_tid][0].R.num_tuples - 1].payload);
    for(int i = 0; i < args->threadrelchunks[args->my_tid][0].R.num_tuples - 1; i++)
    {
      if(rtuples[i].key > rtuples[i+1].key)
       printf("my_tid %d unsorted R elems keys %ld %ld payload %ld %ld at %d %d\n", args->my_tid, rtuples[i].key, rtuples[i+1].key, rtuples[i].payload, rtuples[i+1].payload, i, i+1);
    }
    printf("my_tid %d S num_elems %ld min %ld %ld max %ld %ld\n", args->my_tid, args->threadrelchunks[args->my_tid][0].S.num_tuples, 
        stuples[0].key, stuples[0].payload, stuples[args->threadrelchunks[args->my_tid][0].S.num_tuples - 1].key, stuples[args->threadrelchunks[args->my_tid][0].S.num_tuples - 1].payload);
    for(int i = 0; i < args->threadrelchunks[args->my_tid][0].R.num_tuples - 1; i++)
    {
      if(rtuples[i].key > rtuples[i+1].key)
       printf("my_tid %d unsorted R elems keys %ld %ld payload %ld %ld at %d %d\n", args->my_tid, rtuples[i].key, rtuples[i+1].key, rtuples[i].payload, rtuples[i+1].payload, i, i+1);
    }

    nresults = merge_join(rtuples, stuples, 
                args->threadrelchunks[args->my_tid][0].R.num_tuples, 
                args->threadrelchunks[args->my_tid][0].S.num_tuples, chainedbuf);
    
    args->result = nresults;
    // printf("TID=%d --> #res = %d %d\n", my_tid, args->result, nresults); 
    
    /*
    if(args->nthreads > 1){
        Tuple<KeyType, PayloadType> * rtuples = (Tuple<KeyType, PayloadType> *) mergedRelR->tuples;
        Tuple<KeyType, PayloadType> * stuples = (Tuple<KeyType, PayloadType> *) mergedRelS->tuples;

        nresults = merge_join(rtuples, stuples,
                                   mergedRelR->num_tuples, mergedRelS->num_tuples, chainedbuf);

    } else {
        // single-threaded execution: just join sorted partition-pairs
        for(int i = 0; i < PARTFANOUT; i ++) {
            // evaluate join between each sorted part 
            nresults += merge_join(relRparts[i]->tuples, relSparts[i]->tuples,
                    relRparts[i]->num_tuples, relSparts[i]->num_tuples,
                    chainedbuf);
        }

    }
    */
  }

  void partitioning_cleanup(Relation<KeyType, PayloadType> ** relRparts, Relation<KeyType, PayloadType> ** relSparts) 
  {
    free(relRparts[0]);
    free(relRparts);
    free(relSparts);
  }

  ~LearnedSortMergeMultiwayJoinSteps() {
  }
 
    /**
     * Balanced partition of n sorted lists using equi-depth histograms
     * and combined CDF over them. Returns num X intkey_t, 
     * which is a list of partitioning keys for ranges of each thread.
     */
    void
    balanced_sorted_partition(Relation<KeyType, PayloadType> ** rels, int nrels, int eqhistsamples,
                              KeyType * partkeys /* [out] */, int nparts, int my_tid,
                              Relation<KeyType, PayloadType> ** heavyhitters)
    {
        /* compute equi-depth histogram for each chunk with p=eqhistsamples values */
        Relation<KeyType, PayloadType> eqhist[nrels];
        Tuple<KeyType, PayloadType> tmp[nrels * eqhistsamples];
        int64_t total = 0;
        for(int i = 0; i < nrels; i++){
            eqhist[i].num_tuples = eqhistsamples;
            eqhist[i].tuples     = (Tuple<KeyType, PayloadType>*)(tmp + i * eqhistsamples);
        }

        /* fill in the equi-hist values from each chunk */
        for(int i = 0; i < nrels; i++){
            int eqwidth = rels[i]->num_tuples / eqhistsamples;
            total += rels[i]->num_tuples;
            for(int p = eqhistsamples-1; p >= 0; p--){
                uint64_t  idx               = (p+1) * eqwidth - 1;
                Tuple<KeyType, PayloadType> * tup               = &(rels[i]->tuples[idx]);
                eqhist[i].tuples[p].key     = tup->key;
                eqhist[i].tuples[p].payload = eqwidth;
            }
        }

        /* find heavy hitters */
        Relation<KeyType, PayloadType> * heavyhits = (Relation<KeyType, PayloadType> *) malloc(sizeof(Relation<KeyType, PayloadType>));
        heavyhits->tuples = (Tuple<KeyType, PayloadType>*)malloc(sizeof(Tuple<KeyType, PayloadType>)*SKEW_MAX_HEAVY_HITTERS);
        heavyhits->num_tuples = 0;

        for(int i = 0; i < nrels; i++){
            find_heavy_hitters(&eqhist[i], SKEW_HEAVY_HITTER_THR, heavyhits);
        }
        std::sort((int64_t*)(heavyhits->tuples), (int64_t*)(heavyhits->tuples + heavyhits->num_tuples));
        // for(int h = 0; h < heavyhits->num_tuples; h++){
        //     printf("===> TID = %d, HEAVY-HITTER-%d = %d\n", my_tid, h, heavyhits->tuples[h].key);
        // }
                                    

        /* merge the histograms, new size can be up to (eqhistsamples * nrels) */
        std::map<KeyType, int> cdfhist;
        for(int i = 0; i < nrels; i++){
            //if(my_tid == 3) printf("****************** HIST for %d\n", i);
            for(int p = eqhistsamples-1; p >= 0; p--){
                KeyType key          = eqhist[i].tuples[p].key;
                int      count        = eqhist[i].tuples[p].payload;
                if(cdfhist.find(key) != cdfhist.end()){
                    cdfhist[key]     += count;
                }
                else
                    cdfhist[key] = count;
                //if(my_tid == 3) printf("KEY=%d ; COUNT=%d\n", key,count);
            }
        }

        /* compute a prefix-sum/cdf over the merged histograms */
        Relation<KeyType, PayloadType> prefixcdf;
        prefixcdf.num_tuples = cdfhist.size();
        prefixcdf.tuples     = (Tuple<KeyType, PayloadType> *) malloc(sizeof(Tuple<KeyType, PayloadType>) * cdfhist.size());

        int sum = 0;
        int i = 0;
        for (/*std::map<KeyType,int>::iterator*/ auto it = cdfhist.begin(); it != cdfhist.end(); ++it) 
        {
            sum                         += it->second;
            prefixcdf.tuples[i].key      = it->first;
            prefixcdf.tuples[i].payload  = sum;
            i++;
        }

        // char fn[512];
        // sprintf(fn, "cdf-%d.tbl", my_tid);
        // write_relation(&prefixcdf, fn);
        
        // std::cout << "#KEY; COUNT\n";
        // for (i = 0; i < prefixcdf.num_tuples; i++)
        //     std::cout << prefixcdf.tuples[i].key << " "
        //               << prefixcdf.tuples[i].payload << std::endl;

        /* find partition boundaries */
        unsigned int j = 0;
        int64_t numperpart = total / nparts;

        for(i = 0; i < nparts-1; i++) {
            PayloadType count = (i+1) * numperpart;
            for(; j < prefixcdf.num_tuples; j++) {
                if(prefixcdf.tuples[j].payload == count){
                   partkeys[i] = prefixcdf.tuples[j].key;
                    // fprintf(stderr, "TASK-%d partition-key = %d\n",
                    //         i, partkeys[i]);
                    break;
                }
                else if(prefixcdf.tuples[j].payload > count){
                    int cnt2 = prefixcdf.tuples[j].payload;
                    int key2 = prefixcdf.tuples[j].key;
                    int cnt1 = 0;
                    int key1 = 0;

                    if( j > 0 ) {
                        cnt1 = prefixcdf.tuples[j-1].payload;
                        key1 = prefixcdf.tuples[j-1].key;
                    }

                    double factor = (count - cnt1) / (double)(cnt2 - cnt1);
                    /* find the partitioning key with interpolation */
                    KeyType partkey = (key2 - key1) * factor + key1;
                    partkeys[i] = partkey;
                    // fprintf(stderr, "TASK-%d partition-key = %d\n",
                    //         i, partkey);
                    break;
                }
            }
        }

        free(prefixcdf.tuples);

        /* return heavyhitters */
        *heavyhitters = heavyhits;

    }

    /**
     * Determines and returns at most SKEW_MAX_HEAVY_HITTERS heavy hitter keys
     * from an equi-depth histogram. Heavy hitter is a key that occurs
     * above given threshold (SKEW_HEAVY_HITTER_THR) percent of the buckets in the
     * histogram. 
     */
    void
    find_heavy_hitters(Relation<KeyType, PayloadType> * equidepthhist, double thrpercent, 
                       Relation<KeyType, PayloadType> * heavyhits)
    {
        if(heavyhits->num_tuples > SKEW_MAX_HEAVY_HITTERS)
            return;
        
        /* number of threshold buckets that a value must occur to be a HH */
        int64_t threshold = equidepthhist->num_tuples * thrpercent;
        int64_t count = 0;
        KeyType lastkey = -1;
        
        for(unsigned int i = 0; i < equidepthhist->num_tuples; i++){
            KeyType key = equidepthhist->tuples[i].key;
            
            if(key == lastkey){
                count++;
            }
            else {
                if(lastkey != -1){
                    if(count > threshold){
                        /* search in HH first */
                        int notfound = 1;
                        for(unsigned int h = 0; h < heavyhits->num_tuples; h++){
                            if(heavyhits->tuples[h].key == lastkey){
                                notfound = 0;
                                break;
                            }
                        }

                        if(notfound){
                            heavyhits->tuples[heavyhits->num_tuples].key = lastkey;
                            heavyhits->num_tuples ++;
                            if(heavyhits->num_tuples > SKEW_MAX_HEAVY_HITTERS)
                                break;
                        }
                    }
                }
                lastkey = key;
                count = 1;
            }
        }
    }

     void create_and_add_merge_tasks(taskqueue_t * mergequeue,
                                 MergeTask<KeyType, PayloadType> ** largest_ret,
                                 KeyType * merge_splitters,
                                 Relation<KeyType, PayloadType> ** relSparts,
                                 Tuple<KeyType, PayloadType> * outputptr,
                                 const int ntasks,
                                 const int PARTFANOUT)
    {
        /** keep track of the largest mway task */
        MergeTask<KeyType, PayloadType> * largest = 0;
        //int my_tid = -1;
        int64_t startoffset[PARTFANOUT];
        uint64_t runsizes[PARTFANOUT];
        uint64_t task0tot = 0;
        {
            largest = create_mergetask(PARTFANOUT);
            /* t = 0 --> first task is directly assigned */
            KeyType tk0 = merge_splitters[0];
            for(int p = 0; p < PARTFANOUT; p++){
                if (relSparts[p]->num_tuples > 0 )
                {
                    int64_t stoff = binsearch_lower_bound(relSparts[p]->tuples,
                                                          relSparts[p]->num_tuples,
                                                          tk0);
                    startoffset[p] = stoff;
                    runsizes[p] = relSparts[p]->num_tuples - stoff;
                    task0tot += stoff;

                    largest->runstomerge[p]->tuples = relSparts[p]->tuples;
                    largest->runstomerge[p]->num_tuples = stoff;
                }
                else 
                {
                    startoffset[p] = 0;
                    runsizes[p] = 0;
                    largest->runstomerge[p]->tuples = 0;
                    largest->runstomerge[p]->num_tuples = 0;
                }
            }
            largest->totaltuples = task0tot;
            largest->output = outputptr;
            outputptr += task0tot;
        }
        // printf("[INFO ] TID=%d merge-task-0 total = %d\n", my_tid, task0tot);
        for(int t = 1; t < ntasks-1; t++) {
            KeyType tk = merge_splitters[t];
            MergeTask<KeyType, PayloadType> * mwtask = create_mergetask(PARTFANOUT);
            uint64_t tasktot = 0;
            for(int p = 0; p < PARTFANOUT; p++){
                if ( runsizes[p] > 0 )
                {
                    int64_t stoff = binsearch_lower_bound(relSparts[p]->tuples 
                                                        + startoffset[p],
                                                        runsizes[p],
                                                        tk);
                    mwtask->runstomerge[p]->tuples = relSparts[p]->tuples + startoffset[p];
                    mwtask->runstomerge[p]->num_tuples = stoff;
                    startoffset[p] += stoff;
                    runsizes[p] -= stoff;
                    tasktot += mwtask->runstomerge[p]->num_tuples;
                }
                else 
                {
                    mwtask->runstomerge[p]->tuples = 0;
                    mwtask->runstomerge[p]->num_tuples = 0;
                }
            }
            mwtask->totaltuples = tasktot;
            mwtask->output = outputptr;
            outputptr += tasktot;

            // if(mwtask->totaltuples > DECOMPOSE_THR){
            //     int onefourth = mwtask->runstomerge[0]->num_tuples/4;
            //     int threefourth = 3*onefourth;
            //     printf("[IMPORTANT:TID=%d] HeavyHitter => NUM=%d; at[%d]=%d ; at[%d]=%d\n",
            //            my_tid, mwtask->runstomerge[0]->num_tuples,
            //            onefourth, mwtask->runstomerge[0]->tuples[onefourth].key,
            //            threefourth, mwtask->runstomerge[0]->tuples[threefourth].key);
            // }
                    
            //printf("[INFO ] TID=%d creating merge task with %d total tuples.\n", my_tid, tasktot);
            if(largest == 0){
                if(mwtask->totaltuples > task0tot){
                    largest = mwtask;
                }
                else {
                    //DEBUGMSG(1, "[---- TID=%d adding merge task with %d tuples.\n",
                        //my_tid, mwtask->totaltuples);
                    taskqueue_addtail(mergequeue, mwtask);
                }
            }
            else if(mwtask->totaltuples > largest->totaltuples){
                //DEBUGMSG(1, "[---- TID=%d adding merge task with %d tuples.\n", 
                    //my_tid, largest->totaltuples);
                taskqueue_addtail(mergequeue, largest);
                largest = mwtask;
            }
            else {
                //DEBUGMSG(1, "[---- TID=%d adding merge task with %d tuples.\n", 
                    //my_tid, mwtask->totaltuples);
                taskqueue_addtail(mergequeue, mwtask);
            }
        }
        {
            /* last task until the end of the run */
            MergeTask<KeyType, PayloadType> * mwtask = create_mergetask(PARTFANOUT);
            uint64_t tasktot = 0;
            for(int p = 0; p < PARTFANOUT; p++){
                mwtask->runstomerge[p]->tuples = relSparts[p]->tuples + startoffset[p];
                mwtask->runstomerge[p]->num_tuples = runsizes[p];
                tasktot += mwtask->runstomerge[p]->num_tuples;
                // if(my_tid==3){
                //     printf("TID=3 --> key0=%d key-last=%d\n",
                //            mwtask->runstomerge[p]->tuples[0].key,
                //            mwtask->runstomerge[p]->tuples[mwtask->runstomerge[p]->num_tuples-1].key);
                // }
            }
            mwtask->totaltuples = tasktot;
            mwtask->output = outputptr;
            //printf("[INFO ] TID=%d creating merge task with %d total tuples to merge.\n", my_tid, tasktot);
            if(largest == 0){
                if(mwtask->totaltuples > task0tot){
                    largest = mwtask;
                }
                else {
                    //DEBUGMSG(1, "[---- TID=%d adding merge task with %d tuples.\n", 
                    //    my_tid, mwtask->totaltuples);
                    taskqueue_addtail(mergequeue, mwtask);
                }
            }
            else if(mwtask->totaltuples > largest->totaltuples){
                //DEBUGMSG(1, "[---- TID=%d adding merge task with %d tuples.\n", 
                //    my_tid, largest->totaltuples);
                taskqueue_addtail(mergequeue, largest);
                largest = mwtask;
            }
            else {
                //DEBUGMSG(1, "[---- TID=%d adding merge task with %d tuples.\n", 
                //    my_tid, mwtask->totaltuples);
                taskqueue_addtail(mergequeue, mwtask);
            }
        }

        /* return the largest merge tasks without adding to the task queue. */
        *largest_ret = largest;
    }

    void
    free_mergetask(MergeTask<KeyType, PayloadType> * mt)
    {
        free(mt);
    }

    /**
     * Find the first offset in sortedrun which is >= skey.
     */
    int64_t
    binsearch_lower_bound(Tuple<KeyType, PayloadType> * sortedrun, int64_t ntups, KeyType skey)
    {
        const Tuple<KeyType, PayloadType> * tups = sortedrun;
        int64_t lo, hi, mid;
        lo = 0;
        hi = ntups - 1;
        
        while(lo < hi){
            mid = (lo + hi)/2;

            if(tups[mid].key >= skey)
                hi = mid;
            else
                lo = mid + 1;
        }

        if(tups[lo].key >= skey)
            return lo;
        else {
            //DEBUGMSG(1, "1. TODO FIXME : can be a problem !\n");
            return lo+1;//return -1;/* not found */
        }
    }

    /**
     * Find the maximum offset in sortedrun which is < skey.
     */
    int64_t
    binsearch_upper_bound(Tuple<KeyType, PayloadType> * sortedrun, int64_t ntup,  KeyType skey)
    {
        const Tuple<KeyType, PayloadType> * tups = sortedrun;
        int64_t lo, hi, mid;
        lo = 0;
        hi = ntup - 1;
        
        while(lo < hi){
            mid = (lo + hi + 1)/2;

            if(tups[mid].key < skey)
                lo = mid;
            else
                hi = mid - 1;
        }

        if(tups[lo].key < skey)
            return lo;
        else{
            //DEBUGMSG(1, "2. TODO FIXME : can be a problem !\n");
            return 0; //-1;/* not found */
        }
    }


    MergeTask<KeyType, PayloadType> *
    create_mergetask(int nruns)
    {
        size_t alloc_sz = sizeof(MergeTask<KeyType, PayloadType>) + nruns *
                    (sizeof(Relation<KeyType, PayloadType> *) + sizeof(Relation<KeyType, PayloadType>));
        char * mem = (char *) malloc(alloc_sz);

        MergeTask<KeyType, PayloadType> * mt = (MergeTask<KeyType, PayloadType> *)(mem);
        mem += sizeof(MergeTask<KeyType, PayloadType>);
        mt->runstomerge = (Relation<KeyType, PayloadType> **)(mem);
        mem += nruns * sizeof(Relation<KeyType, PayloadType> *);
        for(int i = 0; i < nruns; i++){
            mt->runstomerge[i] = (Relation<KeyType, PayloadType> *) (mem);
            mem += sizeof(Relation<KeyType, PayloadType>);
        }
        mt->numruns = nruns;
        mt->isheavyhitter = 0;

        return mt;
    }

    void
    detect_heavy_hitters_and_add_tasks(Relation<KeyType, PayloadType> * heavyhits,
                                    taskqueue_t * mergequeue,
                                    MergeTask<KeyType, PayloadType> * largest_task,
                                    int my_tid,
                                    const int PARTFANOUT)
    {
        Relation<KeyType, PayloadType> ** relparts = largest_task->runstomerge;
        /* check if any heavyhits exist in relparts[0...F] then
        separately handle those */
        int heavyhitexist = 0;

        for(unsigned int h = 0; h < heavyhits->num_tuples; h++)
        {
            for(int p = 0; p < PARTFANOUT; p++)
            {
                if( relparts[p]->num_tuples > 0 ) 
                {
                    if(relparts[p]->tuples[0].key <= heavyhits->tuples[h].key &&
                    relparts[p]->tuples[relparts[p]->num_tuples-1].key 
                    >= heavyhits->tuples[h].key)
                    {
                        heavyhitexist = heavyhits->tuples[h].key;
                        break;
                    }
                }
            }

            if( heavyhitexist )
                break;
        }

        if( heavyhitexist )
        {
            /* try to identify heavyhit bounds and directly copy */
            //DEBUGMSG(1, "$$$ HEAVYHIT FOUND = %d\n", heavyhitexist);
            MergeTask<KeyType, PayloadType> * mwtask = create_mergetask(PARTFANOUT);
            MergeTask<KeyType, PayloadType> * mwtask2 = create_mergetask(PARTFANOUT);
            int64_t tasktot = 0;
            mwtask2->totaltuples = 0;
            for(int p = 0; p < PARTFANOUT; p++)
            {
                if ( relparts[p]->num_tuples > 0 ) 
                {
                    /* if(relparts[p]->tuples[0].key <= heavyhitexist && */
                    /*    relparts[p]->tuples[relparts[p]->num_tuples-1].key  */
                    /*    >= heavyhitexist) */
                    /* { */
                        uint64_t stoff = binsearch_lower_bound(relparts[p]->tuples,
                                                            relparts[p]->num_tuples,
                                                            heavyhitexist);
                        uint64_t endoff = binsearch_upper_bound(relparts[p]->tuples,
                                                                relparts[p]->num_tuples,
                                                                heavyhitexist+1);

                        if(endoff > stoff)
                        {
                            // printf("--- HH@[%d]=%d --- HH@[%d]=%d\n",
                            //        stoff, relparts[p]->tuples[stoff].key,
                            //        endoff,
                            //       relparts[p]->tuples[endoff].key);
                            mwtask->runstomerge[p]->tuples = relparts[p]->tuples + stoff;
                            mwtask->runstomerge[p]->num_tuples = endoff - stoff + 1;
                            tasktot += mwtask->runstomerge[p]->num_tuples;
                            if(endoff < (relparts[p]->num_tuples-1)){
                                mwtask2->runstomerge[p]->tuples = relparts[p]->tuples + endoff + 1;
                                mwtask2->runstomerge[p]->num_tuples = relparts[p]->num_tuples - endoff -1;
                                mwtask2->totaltuples += mwtask2->runstomerge[p]->num_tuples;
                            }
                            else {
                                mwtask2->runstomerge[p]->tuples = 0;
                                mwtask2->runstomerge[p]->num_tuples = 0;
                            }
                            /* Update the original task. */
                            relparts[p]->num_tuples = stoff;
                        }
                        /* else if (stoff == relparts[p]->num_tuples)  */
                        /* { */
                        /*     mwtask->runstomerge[p]->tuples = relparts[p]->tuples; */
                        /*     mwtask->runstomerge[p]->num_tuples = stoff; */
                        /*     tasktot += mwtask->runstomerge[p]->num_tuples; */
                        /*     mwtask2->runstomerge[p]->tuples = 0; */
                        /*     mwtask2->runstomerge[p]->num_tuples = 0; */
                        /*     /\* Update the original task. *\/ */
                        /*     relparts[p]->num_tuples = 0; */
                        /* } */
                        else if (endoff == 0) 
                        {
                            mwtask->runstomerge[p]->tuples = 0;
                            mwtask->runstomerge[p]->num_tuples = 0;
                            mwtask2->runstomerge[p]->tuples = relparts[p]->tuples;
                            mwtask2->runstomerge[p]->num_tuples = relparts[p]->num_tuples;
                            mwtask2->totaltuples += mwtask2->runstomerge[p]->num_tuples;
                            /* Update the original task. */
                            relparts[p]->num_tuples = 0;
                        }
                        else {
                            mwtask->runstomerge[p]->tuples = 0;
                            mwtask->runstomerge[p]->num_tuples = 0;
                            mwtask2->runstomerge[p]->tuples = 0;
                            mwtask2->runstomerge[p]->num_tuples = 0;
                        }
                    /* } */
                    /* else { */
                    /*     mwtask->runstomerge[p]->tuples = 0; */
                    /*     mwtask->runstomerge[p]->num_tuples = 0; */
                    /*     mwtask2->runstomerge[p]->tuples = 0; */
                    /*     mwtask2->runstomerge[p]->num_tuples = 0; */
                    /* } */
                }
                else 
                {
                    mwtask->runstomerge[p]->tuples = 0;
                    mwtask->runstomerge[p]->num_tuples = 0;
                    mwtask2->runstomerge[p]->tuples = 0;
                    mwtask2->runstomerge[p]->num_tuples = 0;
                }
            }
            mwtask->totaltuples = tasktot;

            /* Update the original task. */
            largest_task->totaltuples = largest_task->totaltuples 
                                        - mwtask->totaltuples 
                                        - mwtask2->totaltuples;

            /* Update starting output offsets. */
            mwtask->output = largest_task->output + largest_task->totaltuples;
            mwtask->isheavyhitter = 1;
            mwtask2->output = mwtask->output + tasktot;

            /*DEBUGMSG(1, "[---- TID=%d adding merge HH-copy task with %d tuples.\n", 
                    my_tid, mwtask->totaltuples);
            DEBUGMSG(1, "[---- TID=%d adding merge task with %d tuples.\n", 
                    my_tid, mwtask2->totaltuples);
            DEBUGMSG(1, "[---- TID=%d adding merge task with %d tuples.\n",
                    my_tid, largest_task->totaltuples);
            */
            /* Add the largest one of these subtasks to the front of the queue. */
            if ( largest_task->totaltuples > mwtask->totaltuples )
            {
                if ( largest_task->totaltuples > mwtask2->totaltuples )
                {
                    taskqueue_addfront(mergequeue, largest_task);
                    taskqueue_addtail(mergequeue, mwtask);
                    taskqueue_addtail(mergequeue, mwtask2);
                }
                else {
                    taskqueue_addfront(mergequeue, mwtask2);
                    taskqueue_addtail(mergequeue, largest_task);
                    taskqueue_addtail(mergequeue, mwtask);
                }
            }
            else 
            {
                if ( mwtask->totaltuples > mwtask2->totaltuples )
                {
                    taskqueue_addfront(mergequeue, mwtask);
                    taskqueue_addtail(mergequeue, largest_task);
                    taskqueue_addtail(mergequeue, mwtask2);
                }
                else {
                    taskqueue_addfront(mergequeue, mwtask2);
                    taskqueue_addtail(mergequeue, mwtask);
                    taskqueue_addtail(mergequeue, largest_task);
                }
            }

        }
        else 
        {
            //DEBUGMSG(1, "$$$ NO HEAVYHIT FOUND; JUST ADDING TO FRONT OF QUEUE.\n");
            taskqueue_addfront(mergequeue, largest_task);
        }
    }

    /**
     * This method copies given relation tuples to the output array one by
     * one. Used for copying out heavy hitters. There is an optimization
     * opportunity to use AVX/SSE copy.
     */
    void
    do_fast_memcpy(MergeTask<KeyType, PayloadType> * mwtask)
    {
        Tuple<KeyType, PayloadType> * out = mwtask->output;
        for(int i = 0; i < mwtask->numruns; i++){
            int64_t numtuples = mwtask->runstomerge[i]->num_tuples;
            nt_memcpy(out, mwtask->runstomerge[i]->tuples, numtuples * sizeof(Tuple<KeyType, PayloadType>));
            out += numtuples;
        }
    }

    /**
     * Efficient memcpy operation using non-temporal load/stores.
     */
    void
    nt_memcpy(void * dst, void * src, size_t sz)
    {
    // #ifdef __AVX__
    // #define ALIGN_BOUNDARY 64
    //     char * dstptr, * srcptr;
    //     char * endptr = (char *)src + sz;
    //     for(dstptr = (char *)dst, srcptr = (char *)src; srcptr != endptr; ){
    //         if(((uintptr_t)srcptr % ALIGN_BOUNDARY) == 0 &&
    //            ((uintptr_t)dstptr % ALIGN_BOUNDARY) == 0)
    //             break;

    //         *dstptr ++ = *srcptr ++;
    //     }

    //     __m256i * src256 = (__m256i *) srcptr;
    //     __m256i * dst256 = (__m256i *) dstptr;
    //     __m256i * end256 = (__m256i *) ((uintptr_t)endptr & ~(ALIGN_BOUNDARY*2-1));
    //     for( ; src256 != end256; src256 += 2, dst256 += 2) {
    //      __asm volatile(
    //                   "vmovntdqa 0(%0) , %%xmm0;  "
    //                   "vmovntdqa 32(%0), %%xmm1;  " 
    //                   "vmovntdq %%xmm0, 0(%1) ;  "
    //                   "vmovntdq %%xmm1, 32(%1) ;  "
    //                   ::"r"(src256), "r"(dst256));
    //     }

    //     /* handle remainders */
    //     for(dstptr = (char *)dst256, srcptr = (char*)src256; srcptr != endptr; ){
    //         *dstptr ++ = *srcptr ++;
    //     }


    // #elif defined(__SSE4_2__)
    // #define ALIGN_BOUNDARY 16
    //     /* TODO: No AVX but SSE 4.2 or less? */
    // #else
        /* fallback memcpy */
        memcpy(dst, src, sz);
    // #endif
    }

    /**
     * Does merge join on two sorted relations. Just a naive scalar
     * implementation. TODO: consider AVX for this code.
     *
     * @param rtuples sorted relation R
     * @param stuples sorted relation S
     * @param numR number of tuples in R
     * @param numS number of tuples in S
     * @param output join results, if JOIN_MATERIALE defined.
     */
    uint64_t
    merge_join(Tuple<KeyType, PayloadType> * rtuples, Tuple<KeyType, PayloadType> * stuples,
            const uint64_t numR, const uint64_t numS, void * output)
    {
        uint64_t i = 0, j = 0;
        uint64_t matches = 0;

        while( i < numR && j < numS ) {
            if( rtuples[i].key < stuples[j].key )
                i ++;
            else if(rtuples[i].key > stuples[j].key)
                j ++;
            else {
                /* rtuples[i].key is equal to stuples[j].key */
                uint64_t jj;
                do {
                    jj = j;

                    do {
                        matches ++;
                        jj++;

                    } while(jj < numS && rtuples[i].key == stuples[jj].key);

                    i++;

                } while(i < numR && rtuples[i].key == stuples[j].key);

                j = jj;
            }
        }
        /* printf("More-S = %d, More-R = %d, remS=%d\n", (j<numS), (i<numR), numS-j); */
        /* if(rtuples[numR-1].key == stuples[j].key) */
        /*     printf("lastS equal lastR = %d\n", 1); */
        /* matches = merge_join8((int64_t *)rtuples, (int64_t*)stuples, 0, numR); */

        return matches;
    }
};