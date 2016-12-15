#ifndef QPHIX_BARRIER_H_
#define QPHIX_BARRIER_H_

#include <immintrin.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace std;

#ifndef __MIC__
#define _mm_delay_64(x)
#endif

class Barrier {
	typedef struct {
		volatile int counter;
	} atomic_t;

#define atomic_set1(v,i) (((v)->counter) = (i))

	static inline int atomic_dec_and_test( atomic_t *v )
	{
		return !(__sync_sub_and_fetch(&v->counter, 1));
	}


	typedef struct {
		int volatile sense;
	} htbar_t;


	typedef struct {
		int volatile htbar_sense;
		htbar_t htbar[4];
		int volatile *myflags[2];
		int **partnerflags[2];
		int parity;
		int sense;
	} kingbar_t;

	int num_kingbar_threads, num_kingkar_nodes, lg_num_kingbar_nodes, num_threads_per_core;
	atomic_t kingbar_threads_not_alloced;
	volatile int kingbar_alloc_done;
	kingbar_t **bar;
	
	public:
	Barrier(int num_cores, int threads_per_core)
	{
		int num_threads = num_cores * threads_per_core;
		num_threads_per_core = threads_per_core;
		int num_nodes = num_cores;
		bar = (kingbar_t **)_mm_malloc(num_nodes * sizeof(kingbar_t *), 64);

		num_kingbar_threads = num_threads;
		lg_num_kingbar_nodes = (int)ceil(log2(num_nodes));
		atomic_set1(&kingbar_threads_not_alloced, num_kingbar_threads);
		kingbar_alloc_done = 0;
	}

	void init(int tid)
	{
		kingbar_t *node;
		int j, m, n, nid = tid / num_threads_per_core, num_nodes = num_kingbar_threads / num_threads_per_core;

    // only the first of the 4 hardware threads per core sets up a node
		if ((tid % num_threads_per_core) == 0) {
			node = (kingbar_t *)_mm_malloc(sizeof(kingbar_t), 64);

			node->htbar_sense = 1;
			for(int i = 0; i < num_threads_per_core; i++) node->htbar[i].sense = 1;
			if(lg_num_kingbar_nodes > 0) {
				for (int i = 0;  i < 2;  i++) {
					node->myflags[i] = (int *)_mm_malloc(lg_num_kingbar_nodes * sizeof(int), 64);
					node->partnerflags[i] = (int **)MALLOC(lg_num_kingbar_nodes * sizeof(int *), 64);
				}
			}
			node->parity = 0;
			node->sense = 1;

			bar[nid] = node;
		}

    // barrier to let all the nodes get installed into bar
		if (atomic_dec_and_test(&kingbar_threads_not_alloced)) {
			atomic_set1(&kingbar_threads_not_alloced, num_kingbar_threads);
			kingbar_alloc_done = 1;
		} else {
			while(!kingbar_alloc_done);
		}

    // only the first of the 4 hardware threads per core completes setup
		if ((tid % num_threads_per_core) == 0) {
			for (int i = 0;  i < lg_num_kingbar_nodes;  i++) {
				node->myflags[0][i] = node->myflags[1][i] = 0;
				j = (nid + (1 << i)) % num_nodes;
				node->partnerflags[0][i] = (int *)&bar[j]->myflags[0][i];
				node->partnerflags[1][i] = (int *)&bar[j]->myflags[1][i];
			}
		}

    // barrier to let setup finish
		if (atomic_dec_and_test(&kingbar_threads_not_alloced)) {
			atomic_set1(&kingbar_threads_not_alloced, num_kingbar_threads);
			kingbar_alloc_done = 0;
		} else {
			while(kingbar_alloc_done);
		}
		//printf("%s:%d: reached %d\n", __func__, __LINE__, tid);
	}



	~Barrier()
	{
	  // Wot? No cleanup
	}

	void wait(int tid)
	{
		int i, nid = tid / num_threads_per_core, tofs = tid % num_threads_per_core, num_nodes = num_kingbar_threads / num_threads_per_core;
		kingbar_t *node = bar[nid];

    // first sync using a kbar inside the core
		node->htbar[tofs].sense = !node->htbar[tofs].sense;
		if (tofs == 0) {
			for (i = 1;  i < num_threads_per_core;  i++) {
				while (node->htbar[i].sense == node->htbar_sense) _mm_delay_64(100);
			}
		}
		else {
			while (node->htbar_sense != node->htbar[tofs].sense)_mm_delay_64(100);
		}

    // now, the first of the 4 hardware threads per core syncs with the others
		if (tofs == 0) {
			for (i = 0;  i < lg_num_kingbar_nodes;  i++) {
            //printf("thread %d setting %p, waiting on %p (sense %d)\n", tid, node->partnerflags[node->parity][i], &node->myflags[node->parity][i], node->sense);
				*node->partnerflags[node->parity][i] = node->sense;
				while (node->myflags[node->parity][i] != node->sense)_mm_delay_64(100);
			}
			if (node->parity == 1)
				node->sense = !node->sense;
			node->parity = 1 - node->parity;

        // wake up the waiting threads in this core
			node->htbar_sense = node->htbar[tofs].sense;
		}
	}
};

class KBar {
  typedef struct {
    int volatile sense;
  } htbar_t;

  int nThreads;
    int volatile htbar_sense;
    htbar_t htbar[4];

  public:
    KBar(int num_threads)
    {
      nThreads = num_threads;
      htbar_sense = 1;
      for(int i = 0; i < nThreads; i++) htbar[i].sense = 1;
    }

    ~KBar()
    {
    // Wot? No cleanup
    }

    void wait(int tofs)
    {
      int i;

      htbar[tofs].sense = !htbar[tofs].sense;
      if (tofs == 0) {
        for (i = 1;  i < nThreads;  i++) {
          while (htbar[i].sense == htbar_sense) _mm_delay_64(100);
        }
        htbar_sense = htbar[tofs].sense;
      }
      else {
        while (htbar_sense != htbar[tofs].sense)_mm_delay_64(100);
      }
    }
};

#endif
