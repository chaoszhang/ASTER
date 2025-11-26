#include<stdlib.h>

#define BATCH 32
#define BATCH_TYPE unsigned
#define MAX_TAPE_SIZE 4096

#define THREAD_POOL
#ifdef THREAD_POOL
#define TAPE_SIZE MAX_TAPE_SIZE
#else
#define TAPE_SIZE 1
#endif

struct GroupBatch{
	double p[4][BATCH] = {};
	double w[BATCH] = {};
	int start, end;
};

typedef int GroupIBatch[2];
typedef double GroupFBatch[5][BATCH];
typedef short SiteBatch[4][3][BATCH];
typedef BATCH_TYPE HotBatch[4];

void myCudaMalloc(void **p, int size);
void myCudaMemcpyH2D(void *T, void *S, int size);
void myCudaMemcpyD2H(void *T, void *S, int size);
void myCudaFree(void *p);
void cudaScore(GroupIBatch *GI, GroupFBatch *GF, SiteBatch* S, double* BS, int kc, int b);
void cudaUpdate(GroupIBatch *GI, SiteBatch *S, HotBatch *H, int n, int x, int y, int i, int kc, int b);

void myCudaMemcpyH2C(void *S, int size);
void cudaWork(GroupIBatch *GI, GroupFBatch *GF, SiteBatch* S, double* BS, HotBatch *H, int n, int kc, int b, int tapeSize);

//void cudaWork(GroupIBatch *GI, GroupFBatch *GF, SiteBatch* S, double* BS, HotBatch *H, int n, int kc, int b, int tapeSize, short *TAPE);