#ifndef CHACAPK_BATCH_SIZE_H_INCLUDED
#define CHACAPK_BATCH_SIZE_H_INCLUDED
/* Opt-in build setting shared by the C producer and native CUDA consumer.
 * Ordinary CPU builds retain their scalar route and the historical bound. */
#ifndef RAD_HACAPK_BATCH_SIZE
#define RAD_HACAPK_BATCH_SIZE 256
#endif
#if RAD_HACAPK_BATCH_SIZE < 1 || RAD_HACAPK_BATCH_SIZE > 2048
#error Unsupported HACApK entry batch size
#endif
#endif
