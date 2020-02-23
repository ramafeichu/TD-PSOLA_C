#ifndef ARRAYS_PSOLA_H
#define	ARRAYS_PSOLA_H

#include <stdint.h>

typedef struct
{
    int16_t* buffer;
    
    int16_t* p2future;
    int16_t* p2present;
    int16_t* p2past;
    
    uint16_t num_samp_future;
    uint16_t num_samp_present;
    uint16_t num_samp_past;
   
}rt_buffer;

void init_rt_buffer(rt_buffer* rt_buff, uint16_t window_len);
void clear_rt_buffer(rt_buffer* rt_buff);
void push_rt_buffer(rt_buffer* rt_buff, int16_t* p2data, uint16_t number_samples);
void time_swap(rt_buffer* rt_buff);

#endif	/* ARRAYS_PSOLA_H */