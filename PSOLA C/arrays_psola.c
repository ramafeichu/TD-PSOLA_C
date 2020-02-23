#include "arrays_psola.h"
#include <stdlib.h>

void init_rt_buffer(rt_buffer* rt_buff, uint16_t window_len)
{
    rt_buff->buffer = (int16_t *)malloc(3 * sizeof(int16_t) * window_len);
    
    rt_buff->p2future = rt_buff->buffer;
    rt_buff->p2present = rt_buff->buffer + sizeof(int16_t)* window_len;
    rt_buff->p2past = rt_buff->buffer + sizeof(int16_t)* window_len * 2;
}

void clear_rt_buffer(rt_buffer* rt_buff)
{
    free(rt_buff->buffer);
    free(rt_buff); 
}

void push_rt_buffer(rt_buffer* rt_buff, int16_t* p2data, uint16_t number_samples)
{
    // this is the length of the package to send
    rt_buff->num_samp_present = number_samples;
    // this is the package to send
    rt_buff->p2present = p2data;
}

void time_swap(rt_buffer* rt_buff)
{
    rt_buff->p2past = rt_buff->p2present;
    rt_buff->p2present = rt_buff->p2future;
    rt_buff->p2future = rt_buff->buffer;
}
