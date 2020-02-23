#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"

float* linspace(float start, float end, uint16_t division)
{
    float* result = (float*)malloc(sizeof(float) * division);
    for(uint16_t i = 0; i < division; i++)
    {
        result[i] = start + i * (end - start)/(division - 1);
    }
    return result;
}


uint16_t argmin(int16_t * source, uint16_t index_start, uint16_t num_elements)
{
	uint16_t min_value = index_start;
	for(uint32_t i = index_start; i < num_elements; i++)
		if(source[i] < source[min_value])
			min_value = i;
	
	return min_value;
}


uint16_t argmax(int16_t * source, uint16_t index_start, uint16_t num_elements)
{
	uint16_t max_value = index_start;
	for(uint32_t i = index_start; i < num_elements; i++)
		if(source[i] > source[max_value])
			max_value = i;
	
	return max_value;
}


void print_array(int16_t * array, uint16_t len)
{
    printf("[");
    for(uint16_t i = 0; i < len; i++)
    {
        printf("%d,", array[i]);
    }
     printf("] \n");
}


float hamming(int16_t point, uint16_t numSamples)
{
	return 0.54f - 0.46f * cosf ( (2.0f * M_PI_F * point) / ( (float)numSamples - 1.0f) ) ;
}


float mean(int16_t * array, uint32_t len_array) 
{
   float sum = 0;

   for(int i = 0; i < len_array; i++) {
      sum += array[i];
   }
   
   return sum/len_array;
}

/*
//Falta probar
*/
void save_array(char * file_name, int16_t * array, uint32_t length)
{
	FILE * log_file = NULL;
	log_file = fopen("salido_auto_correlation.log", "wb+");

	for(int i = 0; i < 1024; i++)
		fprintf(log_file, "%.4f, ",array[i]);

	fflush(log_file);
}