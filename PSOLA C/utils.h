/* 
 * File:   utils.h
 * Author: Ramita
 */
#ifndef UTILS_H
#define	UTILS_H

#include <stdint.h>

#define M_PI_F 3.14159f

float* linspace(float start, float end, uint16_t division);
uint16_t argmin(int16_t * source, uint16_t index_start, uint16_t num_elements);
uint16_t argmax(int16_t * source, uint16_t index_start, uint16_t num_elements);
void print_array(int16_t * array, uint16_t len);
float hamming(int16_t point, uint16_t numSamples);
float mean(int16_t * array, uint32_t len_array);

#endif	/* UTILS_H */

