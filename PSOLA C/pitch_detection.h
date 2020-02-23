/* 
 * File:   pitch_detection.h
 * Author: Ramita
 *
 * Created on 23 de febrero de 2020, 02:11
 */

#ifndef PITCH_DETECTION_H
#define	PITCH_DETECTION_H

#include <stdint.h>

uint16_t find_peaks(int16_t * signal, uint32_t length, int16_t * sig_peak, float fs);
void auto_correlation(int16_t * signal, uint32_t length , int16_t * result);
int16_t get_fundamental(int16_t * signal, uint32_t length, float fs);

#endif	/* PITCH_DETECTION_H */

