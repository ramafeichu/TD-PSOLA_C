#include <stdint.h>
#include <stdlib.h>
#include "pitch_detection.h"
#include "utils.h"
        
#define MAX_HZ 950
#define MIN_HZ 75
#define MAX_CHANGE	1.3f
#define MIN_CHANGE	0.7f
#define ANALY_WIN_MS	20

//	arm_correlate_fast_q15	(const q15_t * 	pSrcA, uint32_t 	srcALen, 
//				const q15_t * 	pSrcB, uint32_t 	srcBLen, q15_t * 	pDst)	

uint16_t find_peaks(int16_t * signal, uint32_t length, int16_t * sig_peak, float fs)
{
	uint16_t * autocorr = NULL;
	uint16_t * periods = NULL;
	
	uint32_t autocorr_peak_index = 0;
    uint32_t periods_peak_index = 0;
	
    uint32_t idx ,max_change_index, min_change_index;

    uint16_t prev;
    int32_t N = length;
    int16_t min_period = fs / MAX_HZ;
    int16_t max_period = fs / MIN_HZ;
	uint16_t number_peak = 0;
    // compute pitch periodicity
    uint32_t sequence_len = ANALY_WIN_MS / 1000.0f * fs; //sequence length in samples
	
	//Memory allocation
	autocorr = (uint16_t *) malloc(sizeof(uint16_t) * length);
	periods = (uint16_t *) malloc(sizeof(uint16_t) * (length/sequence_len) );
	
	for(uint32_t i = 0; i < N ; i += sequence_len)
	{
		auto_correlation(signal + i, sequence_len, autocorr);
		//Agrega el peak detectado
		periods[autocorr_peak_index++] = min_period + argmax(autocorr, min_period, max_period);
	}

    sig_peak[periods_peak_index++] = argmax(signal, 0, periods[0] * 1.1f);
    
    while(1)
    {
        prev = sig_peak[periods_peak_index-1];//Ultimo elemento 
        idx = prev/sequence_len;// current autocorrelation analysis window
		
        max_change_index = prev + sig_peak[idx] * MAX_CHANGE;
        min_change_index = prev + sig_peak[idx] * MIN_CHANGE; 
        if(prev + max_change_index >= N)
            break;
		
        //find maximum near expected location
        sig_peak[periods_peak_index++] = min_change_index 
				+ argmax(signal,min_change_index, max_change_index);
        number_peak++;
    }
	
	free(autocorr);
	free(periods);
    return number_peak;
}
			
void auto_correlation(int16_t * signal, uint32_t length , int16_t * result)
{
    int16_t sum;
    for(uint32_t i = 0; i < length ; i++) {
        sum=0;
        for (uint32_t j = 0; j < length-i ; j++) {
            sum += signal[j] * signal[j+i];
        }
        result[i]=sum;
    }
}

int16_t get_fundamental(int16_t * signal, uint32_t length, float fs)
{
    uint16_t * sig_peak = NULL;
    uint16_t * periods = NULL;
	uint32_t sequence_len = ANALY_WIN_MS / 1000.0f * fs; //sequence length in samples
    uint16_t period, len_peaks;

	sig_peak = (uint16_t *) malloc(sizeof(uint16_t) * length);
	periods = (uint16_t *) malloc(sizeof(uint16_t) * (length/sequence_len) );

    len_peaks = find_peaks(signal, length, sig_peak, fs);
    
    for(int i = 0; i < len_peaks; i++)
        periods[i] = sig_peak[i + 1] - sig_peak[i];
    
    period = mean(periods, len_peaks);
	free(sig_peak);
	free(periods);
    return fs/period;
}