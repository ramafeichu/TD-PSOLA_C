#include "PSOLA.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "utils.h"



// Some helper functions ==================

// Q15 multiplication
int16_t Q15mult(int16_t x, int16_t y) {
    long temp = (long)x * (long)y;
    temp += Q15_RESOLUTION;
    return temp >> FIXED_FBITS;
}

// Q15 wrapped addition
int16_t Q15addWrap(int16_t x, int16_t y) {
    return x + y;
}

// Q15 saturation addition
int16_t Q15addSat(int16_t x, int16_t y) {
    long temp = (long)x+(long)y;
    if (temp > 0x7FFFFFFF) temp = 0x7FFFFFFF;
    if (temp < -1*0x7FFFFFFF) temp = -1*0x7FFFFFFF;
    return (int16_t)temp;
}
/** ==============================================================================
 * @brief       Function initializes the pitch correction module.
 *
 * @details     The pitch correction module is designed for Q15 data. Pitch correction may not work correctly unless buffer is long enough for the given sampling frequency and pitch.
 *
 * @todo        Improve to TD-PSOLA algorithm to account for phase inconsistencies between buffers
 *
 * @param       bufferLen       Number of data the module expects when pitchCorrect() is called in order to optimize performance.
 *
 * @see
 *              E.moulines and W. Verhelst. Time-domain and frequency-domain techniques for prosodic modifications of speech. In W. Bastiaan Kleijn and K.K. Paliwal, editors, Speech Coding and Synthesis, chapter 15, pages 519-555. Elsevier, 1995.
 * ================================================================================
 */

void initPSOLA(PSOLA* psola, uint16_t bufferLen)
{
    if (bufferLen < 1) 
    {
        psola->bufferLen = DEFAULT_BUFFER_SIZE;
    } 
    else 
    {
        psola->bufferLen = bufferLen;
    }
    psola->workingBuffer = (int16_t *)malloc(2*sizeof(int16_t)*bufferLen);
    psola->storageBuffer = (int16_t *)malloc(2*sizeof(int16_t)*bufferLen);
    psola->window = (int16_t *)malloc(sizeof(int16_t)*bufferLen);
}

void clearPSOLA(PSOLA* psola)
{
    free(psola->workingBuffer);
    free(psola->storageBuffer);
    free(psola->window);
    free(psola);
}


void pitchCorrect(PSOLA* psola, int16_t* input, uint16_t Fs, float inputPitch, float desiredPitch)
{
    // Percent change of frequency
    float scalingFactor = 1 + (inputPitch - desiredPitch)/desiredPitch;
    // PSOLA constants
    int16_t analysisShift = ceil(Fs/inputPitch);
    int16_t analysisShiftHalfed = round(analysisShift/2);
    int16_t synthesisShift = round(analysisShift*scalingFactor);
    int16_t analysisIndex = -1;
    int16_t synthesisIndex = 0;
    int16_t analysisBlockStart;
    int16_t analysisBlockEnd;
    int16_t synthesisBlockEnd;  // ESTE ME DICE DONDE TERMINÓ
    int16_t analysisLimit = psola->bufferLen - analysisShift - 1;
    int16_t inputIndex;
    // Window declaration
    int16_t winLength = analysisShift + analysisShiftHalfed + 1;
    int16_t windowIndex;
	
    for (uint16_t i = 0; i < psola->bufferLen; i++) 
    {
        //slide the past data into the front
        psola->storageBuffer[i] = psola->storageBuffer[i + psola->bufferLen];
        //load up next set of data
        psola->storageBuffer[i + psola->bufferLen] = input[i];
    }
    bartlett(psola->window, winLength);
    // PSOLA Algorithm
    while (analysisIndex < analysisLimit) 
    {
        // Analysis blocks are two pitch periods long
        analysisBlockStart = (analysisIndex + 1) - analysisShiftHalfed;
        if (analysisBlockStart < 0) 
        {
            analysisBlockStart = 0;
        }
        analysisBlockEnd = analysisBlockStart + analysisShift + analysisShiftHalfed;
        if (analysisBlockEnd > psola->bufferLen - 1) 
        {
            analysisBlockEnd = psola->bufferLen - 1;
        }
        // Overlap and add
        synthesisBlockEnd = synthesisIndex + analysisBlockEnd - analysisBlockStart;
        inputIndex = analysisBlockStart;
        windowIndex = 0;
        for (uint16_t j = synthesisIndex; j <= synthesisBlockEnd; j++) 
        {
            psola->workingBuffer[j] = Q15addWrap(psola->workingBuffer[j], Q15mult(input[inputIndex],psola->window[windowIndex]));
            inputIndex++;
            windowIndex++;
        }
        // Update pointers
        analysisIndex += analysisShift;
        synthesisIndex += synthesisShift;
    }
//    psola->num_samp_working = analysisLimit;
//    psola->debugging = (int16_t *)malloc(sizeof(int16_t)*analysisLimit);
//    memcpy(psola->debugging, psola->workingBuffer, sizeof(int16_t)*analysisLimit);
//    print_array(psola->debugging, synthesisBlockEnd);
}

void bartlett(int16_t* window, int16_t length) 
{
    if (length < 1) 
        return;
    if (length == 1) 
    {
        window[0] = 1;
        return;
    }
    
    int16_t N = length - 1;
    int16_t middle = N >> 1;
    int16_t slope = round( ((float)(1<<(FIXED_FBITS - 1)))/ N * 4 );
    if (N%2 == 0) {
        // N even = L odd
        window[0] = 0;
        for (int16_t i = 1; i <= middle; i++) 
        {
            window[i] = window[i - 1] + slope;
        }
        for (int16_t i = middle+1; i <= N; i++)
        {
            window[i] = window[N - i];
        }
        // double check that the middle value is the maximum Q15 number
        window[middle] = LARGEST_Q15_NUM;
    } 
    else 
    {
        // N odd = L even
        window[0] = 0;
        for (int16_t i = 1; i <= middle; i++) 
        {
            window[i] = window[i-1] + slope;
        }
        window[middle + 1] = window[middle];
        for (int16_t i = middle + 1; i <= N; i++) 
        {
            window[i] = window[N - i];
        }
    }
}


int16_t* pitch_correct(int16_t* signal, int16_t* peaks, uint16_t peaks_len , float fs, float f_ratio)
{  
    uint16_t N = 512;
    int16_t * new_signal = (int16_t *)calloc(N, sizeof(int16_t));
    uint16_t len = roundf(peaks_len * f_ratio);
    float* new_peaks_ref = linspace(0, peaks_len - 1, len);
    int16_t * new_peaks = (int16_t *)calloc(len, sizeof(int16_t));
    
    // When creating the new peaks calculates the weighted sum of the original adjacent peaks
    for(uint16_t i = 0; i < len; i++)
    {
        float weight = new_peaks_ref[i] - floor(new_peaks_ref[i]);
        uint16_t left = floor(new_peaks_ref[i]);
        uint16_t right = ceil(new_peaks_ref[i]);
        new_peaks[i] = (int)(peaks[left] * (1 - weight) + peaks[right] * weight);
    }
    
    // Overlap-and-add
    uint16_t P1[2];
    for(uint16_t j; j < len; j++)
    {
        //find the corresponding old peak index
        uint16_t* abs_diff = (uint16_t *)malloc(peaks_len * sizeof(uint16_t));
        
        for(uint16_t g; g < peaks_len; g++)
        {
            abs_diff[g] = abs(peaks[g] - new_peaks[j]);
        }
        uint16_t i = argmin(abs_diff, 0, peaks_len);
        
        P1[0] = (j == 0)? new_peaks[j] : new_peaks[j] - new_peaks[j-1];
        P1[1] = (j == (len - 1))? N - 1 - new_peaks[j] : new_peaks[j+1] - new_peaks[j];
        
        //edge case truncation
        if(peaks[i] - P1[0] < 0)
        {
            P1[0] = peaks[i];
        }
        
        if(peaks[i] + P1[1] > N - 1)
        {
            P1[1] = N - 1 - peaks[i];
        }
           
        //Windowing
        //center window from original signal at the new peak
        for(uint16_t z = 0; z < P1[0] + P1[1] ; z++)
        {
            new_signal[z + new_peaks[j]] += signal[z + peaks[i]] * hamming(z, P1[0] + P1[1]);
        }
        
        free(abs_diff);
    }
    
    free(new_peaks_ref);
    free(new_peaks);
    
    return new_signal;
}


