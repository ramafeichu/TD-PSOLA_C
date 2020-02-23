/**
 *   @mainpage Fast Lifting Wavelet Transform Module (FLWT) for Pitch Detection
 *   @author Terry Kong
 *   @date Mar. 9, 2015
 *
 *   \section desc_sec Description
 *   This is a fast implementation of the FLWT for pitch detection. The implementation
 *      works on integer data, which is the preferred data type on embedded systems. 
 *      The algorithm is discussed in great detail in http://online.physics.uiuc.edu/courses/phys406/NSF_REU_Reports/2005_reu/Real-Time_Time-Domain_Pitch_Tracking_Using_Wavelets.pdf
 *
 *  @n Included in this class are an assortment of wrapper functions that help increase
 *      the consistency between pitch measurements. The most reliable of which involves 
 *      using a median filter to smooth out the fluctuations in frequency.
 *
 *  @n This pitch detection algorithm works best when the audio has a singal strong 
 *      fundamental harmonic. Audio with a mixture of fundamental harmonics, like a 
 *      music track with many instruments, is not likely to be processed well. Also,
 *      The change in pitch with time needs to be relatively slow in order for the
 *      algorithm to estimatethe pitch reliably. If the algorithm decides the windowed 
 *      audio is pitchless, it conservatively returns 0.
 *
 *  \section contents_sec Table of Contents
 *    FLWT.cpp
 *
 *    FLWT.h
 *
 */

/**
 *  @file FLWT.cpp
 *  @brief Source File for FLWT
 *  @file FLWT.h
 *  @brief Header File for FLWT
 */

#include "FLWT.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

// ====================================
#define DEFAULT_WIN_LENGTH  1024
#define MAX_LEVELS          6
#define MAX_INT16           32767
#define MIN_INT16           -32768
#define MAX_INT32           2147483647
#define MIN_INT32           -2147483648
// ====================================

/** Sets the number of differences between peaks to consider as the mode
 *
 * @b Example: 
 *
 * @code
 *    peakIndices = [d0,d1,d2,d3,d4,...]
 *    for index = 0:end
 *      for i = 1:MAX_NUM_OF_PEAKS_BETWEEN_MODE
 *                  possibleMode = abs( peakIndices[j+1] - peakIndices[j] )
 * @endcode
 *
 */
#define MAX_NUM_OF_PEAKS_BETWEEN_MODE   3

/** Parameter used when looking for peaks (0 ≤ GLOBAL_MAX_THRESHOLD ≤ 1) */
#define GLOBAL_MAX_THRESHOLD    0.75

/** Any frequency above max is ignored */
#define MAX_FREQUENCY           3000

/** Tolerance when checking if something is an octave off (in Hz) */
#define OCTAVE_TOLERANCE        3

/** Tolerance when checking if the change in frequency is humanly possible */
#define CHANGE_IN_FREQ_TOLERANCE        0.2

#ifndef max
#define max(x, y) ((x) > (y)) ? (x) : (y)
#endif
#ifndef min
#define min(x, y) ((x) < (y)) ? (x) : (y)
#endif

// abs value
int16_t iabs(int16_t x) 
{
    if (x >= 0) return x;
    return -x;
}

void addToMedianBuffer(FLWT* flwt, float f) 
{
    if (flwt->medianBufferLastIndex == MEDIAN_BUFFER_LENGTH) 
    {
        flwt->medianBufferLastIndex = 0;
    }
    flwt->medianBuffer5[flwt->medianBufferLastIndex] = f;
    flwt->medianBufferLastIndex++;
}

float median5(FLWT* flwt) 
{
    float a = flwt->medianBuffer5[0];
    float b = flwt->medianBuffer5[1];
    float c = flwt->medianBuffer5[2];
    float d = flwt->medianBuffer5[3];
    float e = flwt->medianBuffer5[4];
    // Put the largest value in a
    if (a <= b) { a = a+b; b = a-b; a = a-b; }
    if (a <= c) { a = a+c; c = a-c; a = a-c; }
    if (a <= d) { a = a+d; d = a-d; a = a-d; }
    if (a <= e) { a = a+e; e = a-e; a = a-e; }
    // Put the second largest value in b
    if (b <= c) { b = b+c; c = b-c; b = b-c; }
    if (b <= d) { b = b+d; d = b-d; b = b-d; }
    if (b <= e) { b = b+e; e = b-e; b = b-e; }
    // The 3rd largest is the median
    if (c <= d) { c = c+d; d = c-d; c = c-d; }
    if (c <= e) { c = c+e; e = c-e; c = c-e; }
    return c;
}

/** ==============================================================================
 * @brief       Function initializes fast lifting wavelet transform module for pitch detection.
 *
 * @details     [Verbose description of method
 *              (or function) details.]
 *
 * @todo        [ Finish writing the details]
 *
 * @param       levels       Number of levels for the FLWT algorithm
 * @param       windowLen    Length of the window
 *
 * @see        
 *              http://online.physics.uiuc.edu/courses/phys406/NSF_REU_Reports/2005_reu/Real-Time_Time-Domain_Pitch_Tracking_Using_Wavelets.pdf
 * ================================================================================
 */
void initFLWT(FLWT* flwt, int16_t levels, int16_t windowLen) 
{
    // Error Handle
    if (windowLen < 4 || windowLen > DEFAULT_WIN_LENGTH) 
    {
        flwt->window = (int16_t *)malloc(DEFAULT_WIN_LENGTH * sizeof(int16_t));
    } 
    else 
    {
        flwt->window = (int16_t *)malloc(windowLen * sizeof(int16_t));
    }
    if (levels < 0 || levels > MAX_LEVELS) 
    {
        flwt->levels = MAX_LEVELS;
    } 
    else 
    {
        flwt->levels = levels;
    }
 
    flwt->maxCount = (int16_t *)malloc(sizeof(int16_t) * levels);
    flwt->minCount = (int16_t *)malloc(sizeof(int16_t) * levels);
            
    flwt->maxIndices = (int16_t *)malloc(sizeof(int16_t) * windowLen);
    flwt->minIndices = (int16_t *)malloc(sizeof(int16_t) * windowLen);

    flwt->oldFreq = 0.0;
    flwt->oldMode = 0;
    flwt->mode = (int16_t *)malloc(sizeof(int16_t) * levels);
    flwt->winLength = windowLen;
    flwt->dLength = 0;
    // This buffer can't overflow up unless windowLen < (#peaks + #valleys)*(#peaks + #valleys + 1)/2
    flwt->differs = (int16_t *)malloc(sizeof(int16_t) * windowLen);
    // Median Buffer variables
    flwt->medianBuffer5 = (float *)malloc(sizeof(float) * 5);
    flwt->medianBufferLastIndex = 0;
}

void clearFLWT(FLWT* flwt) 
{
    free(flwt->window);
    free(flwt->maxCount);
    free(flwt->minCount);
    free(flwt->maxIndices);
    free(flwt->minIndices);
    free(flwt->mode);
    free(flwt->medianBuffer5);
    free(flwt);
}

float getPitch(FLWT* flwt, int16_t* data, int16_t datalen, long fs) 
{
    // Calculate Parameters for this window
    int16_t newWidth = (datalen > flwt->winLength) ? flwt->winLength : datalen;
    long average = 0;
    int16_t globalMax = MIN_INT16;
    int16_t globalMin = MAX_INT16;
    int16_t maxThresh;
    int16_t minThresh;
    for (int16_t i = 0; i < datalen; i++) 
    {
        flwt->window[i] = data[i];
        average += data[i];
        if (data[i] > globalMax) 
        {
            globalMax = data[i];
        }
        if (data[i] < globalMin) 
        {
            globalMin = data[i];
        }
    }
    average /= datalen;
    maxThresh = GLOBAL_MAX_THRESHOLD*(globalMax - average) + average;
    minThresh = GLOBAL_MAX_THRESHOLD*(globalMin - average) + average;
    
    // Perform FLWT Algorithm
    int16_t minDist;
    int16_t climber;
    bool isSearching; // flag for whether or not another peak can be found
    int16_t tooClose; // Make sure the peaks aren't too close
    int16_t test;
    for(int16_t lev = 0; lev < flwt->levels; lev++) 
    {
        // Reinitialize level parameters
        flwt->mode[lev] = 0;
        flwt->maxCount[lev] = 0;
        flwt->minCount[lev] = 0;
        isSearching = true;
        tooClose = 0;
        flwt->dLength = 0;
        
        newWidth = newWidth >> 1;
        minDist = max( floor((fs/MAX_FREQUENCY) >> (lev+1)) ,1);
        // First forward difference of new window (a(i,2) - a(i,1) > 0)
        if ((flwt->window[3]+flwt->window[2]-flwt->window[1]-flwt->window[0]) > 0) 
        {
            climber = 1;
        } 
        else 
        {
            climber = -1;
        }
        
        // Calculate the Approximation component (only) inplace
    
        // First and last element of the approximation is calculated separately
        //  to exploit the fact that maxima and minima can be calculated
        //  while next approximation component is being filled
        flwt->window[0] = (flwt->window[1] + flwt->window[0]) >> 1;
        for (int16_t j = 1; j < newWidth; j++) 
        {
            flwt->window[j] = (flwt->window[2*j+1] + flwt->window[2*j]) >> 1;
            
            // While the window is being filled, find max and mins
            test = flwt->window[j] - flwt->window[j-1]; // first backward difference
            
            if (climber >= 0 && test < 0) { // reached a peak
                if (flwt->window[j-1] >= maxThresh && isSearching && !tooClose) 
                {
                    // value is large enough, haven't found peak yet, and not too close
                    flwt->maxIndices[flwt->maxCount[lev]] = j-1;
                    flwt->maxCount[lev]++;
                    isSearching = false;
                    tooClose = minDist;
                }
                climber = -1;
            } 
            else if(climber <= 0 && test > 0) 
            { // reached valley
                if (flwt->window[j-1] <= minThresh && isSearching && !tooClose) 
                {
                    // value is small enough, haven't found peak yet, and not too close
                    flwt->minIndices[flwt->minCount[lev]] = j-1;
                    flwt->minCount[lev]++;
                    isSearching = false;
                    tooClose = minDist;
                }
                climber = 1;
            
            }
            
            // If we reach zero crossing, we can look for another peak
            if ((flwt->window[j] <= average && flwt->window[j-1] > average) ||
                (flwt->window[j] >= average && flwt->window[j-1] < average)) 
            {
                isSearching = true;
            }
            
            if (tooClose) 
            {
                tooClose--;
            }
        }
        
        // Find the mode distance between peaks
        if (flwt->maxCount[lev] >= 2 && flwt->minCount[lev] >= 2) 
        {
            // Find all differences between maxima/minima
            for (int16_t j = 1; j <= MAX_NUM_OF_PEAKS_BETWEEN_MODE; j++) 
            {
                for (int16_t k = 0; k < flwt->maxCount[lev] - j; k++) 
                {
                    flwt->differs[flwt->dLength] = iabs(flwt->maxIndices[k] - flwt->maxIndices[k+j]);
                    flwt->dLength++;
                }
                for (int16_t k = 0; k < flwt->minCount[lev] - j; k++) 
                {
                    flwt->differs[flwt->dLength] = iabs(flwt->minIndices[k] - flwt->minIndices[k+j]);
                    flwt->dLength++;
                }
            }
            
            // Determine the mode
            int16_t numer = 1; // Require at least two agreeing _differs to yield a mode
            int16_t numerJ;
            for (int16_t j = 0; j < flwt->dLength; j++) 
            {
                numerJ = 0;
                
                // Find the number of _differs that are near _differs[j]
                for (int16_t n = 0; n < flwt->dLength; n++) 
                {
                    if (iabs(flwt->differs[j] - flwt->differs[n]) < minDist) 
                    {
                        numerJ++;
                    }
                }
                
                // Check to see if there is a better candidate for the mode
                if (numerJ >= numer && numerJ > floor((newWidth/flwt->differs[j])>>2)) 
                {
                    if (numerJ == numer) 
                    {
                        if (flwt->oldMode && iabs(flwt->differs[j] - (flwt->oldMode >> (lev+1))) < minDist) 
                        {
                            flwt->mode[lev] = flwt->differs[j];
                        } 
                        else if (~flwt->oldMode && (flwt->differs[j] > 1.95*flwt->mode[lev] && flwt->differs[j] < 2.05*flwt->mode[lev])) 
                        {
                            flwt->mode[lev] = flwt->differs[j];
                        }
                    } 
                    else 
                    {
                        numer = numerJ;
                        flwt->mode[lev] = flwt->differs[j];
                    }
                } 
                else if (numerJ == numer-1 && flwt->oldMode && iabs(flwt->differs[j] - (flwt->oldMode >> (lev+1))) < minDist) 
                {
                    flwt->mode[lev] = flwt->differs[j];
                }
            }
            
            // Average to get the mode
            if (flwt->mode[lev]) 
            {
                int16_t numerator = 0;
                int16_t denominator = 0;
                for (int16_t m = 0; m < flwt->dLength; m++) 
                {
                    if (iabs(flwt->mode[lev] - flwt->differs[m]) <= minDist) 
                    {
                        numerator += flwt->differs[m];
                        denominator++;
                    }
                }
                flwt->mode[lev] = numerator/denominator;
            }
            
            // Check if the mode is shared with the previous level
            if (lev == 0) 
            {
                // Do nothing
            } 
            else if (flwt->mode[lev-1] && flwt->maxCount[lev-1] >= 2 && flwt->minCount[lev-1] >= 2) 
            {
                // If the modes are within a sample of one another, return the calculated frequency
                if (iabs(flwt->mode[lev-1] - 2*flwt->mode[lev]) <= minDist) 
                {
                    flwt->oldMode = flwt->mode[lev-1];
                    flwt->oldFreq = ((float)fs)/((float)flwt->mode[lev-1])/((float)(1<<(lev)));
                    // Add the frequency to the median buffer
                    addToMedianBuffer(flwt, flwt->oldFreq);
                    return flwt->oldFreq;
                }
            }
            
        }
    }
    
    // Getting here means the window was pitchless
    
    // Add to this value to the median filter
    addToMedianBuffer(flwt, 0.0);
    return 0.0;
}

/** ====================================================
 * @brief       Calculates the pitch of a set of data using a median filter.
 *
 * @details     Uses a median filter to calculate the median on the past 5
 *              frequencies including the newest one calculated from the new
 *              set of data. Uses the same algorithm as described in 
 *              FLWT::getpitch. This has a slightly larger overhead to return
 *              the pitch.
 *
 * @param       data        Pointer to array of data
 * @param       datalen     Length of the array
 * @param       fs          Sampling frequency of data
 *
 * @return      Pitch of the window passed through a median filter
 *
 * @retval      pitch
 *
 * ======================================================
 */
//float FLWT::getPitchWithMedian5(int* data, int datalen, long fs) {
//    this->getPitch(data,datalen,fs);
//    return this->median5();
//}

float getPitchWithMedian5(FLWT* flwt, int16_t* data, int16_t datalen, long fs)
{
    getPitch(flwt, data, datalen, fs);
    return median5(flwt);
}

/** ====================================================
 * @brief       Calculates the pitch of a set of data.
 *
 * @details     Uses the same algorithm as described in FLWT::getpitch. 
 *              If the current data is pitchless, this returns the last
 *              non-zero result.
 *
 * @param       data        Pointer to array of data
 * @param       datalen     Length of the array
 * @param       fs          Sampling frequency of data
 *
 * @return      Pitch of the window (non-zero)
 *
 * @retval      pitch
 *
 * ======================================================
 */
//float FLWT::getPitchLastReliable(int* data, int datalen, long fs) {
//    this->getPitch(data,datalen,fs);
//    return _oldFreq;
//}

float getPitchLastReliable(FLWT* flwt, int16_t* data, int16_t datalen, long fs) 
{
    getPitch(flwt, data, datalen, fs);
    return flwt->oldFreq;
}

/** ====================================================
 * @brief       Calculates the pitch of a set of data.
 *
 * @details     Uses the same algorithm as described in FLWT::getpitch.
 *              If the new pitch is either twice or half the frequency of 
 *              the previous window, the new pitch is considered an error
 *              and the return value is the pitch of the last window.
 *
 * @param       data        Pointer to array of data
 * @param       datalen     Length of the array
 * @param       fs          Sampling frequency of data
 *
 * @return      Pitch of the window that is octave invariant
 *
 * @retval      pitch
 *
 * ======================================================
 */
//float FLWT::getPitchOctaveInvariant(int* data, int datalen, long fs) {
//    float old = _oldFreq;
//    this->getPitch(data,datalen,fs);
//    if (_oldFreq >= (2*old - OCTAVE_TOLERANCE) && _oldFreq <= (2*old + OCTAVE_TOLERANCE)) {
//        _oldFreq = old;
//    } else if (old >= (2*_oldFreq - OCTAVE_TOLERANCE) && old <= (2*_oldFreq + OCTAVE_TOLERANCE)) {
//        _oldFreq = old;
//    }
//    return _oldFreq;
//}

float getPitchOctaveInvariant(FLWT* flwt, int16_t* data, int16_t datalen, long fs)
{
    float old = flwt->oldFreq;
    getPitch(flwt, data,datalen,fs);
    if (flwt->oldFreq >= (2*old - OCTAVE_TOLERANCE) && flwt->oldFreq <= (2*old + OCTAVE_TOLERANCE))
    {
        flwt->oldFreq = old;
    } 
    else if (old >= (2*flwt->oldFreq - OCTAVE_TOLERANCE) && old <= (2*flwt->oldFreq + OCTAVE_TOLERANCE))
    {
        flwt->oldFreq = old;
    }
    return flwt->oldFreq;
}

/** ====================================================
 * @brief       Calculates the pitch of a set of data (robustly).
 *
 * @details     Uses the same algorithm as described in FLWT::getpitch.
 *              It uses a robust algorithm that is a combination of all the
 *              methods included in this class.
 *
 * @param       data        Pointer to array of data
 * @param       datalen     Length of the array
 * @param       fs          Sampling frequency of data
 *
 * @return      Pitch of the window
 *
 * @retval      pitch
 *
 * @todo        Make this function more reliable and robust
 *
 * ======================================================
 */
float getPitchRobust(FLWT* flwt, int16_t* data, int16_t datalen, long fs) {
    // Calculate Parameters for this window
    float currentFreq;
    int16_t newWidth = (datalen > flwt->winLength) ? flwt->winLength : datalen;
    long average = 0;
    int16_t globalMax = MIN_INT16;
    int16_t globalMin = MAX_INT16;
    int16_t maxThresh;
    int16_t minThresh;
    for (int16_t i = 0; i < datalen; i++) 
    {
        flwt->window[i] = data[i];
        average += data[i];
        if (data[i] > globalMax) 
        {
            globalMax = data[i];
        }
        if (data[i] < globalMin) 
        {
            globalMin = data[i];
        }
    }
    average /= datalen;
    maxThresh = GLOBAL_MAX_THRESHOLD*(globalMax - average) + average;
    minThresh = GLOBAL_MAX_THRESHOLD*(globalMin - average) + average;
    
    // Perform FLWT Algorithm
    int16_t minDist;
    int16_t climber;
    bool isSearching; // flag for whether or not another peak can be found
    int16_t tooClose; // Make sure the peaks aren't too close
    int16_t test;
    for(int16_t lev = 0; lev < flwt->levels; lev++) 
    {
        // Reinitialize level parameters
        flwt->mode[lev] = 0;
        flwt->maxCount[lev] = 0;
        flwt->minCount[lev] = 0;
        isSearching = true;
        tooClose = 0;
        flwt->dLength = 0;
        
        newWidth = newWidth >> 1;
        minDist = max( floor((fs/MAX_FREQUENCY) >> (lev+1)) ,1);
        // First forward difference of new window (a(i,2) - a(i,1) > 0)
        if ((flwt->window[3]+flwt->window[2]-flwt->window[1]-flwt->window[0]) > 0) 
        {
            climber = 1;
        } 
        else 
        {
            climber = -1;
        }
        
        // Calculate the Approximation component (only) inplace
        
        // First and last element of the approximation is calculated separately
        //  to exploit the fact that maxima and minima can be calculated
        //  while next approximation component is being filled
        flwt->window[0] = (flwt->window[1] + flwt->window[0]) >> 1;
        for (int16_t j = 1; j < newWidth; j++) 
        {
            flwt->window[j] = (flwt->window[2*j+1] + flwt->window[2*j]) >> 1;
            
            // While the window is being filled, find max and mins
            test = flwt->window[j] - flwt->window[j-1]; // first backward difference
            
            if (climber >= 0 && test < 0) 
            { // reached a peak
                if (flwt->window[j-1] >= maxThresh && isSearching && !tooClose) 
                {
                    // value is large enough, haven't found peak yet, and not too close
                    flwt->maxIndices[flwt->maxCount[lev]] = j-1;
                    flwt->maxCount[lev]++;
                    isSearching = false;
                    tooClose = minDist;
                }
                climber = -1;
            } else if(climber <= 0 && test > 0) 
            { // reached valley
                if (flwt->window[j-1] <= minThresh && isSearching && !tooClose) 
                {
                    // value is small enough, haven't found peak yet, and not too close
                    flwt->minIndices[flwt->minCount[lev]] = j-1;
                    flwt->minCount[lev]++;
                    isSearching = false;
                    tooClose = minDist;
                }
                climber = 1;
                
            }
            
            // If we reach zero crossing, we can look for another peak
            if ((flwt->window[j] <= average && flwt->window[j-1] > average) ||
                (flwt->window[j] >= average && flwt->window[j-1] < average)) 
            {
                isSearching = true;
            }
            
            if (tooClose) 
            {
                tooClose--;
            }
        }
        
        // Find the mode distance between peaks
        if (flwt->maxCount[lev] >= 2 && flwt->minCount[lev] >= 2) 
        {
            
            // Find all differences between maxima/minima
            for (int16_t j = 1; j <= MAX_NUM_OF_PEAKS_BETWEEN_MODE; j++) 
            {
                for (int16_t k = 0; k < flwt->maxCount[lev] - j; k++) 
                {
                    flwt->differs[flwt->dLength] = iabs(flwt->maxIndices[k] - flwt->maxIndices[k+j]);
                    flwt->dLength++;
                }
                for (int16_t k = 0; k < flwt->minCount[lev] - j; k++)
                {
                    flwt->differs[flwt->dLength] = iabs(flwt->minIndices[k] - flwt->minIndices[k+j]);
                    flwt->dLength++;
                }
            }
            
            // Determine the mode
            int16_t numer = 1; // Require at least two agreeing _differs to yield a mode
            int16_t numerJ;
            for (int16_t j = 0; j < flwt->dLength; j++) 
            {
                numerJ = 0;
                
                // Find the number of _differs that are near _differs[j]
                for (int16_t n = 0; n < flwt->dLength; n++) 
                {
                    if (iabs(flwt->differs[j] - flwt->differs[n]) < minDist) 
                    {
                        numerJ++;
                    }
                }
                
                // Check to see if there is a better candidate for the mode
                if (numerJ >= numer && numerJ > floor((newWidth/flwt->differs[j])>>2)) 
                {
                    if (numerJ == numer) 
                    {
                        if (flwt->oldMode && iabs(flwt->differs[j] - (flwt->oldMode >> (lev+1))) < minDist) 
                        {
                            flwt->mode[lev] = flwt->differs[j];
                        } else if (~flwt->oldMode && (flwt->differs[j] > 1.95*flwt->mode[lev] && flwt->differs[j] < 2.05*flwt->mode[lev])) 
                        {
                            flwt->mode[lev] = flwt->differs[j];
                        }
                    } 
                    else 
                    {
                        numer = numerJ;
                        flwt->mode[lev] = flwt->differs[j];
                    }
                } 
                else if (numerJ == numer-1 && flwt->oldMode && iabs(flwt->differs[j] - (flwt->oldMode >> (lev+1))) < minDist) 
                {
                    flwt->mode[lev] = flwt->differs[j];
                }
            }
            
            // Average to get the mode
            if (flwt->mode[lev]) 
            {
                int16_t numerator = 0;
                int16_t denominator = 0;
                for (int16_t m = 0; m < flwt->dLength; m++) 
                {
                    if (iabs(flwt->mode[lev] - flwt->differs[m]) <= minDist) 
                    {
                        numerator += flwt->differs[m];
                        denominator++;
                    }
                }
                flwt->mode[lev] = numerator/denominator;
            }
            
            // Check if the mode is shared with the previous level
            if (lev == 0) 
            {
                // Do nothing
            } else if (flwt->mode[lev-1] && flwt->maxCount[lev-1] >= 2 && flwt->minCount[lev-1] >= 2) 
            {
                // If the modes are within a sample of one another, return the calculated frequency
                if (iabs(flwt->mode[lev-1] - 2*flwt->mode[lev]) <= minDist) 
                {
                    flwt->oldMode = flwt->mode[lev-1];
                    currentFreq = ((float)fs)/((float)flwt->mode[lev-1])/((float)(1<<(lev)));
                    // Add the frequency to the median buffer
                    //addToMedianBuffer(_oldFreq);
                    //return _oldFreq;
                    goto wrap_up;
                }
            }
            
        }
    }
    
    // Getting here means the window was pitchless
    currentFreq = 0.0;
    // Add to this value to the median filter
    //addToMedianBuffer(0.0);
    //return 0.0;
    
	float currentMedian;

wrap_up:

    if (currentFreq == 0.0) 
    {
        currentMedian = median5(flwt);
        // check if median is zero
        if (currentMedian == 0.0) 
        {
            addToMedianBuffer(flwt, flwt->oldFreq);
            return flwt->oldFreq;
        } 
        else 
        {
            addToMedianBuffer(flwt, currentMedian);
            return currentMedian;
        }
    } 
    else 
    {
        addToMedianBuffer(flwt, currentFreq);
        flwt->oldFreq = median5(flwt);
        return flwt->oldFreq;
    }
}