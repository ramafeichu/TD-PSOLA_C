//
//  FLWT.h
//  Created by Terry Kong on 2/21/15.
//

#ifndef FLWT_H
#define FLWT_H

#include <stdint.h>

#define DEFAULT_WIN_LENGTH  1024
#define MEDIAN_BUFFER_LENGTH 5

/**
 * @brief      Short class description.
 *
 * @details    Verbose description of class details.
 *
 * Note the following example code:
 * @code
 *    Window win = new Window(parent);
 *    win.show();
 * @endcode
 */

typedef struct
{
    int16_t *window;
    int16_t levels;
    int16_t *maxCount;
    int16_t *minCount;
    int16_t *maxIndices;
    int16_t *minIndices;
    float oldFreq;
    int16_t oldMode;
    int16_t *mode;
    int16_t winLength;
    int16_t dLength;
    int16_t *differs;
    float *medianBuffer5;
    int16_t medianBufferLastIndex;
}FLWT;

void initFLWT(FLWT* flwt, int16_t levels, int16_t windowLen);
void clearFLWT(FLWT* flwt);
void addToMedianBuffer(FLWT* flwt, float f);
float median5(FLWT* flwt);
// datalen MUST be divisible by 2^(levels-1)
// Returns 0.0 if it deduces the segment is pitchless
float getPitch(FLWT* flwt, int16_t* data, int16_t datalen, long fs);
float getPitchWithMedian5(FLWT* flwt, int16_t* data, int16_t datalen, long fs);
float getPitchLastReliable(FLWT* flwt, int16_t* data, int16_t datalen, long fs);
float getPitchOctaveInvariant(FLWT* flwt, int16_t* data, int16_t datalen, long fs);
float getPitchRobust(FLWT* flwt, int16_t* data, int16_t datalen, long fs);

#endif /* FLWT_H */
