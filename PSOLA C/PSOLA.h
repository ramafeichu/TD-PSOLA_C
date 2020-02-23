#ifndef ____PSOLA__
#define ____PSOLA__

#include <stdint.h>

#define DEFAULT_BUFFER_SIZE 512
#define FIXED_BITS        16
#define FIXED_WBITS       0
#define FIXED_FBITS       15
#define Q15_RESOLUTION   (1 << (FIXED_FBITS - 1))
#define LARGEST_Q15_NUM   32767


typedef struct
{
    uint16_t bufferLen;
    int16_t* workingBuffer;
    int16_t* storageBuffer;
    int16_t* window;
    
    int16_t num_samp_working;
    int16_t * debugging;
}PSOLA;

void initPSOLA(PSOLA* psola, uint16_t bufferLen);
void clearPSOLA(PSOLA* psola);
void pitchCorrect(PSOLA* psola, int16_t* input, uint16_t Fs, float inputPitch, float desiredPitch);
int16_t* pitch_correct(int16_t* signal, int16_t* peaks, uint16_t peaks_len , uint16_t fs, float f_ratio);
void bartlett(int16_t* window, int16_t length);

#endif /* defined(____PSOLA__) */
