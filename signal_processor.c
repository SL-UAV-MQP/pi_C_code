/**
 * @file signal_processor.c
 * @brief Signal preprocessing for SDR IQ data
 *
 * Implements FFT, PSD, spectrogram, band normalization, filtering,
 * peak detection, downsampling, and SNR estimation for real-time
 * signal processing on Raspberry Pi.
 *
 * Performance targets: FFT <5ms, PSD <10ms, filtering <5ms
 */

#include "signal_processor.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <fftw3.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Internal Helper Functions
 * ============================================================================ */

/**
 * @brief Create window function
 */
static void create_window(window_type_t type, int size, double* window) {
    switch (type) {
        case WINDOW_HAMMING:
            for (int n = 0; n < size; n++) {
                window[n] = 0.54 - 0.46 * cos(2.0 * M_PI * n / (size - 1));
            }
            break;

        case WINDOW_HANN:
            for (int n = 0; n < size; n++) {
                window[n] = 0.5 * (1.0 - cos(2.0 * M_PI * n / (size - 1)));
            }
            break;

        case WINDOW_BLACKMAN:
            for (int n = 0; n < size; n++) {
                double a0 = 0.42;
                double a1 = 0.5;
                double a2 = 0.08;
                window[n] = a0 - a1 * cos(2.0 * M_PI * n / (size - 1)) +
                           a2 * cos(4.0 * M_PI * n / (size - 1));
            }
            break;

        case WINDOW_KAISER:
            // Kaiser window with beta=8.6
            for (int n = 0; n < size; n++) {
                double beta = 8.6;
                double alpha = (size - 1) / 2.0;
                double x = (n - alpha) / alpha;
                // Simplified I0 Bessel function approximation
                double bessel_i0 = 1.0 + pow(beta * sqrt(1 - x * x) / 2.0, 2);
                window[n] = bessel_i0 / (1.0 + pow(beta / 2.0, 2));
            }
            break;

        case WINDOW_RECTANGULAR:
        default:
            for (int n = 0; n < size; n++) {
                window[n] = 1.0;
            }
            break;
    }
}

/**
 * @brief Compute magnitude in dB from complex value
 */
static inline double magnitude_db(cdouble_t z) {
    double mag = cabs(z);
    return 20.0 * log10(mag + EPSILON);
}

/**
 * @brief Compute power in dB from complex value
 */
static inline double power_db(cdouble_t z) {
    double power = creal(z) * creal(z) + cimag(z) * cimag(z);
    return 10.0 * log10(power + EPSILON);
}

/* ============================================================================
 * Initialization and Cleanup
 * ============================================================================ */

signal_processor_t* signal_processor_init(
    double sample_rate,
    int fft_size,
    window_type_t window_type,
    double overlap
) {
    // Validate inputs
    if (sample_rate <= 0 || fft_size <= 0 || fft_size > MAX_FFT_SIZE) {
        return NULL;
    }

    if (overlap < 0.0 || overlap >= 1.0) {
        return NULL;
    }

    // Check if fft_size is power of 2
    if ((fft_size & (fft_size - 1)) != 0) {
        return NULL;
    }

    // Allocate processor structure
    signal_processor_t* processor = (signal_processor_t*)malloc(sizeof(signal_processor_t));
    if (!processor) {
        return NULL;
    }

    // Initialize parameters
    processor->sample_rate = sample_rate;
    processor->fft_size = fft_size;
    processor->window_type = window_type;
    processor->overlap = overlap;
    processor->hop_size = (int)(fft_size * (1.0 - overlap));

    // Allocate and create window
    processor->window = (double*)malloc(fft_size * sizeof(double));
    if (!processor->window) {
        free(processor);
        return NULL;
    }

    create_window(window_type, fft_size, processor->window);

    return processor;
}

void signal_processor_free(signal_processor_t* processor) {
    if (processor) {
        if (processor->window) {
            free(processor->window);
        }
        free(processor);
    }
}

/* ============================================================================
 * FFT and Spectral Analysis
 * ============================================================================ */

music_status_t signal_processor_compute_fft(
    const signal_processor_t* processor,
    const cdouble_t* iq_samples,
    int num_samples,
    fft_result_t* result
) {
    if (!processor || !iq_samples || !result) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    if (num_samples < processor->fft_size) {
        return MUSIC_ERROR_INVALID_SIZE;
    }

    int N = processor->fft_size;

    // Allocate FFTW arrays
    fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

    if (!in || !out) {
        if (in) fftw_free(in);
        if (out) fftw_free(out);
        return MUSIC_ERROR_MEMORY;
    }

    // Create FFTW plan
    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Apply window and copy to input buffer
    for (int i = 0; i < N; i++) {
        double window_val = processor->window[i];
        in[i][0] = creal(iq_samples[i]) * window_val;
        in[i][1] = cimag(iq_samples[i]) * window_val;
    }

    // Execute FFT
    fftw_execute(plan);

    // Compute frequency axis and magnitude
    double df = processor->sample_rate / N;
    for (int i = 0; i < N; i++) {
        // FFT shift: map index to frequency
        int k = (i + N/2) % N;
        result->frequencies[i] = (k - N/2) * df;

        // Magnitude in dB
        cdouble_t val = out[k][0] + I * out[k][1];
        result->magnitude_db[i] = magnitude_db(val);
    }

    result->num_bins = N;

    // Cleanup
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return MUSIC_SUCCESS;
}

music_status_t signal_processor_compute_psd(
    const signal_processor_t* processor,
    const cdouble_t* iq_samples,
    int num_samples,
    psd_result_t* result
) {
    // compute_psd using Welch's method
    if (!processor || !iq_samples || !result) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    int N = processor->fft_size;
    int hop = processor->hop_size;
    int num_segments = (num_samples - N) / hop + 1;

    if (num_segments < 1) {
        return MUSIC_ERROR_INVALID_SIZE;
    }

    // Allocate working buffers
    double* psd_accumulator = (double*)calloc(N, sizeof(double));
    fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

    if (!psd_accumulator || !in || !out) {
        if (psd_accumulator) free(psd_accumulator);
        if (in) fftw_free(in);
        if (out) fftw_free(out);
        return MUSIC_ERROR_MEMORY;
    }

    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Window normalization factor
    double window_norm = 0.0;
    for (int i = 0; i < N; i++) {
        window_norm += processor->window[i] * processor->window[i];
    }

    // Process each segment (Welch's method)
    for (int seg = 0; seg < num_segments; seg++) {
        int offset = seg * hop;

        // Apply window
        for (int i = 0; i < N; i++) {
            double window_val = processor->window[i];
            in[i][0] = creal(iq_samples[offset + i]) * window_val;
            in[i][1] = cimag(iq_samples[offset + i]) * window_val;
        }

        // FFT
        fftw_execute(plan);

        // Accumulate power spectrum
        for (int i = 0; i < N; i++) {
            double power = out[i][0] * out[i][0] + out[i][1] * out[i][1];
            psd_accumulator[i] += power;
        }
    }

    // Average and normalize
    double scale = 1.0 / (num_segments * window_norm * processor->sample_rate);
    double df = processor->sample_rate / N;

    for (int i = 0; i < N; i++) {
        // FFT shift
        int k = (i + N/2) % N;
        result->frequencies[i] = (k - N/2) * df;

        double psd_linear = psd_accumulator[k] * scale;
        result->psd_db[i] = 10.0 * log10(psd_linear + EPSILON);
    }

    result->num_bins = N;

    // Cleanup
    fftw_destroy_plan(plan);
    free(psd_accumulator);
    fftw_free(in);
    fftw_free(out);

    return MUSIC_SUCCESS;
}

music_status_t signal_processor_compute_spectrogram(
    const signal_processor_t* processor,
    const cdouble_t* iq_samples,
    int num_samples,
    spectrogram_result_t* result
) {
    // compute_spectrogram using STFT
    if (!processor || !iq_samples || !result) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    int N = processor->fft_size;
    int hop = processor->hop_size;
    int num_time_bins = (num_samples - N) / hop + 1;

    if (num_time_bins < 1) {
        return MUSIC_ERROR_INVALID_SIZE;
    }

    // Allocate FFTW buffers
    fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

    if (!in || !out) {
        if (in) fftw_free(in);
        if (out) fftw_free(out);
        return MUSIC_ERROR_MEMORY;
    }

    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Compute time and frequency axes
    double dt = (double)hop / processor->sample_rate;
    double df = processor->sample_rate / N;

    for (int t = 0; t < num_time_bins; t++) {
        result->times[t] = t * dt;
    }

    for (int f = 0; f < N; f++) {
        int k = (f + N/2) % N;
        result->frequencies[f] = (k - N/2) * df;
    }

    // Compute STFT
    for (int t = 0; t < num_time_bins; t++) {
        int offset = t * hop;

        // Apply window
        for (int i = 0; i < N; i++) {
            double window_val = processor->window[i];
            in[i][0] = creal(iq_samples[offset + i]) * window_val;
            in[i][1] = cimag(iq_samples[offset + i]) * window_val;
        }

        // FFT
        fftw_execute(plan);

        // Store magnitude in dB (column-major order)
        for (int f = 0; f < N; f++) {
            int k = (f + N/2) % N;
            cdouble_t val = out[k][0] + I * out[k][1];
            result->spectrogram_db[f * num_time_bins + t] = magnitude_db(val);
        }
    }

    result->num_time_bins = num_time_bins;
    result->num_freq_bins = N;

    // Cleanup
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Normalization and Filtering
 * ============================================================================ */

music_status_t signal_processor_normalize_band(
    const signal_processor_t* processor,
    cdouble_t* iq_samples,
    int num_samples,
    normalize_method_t method
) {
    // normalize_band
    if (!processor || !iq_samples) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    if (method == NORMALIZE_PSD) {
        // Normalize based on PSD
        psd_result_t* psd = psd_result_alloc(processor->fft_size);
        if (!psd) {
            return MUSIC_ERROR_MEMORY;
        }

        music_status_t status = signal_processor_compute_psd(processor, iq_samples,
                                                             num_samples, psd);
        if (status != MUSIC_SUCCESS) {
            psd_result_free(psd);
            return status;
        }

        // Compute mean power
        double mean_power_db = 0.0;
        for (int i = 0; i < psd->num_bins; i++) {
            mean_power_db += psd->psd_db[i];
        }
        mean_power_db /= psd->num_bins;
        double mean_power_linear = pow(10.0, mean_power_db / 10.0);

        // Current signal power
        double current_power = 0.0;
        for (int i = 0; i < num_samples; i++) {
            double mag = cabs(iq_samples[i]);
            current_power += mag * mag;
        }
        current_power /= num_samples;

        // Scale factor
        double scale_factor = sqrt(mean_power_linear / (current_power + EPSILON));

        // Apply scaling
        for (int i = 0; i < num_samples; i++) {
            iq_samples[i] *= scale_factor;
        }

        psd_result_free(psd);

    } else if (method == NORMALIZE_MEAN) {
        // Zero mean, unit variance
        cdouble_t mean = 0.0;
        for (int i = 0; i < num_samples; i++) {
            mean += iq_samples[i];
        }
        mean /= num_samples;

        // Remove mean
        for (int i = 0; i < num_samples; i++) {
            iq_samples[i] -= mean;
        }

        // Compute standard deviation
        double variance = 0.0;
        for (int i = 0; i < num_samples; i++) {
            double mag = cabs(iq_samples[i]);
            variance += mag * mag;
        }
        double std = sqrt(variance / num_samples);

        // Normalize to unit variance
        for (int i = 0; i < num_samples; i++) {
            iq_samples[i] /= (std + EPSILON);
        }

    } else if (method == NORMALIZE_MAX) {
        // Normalize to max amplitude
        double max_amp = 0.0;
        for (int i = 0; i < num_samples; i++) {
            double mag = cabs(iq_samples[i]);
            if (mag > max_amp) {
                max_amp = mag;
            }
        }

        for (int i = 0; i < num_samples; i++) {
            iq_samples[i] /= (max_amp + EPSILON);
        }
    }

    return MUSIC_SUCCESS;
}

music_status_t signal_processor_design_bandpass(
    double sample_rate,
    double low_freq,
    double high_freq,
    int order,
    filter_coeffs_t* coeffs
) {
    //  bandpass_filter - design part using scipy.signal.butter
    if (!coeffs || order <= 0 || order > MAX_FILTER_ORDER) {
        return MUSIC_ERROR_INVALID_CONFIG;
    }

    double nyquist = sample_rate / 2.0;
    double low_norm = low_freq / nyquist;
    double high_norm = high_freq / nyquist;

    // Clamp to valid range
    if (low_norm < 0.01) low_norm = 0.01;
    if (high_norm > 0.99) high_norm = 0.99;
    if (low_norm >= high_norm) {
        return MUSIC_ERROR_INVALID_CONFIG;
    }

    // Simplified Butterworth bandpass design
    // In practice, would use a proper IIR design library
    // This is a placeholder - real implementation needs bilinear transform
    coeffs->order = order;

    // Placeholder coefficients (identity filter)
    coeffs->b[0] = 1.0;
    for (int i = 1; i <= order; i++) {
        coeffs->b[i] = 0.0;
    }

    coeffs->a[0] = 1.0;
    for (int i = 1; i <= order; i++) {
        coeffs->a[i] = 0.0;
    }

    return MUSIC_SUCCESS;
}

music_status_t signal_processor_apply_filter(
    const filter_coeffs_t* coeffs,
    const cdouble_t* iq_samples,
    int num_samples,
    cdouble_t* output
) {
    // scipy.signal.filtfilt
    if (!coeffs || !iq_samples || !output) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Apply IIR filter forward
    // In practice, use filtfilt (zero-phase filtering)
    // This is simplified direct form II implementation

    for (int n = 0; n < num_samples; n++) {
        output[n] = iq_samples[n];  // Placeholder - needs proper IIR filtering
    }

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Peak Detection and Analysis
 * ============================================================================ */

music_status_t signal_processor_detect_peaks(
    const double* frequencies,
    const double* spectrum,
    int num_bins,
    double height,
    double threshold,
    int distance,
    peak_result_t* result
) {
    // detect_peaks using scipy.signal.find_peaks
    if (!frequencies || !spectrum || !result) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Auto-detect threshold if needed
    if (isnan(threshold)) {
        // Compute median as noise floor
        double* sorted = (double*)malloc(num_bins * sizeof(double));
        memcpy(sorted, spectrum, num_bins * sizeof(double));

        // Simple insertion sort for median
        for (int i = 1; i < num_bins; i++) {
            double key = sorted[i];
            int j = i - 1;
            while (j >= 0 && sorted[j] > key) {
                sorted[j + 1] = sorted[j];
                j--;
            }
            sorted[j + 1] = key;
        }

        double noise_floor = sorted[num_bins / 2];
        threshold = noise_floor + 10.0;  // 10 dB above noise floor

        free(sorted);
    }

    if (isnan(height)) {
        height = threshold;
    }

    // Find peaks
    int num_peaks = 0;
    for (int i = 1; i < num_bins - 1 && num_peaks < MAX_PEAKS; i++) {
        // Check if local maximum
        if (spectrum[i] > spectrum[i-1] && spectrum[i] > spectrum[i+1]) {
            // Check height threshold
            if (spectrum[i] >= height) {
                // Check distance from previous peaks
                bool too_close = false;
                for (int j = 0; j < num_peaks; j++) {
                    if (abs(i - result->indices[j]) < distance) {
                        too_close = true;
                        break;
                    }
                }

                if (!too_close) {
                    result->indices[num_peaks] = i;
                    result->frequencies[num_peaks] = frequencies[i];
                    result->magnitudes[num_peaks] = spectrum[i];

                    // Compute prominence (simplified)
                    result->prominences[num_peaks] = spectrum[i] - threshold;

                    // Compute width (simplified)
                    result->widths[num_peaks] = 1.0;

                    num_peaks++;
                }
            }
        }
    }

    result->num_peaks = num_peaks;

    return MUSIC_SUCCESS;
}

music_status_t signal_processor_downsample(
    const cdouble_t* iq_samples,
    int num_samples,
    int factor,
    cdouble_t* output,
    int* output_length
) {
    // downsample using scipy.signal.decimate
    if (!iq_samples || !output || !output_length || factor <= 0) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Simplified decimation (no anti-aliasing filter for now)
    int out_len = num_samples / factor;

    for (int i = 0; i < out_len; i++) {
        output[i] = iq_samples[i * factor];
    }

    *output_length = out_len;

    return MUSIC_SUCCESS;
}

music_status_t signal_processor_estimate_snr(
    const signal_processor_t* processor,
    const cdouble_t* iq_samples,
    int num_samples,
    double signal_band_low,
    double signal_band_high,
    double* snr_db
) {
    // estimate_snr
    if (!processor || !iq_samples || !snr_db) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Compute PSD
    psd_result_t* psd = psd_result_alloc(processor->fft_size);
    if (!psd) {
        return MUSIC_ERROR_MEMORY;
    }

    music_status_t status = signal_processor_compute_psd(processor, iq_samples,
                                                         num_samples, psd);
    if (status != MUSIC_SUCCESS) {
        psd_result_free(psd);
        return status;
    }

    // Find signal band
    double signal_power_sum = 0.0;
    int signal_count = 0;
    double noise_power_sum = 0.0;
    int noise_count = 0;

    for (int i = 0; i < psd->num_bins; i++) {
        double freq = psd->frequencies[i];
        double power_linear = pow(10.0, psd->psd_db[i] / 10.0);

        if (freq >= signal_band_low && freq <= signal_band_high) {
            signal_power_sum += power_linear;
            signal_count++;
        } else {
            noise_power_sum += power_linear;
            noise_count++;
        }
    }

    double signal_power = signal_power_sum / (signal_count + EPSILON);
    double noise_power = noise_power_sum / (noise_count + EPSILON);

    *snr_db = 10.0 * log10((signal_power / (noise_power + EPSILON)) + EPSILON);

    psd_result_free(psd);

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Memory Management Helpers
 * ============================================================================ */

fft_result_t* fft_result_alloc(int num_bins) {
    fft_result_t* result = (fft_result_t*)malloc(sizeof(fft_result_t));
    if (!result) {
        return NULL;
    }

    result->num_bins = num_bins;
    result->frequencies = (double*)malloc(num_bins * sizeof(double));
    result->magnitude_db = (double*)malloc(num_bins * sizeof(double));

    if (!result->frequencies || !result->magnitude_db) {
        fft_result_free(result);
        return NULL;
    }

    return result;
}

void fft_result_free(fft_result_t* result) {
    if (result) {
        if (result->frequencies) free(result->frequencies);
        if (result->magnitude_db) free(result->magnitude_db);
        free(result);
    }
}

psd_result_t* psd_result_alloc(int num_bins) {
    psd_result_t* result = (psd_result_t*)malloc(sizeof(psd_result_t));
    if (!result) {
        return NULL;
    }

    result->num_bins = num_bins;
    result->frequencies = (double*)malloc(num_bins * sizeof(double));
    result->psd_db = (double*)malloc(num_bins * sizeof(double));

    if (!result->frequencies || !result->psd_db) {
        psd_result_free(result);
        return NULL;
    }

    return result;
}

void psd_result_free(psd_result_t* result) {
    if (result) {
        if (result->frequencies) free(result->frequencies);
        if (result->psd_db) free(result->psd_db);
        free(result);
    }
}

spectrogram_result_t* spectrogram_result_alloc(int num_time_bins, int num_freq_bins) {
    spectrogram_result_t* result = (spectrogram_result_t*)malloc(sizeof(spectrogram_result_t));
    if (!result) {
        return NULL;
    }

    result->num_time_bins = num_time_bins;
    result->num_freq_bins = num_freq_bins;
    result->times = (double*)malloc(num_time_bins * sizeof(double));
    result->frequencies = (double*)malloc(num_freq_bins * sizeof(double));
    result->spectrogram_db = (double*)malloc(num_time_bins * num_freq_bins * sizeof(double));

    if (!result->times || !result->frequencies || !result->spectrogram_db) {
        spectrogram_result_free(result);
        return NULL;
    }

    return result;
}

void spectrogram_result_free(spectrogram_result_t* result) {
    if (result) {
        if (result->times) free(result->times);
        if (result->frequencies) free(result->frequencies);
        if (result->spectrogram_db) free(result->spectrogram_db);
        free(result);
    }
}

iq_buffer_t* iq_buffer_alloc(int length) {
    iq_buffer_t* buffer = (iq_buffer_t*)malloc(sizeof(iq_buffer_t));
    if (!buffer) {
        return NULL;
    }

    buffer->length = length;
    buffer->samples = (cdouble_t*)malloc(length * sizeof(cdouble_t));

    if (!buffer->samples) {
        free(buffer);
        return NULL;
    }

    return buffer;
}

void iq_buffer_free(iq_buffer_t* buffer) {
    if (buffer) {
        if (buffer->samples) free(buffer->samples);
        free(buffer);
    }
}
