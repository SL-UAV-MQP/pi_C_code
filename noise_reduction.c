/**
 * @file noise_reduction.c
 * @brief Advanced Noise Reduction Algorithms
 *
 * Implements spectral and spatial noise reduction for urban RF environment.
 *
 * Performance targets: Single-channel <50ms, SNR improvement 5-15 dB
 */

#include "noise_reduction.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define SPEED_OF_LIGHT 299792458.0

/* ============================================================================
 * Configuration Helpers
 * ============================================================================ */

void noise_reduction_config_default(noise_reduction_config_t* config) {
    config->alpha = NR_ALPHA;
    config->beta = NR_SPECTRAL_FLOOR;
    config->noise_estimation_frames = 10;
    config->wiener_gain_min = NR_WIENER_GAIN_MIN;
    config->num_elements = 6;
    config->null_depth = 30.0;
    config->adaptation_rate = 0.1;
    config->noise_percentile = 10.0;
    config->fft_size = NR_FFT_SIZE;
    config->overlap = NR_OVERLAP;
    config->sample_rate = 61.44e6;
}

/* ============================================================================
 * Noise Floor Estimator
 * ============================================================================ */

music_status_t noise_floor_estimator_init(
    noise_floor_estimator_t* estimator,
    const noise_reduction_config_t* config
) {
    if (!estimator || !config) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    estimator->config = *config;
    estimator->num_bins = config->fft_size / 2;
    estimator->initialized = false;

    estimator->noise_floor = (double*)calloc(estimator->num_bins, sizeof(double));
    if (!estimator->noise_floor) {
        return MUSIC_ERROR_MEMORY;
    }

    return MUSIC_SUCCESS;
}

void noise_floor_estimator_free(noise_floor_estimator_t* estimator) {
    if (estimator && estimator->noise_floor) {
        free(estimator->noise_floor);
        estimator->noise_floor = NULL;
    }
}

music_status_t noise_floor_estimate(
    noise_floor_estimator_t* estimator,
    const double* spectrum,
    int num_bins,
    double* noise_floor
) {
    // estimate_noise_floor
    if (!estimator || !spectrum || !noise_floor) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Use percentile-based estimation
    double* sorted = (double*)malloc(num_bins * sizeof(double));
    if (!sorted) {
        return MUSIC_ERROR_MEMORY;
    }

    memcpy(sorted, spectrum, num_bins * sizeof(double));

    // Sort for percentile calculation
    for (int i = 1; i < num_bins; i++) {
        double key = sorted[i];
        int j = i - 1;
        while (j >= 0 && sorted[j] > key) {
            sorted[j + 1] = sorted[j];
            j--;
        }
        sorted[j + 1] = key;
    }

    int percentile_idx = (int)(num_bins * estimator->config.noise_percentile / 100.0);
    double percentile_value = sorted[percentile_idx];

    // Fill noise floor with minimum of spectrum and percentile
    for (int i = 0; i < num_bins; i++) {
        noise_floor[i] = (spectrum[i] < percentile_value) ? spectrum[i] : percentile_value;
    }

    free(sorted);
    return MUSIC_SUCCESS;
}

music_status_t noise_floor_update(
    noise_floor_estimator_t* estimator,
    const double* spectrum,
    int num_bins
) {
    if (!estimator || !spectrum) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    double* current_estimate = (double*)malloc(num_bins * sizeof(double));
    if (!current_estimate) {
        return MUSIC_ERROR_MEMORY;
    }

    music_status_t status = noise_floor_estimate(estimator, spectrum, num_bins,
                                                  current_estimate);
    if (status != MUSIC_SUCCESS) {
        free(current_estimate);
        return status;
    }

    if (!estimator->initialized) {
        memcpy(estimator->noise_floor, current_estimate, num_bins * sizeof(double));
        estimator->initialized = true;
    } else {
        // Exponential moving average
        double alpha = estimator->config.adaptation_rate;
        for (int i = 0; i < num_bins; i++) {
            estimator->noise_floor[i] = (1.0 - alpha) * estimator->noise_floor[i] +
                                        alpha * current_estimate[i];
        }
    }

    free(current_estimate);
    return MUSIC_SUCCESS;
}

music_status_t noise_floor_get_snr(
    const noise_floor_estimator_t* estimator,
    const double* spectrum,
    int num_bins,
    double* snr_db
) {
    // get_snr_estimate
    if (!estimator || !spectrum || !snr_db) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    if (!estimator->initialized) {
        *snr_db = 0.0;
        return MUSIC_SUCCESS;
    }

    double signal_power = 0.0;
    double noise_power = 0.0;

    for (int i = 0; i < num_bins; i++) {
        signal_power += spectrum[i];
        noise_power += estimator->noise_floor[i];
    }

    signal_power /= num_bins;
    noise_power /= num_bins;

    *snr_db = 10.0 * log10((signal_power / (noise_power + 1e-12)) + 1e-12);

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Spectral Subtraction
 * ============================================================================ */

music_status_t spectral_subtractor_init(
    spectral_subtractor_t* subtractor,
    const noise_reduction_config_t* config
) {
    if (!subtractor || !config) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    subtractor->config = *config;

    music_status_t status = noise_floor_estimator_init(&subtractor->noise_estimator, config);
    if (status != MUSIC_SUCCESS) {
        return status;
    }

    // Create Hann window
    subtractor->window = (double*)malloc(config->fft_size * sizeof(double));
    if (!subtractor->window) {
        noise_floor_estimator_free(&subtractor->noise_estimator);
        return MUSIC_ERROR_MEMORY;
    }

    for (int i = 0; i < config->fft_size; i++) {
        subtractor->window[i] = 0.5 * (1.0 - cos(2.0 * M_PI * i / (config->fft_size - 1)));
    }

    return MUSIC_SUCCESS;
}

void spectral_subtractor_free(spectral_subtractor_t* subtractor) {
    if (subtractor) {
        noise_floor_estimator_free(&subtractor->noise_estimator);
        if (subtractor->window) {
            free(subtractor->window);
        }
    }
}

music_status_t spectral_subtractor_process(
    spectral_subtractor_t* subtractor,
    const cdouble_t* iq_samples,
    int num_samples,
    cdouble_t* output
) {
    // process
    // Simplified implementation - full STFT processing would go here

    if (!subtractor || !iq_samples || !output) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // For now, pass through (placeholder for full STFT implementation)
    memcpy(output, iq_samples, num_samples * sizeof(cdouble_t));

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Wiener Filter
 * ============================================================================ */

music_status_t wiener_filter_init(
    wiener_filter_t* filter,
    const noise_reduction_config_t* config
) {
    if (!filter || !config) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    filter->config = *config;

    music_status_t status = noise_floor_estimator_init(&filter->noise_estimator, config);
    if (status != MUSIC_SUCCESS) {
        return status;
    }

    filter->window = (double*)malloc(config->fft_size * sizeof(double));
    if (!filter->window) {
        noise_floor_estimator_free(&filter->noise_estimator);
        return MUSIC_ERROR_MEMORY;
    }

    for (int i = 0; i < config->fft_size; i++) {
        filter->window[i] = 0.5 * (1.0 - cos(2.0 * M_PI * i / (config->fft_size - 1)));
    }

    return MUSIC_SUCCESS;
}

void wiener_filter_free(wiener_filter_t* filter) {
    if (filter) {
        noise_floor_estimator_free(&filter->noise_estimator);
        if (filter->window) {
            free(filter->window);
        }
    }
}

music_status_t wiener_compute_gain(
    const double* signal_power,
    const double* noise_power,
    int num_bins,
    double gain_min,
    double* wiener_gain
) {
    // compute_wiener_gain
    if (!signal_power || !noise_power || !wiener_gain) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    for (int i = 0; i < num_bins; i++) {
        double snr = signal_power[i] / (noise_power[i] + 1e-12);
        double gain = snr / (snr + 1.0);
        wiener_gain[i] = (gain > gain_min) ? gain : gain_min;
    }

    return MUSIC_SUCCESS;
}

music_status_t wiener_filter_process(
    wiener_filter_t* filter,
    const cdouble_t* iq_samples,
    int num_samples,
    cdouble_t* output
) {
    // process
    // Simplified implementation - full STFT processing would go here

    if (!filter || !iq_samples || !output) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Placeholder - copy input to output
    memcpy(output, iq_samples, num_samples * sizeof(cdouble_t));

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Spatial Noise Canceller
 * ============================================================================ */

music_status_t spatial_canceller_init(
    spatial_noise_canceller_t* canceller,
    const noise_reduction_config_t* config,
    const double* antenna_positions
) {
    if (!canceller || !config) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    canceller->config = *config;
    canceller->num_elements = config->num_elements;

    canceller->antenna_positions = (double*)malloc(config->num_elements * 3 * sizeof(double));
    if (!canceller->antenna_positions) {
        return MUSIC_ERROR_MEMORY;
    }

    if (antenna_positions) {
        memcpy(canceller->antenna_positions, antenna_positions,
               config->num_elements * 3 * sizeof(double));
    } else {
        // Default UCA positions
        double radius = 0.176;  // meters
        for (int i = 0; i < config->num_elements; i++) {
            double phi = 2.0 * M_PI * i / config->num_elements;
            canceller->antenna_positions[i * 3 + 0] = radius * cos(phi);  // x
            canceller->antenna_positions[i * 3 + 1] = radius * sin(phi);  // y
            canceller->antenna_positions[i * 3 + 2] = 0.0;                // z
        }
    }

    return MUSIC_SUCCESS;
}

void spatial_canceller_free(spatial_noise_canceller_t* canceller) {
    if (canceller && canceller->antenna_positions) {
        free(canceller->antenna_positions);
        canceller->antenna_positions = NULL;
    }
}

music_status_t spatial_compute_null_weights(
    const spatial_noise_canceller_t* canceller,
    const double* interference_azimuths,
    int num_interferers,
    double frequency,
    double desired_azimuth,
    cdouble_t* weights
) {
    // compute_null_steering_weights
    // Simplified LCMV implementation

    if (!canceller || !weights) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Placeholder - set uniform weights
    for (int i = 0; i < canceller->num_elements; i++) {
        weights[i] = 1.0 / canceller->num_elements;
    }

    return MUSIC_SUCCESS;
}

music_status_t spatial_canceller_process(
    spatial_noise_canceller_t* canceller,
    const cdouble_t* multichannel_samples,
    int num_elements,
    int num_samples,
    const double* interference_azimuths,
    int num_interferers,
    double frequency,
    double desired_azimuth,
    cdouble_t* output
) {
    // process_multichannel
    if (!canceller || !multichannel_samples || !output) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Compute weights
    cdouble_t* weights = (cdouble_t*)malloc(num_elements * sizeof(cdouble_t));
    if (!weights) {
        return MUSIC_ERROR_MEMORY;
    }

    music_status_t status = spatial_compute_null_weights(canceller, interference_azimuths,
                                                          num_interferers, frequency,
                                                          desired_azimuth, weights);
    if (status != MUSIC_SUCCESS) {
        free(weights);
        return status;
    }

    // Apply beamforming: y = w^H * x
    for (int n = 0; n < num_samples; n++) {
        output[n] = 0.0;
        for (int m = 0; m < num_elements; m++) {
            output[n] += conj(weights[m]) * multichannel_samples[m + n * num_elements];
        }
    }

    free(weights);
    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Unified Noise Reducer
 * ============================================================================ */

music_status_t noise_reducer_init(
    noise_reducer_t* reducer,
    const noise_reduction_config_t* config
) {
    if (!reducer) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    noise_reduction_config_t default_config;
    if (!config) {
        noise_reduction_config_default(&default_config);
        config = &default_config;
    }

    reducer->config = *config;

    music_status_t status = spectral_subtractor_init(&reducer->spectral_subtractor, config);
    if (status != MUSIC_SUCCESS) {
        return status;
    }

    status = wiener_filter_init(&reducer->wiener_filter, config);
    if (status != MUSIC_SUCCESS) {
        spectral_subtractor_free(&reducer->spectral_subtractor);
        return status;
    }

    status = spatial_canceller_init(&reducer->spatial_canceller, config, NULL);
    if (status != MUSIC_SUCCESS) {
        spectral_subtractor_free(&reducer->spectral_subtractor);
        wiener_filter_free(&reducer->wiener_filter);
        return status;
    }

    return MUSIC_SUCCESS;
}

void noise_reducer_free(noise_reducer_t* reducer) {
    if (reducer) {
        spectral_subtractor_free(&reducer->spectral_subtractor);
        wiener_filter_free(&reducer->wiener_filter);
        spatial_canceller_free(&reducer->spatial_canceller);
    }
}

music_status_t noise_reducer_process_single(
    noise_reducer_t* reducer,
    const cdouble_t* iq_samples,
    int num_samples,
    noise_reduction_method_t method,
    cdouble_t* output,
    noise_reduction_metrics_t* metrics
) {
    // process_single_channel
    if (!reducer || !iq_samples || !output || !metrics) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // Measure input SNR
    double input_snr;
    noise_reducer_estimate_snr(reducer, iq_samples, num_samples, &input_snr);

    music_status_t status = MUSIC_SUCCESS;

    if (method == NR_SPECTRAL_SUBTRACTION) {
        status = spectral_subtractor_process(&reducer->spectral_subtractor,
                                             iq_samples, num_samples, output);
        strcpy(metrics->method, "SpectralSubtraction");
    } else if (method == NR_WIENER_FILTER) {
        status = wiener_filter_process(&reducer->wiener_filter,
                                       iq_samples, num_samples, output);
        strcpy(metrics->method, "WienerFilter");
    } else {
        memcpy(output, iq_samples, num_samples * sizeof(cdouble_t));
        strcpy(metrics->method, "None");
    }

    // Measure output SNR
    double output_snr;
    noise_reducer_estimate_snr(reducer, output, num_samples, &output_snr);

    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed_ms = (end.tv_sec - start.tv_sec) * 1000.0 +
                       (end.tv_nsec - start.tv_nsec) / 1000000.0;

    // Fill metrics
    metrics->input_snr_db = input_snr;
    metrics->output_snr_db = output_snr;
    metrics->snr_improvement_db = output_snr - input_snr;
    metrics->processing_time_ms = elapsed_ms;
    metrics->samples_processed = num_samples;

    return status;
}

music_status_t noise_reducer_estimate_snr(
    const noise_reducer_t* reducer,
    const cdouble_t* iq_samples,
    int num_samples,
    double* snr_db
) {
    // estimate_snr
    if (!reducer || !iq_samples || !snr_db) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Simple FFT-based SNR estimation
    int fft_size = (num_samples < reducer->config.fft_size) ?
                   num_samples : reducer->config.fft_size;

    fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fft_size);
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fft_size);

    if (!in || !out) {
        if (in) fftw_free(in);
        if (out) fftw_free(out);
        return MUSIC_ERROR_MEMORY;
    }

    // Copy samples
    for (int i = 0; i < fft_size; i++) {
        in[i][0] = creal(iq_samples[i]);
        in[i][1] = cimag(iq_samples[i]);
    }

    // FFT
    fftw_plan plan = fftw_plan_dft_1d(fft_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    // Compute power spectrum
    double max_power = 0.0;
    double* power = (double*)malloc(fft_size * sizeof(double));
    for (int i = 0; i < fft_size; i++) {
        power[i] = out[i][0] * out[i][0] + out[i][1] * out[i][1];
        if (power[i] > max_power) {
            max_power = power[i];
        }
    }

    // Estimate noise as median
    double* sorted = (double*)malloc(fft_size * sizeof(double));
    memcpy(sorted, power, fft_size * sizeof(double));

    for (int i = 1; i < fft_size; i++) {
        double key = sorted[i];
        int j = i - 1;
        while (j >= 0 && sorted[j] > key) {
            sorted[j + 1] = sorted[j];
            j--;
        }
        sorted[j + 1] = key;
    }

    double noise_power = sorted[fft_size / 2];

    *snr_db = 10.0 * log10((max_power / (noise_power + 1e-12)) + 1e-12);

    // Cleanup
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
    free(power);
    free(sorted);

    return MUSIC_SUCCESS;
}

music_status_t create_noise_reducer(
    noise_reducer_t* reducer,
    double sample_rate,
    int fft_size,
    double alpha,
    double adaptation_rate
) {
    // create_noise_reducer
    noise_reduction_config_t config;
    noise_reduction_config_default(&config);

    config.sample_rate = sample_rate;
    config.fft_size = fft_size;
    config.alpha = alpha;
    config.adaptation_rate = adaptation_rate;

    return noise_reducer_init(reducer, &config);
}
