/**
 * Copyright (C) 2019, Sorbonne Universite, LIP6
 * This file is part of the FastCPA project, under the GPL v3.0 license
 * See https://www.gnu.org/licenses/gpl-3.0.en.html for license information
 * SPDX-License-Identifier: GPL-3.0-only
 * Author(s): Quentin L. Meunier
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>

#include <math.h>

#include "common.h"

#define TRACES_PREFIX "trace_example/"
#include              "trace_example/config.h"


void init_prediction_table(float_f t[256][256]) {
    for (int32_t pt = 0; pt < 256; pt += 1) {
        for (int32_t k = 0; k < 256; k += 1) {
            // a row is a key hypothesis for all PT
            t[k][pt] = (float_f) hamming_weight(sbox[k ^ pt]);
        }
    }
}


float_f compute_correlation_coeff_opt(float_f mu_0, float_f var_0, float_f mu_1, float_f var_1, float_f t0[256], float_f t1[256]) {
    float_f mu_01 = 0;
    for (int32_t i = 0; i < 256; i += 1) {
        mu_01 += (t0[i] * t1[i] - mu_01) / (i + 1);
    }

    return (mu_01 - mu_0 * mu_1) / sqrt(var_0 * var_1);
}


float_f compute_mu(float_f t[256]) {
    float_f mu = 0;
    for (int32_t i = 0; i < 256; i += 1) {
        mu += (t[i] - mu) / (i + 1);
    }
    return mu;
}


float_f compute_sigma(float_f t[256], float_f mu) {
    float_f sigma = 0;
    for (int32_t i = 0; i < 256; i += 1) {
        float_f v = square(t[i] - mu);
        sigma += (v - sigma) / (i + 1);
    }
    return sigma;
}


void init_prediction_table_mu_sigma(float_f prediction_table[256][256], float_f prediction_table_mu[256], float_f prediction_table_sigma[256]) {
    for (int32_t k = 0; k < 256; k += 1) {
        prediction_table_mu[k] = 0;
        for (int32_t pt = 0; pt < 256; pt += 1) {
            prediction_table_mu[k] += (prediction_table[k][pt] - prediction_table_mu[k]) / (pt + 1);
        }

        prediction_table_sigma[k] = 0;
        for (int32_t pt = 0; pt < 256; pt += 1) {
            float_f v = square(prediction_table[k][pt] - prediction_table_mu[k]);
            prediction_table_sigma[k] += (v - prediction_table_sigma[k]) / (pt + 1);
        }
    }
}



int main(int32_t argc, char ** argv) {
    const int32_t nb_bytes = 16;
    const int32_t start_sample = 0;
    const int32_t end_sample = FILE_SAMPLES;
    const int32_t nb_samples = end_sample - start_sample;
    
    float_f prediction_table[256][256];
    float_f prediction_table_mu[256];
    float_f prediction_table_sigma[256];

    if (end_sample > FILE_SAMPLES) {
        fprintf(stderr, "*** Error: end_sample index greater than the number of points in a trace.\n");
    }
    if (start_sample > end_sample) {
        fprintf(stderr, "*** Error: start_sample index greater than end_sample index.\n");
    }

    print_setup(start_sample, end_sample, NB_TRACES, TRACES_PREFIX, "FastCPA");

    float_f ** traces = malloc_err(sizeof(float_f *) * NB_TRACES);
    for (int32_t i = 0; i < NB_TRACES; i += 1) {
        traces[i] = malloc_err(sizeof(float_f) * nb_samples);
    }

    uint8_t ** textin = malloc_err(sizeof(uint8_t *) * NB_TRACES);
    for (int32_t i = 0; i < NB_TRACES; i += 1) {
        textin[i] = malloc_err(sizeof(uint8_t) * nb_bytes);
    }

    uint8_t knownkey[16];

    if (read_files(traces, textin, knownkey, start_sample, end_sample, FILE_SAMPLES, NB_TRACES, TRACES_PREFIX) == -1) {
        exit(-1);
    }

    print_key(knownkey);


    init_prediction_table(prediction_table);
    init_prediction_table_mu_sigma(prediction_table, prediction_table_mu, prediction_table_sigma);

    // mean power consumption
    // indexed by PT
    float_f * mean_power[256];
    // number of power values used for computing the mean in mean_power
    int32_t * num_values[256];
    for (int32_t pt = 0; pt < 256; pt += 1) {
        mean_power[pt] = malloc_err(sizeof(float_f) * nb_samples);
        num_values[pt] = malloc_err(sizeof(int32_t) * nb_samples);
    }

    uint64_t start_time = clock();
    int32_t nb_bytes_ok = 0;

    for (int32_t byte = 0; byte < nb_bytes; byte += 1) {
        // time index at which to look at

        for (int32_t pt = 0; pt < 256; pt += 1) {
            for (int32_t idx = 0; idx < nb_samples; idx += 1) {
                mean_power[pt][idx] = 0;
                num_values[pt][idx] = 0;
            }
        }

        for (int32_t i = 0; i < NB_TRACES; i += 1) {
            uint8_t pt = textin[i][byte];
            for (int32_t idx = 0; idx < nb_samples; idx += 1) {
                num_values[pt][idx] += 1;
                // Adding to the mean consumption for this PT value
                mean_power[pt][idx] += (traces[i][idx] - mean_power[pt][idx]) / num_values[pt][idx];
            }
        }

        float_f max_corr_coeff = 0;
        int32_t max_k_hyp = 0;
        int32_t max_idx = 0;
        for (int32_t idx = 0; idx < nb_samples; idx += 1) {
            float_f mp[256];
            for (int pt = 0; pt < 256; pt += 1) {
                mp[pt] = mean_power[pt][idx];
            }

            float_f mu_0 = compute_mu(mp);
            float_f sigma_0 = compute_sigma(mp, mu_0);

            for (int32_t k_hyp = 0; k_hyp < 256; k_hyp += 1) {
                float_f mu_1 = prediction_table_mu[k_hyp];
                float_f sigma_1 = prediction_table_sigma[k_hyp];
                float_f corr_coeff = compute_correlation_coeff_opt(mu_0, sigma_0, mu_1, sigma_1, mp, &prediction_table[k_hyp][0]);
                if (mabs(corr_coeff) > max_corr_coeff) {
                    max_corr_coeff = mabs(corr_coeff);
                    max_k_hyp = k_hyp;
                    max_idx = idx;
                }
            }
        }

        if (knownkey[byte] == max_k_hyp) {
          nb_bytes_ok += 1;
        }
        printf("key byte %02d : [expected 0x%.2x - found 0x%.2x] / corr_coeff = %f at sample [%04d]\n", byte, knownkey[byte], max_k_hyp, max_corr_coeff, max_idx + start_sample);
    }

    uint64_t end_time = clock();
    printf("[EXEC_TIME]: %f\n", (float_f) (end_time - start_time) / 1000000);
    printf("[NB_BYTES_OK]: %d\n", nb_bytes_ok);

    for (int32_t pt = 0; pt < 256; pt += 1) {
        free(mean_power[pt]);
        free(num_values[pt]);
    }

    for (int32_t i = 0; i < NB_TRACES; i += 1) {
        free(textin[i]);
        free(traces[i]);
    }
    free(textin);
    free(traces);

    return 0;
}

