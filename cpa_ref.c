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



int main(int32_t argc, char ** argv) {
    const int32_t nb_bytes = 16;

    const int32_t start_sample = 0;
    const int32_t end_sample = FILE_SAMPLES;
    const int32_t nb_samples = end_sample - start_sample;
    
    float_f prediction_table[256][256];

    if (end_sample > FILE_SAMPLES) {
        fprintf(stderr, "*** Error: end_sample index greater than the number of points in a trace.\n");
    }
    if (start_sample > end_sample) {
        fprintf(stderr, "*** Error: start_sample index greater than end_sample index.\n");
    }

    print_setup(start_sample, end_sample, NB_TRACES, TRACES_PREFIX, "Ref. CPA");

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

    float_f ** sum_ht;
    float_f * sum_t;
    float_f * sum_h;
    float_f * sum_h_sq;
    float_f * sum_t_sq;

    sum_ht = malloc_err(sizeof(float_f *) * nb_samples);
    sum_h = malloc_err(sizeof(float_f) * 256);
    sum_h_sq = malloc_err(sizeof(float_f) * 256);
    sum_t = malloc_err(sizeof(float_f) * nb_samples);
    sum_t_sq = malloc_err(sizeof(float_f) * nb_samples);

    for (int32_t i = 0; i < nb_samples; i += 1) {
        sum_ht[i] = malloc(sizeof(float_f) * 256);
    }

    uint64_t start_time = clock();
    int32_t nb_bytes_ok = 0;

    for (int32_t byte = 0; byte < nb_bytes; byte += 1) {
        printf("### Processing byte %d\n", byte);

        for (int32_t i = 0; i < nb_samples; i += 1) {
            sum_t[i] = 0;
            sum_t_sq[i] = 0;
        }

        for (int32_t k_hyp = 0; k_hyp < 256; k_hyp += 1) {
            sum_h[k_hyp] = 0;
            sum_h_sq[k_hyp] = 0;
            for (int32_t i = 0; i < nb_samples; i += 1) {
                sum_ht[i][k_hyp] = 0;
            }
        }


        float_f max_corr_coeff = 0;
        int32_t max_k_hyp = 0;
        int32_t max_idx = 0;

        for (int32_t n = 0; n < NB_TRACES; n += 1) {
            uint8_t pt = textin[n][byte];

            // See: http://wiki.newae.com/Correlation_Power_Analysis
            for (int32_t i = 0; i < nb_samples; i += 1) {
                float_f t = traces[n][i];
                sum_t[i] += t;
                sum_t_sq[i] += square(t);
            }
            for (int32_t k_hyp = 0; k_hyp < 256; k_hyp += 1) {
                float_f h = prediction_table[k_hyp][pt];
                
                sum_h[k_hyp] += h;
                sum_h_sq[k_hyp] += square(h);
                for (int32_t i = 0; i < nb_samples; i += 1) {
                    float_f t = traces[n][i];
                    sum_ht[i][k_hyp] += t * h;
                }
            }
        }


        for (int32_t i = 0; i < nb_samples; i += 1) {
            for (int32_t k_hyp = 0; k_hyp < 256; k_hyp += 1) {
                float_f denom = sqrt((square(sum_h[k_hyp]) - (NB_TRACES) * sum_h_sq[k_hyp]) * (square(sum_t[i]) - (NB_TRACES) * sum_t_sq[i]));
                float_f corr_coeff = ((NB_TRACES) * sum_ht[i][k_hyp] - sum_h[k_hyp] * sum_t[i]) / denom;
                if (mabs(corr_coeff) > max_corr_coeff) {
                    max_corr_coeff = mabs(corr_coeff);
                    max_k_hyp = k_hyp;
                    max_idx = i;
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

    for (int32_t i = 0; i < nb_samples; i += 1) {
        free(sum_ht[i]);
    }
    free(sum_ht);
    free(sum_h);
    free(sum_h_sq);
    free(sum_t);
    free(sum_t_sq);


    for (int32_t i = 0; i < NB_TRACES; i += 1) {
        free(textin[i]);
        free(traces[i]);
    }
    free(textin);
    free(traces);

    return 0;
}

