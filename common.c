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
#include <stdbool.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>

#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include <math.h>

#include "common.h"



uint64_t clock() {
    struct timeval full_time;
    gettimeofday(&full_time, NULL);
    return (full_time.tv_usec + full_time.tv_sec * 1000000);
}


int32_t hamming_weight(int32_t n) {
    return __builtin_popcountll(n);
}


void * malloc_err(int32_t size) {
    void * ret = malloc(size);
    if (ret == NULL) {
        perror("malloc");
        exit(-1);
    }

    return ret;
}


void print_pt(uint8_t text[16]) {
    printf("[ ");
    for (int32_t i = 0; i < 15; i++) {
        printf("0x%.2x, ", text[i]);
    }
    printf("0x%.2x ]\n", text[15]);
}


void print_key(uint8_t knownkey[16]) {
    printf("Key : [ ");
    for (int32_t i = 0; i < 15; i++) {
        printf("0x%.2x, ", knownkey[i]);
    }
    printf("0x%.2x ]\n", knownkey[15]);
}


void print_setup(int32_t start_sample, int32_t end_sample, int32_t nb_traces, char * traces_prefix, char * desc) {
    int32_t nb_samples = end_sample - start_sample;
    printf("*** %s ***\n", desc);
    printf(" * Interval [%d, %d] (%d samples)\n", start_sample, end_sample - 1, nb_samples);
    printf(" * %d traces\n", nb_traces);
    printf(" * traces prefix: %s\n", traces_prefix);
}


int32_t read_files(float_f ** traces, uint8_t ** textin, uint8_t * knownkey, int32_t start_sample, int32_t end_sample, int32_t file_samples, int32_t nb_traces, char * traces_prefix) {
    int32_t fd;
    int32_t nb_samples = end_sample - start_sample;
    #define BUF_SIZE 256
    char path[BUF_SIZE];
    int32_t len = strlen(traces_prefix) + 1;

    strncpy(&path[0], traces_prefix, BUF_SIZE);
    strncat(&path[0], TEXTIN, BUF_SIZE - len);

    printf("Reading file %s\n", path);
    if ((fd = open(path, O_RDONLY)) == -1) {
        perror("open");
        return -1;
    }

    for (int32_t i = 0; i < nb_traces; i++) {
        if (read(fd, textin[i], sizeof(uint8_t) * 16) == -1) {
            perror("read");
            return -1;
        }
    }
    close(fd);

    strncpy(&path[0], traces_prefix, BUF_SIZE);
    strncat(&path[0], TRACES, BUF_SIZE - len);

    printf("Reading file %s\n", path);
    if ((fd = open(path, O_RDONLY)) == -1) {
        perror("open");
        return -1;
    }

    for (int32_t i = 0; i < nb_traces; i++) {
        lseek(fd, sizeof(float_f) * (file_samples * i + start_sample) , SEEK_SET);
        if (read(fd, traces[i], sizeof(float_f) * nb_samples) == -1) {
            perror("read");
            return -1;
        }
    }
    close(fd);

    strncpy(&path[0], traces_prefix, BUF_SIZE);
    strncat(&path[0], KEY, BUF_SIZE - len);

    printf("Reading file %s\n", path);
    if ((fd = open(path, O_RDONLY)) == -1) {
        perror("open");
        return -1;
    }

    if (read(fd, knownkey, sizeof(uint8_t) * 16) == -1) {
        perror("read");
        return -1;
    }
    close(fd);
    #undef BUF_SIZE


    return 0;
}


