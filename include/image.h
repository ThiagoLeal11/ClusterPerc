#ifndef CONNECTED_COMPONENTS_LABELING_IMAGE_H
#define CONNECTED_COMPONENTS_LABELING_IMAGE_H

#pragma once

//#define STB_IMAGE_IMPLEMENTATION
//#include "stb_image.h"

#include "buffer.h"

struct Image {
    // Image encoded with RGBRGBRGBRGB
    Image() = default;
    Image(int width, int height, int channels, unsigned char * data) :
        width(width), height(height), channels(channels), dataStream(data),
        rowLen(width*channels), columnLen(channels) {};

    int width;
    int height;
    int channels;
    unsigned char *dataStream;

    int rowLen; // width * channels
    int columnLen; // channels

    unsigned char * getPixel(int x, int y, int c) {
        return &dataStream[y*rowLen + x*columnLen + c];
    }

    unsigned char pixel(int x, int y, int c) {
        return dataStream[y*rowLen + x*columnLen + c];
    }

    vec3i pixel(int x, int y) {
        unsigned char * pixel = getPixel(x, y, 0);
        return {pixel[0], pixel[1], pixel[2]};
    }
};

struct Image loadImage(const char * filepath) {
    int width, height, channels;
    unsigned char *image = stbi_load(filepath, &width, &height, &channels, STBI_rgb);

    return {width, height, 3, image};
}


#endif //CONNECTED_COMPONENTS_LABELING_IMAGE_H
