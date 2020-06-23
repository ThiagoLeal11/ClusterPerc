#include <iostream>
#include <cmath>
#include <cstdio>
#include <thread>
#include <chrono>
#include "omp.h"

#define STB_IMAGE_IMPLEMENTATION

#include "include/stb_image.h"
#include "include/image.h"
#include "include/buffer.h"

#include <iostream>
#include <chrono>

typedef unsigned short int vt;

class Timer
{
public:
    Timer() : beg_(clock_::now()) {}
    void reset() { beg_ = clock_::now(); }
    double elapsed() const {
        return std::chrono::duration_cast<second_>
                (clock_::now() - beg_).count(); }

private:
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > second_;
    std::chrono::time_point<clock_> beg_;
};

enum {
    LABEL = 0,
    MAX_CLUSTER_SIZE = 1,
    PERC_COUNTER = 2
};

using namespace std;

int calcRadius(const int radius) {
    return 3 + radius*2;
}

int PercReduceChannels(vec3i pixel, vec3i middle, int radius) {
    return abs(pixel - middle) <= radius;
}

void updateLabel(const int label1, const int label2, vt *labels) {
    if (label1 > label2) {
        if (labels[label1] != label2) { // If need a update.
            if (labels[label1] != label1) { // IF need update point.
                updateLabel(labels[label1], label2, labels);
            }
            labels[label1] = label2;
        }
    } else {
        if (labels[label2] != label1) {
            if (labels[label2] != label1) {
                updateLabel(labels[label2], label1, labels);
            }
            labels[label2] = label1;
        }
    }
}


int connectedComponentsIt(int actualLabel, int upLabel, int backLabel, vt * labels, vt * counter, int *newLabel) {
    if (!actualLabel) {
        return 0;
    } else if (upLabel) {
        if (backLabel && backLabel != upLabel) updateLabel(backLabel, upLabel, labels);
        counter[upLabel]++;
        return upLabel;
    } else if (backLabel) {
        counter[backLabel]++;
        return backLabel;
    } else {
        counter[*newLabel]++;
        labels[*newLabel] = *newLabel;
        return (*newLabel)++;
    }
}

int searchIdentityRelation(const vt * labels, int i) {
    while (labels[i] != i) i = labels[i];
    return i;
}

Vector2<int> compressVector(vt * vector, const vt * labels, const int size) {
    int labelsNumber = size;
    int maxClusterSize = 0;

    // Sum compress counter with mask
    for (auto k = 0; k < size; k++) {
        if (labels[k] != k) {
            vector[searchIdentityRelation(labels, k)] += vector[k];
            vector[k] = 0;
            labelsNumber--;
        }
    }

    // get the max value.
    for (auto k = 1, l = 0; k < labelsNumber; l++) {
        if (vector[l] > maxClusterSize) maxClusterSize = vector[l];
        if (vector[l]) k++;
    }

    return {labelsNumber, maxClusterSize};
}

void flipLines(vt *&lineUp, vt *&lineDown) {
    vt *temp = lineUp;
    lineUp = lineDown;
    lineDown = temp;
}

vec3i boxReduce(struct Image input, const int x, const int y, const int radius,
                vt * lineUp, vt * lineDown, vt * labels, vt * counter) {

//vec3i boxReduce(struct Image input, const int x, const int y, const int radius) {
    int bufferSize = radius + 1;
    int radiusSquare = radius * radius;
    int maxNumberLabels = (radiusSquare / 2) + 2;

//    int * lineUp = allocIntVector(bufferSize);
//    int * lineDown = allocIntVector(bufferSize);
//    int * labels = allocIntVector(maxNumberLabels);
//    int * counter = allocIntVector(maxNumberLabels);

    // Create the buffer to allow connected components cluster count.
    clearBuffer(lineUp, bufferSize);
    lineDown[0] = 0;

    clearBuffer(labels, maxNumberLabels);
    clearBuffer(counter, maxNumberLabels);
    int newLabel = 1;

    int percCounter = 0;
    int lim = radius / 2;

    vec3i middle = input.pixel(x, y);
//    printf("(%d %d)\n", x, y);

    for (auto j=-lim; j <= lim; j++) {
        for (auto i=-lim; i <= lim; i++) {
            int actualIdx = i + lim + 1; // + 1 for considering the buffer padding;
            int upValue = lineUp[actualIdx];
            int backValue = lineDown[actualIdx-1];

            vec3i pixel = input.pixel(x+i, y+j);
            int actualValue = PercReduceChannels(pixel, middle, radius);
            percCounter += actualValue;
            lineDown[actualIdx] = connectedComponentsIt(actualValue, upValue, backValue, labels, counter, &newLabel);
        }
//        if (x==55 && y==11) {
//            printBuffer(lineDown, bufferSize);
//        }
        swap(lineUp, lineDown);
//        flipLines(lineUp, lineDown);
    }

//    if (x==55 && y==11) {
//        printBuffer(labels, maxNumberLabels);
//        printBuffer(counter, maxNumberLabels);
//    }


    vec2i cv = compressVector(counter, labels, newLabel);

//    printf("<%d %d %d>\n", cv[LABEL]-1, cv[MAX_CLUSTER_SIZE], percCounter);
//    std::this_thread::sleep_for(std::chrono::milliseconds(3000));

    // Free the buffers vectors.
//    free(labels);
//    free(counter);
//    free(lineUp);
//    free(lineDown);

    return {cv[LABEL], cv[MAX_CLUSTER_SIZE], percCounter};
}

int clusterPerc(struct Image input, const int maxRadius) {
    int numberOfRadios = (maxRadius - 1) / 2;
//    double * p = allocDoubleVector(numberOfRadios);
//    double * g = allocDoubleVector(numberOfRadios);
//    double * h = allocDoubleVector(numberOfRadios);

    double * p = allocDoubleVector(numberOfRadios);
    double * g = allocDoubleVector(numberOfRadios);
    double * h = allocDoubleVector(numberOfRadios);

    for (auto r = 0; r < numberOfRadios; r++) {
        int actualRadius = calcRadius(r);
        int actualRadiusSquare = actualRadius * actualRadius;

        int width = input.width;
        int height = input.height;

        printf("Kernel (%d)\n", actualRadius);

        int pTemp = 0;
        int gTemp = 0;

        int numberBoxes = (width - actualRadius + 1) * (height - actualRadius + 1);
        double *bigClusters = allocDoubleVector(numberBoxes);

        int lim = actualRadius / 2;

        // Allocating the buffer vectors.
        int bufferSize = actualRadius+1;
        int maxNumberLabels = (int) (actualRadiusSquare / 2) + 2;
//        int comp = 0;

        #pragma omp parallel for schedule (dynamic, 16) default(none) shared(bufferSize, maxNumberLabels, actualRadius, actualRadiusSquare, width, height, lim, input, bigClusters) reduction(+:pTemp) reduction(+:gTemp)
        for (auto y=lim; y < height - lim; y++ ) {

            auto * lineUp = allocShortVector(bufferSize);
            auto * lineDown = allocShortVector(bufferSize);
            auto * labels = allocShortVector(maxNumberLabels);
            auto * counter = allocShortVector(maxNumberLabels);

            for (auto x=lim; x < input.width - lim; x++) {
//                if (x==55 && y==11) {
//                    printf("");
//                }
                vec3i ccl = boxReduce(input, x, y, actualRadius, lineUp, lineDown, labels, counter);
//                vec3i ccl = boxReduce(input, x, y, actualRadius);

                int boxIdx = (y - lim) * (width-lim-lim) + x-lim;
                bigClusters[boxIdx] = (double) ccl[MAX_CLUSTER_SIZE] / (double) actualRadiusSquare;
//                if (bigClusters[boxIdx] > 1) {
//                    printf("erro (%d %d) %d", x, y, bigClusters[boxIdx]);
//                }

                pTemp += ccl[LABEL] - 1;

                if (((double) ccl[PERC_COUNTER] / (double) actualRadiusSquare) >= 0.59275) {
                    gTemp++;
                }
            }

            free(labels);
            free(counter);
            free(lineUp);
            free(lineDown);
        }

//        printf("p: %d, g: %d\n", pTemp, gTemp);
        p[r] = (double) pTemp / (double) numberBoxes;
        g[r] = (double) gTemp / (double) numberBoxes;
        h[r] = mean(bigClusters, numberBoxes);
        free(bigClusters);
    }

    double areaCluster = trapz(p, 0, numberOfRadios);
    double areaPerc = trapz(g, 0, numberOfRadios);
    double areaMaxCluster = trapz(h, 0, numberOfRadios);

    double skewCluster = skewness(p, numberOfRadios);
    double skewPerc = skewness(g, numberOfRadios);
    double skewMaxCluster = skewness(h, numberOfRadios);

    int maxClusterIndex = argmax(p, numberOfRadios);
    int maxPercIndex = argmax(g, numberOfRadios);
    int maxMaxClusterIndex = argmax(h, numberOfRadios);
    double maxCluster = p[maxClusterIndex];
    double maxPerc = g[maxPercIndex];
    double maxMaxCluster = h[maxMaxClusterIndex];

    int half = ceil(numberOfRadios / 2.0) ;

    double areaRatioCluster = trapz(p, half, numberOfRadios) / trapz(p, 0, half);
    double areaRatioPerc = trapz(g, half, numberOfRadios) / trapz(g, 0, half);
    double areaRatioMaxCluster = trapz(h, half, numberOfRadios) / trapz(h, 0, half);

    printf("%d %d %d \n%.10f %.10f %.10f \n%.10f %.10f %.10f \n%.10f %.10f %.10f \n%.10f %.10f %.10f \n\n",
           maxClusterIndex, maxPercIndex, maxMaxClusterIndex,
           areaRatioMaxCluster, maxMaxCluster, skewMaxCluster,
           areaMaxCluster, areaRatioCluster, areaRatioPerc,
           maxCluster, maxPerc, skewCluster,
           skewPerc, areaPerc, areaCluster);

    printBuffer(p, numberOfRadios);
    printBuffer(g, numberOfRadios);
    printBuffer(h, numberOfRadios);

    free(p);
    free(g);
    free(h);

    return 0;
}

int main() {
    Timer tmr;
    struct Image image = loadImage(R"(C:\Users\thiag\CLionProjects\connected_components_labeling\images\sj-03-476_001.png)");

    tmr.reset();
    clusterPerc(image, 45);
    double t = tmr.elapsed();
    std::cout << "time:" << t << std::endl;

    return 0;
}
