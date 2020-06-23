#ifndef CONNECTED_COMPONENTS_LABELING_CLUSTER_PERC_H
#define CONNECTED_COMPONENTS_LABELING_CLUSTER_PERC_H

enum {
    LABEL = 0,
    MAX_CLUSTER_SIZE = 1,
    PERC_COUNTER = 2
};

int calcRadius(const int radius) {
    return 3 + radius*2;
}

int PercReduceChannels(vec3i pixel, vec3i middle, int radius) {
    return (abs(pixel[0]- middle[0]) <= radius &&
            abs(pixel[1] - middle[1]) <= radius &&
            abs(pixel[2] - middle[2]) <= radius) ? 1 : 0;
}

void updateLabel(const int label1, const int label2, Buffer<int> labels) {
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


int connectedComponentsIt(int actualLabel, int upLabel, int backLabel, Buffer<int> labels, Buffer<int> counter, int *newLabel) {
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

int searchIdentityRelation(Buffer<int> labels, int i) {
    while (labels[i] != i) i = labels[i];
    return i;
}

Vector2<int> compressVector(Buffer<int> vector, Buffer<int> labels, const int size) {
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

void flipLines(Buffer<int> lineUp, Buffer<int> lineDown) {
    Buffer<int> temp = lineUp;
    lineUp = lineDown;
    lineDown = temp;
}

//vec3i boxReduce(struct Image input, const int x, const int y, const int radius,
//        int * lineUp, int * lineDown, int * labels, int * counter) {

vec3i boxReduce(struct Image input, const int x, const int y, const int radius) {
    int bufferSize = radius + 1;
    int radiusSquare = radius * radius;
    int maxNumberLabels = (radiusSquare / 2) + 1;

    Buffer<int> lineUp(bufferSize);
    Buffer<int> lineDown(bufferSize);
    Buffer<int> labels(maxNumberLabels);
    Buffer<int> counter(maxNumberLabels);

    lineUp.clear();
    lineDown[0] = 0;
    counter.clear();
    labels.clear();

    int newLabel = 1;
    int percCounter = 0;
    int lim = radius / 2;

    vec3i middle = input.pixel(x, y);

    for (auto j=-lim; j <= lim; j++) {
        for (auto i=-lim; i <= lim; i++) {
            int actualIdx = i + lim + 1; // + 1 for considering the buffer padding;
            int upValue = lineUp[actualIdx];
            int backValue = lineDown[actualIdx-1];

            int actualValue = PercReduceChannels(input.pixel(x+i, y+i), middle, radius);
            percCounter += actualValue;
            lineDown[actualIdx] = connectedComponentsIt(actualValue, upValue, backValue, labels, counter, &newLabel);
        }
        flipLines(lineUp, lineDown);
    }

    vec2i cv = compressVector(counter, labels, newLabel);

    return {cv[LABEL], cv[MAX_CLUSTER_SIZE], percCounter};
}

int clusterPerc(struct Image input, const int maxRadius) {
    int numberOfRadios = (maxRadius - 1) / 2;
    Buffer<double> p(numberOfRadios);
    Buffer<double> g(numberOfRadios);
    Buffer<double> h(numberOfRadios);

    for (auto r = 0; r < numberOfRadios; r++) {
        int actualRadius = calcRadius(r);
        int actualRadiusSquare = actualRadius * actualRadius;

        printf("Kernel (%d)\n", actualRadius);

        int pTemp = 0;
        int gTemp = 0;

        int numberBoxes = (input.width - actualRadius + 1) * (input.height - actualRadius + 1);
        double *bigClusters = allocDoubleVector(numberBoxes);
        int boxIdx = 0;

        int lim = actualRadius / 2;

        // Allocating the buffer vectors.
        int bufferSize = actualRadius+1;
//        int * lineUp = allocIntVector(bufferSize);
//        int * lineDown = allocIntVector(bufferSize);
//        int maxNumberLabels = (int) (actualRadiusSquare / 2) + 1;
//        int * labels = new int[maxNumberLabels];
//        int * counter = new int[maxNumberLabels];

        for (auto y=lim; y < input.height - lim; y++ ) {
            for (auto x=lim; x < input.width - lim; x++) {
//                vec3i ccl = boxReduce(input, x, y, actualRadius, lineUp, lineDown, labels, counter);
                vec3i ccl = boxReduce(input, x, y, actualRadius);

                bigClusters[boxIdx++] = (double) ccl[MAX_CLUSTER_SIZE] / (double) actualRadiusSquare;
                pTemp += ccl[LABEL] - 1;

                if (((double) ccl[PERC_COUNTER] / (double) actualRadiusSquare) >= 0.59275) {
                    gTemp++;
                }
            }
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

    p.print();
    g.print();
    h.print();

    return 0;
}


#endif //CONNECTED_COMPONENTS_LABELING_CLUSTER_PERC_H
