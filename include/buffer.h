#ifndef CONNECTED_COMPONENTS_LABELING_BUFFER_H
#define CONNECTED_COMPONENTS_LABELING_BUFFER_H

#include <cstdlib>
#include <cstring>
#include <memory>

template <typename T>
struct Vector2 {
    Vector2() = default;
    Vector2(T x, T y) : x(x), y(y) {}
    T x;
    T y;

    T operator[] (int index) const {
        switch(index) {
            case 0: return x;
            case 1: return y;
            default: exit(0);
        }
    }
};

template <typename T>
struct Vector3 {
    // Constructors
    Vector3() = default;
    explicit Vector3(T v) : x(v), y(v), z(v) {};
    Vector3(T x, T y, T z) : x(x), y(y), z(z) {};

    // Values
    T x;
    T y;
    T z;

    T operator[](int index) const {
        switch (index) {
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
            default:
                exit(0);
        }
    }

    Vector3 operator - (const Vector3& v) const {
        return Vector3(x - v.x, y - v.y, z - v.z);
    }

    bool operator <= (const T v) const {
        return (x <= v && y <= v && z <=v);
    }
};

template <typename T>
Vector3<T> abs(Vector3<T> v) {
    return {abs(v.x), abs(v.y), abs(v.z)};
}

typedef Vector3<int> vec3i;
typedef Vector2<int> vec2i;
typedef Vector2<float> vec2;


template <typename T>
class Buffer {
public:
    Buffer() : _ptr(nullptr), _size(0) {};
    explicit Buffer(std::size_t size) : _ptr( (T *) malloc(size * sizeof(T))), _size(size) {}

    Buffer(const Buffer& other) : _ptr(other._ptr), _size(other._size) {}
//    explicit Buffer(Buffer<int> other) : _ptr(other._ptr), _size(other._size) {}

//    ~Buffer() { _free(); }

    T& operator[](int index) {
        return _ptr[index];
    }

    const T& operator[](int index) const {
        return _ptr[index];
    }

    // Copy assigment
//    Buffer& operator=(const Buffer& other) {
//        if (this == &other) return *this;
//
//        _free();
//        _ptr = other._ptr;
//        _size = other._size;
//
//        return *this;
//    }

    // Move assigment
//    Buffer& operator=(Buffer&& other) noexcept {
//        if (this != &other) {
//            _free();
//            _ptr = std::exchange(other._ptr, nullptr);
//            _size = std::exchange(other.size, 0);
//        }
//        return *this;
//    }

//    // Copy/move contructor
//    Buffer& operator=(Buffer other) noexcept {
//        std::swap(_size, other._size);
//        std::swap(_ptr, other._ptr);
//        return *this;
//    }

    void clear() {
        for (auto i=0; i<_size; i++) {
            _ptr[i] = 0;
        }
    }

    void print() const {
        for (auto i=0; i<_size; i++) {
            if (*(typeid(T).name()) == 'i') {
                printf("%d ", _ptr[i]);
            } else {
                printf("%ld ", _ptr[i]);
            }
        }
        printf("\n");
    }

    std::size_t size() {
        return _size;
    }
private:
    std::shared_ptr<T[]> _ptr;
    std::size_t _size = 0;

//    void _free() {
//        if (_size > 0 && _ptr != nullptr)
//            free(_ptr);
//        _size = 0;
//        _ptr = nullptr;
//    }
};



int * allocIntVector(const int size) {
    return (int *) malloc(size * sizeof(int));
}

unsigned short int * allocShortVector(const int size) {
    return (unsigned short int *) malloc(size * sizeof(unsigned short int));
}

double * allocDoubleVector(const int size) {
    return (double *) malloc(size * sizeof(double));
}

void freeBuffer(int * ptr) {
    free(ptr);
}

void clearBuffer(int * ptr, const int size) {
    memset(ptr, 0, size * sizeof(int));
}

void clearBuffer(unsigned short int * ptr, const int size) {
    memset(ptr, 0, size * sizeof(unsigned short int));
}

void printBuffer(const int * m, const int size) {
    for (auto i = 1; i < size; i++) {
        std::cout << m[i] << ", ";
    }
    std::cout << "\n";
}

void printBuffer(const unsigned short int * m, const int size) {
    for (auto i = 1; i < size; i++) {
        std::cout << m[i] << ", ";
    }
    std::cout << "\n";
}

void printBuffer(const double * m, const int size) {
    std::cout << "[";
    for (auto i = 0; i < size; i++) {
        std::cout << m[i] << ", ";
    }
    std::cout << "]\n";
}

double mean(const double * m, const int size) {
    double sum = 0;

    for (auto i = 0; i < size; i++) {
        sum += m[i];
    }

    auto n = (double) size;

    return sum / n;
}


double skewness(const double * m, const int size) {
    double topSum = 0;
    double bottomSum = 0;

    double median = mean(m, size);

    for (auto i = 0; i < size; i++) {
        double x = m[i] - median;
        topSum += x * x * x;
        bottomSum += x * x;
    }

    auto n = (double) size;

    double bottom = sqrt(bottomSum / n);

    return (topSum / n) / (bottom * bottom * bottom);
}

double trapz(const double * m, const int start, const int end) {
    long double sum = 0;

    for (auto i = start+1; i < end; i++) {
        sum += m[i-1] + m[i];
    }

    return (double) sum / 2;
}

int argmax(const double * m, const int size) {
    int maxIdx = 0;
    double maxValue = m[0];

    for (auto i = 1; i < size; i++) {
        if (m[i] > maxValue) {
            maxIdx = i;
            maxValue = m[i];
        }
    }

    return maxIdx;
}



#endif //CONNECTED_COMPONENTS_LABELING_BUFFER_H
