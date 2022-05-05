#include <iostream>
#include <fstream>

template<class T>
void csvPrint(std::ofstream& out, std::vector<T> v) {
    if (v.empty()) {
        out << std::endl;
        std::cout << std::endl;
        return;
    }

    for (std::size_t i = 0; i < v.size() - 1; ++i) {
        out << v[i] << ",";
        std::cout << v[i] << ",";
    }

    out << v[v.size() - 1] << std::endl;
    std::cout << v[v.size() - 1] << std::endl;
}

template<class T>
void csvPrint(std::ofstream& out, T t) {
    out << t << std::endl;
    std::cout << t << std::endl;
}

template<class T, class... Args>
void csvPrint(std::ofstream& out, T t, Args... args) {
    out << t << ",";
    std::cout << t << ",";

    csvPrint(out, args...);
}