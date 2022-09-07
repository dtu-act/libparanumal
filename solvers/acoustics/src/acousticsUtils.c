#include <random>
#include <iostream>
#include <sstream>

std::string generateUUID(int length)
{
    std::ostringstream out;

    for (int i=0; i<length; i++) {
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<std::mt19937::result_type> dist6(0,9);

        out << dist6(rng);
    }
    return out.str();
}