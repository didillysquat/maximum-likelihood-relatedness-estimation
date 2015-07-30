#include <iostream>
#include <sstream>
#include <algorithm>

#include "utils.hpp"

//Split a string with the given delimiter
std::vector<std::string> split(std::string &s, char delim) {

    /*
    std::vector<std::string> elements;

    std::stringstream stream(s);
    std::string element;
    while (std::getline(stream, element, delim)) {
        elements.push_back(element);
    }

    return elements;
    */

    // implementation 2.0
    std::vector<std::string> elements2;
    //auto it = s.begin();
    size_t start = 0;
    size_t next = 0;
    while ((next = s.find_first_of(delim, start)) != std::string::npos) {
        elements2.emplace_back(s, start, (next-start));
        start = next + 1;
    }
    elements2.emplace_back(s, start, (s.length()-start));
    /*
    if ( elements != elements2 ) {
        std::cerr << "elements vectors differ\n";
        std::cerr << "sizes = " << elements.size() << ", " << elements2.size() << "\n";
    }
    */
    return elements2;
}

std::vector<StringItPair> split_it(std::string& s, char delim){
    // implementation 2.0
    std::vector<StringItPair> elements;
    auto front = s.begin();
    size_t start = 0;
    size_t next = 0;
    while ((next = s.find_first_of(delim, start)) != std::string::npos) {
        elements.emplace_back(front + start, front + next);
        start = next + 1;
    }
    elements.emplace_back(front + start, front + start + (s.length() - start));

    // Check consistency with "simple" split function
    /*
    std::vector<std::string> elements2 = split(s, delim);
    std::vector<std::string> elements3;
    for (auto e : elements) {
        elements3.emplace_back(e.begin, e.end);
    }
    if ( elements2 != elements3 ) {
        std::cerr << "elements vectors differ\n";
        std::cerr << "sizes = " << elements.size() << ", " << elements2.size() << "\n";
    }
    */
    return elements;
}

std::vector<StringItPair> split_it(StringItPair& p, char delim) {
    // implementation 2.0
    std::vector<StringItPair> elements;
    auto front = p.begin;
    auto end = p.end;

    std::string::iterator start = front;
    std::string::iterator next;
    while ((next = std::find(start, end, delim)) != end) {
        elements.emplace_back(start, next);
        start = next + 1;
    }
    elements.emplace_back(start, start + (end - start));

    // Check consistency with "simple" split function
    /*
    std::string s(p.begin, p.end);
    std::vector<std::string> elements2 = split(s, delim);
    std::vector<std::string> elements3;
    for (auto e : elements) {
        elements3.emplace_back(e.begin, e.end);
    }
    if ( elements2 != elements3 ) {
        std::cerr << "elements vectors differ\n";
        std::cerr << "sizes = " << elements.size() << ", " << elements2.size() << "\n";
    }
    */

    return elements;
}

//Print time elapsed in seconds
void print_time_elapsed(std::string desc, struct timeval* start, struct timeval* end) {

    /*
    struct timeval {
        time_t      tv_sec;
        suseconds_t tv_usec;
    }*/
    struct timeval elapsed;

    if(start->tv_usec > end->tv_usec) {
        end->tv_usec += 1000000;
        end->tv_sec--;
    }
    elapsed.tv_usec = end->tv_usec - start->tv_usec;
    elapsed.tv_sec  = end->tv_sec  - start->tv_sec;
    float time_elapsed = (elapsed.tv_sec*1000000 + elapsed.tv_usec)/1000000.f;
    std::cout << desc << " Total Time Elapsed = " << time_elapsed << std::endl;

    return;
}
