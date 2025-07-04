#ifndef NEWSCF_COMMONS_HPP
#define NEWSCF_COMMONS_HPP
#include <cstddef>
#include <cstdlib>
#include <cstring>

namespace newscf {

    inline void* newscf_malloc (unsigned long bytes) {
        return malloc(bytes);
    }

    inline void newscf_free (void* ptr) {
        free(ptr);
    }

    inline void* newscf_memset(void* ptr, int c, unsigned long n) {
        return memset(ptr, c, n);
    }

}; // namespace newscf

#endif // NEWSCF_COMMONS_HPP