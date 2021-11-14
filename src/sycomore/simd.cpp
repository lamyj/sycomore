#include "simd.h"

#include <vector>

#ifdef _WIN32
#include <intrin.h>
#elif defined __GNUC__ || defined __clang__
#include <cpuid.h>
#endif

#include <xsimd/xsimd.hpp>

namespace sycomore
{

namespace simd
{

std::vector<unsigned int> get_cpu_info(unsigned int leaf, unsigned int subleaf)
{
    std::vector<unsigned int> info(4, 0);
#ifdef _WIN32
    __cpuidex(
        reinterpret_cast<int*>(info.data()), 
        static_cast<int>(leaf), static_cast<int>(subleaf));
#elif defined __GNUC__ || defined __clang__
    unsigned int * eax = info.data()+0;
    unsigned int * ebx = info.data()+1;
    unsigned int * ecx = info.data()+2;
    unsigned int * edx = info.data()+3;
    __get_cpuid_count(leaf, subleaf, eax, ebx, ecx, edx);
#else
    // FIXME: icc
#endif
    return info;
}

int get_instruction_set()
{
#if XSIMD_VERSION_MAJOR >= 8
    return xsimd::available_architectures().best;
#else
    auto info = get_cpu_info(1);
    auto const ecx = info[2];
    auto const edx = info[3];
    
    info = get_cpu_info(7, 0);
    auto const ebx = info[1];
    
    int instruction_set = 0;
    if((edx & 1<<25) != 0) { instruction_set = XSIMD_X86_SSE_VERSION; }
    if((edx & 1<<26) != 0) { instruction_set = XSIMD_X86_SSE2_VERSION; }
    if((ecx & 1<< 0) != 0) { instruction_set = XSIMD_X86_SSE2_VERSION; }
    if((ecx & 1<< 9) != 0) { instruction_set = XSIMD_X86_SSSE3_VERSION; }
    if((ecx & 1<<19) != 0) { instruction_set = XSIMD_X86_SSE4_1_VERSION; }
    if((ecx & 1<<20) != 0) { instruction_set = XSIMD_X86_SSE4_2_VERSION; }
    if((ecx & 1<<28) != 0) { instruction_set = XSIMD_X86_AVX_VERSION; }
    if((ebx & 1<< 5) != 0) { instruction_set = XSIMD_X86_AVX2_VERSION; }
    if((ebx & 1<<16) != 0) { instruction_set = XSIMD_X86_AVX512_VERSION; }
    
    return instruction_set;
#endif
}

}

}
