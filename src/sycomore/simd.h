#ifndef _aec30e56_9250_476a_8b0d_0981a035c57b
#define _aec30e56_9250_476a_8b0d_0981a035c57b

// NOTE: include required by one of the trigonometric function in xsimd
#include <array>
#include <vector>
#include <xsimd/xsimd.hpp>

#include "sycomore/sycomore_api.h"

namespace sycomore
{

#if XSIMD_VERSION_MAJOR >= 8

#define INSTRUCTION_SET_TYPE typename

using XSIMD_X86_SSE2_VERSION = xsimd::sse2;
using XSIMD_X86_AVX_VERSION = xsimd::avx;
using XSIMD_X86_AVX512_VERSION = xsimd::avx512f;

template<typename T>
using is_batch = xsimd::is_batch<T>;

using unsupported = xsimd::unsupported;

#else

#define INSTRUCTION_SET_TYPE int

template<typename T>
using is_batch = xsimd::detail::is_simd_type<T>;

int const unsupported = 0;

#endif

namespace simd
{

template<INSTRUCTION_SET_TYPE InstructionSet, typename T> constexpr std::size_t width();

template<> 
constexpr std::size_t width<XSIMD_X86_SSE2_VERSION, double>()
{
    return 2;
}

template<> 
constexpr std::size_t width<XSIMD_X86_SSE2_VERSION, std::complex<double>>()
{
    return 2;
}

template<> 
constexpr std::size_t width<XSIMD_X86_AVX_VERSION, double>()
{
    return 4;
}

template<> 
constexpr std::size_t width<XSIMD_X86_AVX_VERSION, std::complex<double>>()
{
    return 4;
}

template<> 
constexpr std::size_t width<XSIMD_X86_AVX512_VERSION, double>()
{
    return 8;
}

template<> 
constexpr std::size_t width<XSIMD_X86_AVX512_VERSION, std::complex<double>>()
{
    return 8;
}

#if XSIMD_VERSION_MAJOR >= 8
template<typename T, typename InstructionSet>
using Batch = xsimd::batch<T, InstructionSet>;
#else
template<typename T, int InstructionSet>
using Batch = xsimd::batch<T, simd::width<InstructionSet, T>()>;
#endif

/// @brief Return the CPU information given by CPUID (x86/x64 only).
SYCOMORE_API std::vector<unsigned int> get_cpu_info(
    unsigned int leaf, unsigned int subleaf=0);

/// @brief Return the instruction set (x86/x64 only).
int get_instruction_set();

template<typename T1, typename T2>
typename std::enable_if<is_batch<T2>::value, void>::type
load_aligned(T1 const * source, T2 & destination)
{
#if XSIMD_VERSION_MAJOR >= 8
    destination = destination.load_aligned(source);
#else
    destination.load_aligned(source);
#endif
}

template<typename T>
typename std::enable_if<!is_batch<T>::value, void>::type
load_aligned(T const * source, T & destination)
{
    destination = *source;
}

template<typename T1, typename T2>
typename std::enable_if<is_batch<T1>::value, void>::type
store_aligned(T1 const & source, T2 * destination)
{
    source.store_aligned(destination);
}

template<typename T>
typename std::enable_if<!is_batch<T>::value, void>::type
store_aligned(T const & source, T * destination)
{
    *destination = source;
}

template<typename T>
typename std::enable_if<is_batch<T>::value, T>::type
pow(T const & base, int exp)
{
    return xsimd::pow(base, exp);
}

template<typename T>
typename std::enable_if<!is_batch<T>::value, T>::type
pow(T base, int exp)
{
    return std::pow(base, exp);
}

template<typename T>
typename std::enable_if<is_batch<T>::value, T>::type
exp(T const & arg)
{
    return xsimd::exp(arg);
}

template<typename T>
typename std::enable_if<!is_batch<T>::value, T>::type
exp(T arg)
{
    return std::exp(arg);
}

template<typename T>
typename std::enable_if<is_batch<T>::value, T>::type
conj(T const & arg)
{
    return xsimd::conj(arg);
}

template<typename T>
typename std::enable_if<!is_batch<T>::value, T>::type
conj(T arg)
{
    return std::conj(arg);
}

}

}

#if XSIMD_VERSION_MAJOR >= 8
#define SYCOMORE_SET_API_FUNCTION(name) \
    name = &name##_d<xsimd::unsupported>; \
    if(instruction_set >= XSIMD_X86_SSE2_VERSION::version()) \
    { \
        name = &name##_d<XSIMD_X86_SSE2_VERSION>; \
    } \
    if(instruction_set >= XSIMD_X86_AVX_VERSION::version()) \
    { \
        name = &name##_d<XSIMD_X86_AVX_VERSION>; \
    } \
    if(instruction_set >= XSIMD_X86_AVX512_VERSION::version()) \
    { \
        name = &name##_d<XSIMD_X86_AVX512_VERSION>; \
    }
#else
#define SYCOMORE_SET_API_FUNCTION(name) \
    name = &name##_d<0>; \
    if(instruction_set >= XSIMD_X86_SSE2_VERSION) \
    { \
        name = &name##_d<XSIMD_X86_SSE2_VERSION>; \
    } \
    if(instruction_set >= XSIMD_X86_AVX_VERSION) \
    { \
        name = &name##_d<XSIMD_X86_AVX_VERSION>; \
    } \
    if(instruction_set >= XSIMD_X86_AVX512_VERSION) \
    { \
        name = &name##_d<XSIMD_X86_AVX512_VERSION>; \
    }
#endif

#endif // _aec30e56_9250_476a_8b0d_0981a035c57b
