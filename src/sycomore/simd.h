#ifndef _aec30e56_9250_476a_8b0d_0981a035c57b
#define _aec30e56_9250_476a_8b0d_0981a035c57b

#include <vector>
#include <xsimd/xsimd.hpp>

#include "sycomore/sycomore_api.h"

namespace sycomore
{

namespace simd
{

template<int InstructionSet, typename T> constexpr std::size_t width();

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

template<typename T, int InstructionSet>
using Batch = xsimd::batch<T, simd::width<InstructionSet, T>()>;

/// @brief Return the CPU information given by CPUID (x86/x64 only).
SYCOMORE_API std::vector<unsigned int> get_cpu_info(
    unsigned int leaf, unsigned int subleaf=0);

/// @brief Return the instruction set (x86/x64 only).
int get_instruction_set();

template<typename T1, typename T2>
typename std::enable_if<xsimd::detail::is_simd_type<T2>::value, void>::type
load_aligned(T1 const * source, T2 & destination)
{
    destination.load_aligned(source);
}

template<typename T>
typename std::enable_if<!xsimd::detail::is_simd_type<T>::value, void>::type
load_aligned(T const * source, T & destination)
{
    destination = *source;
}

template<typename T1, typename T2>
typename std::enable_if<xsimd::detail::is_simd_type<T1>::value, void>::type
store_aligned(T1 const & source, T2 * destination)
{
    source.store_aligned(destination);
}

template<typename T>
typename std::enable_if<!xsimd::detail::is_simd_type<T>::value, void>::type
store_aligned(T const & source, T * destination)
{
    *destination = source;
}

template<typename T>
typename std::enable_if<xsimd::detail::is_simd_type<T>::value, T>::type
pow(T const & base, int exp)
{
    return xsimd::pow(base, exp);
}

template<typename T>
typename std::enable_if<!xsimd::detail::is_simd_type<T>::value, T>::type
pow(T base, int exp)
{
    return std::pow(base, exp);
}

template<typename T>
typename std::enable_if<xsimd::detail::is_simd_type<T>::value, T>::type
exp(T const & arg)
{
    return xsimd::exp(arg);
}

template<typename T>
typename std::enable_if<!xsimd::detail::is_simd_type<T>::value, T>::type
exp(T arg)
{
    return std::exp(arg);
}

template<typename T>
typename std::enable_if<xsimd::detail::is_simd_type<T>::value, T>::type
conj(T const & arg)
{
    return xsimd::conj(arg);
}

template<typename T>
typename std::enable_if<!xsimd::detail::is_simd_type<T>::value, T>::type
conj(T arg)
{
    return std::conj(arg);
}

}

}

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

#endif // _aec30e56_9250_476a_8b0d_0981a035c57b
