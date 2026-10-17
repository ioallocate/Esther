#pragma once

#include <cfloat>
#include <cstddef>
#include <cstdint>
#include <intrin.h>
#include <limits>

#if defined( __FAST_MATH__ ) || defined( _M_FP_FAST )
#error "esther: fast floating-point modes break error-free transformations"
#endif

#if !defined( _M_X64 ) && !defined( _M_ARM64 )
#error "esther: supported targets are x64 and ARM64"
#endif

#if defined( FLT_EVAL_METHOD ) && FLT_EVAL_METHOD != 0
#error "esther: FLT_EVAL_METHOD must be 0"
#endif

#if defined( __clang__ )
#pragma clang fp contract( off )
#elif defined( _MSC_VER )
#pragma fp_contract( off )
#pragma float_control( precise, on )
#endif

static_assert( std::numeric_limits< double >::is_iec559, "esther: IEEE 754 binary64 required" );

extern "C" __declspec( dllimport ) void* __stdcall VirtualAlloc( void* address, std::size_t size, unsigned long type, unsigned long protect );
extern "C" __declspec( dllimport ) int __stdcall VirtualFree( void* address, std::size_t size, unsigned long type );
extern "C" __declspec( dllimport ) void* __stdcall LoadLibraryA( const char* name );
extern "C" __declspec( dllimport ) void* __stdcall GetProcAddress( void* module, const char* name );
extern "C" __declspec( dllimport ) int __stdcall FreeLibrary( void* module );

namespace esther
{
    inline constexpr unsigned long memory_commit_reserve = 0x3000ul;
    inline constexpr unsigned long memory_read_write = 0x04ul;
    inline constexpr unsigned long memory_release = 0x8000ul;

    inline void* alloc_memory( std::size_t size ) noexcept
    {
        return VirtualAlloc( nullptr, size, memory_commit_reserve, memory_read_write );
    }

    inline void free_memory( void* memory ) noexcept
    {
        VirtualFree( memory, 0, memory_release );
    }

    inline void* load_module( const char* name ) noexcept
    {
        return LoadLibraryA( name );
    }

    inline void* module_symbol( void* module, const char* name ) noexcept
    {
        return GetProcAddress( module, name );
    }

    inline void release_module( void* module ) noexcept
    {
        FreeLibrary( module );
    }

#if defined( _M_X64 )
    using fp_word = unsigned int;
    inline constexpr fp_word fp_guard_mask = 0x8000u | 0x0040u | 0x6000u | 0x1F80u;
    inline constexpr fp_word fp_ieee_word = 0x1F80u;

    inline fp_word fp_load( ) noexcept
    {
        return _mm_getcsr( );
    }

    inline void fp_store( fp_word word ) noexcept
    {
        _mm_setcsr( word );
    }
#else
    using fp_word = std::uint64_t;
    inline constexpr fp_word fp_guard_mask = 0x00C00000ull | 0x01000000ull | 0x00001F00ull;
    inline constexpr fp_word fp_ieee_word = 0;
    inline constexpr int fpcr_sysreg = 0x5A20;

    inline fp_word fp_load( ) noexcept
    {
        return static_cast< fp_word >( _ReadStatusReg( fpcr_sysreg ) );
    }

    inline void fp_store( fp_word word ) noexcept
    {
        _WriteStatusReg( fpcr_sysreg, static_cast< __int64 >( word ) );
    }
#endif
}
