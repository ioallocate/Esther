#pragma once

#include "esther_platform.hpp"

#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace esther
{
    inline constexpr std::uint64_t sign_mask = 0x8000000000000000ull;
    inline constexpr std::uint64_t exponent_mask = 0x7FF0000000000000ull;
    inline constexpr std::uint64_t fraction_mask = 0x000FFFFFFFFFFFFFull;
    inline constexpr std::uint64_t implicit_bit = 0x0010000000000000ull;
    inline constexpr std::uint64_t magnitude_mask = 0x7FFFFFFFFFFFFFFFull;
    inline constexpr std::uint64_t quiet_nan_pattern = 0x7FF8000000000000ull;
    inline constexpr double fixed_scale = 4294967296.0;
    inline constexpr std::size_t max_reproducible_count = std::size_t( 1 ) << 40;
    inline constexpr std::size_t max_folds = 16;
    inline constexpr std::size_t pairwise_leaf = 16;
    inline constexpr std::size_t sample_count = 4096;

    struct dd_t
    {
        double hi;
        double lo;
    };

    inline std::uint64_t to_bits( double x ) noexcept
    {
        return std::bit_cast< std::uint64_t >( x );
    }

    inline double from_bits( std::uint64_t bits ) noexcept
    {
        return std::bit_cast< double >( bits );
    }

    inline double absolute( double x ) noexcept
    {
        return from_bits( to_bits( x ) & magnitude_mask );
    }

    inline bool is_nan( double x ) noexcept
    {
        return ( to_bits( x ) & magnitude_mask ) > exponent_mask;
    }

    inline bool is_inf( double x ) noexcept
    {
        return ( to_bits( x ) & magnitude_mask ) == exponent_mask;
    }

    inline bool is_finite( double x ) noexcept
    {
        return ( to_bits( x ) & exponent_mask ) != exponent_mask;
    }

    inline double quiet_nan( ) noexcept
    {
        return from_bits( quiet_nan_pattern );
    }

    inline double infinity( bool negative ) noexcept
    {
        return from_bits( negative ? ( sign_mask | exponent_mask ) : exponent_mask );
    }

    inline double pow2( int e ) noexcept
    {
        return from_bits( static_cast< std::uint64_t >( e + 1023 ) << 52 );
    }

    inline double scale_pow2( double x, int e ) noexcept
    {
        const int half = e / 2;
        return ( x * pow2( half ) ) * pow2( e - half );
    }

    inline int binary_exponent( double x ) noexcept
    {
        const std::uint64_t bits = to_bits( x );
        const int biased = static_cast< int >( ( bits >> 52 ) & 0x7FF );
        if ( biased != 0 )
            return biased - 1023;

        std::uint64_t mantissa = bits & fraction_mask;
        int e = -1022;
        while ( ( mantissa & implicit_bit ) == 0 )
        {
            mantissa <<= 1;
            e--;
        }
        return e;
    }

    inline bool fp_env_is_ieee( ) noexcept
    {
        return ( fp_load( ) & fp_guard_mask ) == fp_ieee_word;
    }

    inline void fp_env_make_ieee( ) noexcept
    {
        fp_store( ( fp_load( ) & ~fp_guard_mask ) | fp_ieee_word );
    }

    class fp_env_guard_t
    {
    public:
        fp_env_guard_t( ) noexcept : saved( fp_load( ) )
        {
            fp_store( ( saved & ~fp_guard_mask ) | fp_ieee_word );
        }

        ~fp_env_guard_t( ) noexcept
        {
            fp_store( saved );
        }

        fp_env_guard_t( const fp_env_guard_t& ) = delete;
        fp_env_guard_t& operator=( const fp_env_guard_t& ) = delete;

    private:
        fp_word saved;
    };

    inline dd_t two_sum( double a, double b ) noexcept
    {
        const double s = a + b;
        const double b_virtual = s - a;
        const double a_virtual = s - b_virtual;
        return { s, ( a - a_virtual ) + ( b - b_virtual ) };
    }

    inline dd_t fast_two_sum( double a, double b ) noexcept
    {
        const double s = a + b;
        return { s, b - ( s - a ) };
    }

    inline dd_t split( double a ) noexcept
    {
        const double t = 134217729.0 * a;
        const double hi = t - ( t - a );
        return { hi, a - hi };
    }

    inline dd_t two_prod( double a, double b ) noexcept
    {
        const double p = a * b;
        const dd_t x = split( a );
        const dd_t y = split( b );
        return { p, ( ( x.hi * y.hi - p ) + x.hi * y.lo + x.lo * y.hi ) + x.lo * y.lo };
    }

    inline dd_t dd_add( dd_t a, dd_t b ) noexcept
    {
        dd_t s = two_sum( a.hi, b.hi );
        const dd_t t = two_sum( a.lo, b.lo );
        s.lo += t.hi;
        s = fast_two_sum( s.hi, s.lo );
        s.lo += t.lo;
        return fast_two_sum( s.hi, s.lo );
    }

    inline dd_t dd_sub( dd_t a, dd_t b ) noexcept
    {
        return dd_add( a, { -b.hi, -b.lo } );
    }

    inline dd_t dd_mul( dd_t a, dd_t b ) noexcept
    {
        dd_t p = two_prod( a.hi, b.hi );
        p.lo += a.hi * b.lo + a.lo * b.hi;
        return fast_two_sum( p.hi, p.lo );
    }

    inline double sum2( const double* x, std::size_t n ) noexcept
    {
        dd_t acc{ 0.0, 0.0 };
        for ( std::size_t i = 0; i < n; i++ )
        {
            const dd_t t = two_sum( acc.hi, x[ i ] );
            acc.hi = t.hi;
            acc.lo += t.lo;
        }
        return acc.hi + acc.lo;
    }

    inline double dot2( const double* x, const double* y, std::size_t n ) noexcept
    {
        dd_t acc{ 0.0, 0.0 };
        for ( std::size_t i = 0; i < n; i++ )
        {
            const dd_t p = two_prod( x[ i ], y[ i ] );
            const dd_t t = two_sum( acc.hi, p.hi );
            acc.hi = t.hi;
            acc.lo += t.lo + p.lo;
        }
        return acc.hi + acc.lo;
    }

    inline double neumaier_sum( const double* x, std::size_t n ) noexcept
    {
        double sum = 0.0;
        double carry = 0.0;
        for ( std::size_t i = 0; i < n; i++ )
        {
            const double t = sum + x[ i ];
            if ( absolute( sum ) >= absolute( x[ i ] ) )
                carry += ( sum - t ) + x[ i ];
            else
                carry += ( x[ i ] - t ) + sum;
            sum = t;
        }
        return sum + carry;
    }

    inline double pairwise_sum( const double* x, std::size_t n ) noexcept
    {
        if ( n <= pairwise_leaf )
        {
            double sum = 0.0;
            for ( std::size_t i = 0; i < n; i++ )
                sum += x[ i ];
            return sum;
        }

        const std::size_t half = n / 2;
        return pairwise_sum( x, half ) + pairwise_sum( x + half, n - half );
    }

    struct fixed_sum_t
    {
        std::uint64_t units = 0;

        void add( double x ) noexcept
        {
            units += static_cast< std::uint64_t >( static_cast< std::int64_t >( x * fixed_scale ) );
        }

        void merge( const fixed_sum_t& other ) noexcept
        {
            units += other.units;
        }

        double value( ) const noexcept
        {
            return static_cast< double >( static_cast< std::int64_t >( units ) ) / fixed_scale;
        }
    };

    struct magnitude_scan_t
    {
        std::uint64_t largest = 0;
        bool nan = false;
        bool pos_inf = false;
        bool neg_inf = false;
    };

    inline magnitude_scan_t scan_magnitudes( const double* x, std::size_t n ) noexcept
    {
        magnitude_scan_t scan;
        for ( std::size_t i = 0; i < n; i++ )
        {
            const std::uint64_t bits = to_bits( x[ i ] );
            const std::uint64_t magnitude = bits & magnitude_mask;

            if ( magnitude < exponent_mask )
            {
                if ( magnitude > scan.largest )
                    scan.largest = magnitude;
            }
            else if ( magnitude > exponent_mask )
                scan.nan = true;
            else if ( bits & sign_mask )
                scan.neg_inf = true;
            else
                scan.pos_inf = true;
        }
        return scan;
    }

    inline int fold_headroom( std::size_t n ) noexcept
    {
        const std::uint64_t target = 2 * static_cast< std::uint64_t >( n );
        int k = 1;
        while ( ( std::uint64_t( 1 ) << k ) < target )
            k++;
        return k;
    }

    inline double reproducible_sum( const double* x, std::size_t n ) noexcept
    {
        const magnitude_scan_t scan = scan_magnitudes( x, n );
        if ( scan.nan || ( scan.pos_inf && scan.neg_inf ) )
            return quiet_nan( );
        if ( scan.pos_inf )
            return infinity( false );
        if ( scan.neg_inf )
            return infinity( true );
        if ( n > max_reproducible_count )
            return quiet_nan( );
        if ( scan.largest == 0 )
            return 0.0;

        const int shift = -binary_exponent( from_bits( scan.largest ) );
        const int headroom = fold_headroom( n );
        const int width = 52 - headroom;
        const int folds = ( 106 + width - 1 ) / width;

        double sigma[ max_folds ];
        for ( int f = 0; f < folds; f++ )
            sigma[ f ] = 1.5 * pow2( 1 + headroom - f * width );

        double acc[ max_folds ];
        for ( std::size_t i = 0; i < n; i++ )
        {
            double r = scale_pow2( x[ i ], shift );
            for ( int f = 0; f < folds; f++ )
            {
                const double q = ( sigma[ f ] + r ) - sigma[ f ];
                r -= q;
                acc[ f ] = i == 0 ? q : acc[ f ] + q;
            }
        }

        double total = 0.0;
        for ( int f = folds - 1; f >= 0; f-- )
            total += acc[ f ];
        return scale_pow2( total, -shift );
    }

    inline std::uint64_t canonical_bits( double x ) noexcept
    {
        return is_nan( x ) ? quiet_nan_pattern : to_bits( x );
    }

    struct hasher_t
    {
        std::uint64_t state = 1469598103934665603ull;

        void feed_word( std::uint64_t word ) noexcept
        {
            for ( int i = 0; i < 8; i++ )
            {
                state ^= word & 0xFF;
                state *= 1099511628211ull;
                word >>= 8;
            }
        }

        void feed( double x ) noexcept
        {
            feed_word( canonical_bits( x ) );
        }

        void feed_array( const double* x, std::size_t n ) noexcept
        {
            for ( std::size_t i = 0; i < n; i++ )
                feed( x[ i ] );
        }
    };

    inline std::int64_t ordered_key( double x ) noexcept
    {
        const std::int64_t i = std::bit_cast< std::int64_t >( x );
        return i >= 0 ? i : std::numeric_limits< std::int64_t >::min( ) - i;
    }

    inline std::uint64_t ulp_distance( double a, double b ) noexcept
    {
        if ( is_nan( a ) || is_nan( b ) )
            return ( is_nan( a ) && is_nan( b ) ) ? 0 : std::numeric_limits< std::uint64_t >::max( );

        const std::int64_t ka = ordered_key( a );
        const std::int64_t kb = ordered_key( b );
        return ka >= kb ? static_cast< std::uint64_t >( ka ) - static_cast< std::uint64_t >( kb )
                        : static_cast< std::uint64_t >( kb ) - static_cast< std::uint64_t >( ka );
    }

    inline std::uint64_t splitmix64( std::uint64_t& state ) noexcept
    {
        state += 0x9E3779B97F4A7C15ull;
        std::uint64_t z = state;
        z = ( z ^ ( z >> 30 ) ) * 0xBF58476D1CE4E5B9ull;
        z = ( z ^ ( z >> 27 ) ) * 0x94D049BB133111EBull;
        return z ^ ( z >> 31 );
    }

    inline double next_sample( std::uint64_t& state ) noexcept
    {
        const std::uint64_t r = splitmix64( state );
        const std::uint64_t sign = ( r & 1 ) ? sign_mask : 0;
        const std::uint64_t exponent = 1000 + r % 41;
        return from_bits( sign | ( exponent << 52 ) | ( splitmix64( state ) & fraction_mask ) );
    }

    inline void hash_arithmetic( hasher_t& h, const double* x, std::size_t n ) noexcept
    {
        const std::size_t half = n / 2;

        fp_env_make_ieee( );
        fp_env_guard_t guard;
        fp_store( fp_load( ) );
        h.feed_word( fp_env_is_ieee( ) );

        h.feed( reproducible_sum( x, n ) );
        h.feed( reproducible_sum( x, half ) );
        h.feed( sum2( x, n ) );
        h.feed( neumaier_sum( x, n ) );
        h.feed( pairwise_sum( x, n ) );
        h.feed( dot2( x, x + half, half ) );

        fixed_sum_t first;
        fixed_sum_t second;
        for ( std::size_t i = 0; i < half; i++ )
            first.add( x[ i ] );
        for ( std::size_t i = half; i < n; i++ )
            second.add( x[ i ] );
        first.merge( second );
        h.feed( first.value( ) );

        const magnitude_scan_t scan = scan_magnitudes( x, n );
        h.feed_word( scan.largest );
        h.feed_word( scan.nan );
        h.feed_word( scan.pos_inf );
        h.feed_word( scan.neg_inf );
        h.feed_word( static_cast< std::uint64_t >( fold_headroom( n ) ) );

        h.feed_word( to_bits( x[ 0 ] ) );
        h.feed( from_bits( to_bits( x[ 1 ] ) ^ 1 ) );
        h.feed( absolute( x[ 2 ] ) );
        h.feed_word( is_nan( x[ 3 ] ) );
        h.feed_word( is_inf( x[ 3 ] ) );
        h.feed_word( is_finite( x[ 3 ] ) );
        h.feed( quiet_nan( ) );
        h.feed( infinity( true ) );
        h.feed( pow2( binary_exponent( x[ 4 ] ) ) );
        h.feed( scale_pow2( x[ 5 ], 12 ) );
        h.feed_word( static_cast< std::uint64_t >( binary_exponent( x[ 6 ] ) ) );
        h.feed_word( static_cast< std::uint64_t >( ordered_key( x[ 7 ] ) ) );
        h.feed_word( ulp_distance( x[ 8 ], x[ 9 ] ) );
        h.feed_word( canonical_bits( x[ 10 ] ) );

        const dd_t ts = two_sum( x[ 11 ], x[ 12 ] );
        const dd_t fts = fast_two_sum( x[ 13 ], x[ 14 ] );
        const dd_t sp = split( x[ 15 ] );
        const dd_t tp = two_prod( x[ 16 ], x[ 17 ] );
        const dd_t sum = dd_add( ts, fts );
        const dd_t diff = dd_sub( tp, sp );
        const dd_t prod = dd_mul( sum, diff );

        const double pieces[ 6 ] = { ts.lo, fts.lo, sp.hi, tp.lo, sum.hi, prod.lo };
        h.feed_array( pieces, 6 );
        h.feed( diff.hi );
    }

    inline std::uint64_t arithmetic_fingerprint( ) noexcept
    {
        double* samples = static_cast< double* >( alloc_memory( sample_count * sizeof( double ) ) );
        if ( samples == nullptr )
            return 0;

        std::uint64_t state = 0x5E57E4A11CEull;
        for ( std::size_t i = 0; i < sample_count; i++ )
            samples[ i ] = next_sample( state );

        hasher_t h;
        hash_arithmetic( h, samples, sample_count );
        free_memory( samples );
        return h.state;
    }

    inline double sample_in( std::uint64_t& state, std::uint64_t base_exponent, std::uint64_t span, bool negative ) noexcept
    {
        const std::uint64_t r = splitmix64( state );
        const std::uint64_t sign = ( negative && ( r & 1 ) ) ? sign_mask : 0;
        const std::uint64_t exponent = base_exponent + ( r >> 1 ) % span;
        return from_bits( sign | ( exponent << 52 ) | ( splitmix64( state ) & fraction_mask ) );
    }

    inline constexpr std::size_t ucrt_count = 64;

    using unary_fn = double ( * )( double );
    using binary_fn = double ( * )( double, double );

    inline void hash_ucrt_math( hasher_t& h ) noexcept
    {
        void* module = load_module( "ucrtbase.dll" );
        if ( module == nullptr )
            return;

        const unary_fn exp_fn = reinterpret_cast< unary_fn >( module_symbol( module, "exp" ) );
        const unary_fn log_fn = reinterpret_cast< unary_fn >( module_symbol( module, "log" ) );
        const unary_fn sin_fn = reinterpret_cast< unary_fn >( module_symbol( module, "sin" ) );
        const binary_fn pow_fn = reinterpret_cast< binary_fn >( module_symbol( module, "pow" ) );
        if ( exp_fn == nullptr || log_fn == nullptr || sin_fn == nullptr || pow_fn == nullptr )
        {
            release_module( module );
            return;
        }

        fp_env_make_ieee( );
        fp_env_guard_t guard;

        h.feed( exp_fn( -0x1.ffffffa000001p-27 ) );
        h.feed( exp_fn( -0x1.7411ffffffff6p+9 ) );

        std::uint64_t state = 0xC0FFEE15DEADull;
        for ( std::size_t i = 0; i < ucrt_count; i++ )
        {
            const double e = sample_in( state, 1020, 8, true );
            h.feed( exp_fn( e ) );

            const double l = sample_in( state, 1000, 40, false );
            h.feed( log_fn( l ) );

            const double s = sample_in( state, 1020, 20, true );
            h.feed( sin_fn( s ) );

            const double px = sample_in( state, 1020, 8, false );
            const double py = sample_in( state, 1020, 8, true );
            h.feed( pow_fn( px, py ) );
        }

        release_module( module );
    }

    inline std::uint64_t hardware_identification( ) noexcept
    {
        hasher_t h;
        h.feed_word( arithmetic_fingerprint( ) );
        hash_ucrt_math( h );
        return h.state;
    }
}
