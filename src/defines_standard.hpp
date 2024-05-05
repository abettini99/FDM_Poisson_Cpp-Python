#pragma once

// ------------------ //
// unsigned int types //
// ------------------ // 
/** 8-bit unsigned integer */
typedef unsigned char       u8;
/** 16-bit unsigned integer */
typedef unsigned short     u16;
/** 32-bit unsigned integer */
typedef unsigned int       u32;
/** 64-bit unsigned integer */
typedef unsigned long long u64;

// ---------------- //
// signed int types //
// ---------------- // 
/** 8-bit signed integer */
typedef signed char         i8;
/** 16-bit signed integer */
typedef signed short       i16;
/** 32-bit signed integer */
typedef signed int         i32;
/** 64-bit signed integer */
typedef signed long long   i64;


// -------------------- //
// floating point types //
// -------------------- // 
/** 32-bit floating point */
typedef float              f32;
/** 64-bit floating point */
typedef double             f64;

// ------------- //
// boolean types //
// ------------- // 
/** 32-bit boolean */
typedef int                b32;
/** 8-bit boolean */
typedef char                b8;

// ------- //
// boolean //
// ------- // 
/** true boolean (in text form) */
#define TRUE         1
#define EXIT_SUCCESS 1
/** false boolean (in text form) */
#define FALSE        0
#define EXIT_FAILURE 0

// define static assertions
#if defined(__clang__) || defined(__gcc__)
/** compile time assertion (clang, gcc) */
#define STATIC_ASSERT _Static_assert
#else
/** compile time assertion */
#define STATIC_ASSERT static_assert
#endif

// check all types are correctly size
STATIC_ASSERT(sizeof(u8)  == 1,  "u8 should be 1 byte.");
STATIC_ASSERT(sizeof(u16) == 2, "u16 should be 2 byte.");
STATIC_ASSERT(sizeof(u32) == 4, "u32 should be 4 byte.");
STATIC_ASSERT(sizeof(u64) == 8, "u64 should be 8 byte.");

STATIC_ASSERT(sizeof(i8)  == 1,  "i8 should be 1 byte.");
STATIC_ASSERT(sizeof(i16) == 2, "i16 should be 2 byte.");
STATIC_ASSERT(sizeof(i32) == 4, "i32 should be 4 byte.");
STATIC_ASSERT(sizeof(i64) == 8, "i64 should be 8 byte.");

STATIC_ASSERT(sizeof(f32) == 4, "f32 should be 4 byte.");
STATIC_ASSERT(sizeof(f64) == 8, "f64 should be 8 byte.");
